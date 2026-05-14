#Check if ort is open source

#Imports necessary for this code, make sure you are running the code on an environment that has these libraries installed.
import csv
import os
import json
import datetime
import math
import numpy as np
from ortools.sat.python import cp_model
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, FK5
import astropy.units as u
from astropy.utils.iers import conf
conf.auto_max_age = None

# --- OBSERVATORY CONFIGURATION ---
OBS_LAT = 11.3850 * u.deg
OBS_LON = 76.6662 * u.deg
OBS_HEIGHT = 2200 * u.m
obs_location = EarthLocation(lat=OBS_LAT, lon=OBS_LON, height=OBS_HEIGHT)

#The function to read the observatory rules from a JSON file.
def load_observatory_rules(filename="observatory_rules.json"):
    if not os.path.exists(filename): return {}
    with open(filename, 'r') as f: return json.load(f)

#This function converts date object in IST to LST
def get_lst(dt_obj):
    t = Time(dt_obj - datetime.timedelta(hours=5, minutes=30))
    return t.sidereal_time('apparent', longitude=OBS_LON).hour

#This function calculates the LST range for a calibrator source based on its RA/Dec and the desired hour angle range, for a given start date.
def get_calibrator_lst_range(start_date_str, target_ra, target_dec, ha_start, ha_end):
    source = SkyCoord(ra=target_ra, dec=target_dec, unit=(u.hourangle, u.deg), frame='icrs')
    t_start = Time(datetime.datetime.strptime(start_date_str, "%Y-%m-%d") - datetime.timedelta(hours=5, minutes=30))
    source_now = source.transform_to(FK5(equinox=t_start))
    return [((source_now.ra + ha_start * u.hourangle).wrap_at(24 * u.hourangle).hour, 
             (source_now.ra + ha_end * u.hourangle).wrap_at(24 * u.hourangle).hour)]

#This function checks if a given LST value falls within any of the specified valid LST ranges, correctly handling wrap-around at 24 hours.
def check_lst_range(lst, valid_ranges):
    for (start, end) in valid_ranges:
        if start <= end:
            if start <= lst <= end: return True
        else: 
            if lst >= start or lst <= end: return True
    return False

#This function parses a strict IST start time string and converts it to the corresponding hour index relative to the cycle start date, which is used for scheduling fixed-time tasks.
def parse_ist_start(time_str, cycle_start_dt):
    fmt = "%d %m %Y %H %M %S"
    s_dt = datetime.datetime.strptime(time_str, fmt)
    return int((s_dt - cycle_start_dt).total_seconds() / 3600)

# --- CSV INPUT / OUTPUT FUNCTIONS ---
def load_projects_from_csv(filepath):
    projects = []
    with open(filepath, mode='r', encoding='utf-8-sig') as file:
        reader = csv.DictReader(file)
        for raw_row in reader:
            row = {}
            for key, value in raw_row.items():
                if key is not None:
                    row[key.strip().lower()] = value.strip() if value is not None else ""

            if not row or not row.get('id') or row['id'] == '':
                continue

            # PARSE LST RANGES
            lst_str = row.get('lst_ranges', '')
            lst_ranges = []
            if lst_str and lst_str.lower() != 'none':
                for r in lst_str.split(';'):
                    if '-' in r:
                        s, e = r.split('-')
                        lst_ranges.append((float(s.strip()), float(e.strip())))
            
            # PARSE STRICT TIMES (Updated for the new format without dashes)
            strict_str = row.get('strict_time', '')
            strict_time = []
            if strict_str and strict_str.lower() != 'none':
                for window in strict_str.split(';'):
                    # We no longer check for a dash '-', we just grab the start time string
                    strict_time.append(window.strip())
            else:
                strict_time = -1 # Set to -1 for the OR-Tools algorithm if no strict time exists

            proj = {
                "id": row["id"],
                "time_per_rep": int(row.get("time_per_rep", 0) or 0),
                "repetitions": int(row.get("repetitions", 0) or 0),
                "lst_ranges": lst_ranges,
                "mode": row.get("mode", "N/A"),
                "backend": row.get("backend", "N/A"),
                "contact": row.get("contact", "N/A"),
                "time_remaining": float(row.get("time_remaining", 0) or 0),
                "priority": int(row.get("priority", 3) or 3),
                "strict_time": strict_time,
                
                # --- NEW COLUMNS ADDED HERE ---
                "min_gap_days": int(row.get("min_gap_days", 0) or 0), 
                "pi_name": row.get("pi_name", "N/A") 
            }
            projects.append(proj)
    return projects

def save_schedule_to_csv(schedule, start_date_str, filepath, days=14):
    start_date = datetime.datetime.strptime(start_date_str, "%Y-%m-%d")
    with open(filepath, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        ist_row = ["Date / IST Hour"] + [f"{h:02d}:00" for h in range(24)]
        writer.writerow(ist_row)
        
        lst_row = ["LST Approx"]
        for h in range(24):
            dt = start_date + datetime.timedelta(hours=h)
            lst_row.append(f"{int(get_lst(dt)):02d}:00")
        writer.writerow(lst_row)
        
        for day in range(days):
            current_date = start_date + datetime.timedelta(days=day)
            date_str = current_date.strftime("%a %d%b%Y")
            start_idx = day * 24
            end_idx = start_idx + 24
            writer.writerow([date_str] + schedule[start_idx:end_idx])


# --- THE GLOBAL OR-TOOLS OPTIMIZER ---
def generate_schedule(start_date_str, projects, days=14, slot_lengths=1, obs_string="observatory_rules_final.json", time_limit=600, relative_gap=1e-6):
    rules = load_observatory_rules(obs_string) 
    #Calculate total number of slots based on the number of days and slot length, and initialize the schedule with 'w' for white (available) slots. Also, precompute the LST for each time slot to speed up later checks.
    start_date = datetime.datetime.strptime(start_date_str, "%Y-%m-%d")
    total_slots = int(days * 24 / slot_lengths)
    schedule = ['w'] * total_slots
    times = [start_date + datetime.timedelta(hours=i*slot_lengths) for i in range(total_slots)]
    lst_array = [get_lst(t) for t in times]
    slots_per_day = int(24 / slot_lengths)
    
    # MODULAR RULE 1: Paint Fixed Downtime Blocks
    for i, t in enumerate(times):
        for rule in rules.get("fixed_blocks", []):
            if t.month in rule.get("valid_months", []) and t.weekday() in rule.get("valid_days_of_week", []):
                if rule.get("start_hour_ist", 0) <= t.hour < rule.get("end_hour_ist", 24):
                    schedule[i] = rule.get("tag", "b")

    print("\n--- INITIATING GLOBAL OR-TOOLS OPTIMIZER ---")
    model = cp_model.CpModel()
    pacing_penalties = [] 
    
    # ---------------------------------------------------------
    # PART A: USER PROJECTS
    # ---------------------------------------------------------
    X_proj = {}
    parsed_strict_windows = {} 

    for proj in projects:
        pid, dur, reps = proj['id'], proj['time_per_rep'], proj['repetitions']
        strict = proj.get('strict_time', -1)
        X_proj[pid] = {}
        
        # 1. Parse Single-String Strict Times
        if strict != -1:
            if len(strict) != reps:
                strict = -1
            else:
                parsed_strict_windows[pid] = []
                for w in strict:
                    s_hr = parse_ist_start(w, start_date)
                    e_hr = s_hr + dur # End hour dynamically calculated!
                    parsed_strict_windows[pid].append((s_hr, e_hr))
                if not parsed_strict_windows[pid]: strict = -1

        # 2. Variable Generation
        for i in range(total_slots - dur + 1):
            if not all(schedule[i+j] == 'w' for j in range(dur)): continue
            current_start_hr, current_end_hr = i * slot_lengths, (i + dur) * slot_lengths
            
            if strict != -1:
                fits = any(s_hr <= current_start_hr and current_end_hr <= e_hr for s_hr, e_hr in parsed_strict_windows[pid])
                if fits: X_proj[pid][i] = model.NewBoolVar(f"x_{pid}_start{i}")
            else:
                if all(check_lst_range(lst_array[i+j], proj['lst_ranges']) for j in range(dur)):
                    X_proj[pid][i] = model.NewBoolVar(f"x_{pid}_start{i}")

        # 3. Project Constraints
        if X_proj[pid]: 
            if strict != -1:
                for (s_hr, e_hr) in parsed_strict_windows[pid]:
                    w_vars = [X_proj[pid][s] for s in X_proj[pid] if s * slot_lengths >= s_hr and (s + dur) * slot_lengths <= e_hr]
                    if w_vars: model.AddAtMostOne(w_vars)
            else:
                model.Add(sum(X_proj[pid].values()) <= reps)
                
                # NEW: HARD Minimum Gap Constraint
                min_gap_days = proj.get("min_gap_days", 0)
                if min_gap_days > 0:
                    min_gap_hours = int(min_gap_days * 24 / slot_lengths)
                    start_times = sorted(list(X_proj[pid].keys()))
                    
                    for i in range(len(start_times)):
                        t1 = start_times[i]
                        for j in range(i + 1, len(start_times)):
                            t2 = start_times[j]
                            # If t2 is closer than the gap allows, they cannot both be true
                            if t2 - t1 < min_gap_hours:
                                model.AddImplication(X_proj[pid][t1], X_proj[pid][t2].Not())
                            else:
                                break # Time is sorted, skip checking later slots
                
                # Soft constraint of spacing
                if reps > 1:
                    ideal_gap = int(total_slots / reps)
                    start_times = sorted(list(X_proj[pid].keys()))
                    
                    for i in range(len(start_times)):
                        t1 = start_times[i]
                        for j in range(i + 1, len(start_times)):
                            t2 = start_times[j]
                            gap = t2 - t1
                            if gap < ideal_gap:
                                overlap = model.NewBoolVar(f"ov_proj_{pid}_{t1}_{t2}")
                                model.Add(overlap >= X_proj[pid][t1] + X_proj[pid][t2] - 1)
                                penalty_weight = ideal_gap - gap
                                pacing_penalties.append(overlap * penalty_weight
                                                        )
                            else:
                                break

    # ---------------------------------------------------------
    # PART B: Adding the observatory tasks to the model
    # ---------------------------------------------------------
    X_tasks = {}
    task_durations = {}

    for task in rules.get("solver_tasks", []):
        tid = task["id"]
        dur = int(task["duration_hours"] / slot_lengths)
        task_durations[tid] = dur
        X_tasks[tid] = {}
        
        timing, constraints = task.get("timing", {}), task.get("constraints", {})
        active_days = []
        
        for day in range(int(days)):
            day_t = start_date + datetime.timedelta(days=day)
            if "valid_months" in task and day_t.month not in task["valid_months"]: continue
            if "valid_days_of_week" in task and day_t.weekday() not in task["valid_days_of_week"]: continue
            
            active_days.append(day)
            X_tasks[tid][day] = {}
            day_start = int(day * slots_per_day)
            
            if timing.get("mode") == "lst_tracking":
                lst_range = get_calibrator_lst_range(start_date_str, timing["target_ra"], timing["target_dec"], timing["ha_range_start"], timing["ha_range_end"])
                for i in range(day_start, day_start + slots_per_day - dur + 1):
                    if all(schedule[i+j] == 'w' for j in range(dur)):
                        if all(check_lst_range(lst_array[i+j], lst_range) for j in range(dur)):
                            X_tasks[tid][day][i] = model.NewBoolVar(f"{tid}_d{day}_s{i}")
                            
            elif timing.get("mode") == "ist_start_times":
                for h in timing.get("allowed_start_hours", []):
                    start_idx = int(day_start + h / slot_lengths)
                    if all(schedule[start_idx+j] == 'w' for j in range(dur)):
                        X_tasks[tid][day][start_idx] = model.NewBoolVar(f"{tid}_d{day}_hr{h}")

        q_type = constraints.get("quota_type")
        target_instances = 0
        
        if q_type == "exact_per_valid_day":
            target_instances = len(active_days)
            for day in active_days:
                if day in X_tasks[tid] and X_tasks[tid][day]:
                    model.AddExactlyOne(list(X_tasks[tid][day].values()))
                    
        elif q_type == "exact_per_allowed_hour" and len(active_days) > 0:
            reps_per_slot = constraints.get("repetitions_per_slot", 0)
            allowed_hours = timing.get("allowed_start_hours", [])
            
            prop_target = int(round(reps_per_slot * (len(active_days) / days)))
            target_instances = prop_target * len(allowed_hours)
            
            for h in allowed_hours:
                h_vars = []
                for day in active_days:
                    if day in X_tasks[tid]:
                        idx = int(day * slots_per_day + h / slot_lengths)
                        if idx in X_tasks[tid][day]: h_vars.append(X_tasks[tid][day][idx])
                
                if prop_target > 0 and h_vars:
                    model.Add(sum(h_vars) == min(prop_target, len(h_vars)))

        max_per_day = constraints.get("max_per_day")
        if max_per_day is not None:
            for day in active_days:
                if day in X_tasks[tid]:
                    day_vars = list(X_tasks[tid][day].values())
                    if day_vars: model.Add(sum(day_vars) <= max_per_day)

        if constraints.get("uniformity_spacing") == True and target_instances > 1:
            ideal_gap = int(total_slots / target_instances)
            all_starts = []
            for d in X_tasks[tid]:
                for s, var in X_tasks[tid][d].items():
                    all_starts.append((s, var))
            
            all_starts.sort(key=lambda x: x[0])
            
            for i in range(len(all_starts)):
                s1, var1 = all_starts[i]
                for j in range(i + 1, len(all_starts)):
                    s2, var2 = all_starts[j]
                    gap = s2 - s1
                    if gap < ideal_gap:
                        overlap = model.NewBoolVar(f"ov_task_{tid}_{s1}_{s2}")
                        model.Add(overlap >= var1 + var2 - 1)
                        penalty_weight = ideal_gap - gap
                        pacing_penalties.append(overlap * penalty_weight)
                    else:
                        break

    # ---------------------------------------------------------
    # PART C: GLOBAL OVERLAP RESOLUTION
    # ---------------------------------------------------------
    for k in range(total_slots):
        active_in_k = []
        for pid in X_proj:
            dur = next(p['time_per_rep'] for p in projects if p['id'] == pid)
            for start in X_proj[pid]:
                if start <= k < start + dur: active_in_k.append(X_proj[pid][start])
        
        for tid in X_tasks:
            dur = task_durations[tid]
            for day in X_tasks[tid]:
                for start in X_tasks[tid][day]:
                    if start <= k < start + dur: active_in_k.append(X_tasks[tid][day][start])
                    
        if active_in_k: model.AddAtMostOne(active_in_k)

    # ---------------------------------------------------------
    # OBJECTIVE FUNCTION & SOLVER
    # ---------------------------------------------------------
    objective_terms = []
    weight_map = {1: int(1e8), 2: int(1e7), 3: int(1e6)} 
    duration_bonus, pacing_weight = int(1e3), 30

    for pid in X_proj:
        priority, dur = next(p['priority'] for p in projects if p['id'] == pid), next(p['time_per_rep'] for p in projects if p['id'] == pid)
        true_value = weight_map.get(priority, 3) + (dur * duration_bonus)
        for start in X_proj[pid]: objective_terms.append(true_value * X_proj[pid][start])
                
    model.Maximize(sum(objective_terms) - (sum(pacing_penalties) * pacing_weight))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds, solver.parameters.relative_gap_limit = time_limit, relative_gap
    
    print("Solving Matrix using CP-SAT... (600s time limit)")
    status = solver.Solve(model)

    if status not in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        print("⚠️ CRITICAL: Constraints mathematically impossible.")
        return schedule, 0

    task_hours_used = 0
    for tid in X_tasks:
        dur = task_durations[tid]
        for day in X_tasks[tid]:
            for start in X_tasks[tid][day]:
                if solver.BooleanValue(X_tasks[tid][day][start]):
                    for j in range(dur): schedule[start + j] = tid 
                    task_hours_used += dur

    unscheduled_summary = []
    for proj in projects:
        pid, dur, reps = proj['id'], proj['time_per_rep'], proj['repetitions']
        scheduled_reps = 0
        if pid in X_proj:
            for start in X_proj[pid]:
                if solver.BooleanValue(X_proj[pid][start]):
                    for j in range(dur): schedule[start + j] = pid 
                    scheduled_reps += 1
        
        missed = reps - scheduled_reps
        if missed > 0: unscheduled_summary.append({'id': pid, 'priority': proj['priority'], 'missed': missed, 'req': reps, 'sch': scheduled_reps})

    print("\n" + "="*50 + "\n UNSCHEDULED SUMMARY\n" + "="*50)
    if not unscheduled_summary: print("✅ SUCCESS: All projects scheduled!")
    else:
        for m in sorted(unscheduled_summary, key=lambda x: x['priority']): print(f"[{m['id']}] (Priority {m['priority']}) -> {m['sch']} / {m['req']} reps")

    schedule = ['white_res' if slot == 'WHITE_SLOT_SUMMER' else slot for slot in schedule]
    schedule = ['white_res' if slot == 'WHITE_SLOT_WINTER' else slot for slot in schedule]
    
    return schedule, task_hours_used


# --- EXECUTION ---
if __name__ == "__main__":
    cycle_start_date = "2026-05-14"
    days = 14
    
    print("Loading projects from 'pro_info_v5.csv'...")
    try:
        active_projects = load_projects_from_csv("pro_info_v5.csv")
    except FileNotFoundError:
        print("ERROR: Could not find 'pro_info_v5.csv'. Please check your filename.")
        exit()

    # FIX 1: Explicitly pass the correct JSON filename! We use "_" to ignore the combined task hours variable.
    final_schedule, _ = generate_schedule(cycle_start_date, active_projects, days, obs_string="observatory_rules_final.json")

    # FIX 2: Count everything directly from the final array to prevent double-counting bugs
    total_slots = len(final_schedule)
    maintenance_slots = final_schedule.count('b')
    crab_hours = final_schedule.count('CRAB')
    white_slots_secured = final_schedule.count('white_res')
    free_slots = final_schedule.count('w')
    
    scheduled_proj_slots = total_slots - maintenance_slots - free_slots - white_slots_secured - crab_hours

    print("\n--- SCHEDULING COMPLETE ---")
    print(f"Total Calendar Hours: {total_slots}")
    print(f"Maintenance ('b'):    {maintenance_slots} hours")
    print(f"Crab Pulsar:          {crab_hours} hours")
    print(f"Project Time:         {scheduled_proj_slots} hours")
    print(f"Secured White Slots:  {white_slots_secured} hours")
    
    # Save directly to the CSV that your web app reads
    save_schedule_to_csv(final_schedule, cycle_start_date, "schedule_output.csv", days)
    print("Schedule successfully saved to 'schedule_output.csv' for the web interface.")
    print(f"Total Requested Observing Time: {sum(proj['time_per_rep'] * proj['repetitions'] for proj in active_projects)} hours")