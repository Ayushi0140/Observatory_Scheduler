import csv
import json
import ast
import argparse
from flask import Flask, render_template
from datetime import datetime, timedelta
from astropy.time import Time
import astropy.units as u

app = Flask(__name__)

#Longitude for Ooty Radio Telescope
OBS_LON = 76.6667 * u.deg

def load_project_metadata(filename='projects.csv'):
    projects_dict = {}
    try:
        with open(filename, mode='r', encoding='utf-8-sig') as file:
            reader = csv.DictReader(file)
            reader.fieldnames = [field.strip() for field in reader.fieldnames]

            for raw_row in reader:
                row = {key.strip().lower(): (val.strip() if val else "") for key, val in raw_row.items() if key}
                pid = row.get('id')
                if not pid:
                    continue
                
                # Parse LST ranges
                lst_str = row.get('lst_ranges', '')
                lst_ranges = []
                if lst_str and lst_str.lower() != 'none':
                    for r in lst_str.split(';'):
                        if '-' in r:
                            s, e = r.split('-')
                            lst_ranges.append((float(s.strip()), float(e.strip())))

                # Parse Strict Times for the frontend to display
                strict_str = row.get('strict_time', '')
                strict_time_display = []
                if strict_str and strict_str.lower() != 'none':
                    for window in strict_str.split(';'):
                        if '-' in window:
                            s, e = window.split('-')
                            strict_time_display.append(f"{s.strip()} to {e.strip()}")

                projects_dict[pid] = {
                    "id": pid,
                    "time_per_rep": float(row.get('time_per_rep', 0) or 0),
                    "repetitions": int(row.get('repetitions', 0) or 0),
                    "lst_ranges": lst_ranges,
                    "Mode": row.get('mode', 'N/A'),
                    "backend": row.get('backend', 'N/A'),
                    "contact": row.get('contact', 'N/A'),
                    "time_remaining": float(row.get('time_remaining', 0) or 0),
                    "priority": int(row.get('priority', 3) or 3),
                    "strict_time": strict_time_display,
                    "pi_name": row.get('pi_name', 'N/A')
                }
    except FileNotFoundError:
        print(f"Warning: '{filename}' not found. Metadata popups will be empty.")
    return projects_dict

def get_lst(dt_obj):
    """Converts a standard datetime object (IST) to Apparent LST in hours for ORT."""
    utc_time = dt_obj - timedelta(hours=5, minutes=30)
    t = Time(utc_time)
    return t.sidereal_time('apparent', longitude=OBS_LON).hour

def get_lst_str(lst_float):
    """Converts float LST to HH:MM format."""
    lst_float = lst_float % 24
    h = int(lst_float)
    m = int(round((lst_float - h) * 60))
    if m == 60:
        h = (h + 1) % 24
        m = 0
    return f"{h:02d}:{m:02d}"

@app.route('/')
def schedule():
    # Retrieve the filenames from Flask's config
    schedule_file = app.config.get('SCHEDULE_FILE', 'schedule_output.csv')
    metadata_file = app.config.get('METADATA_FILE', 'projects.csv')

    # --- LOAD SCHEDULE DATA ---
    observation_data = []
    start_date_str = None  
    
    try:
        with open(schedule_file, mode='r', encoding='utf-8') as file:
            reader = csv.reader(file)
            next(reader) # Skip IST Headers
            next(reader) # Skip LST Headers
            for row in reader:
                if row:
                    if start_date_str is None:
                        start_date_str = row[0]
                    observation_data.append(row[1:25])
    except FileNotFoundError:
        return f"Error: '{schedule_file}' not found." # Updated error message

    if not start_date_str:
        return f"Error: '{schedule_file}' contains no data rows." # Updated error message

    # --- CALCULATE STATS FOR DASHBOARD ---
    schedule_flat = [slot for day_row in observation_data for slot in day_row]
    
    total_slots = len(schedule_flat)
    maintenance_hours = schedule_flat.count('b')
    free_hours = schedule_flat.count('w')
    white_hours = schedule_flat.count('white_res')
    crab_hours = schedule_flat.count('CRAB')
    
    project_hours = total_slots - maintenance_hours - free_hours - white_hours - crab_hours

    # --- LOAD PROJECT METADATA FROM CSV ---
    project_metadata_dict = load_project_metadata(metadata_file)


    # --- TIME & DATE SETUP (DYNAMIC) ---
    # Convert the extracted string into a Python datetime object
    start_ist = datetime.strptime(start_date_str, "%a %d%b%Y")
    start_date_iso = start_ist.strftime("%Y-%m-%dT%H:%M:%S+05:30")
    
    # Calculate total days dynamically based on rows in CSV
    total_days = len(observation_data)
    
    days = [(start_ist + timedelta(days=i)).strftime('%a %d%b%Y') for i in range(total_days)]
    ist_hours = [f"{h:02d}" for h in range(24)]
    
    base_lst_float = get_lst(start_ist)
    base_lst = int(round(base_lst_float))
    lst_hours = [f"{(base_lst + h) % 24:02d}" for h in range(24)]
        
    exact_lst_ranges = []
    for i in range(total_days):
        day_lsts = []
        for j in range(24):
            current_dt = start_ist + timedelta(days=i, hours=j)
            lst_s = get_lst(current_dt)
            lst_e = get_lst(current_dt + timedelta(hours=1))
            day_lsts.append(f"{get_lst_str(lst_s)} - {get_lst_str(lst_e)}")
        exact_lst_ranges.append(day_lsts)
    
    return render_template('calendar.html', 
                           data=observation_data, 
                           days=days, 
                           ist_hours=ist_hours,
                           lst_hours=lst_hours,
                           exact_lst_ranges=exact_lst_ranges,
                           start_date_iso=start_date_iso,
                           project_metadata=json.dumps(project_metadata_dict),
                           total_hours=total_slots,
                           maint_hours=maintenance_hours,
                           free_hours=free_hours,
                           white_hours=white_hours,
                           crab_hours=crab_hours,
                           proj_hours=project_hours)

if __name__ == '__main__':
    # 1. Set up the argument parser
    parser = argparse.ArgumentParser(description="Run the ORT Schedule Web Dashboard")
    
    # Optional arguments with your original hardcoded files as defaults
    parser.add_argument('--schedule', type=str, default='schedule_output.csv', help='Path to the scheduled CSV file')
    parser.add_argument('--projects', type=str, default='projects.csv', help='Path to the projects metadata CSV file')
    parser.add_argument('--port', type=int, default=5001, help='Port to run the Flask app on')
    
    args = parser.parse_args()

    # 2. Store the parsed arguments inside Flask's global config dictionary
    app.config['SCHEDULE_FILE'] = args.schedule
    app.config['METADATA_FILE'] = args.projects
    
    print(f"Starting server reading from '{args.schedule}' and '{args.projects}'...")

    # 3. Run the app
    app.run(debug=True, port=args.port)