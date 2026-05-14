# ORT Automated Scheduling System

An automated, constraint-based observation scheduling system and interactive web dashboard for the Ooty Radio Telescope (ORT). 

This project utilizes Google's OR-Tools (CP-SAT solver) to generate optimal observing schedules based on dynamic project requests, specific LST (Local Sidereal Time) ranges, and observatory downtime/maintenance rules. The generated schedule is then visualized using a Flask-based interactive web calendar.

## Features
* **Algorithmic Scheduling:** Uses integer programming to resolve scheduling conflicts, respect minimum gap days, and maximize priority observation time.
* **Dynamic Constraints:** Imports constraints such as summer/winter maintenance blocks and calibrator quotas (e.g., Crab pulsar) from a JSON configuration.
* **LST Tracking:** Automatically calculates Apparent LST based on the observatory's longitude and time slots.
* **Interactive Dashboard:** A Flask-powered HTML frontend displaying a 14-day observation matrix with live tracking, tooltips, and detailed project pop-ups.

## Requirements
* Python 3.8+
* See `requirements.txt` for package dependencies.

## Installation

1. Clone the repository:
   ```bash
   git clone [https://github.com/yourusername/ort-scheduling-system.git](https://github.com/yourusername/ort-scheduling-system.git)
   cd ort-scheduling-system
