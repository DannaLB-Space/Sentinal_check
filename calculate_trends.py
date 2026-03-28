import json
import math

# Constants
mu = 398600.4418
earth_radius = 6378.137 # WGS84 Semi-major axis
omega_earth = 7.292115e-5 # rad/s

input_file = "ssas6a30.b24011.e24021.DG_.sp3.001"
output_file = "trend_data.js"
threshold_km = 0.200 # 200 meters

print(f"Analyzing {input_file} for SMA jumps (> {threshold_km*1000}m)...")

try:
    with open(input_file, "r") as f:
        lines = f.readlines()
except Exception as e:
    print(f"Error reading file: {e}")
    exit(1)

raw_records = []
current_time = ""
orbit_id = 0
prev_pz = None

# Global stats tracking
min_alt_val = float('inf')
min_alt_time = ""

# Pass 1: Parsing and Orbit Segmentation
for i in range(len(lines)):
    line = lines[i].strip()
    if line.startswith("*"):
        parts = line.split()
        if len(parts) >= 6:
            current_time = f"{parts[1]}-{parts[2].zfill(2)}-{parts[3].zfill(2)} {parts[4].zfill(2)}:{parts[5].zfill(2)}"
    elif line.startswith("PL40") and i + 1 < len(lines) and lines[i+1].startswith("VL40"):
        v_line = lines[i+1]
        try:
            px, py, pz = float(line[4:18]), float(line[18:32]), float(line[32:46])
            vx_ecef, vy_ecef, vz_ecef = float(v_line[4:18])*1e-4, float(v_line[18:32])*1e-4, float(v_line[32:46])*1e-4
            
            # Identify Ascending Node crossings
            if prev_pz is not None and prev_pz < 0 and pz >= 0:
                orbit_id += 1
            prev_pz = pz
            
            # ECEF to ECI (simplified)
            vx_eci = vx_ecef - (omega_earth * py)
            vy_eci = vy_ecef + (omega_earth * px)
            vz_eci = vz_ecef
            
            r_vec = [px, py, pz]
            v_vec = [vx_eci, vy_eci, vz_eci]
            r = math.sqrt(px**2 + py**2 + pz**2)
            v = math.sqrt(vx_eci**2 + vy_eci**2 + vz_eci**2)
            
            # Osculating Semi-major axis
            energy = (v**2 / 2.0) - (mu / r)
            sma = -mu / (2.0 * energy)
            
            # Eccentricity vector e = ((v^2 - mu/r)r - (r.v)v) / mu
            rdotv = px*vx_eci + py*vy_eci + pz*vz_eci
            ex = ((v**2 - mu/r)*px - rdotv*vx_eci)/mu
            ey = ((v**2 - mu/r)*py - rdotv*vy_eci)/mu
            ez = ((v**2 - mu/r)*pz - rdotv*vz_eci)/mu
            ecc = math.sqrt(ex**2 + ey**2 + ez**2)
            
            alt = r - earth_radius
            if alt < min_alt_val:
                min_alt_val = alt
                min_alt_time = current_time
            
            raw_records.append({
                "time": current_time,
                "sma": sma,
                "ecc": ecc,
                "orbit": orbit_id,
                "alt": alt
            })
        except: continue

if not raw_records:
    print("No valid data found.")
    exit(1)

# Pass 2: Calculate Orbit Averages
orbit_groups = {}
orbit_start_times = {}
for rec in raw_records:
    oid = rec["orbit"]
    if oid not in orbit_groups:
        orbit_groups[oid] = []
        orbit_start_times[oid] = rec["time"]
    orbit_groups[oid].append(rec["sma"])

orbit_means = {oid: sum(s)/len(s) for oid, s in orbit_groups.items()}

# Pass 3: Detect Maneuvers (Mean SMA Jumps)
print("\n[MANEUVER CANDIDATES - Mean SMA Jumps > 200m]")
print(f"{'Orbit Range':<15} | {'Timestamp':<17} | {'Mean SMA Jump (m)':<20}")
print("-" * 60)

maneuver_orbits = []
for oid in range(1, orbit_id + 1):
    m1 = orbit_means.get(oid - 1)
    m2 = orbit_means.get(oid)
    if m1 and m2:
        diff_m = (m2 - m1) * 1000
        if abs(diff_m) > (threshold_km * 1000):
            print(f"Orb {oid-1} -> {oid:<4} | {orbit_start_times[oid]:<17} | {diff_m:>+18.2f} m")
            maneuver_orbits.append(oid)

# Export data
processed = []
for rec in raw_records:
    oid = rec["orbit"]
    processed.append({
        "time": rec["time"],
        "attr": {
            "sma": round(rec["sma"], 4),
            "mean_sma": round(orbit_means[oid], 4),
            "ecc": round(rec["ecc"], 8),
            "alt": round(rec["alt"], 2),
            "orbit": oid,
            "is_maneuver_orbit": oid in maneuver_orbits
        }
    })

# Downsample for UI
step = max(1, len(processed) // 1000)
output = {
    "trends": processed[::step], 
    "stats": {
        "total_orbits": orbit_id + 1, 
        "jumps_found": len(maneuver_orbits),
        "min_alt_val": round(min_alt_val, 2),
        "min_alt_time": min_alt_time
    }
}
with open(output_file, "w") as f:
    f.write(f"const trendData = {json.dumps(output)};")

print(f"\nAnalysis complete. Found {len(maneuver_orbits)} significant mean SMA shifts.")
print(f"Minimum Altitude: {min_alt_val:.2f} km at {min_alt_time}")
