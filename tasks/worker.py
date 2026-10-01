from celery import Celery
import subprocess
from pathlib import Path
import os
import zipfile
import re
import shutil
import glob

# Configure Redis connection with dynamic fallback to 127.0.0.1 for single-container hosting
REDIS_URL = os.getenv("REDIS_URL", "redis://127.0.0.1:6379/0")

celery_app = Celery(
    "colabmd_tasks",
    broker=REDIS_URL,
    backend=REDIS_URL
)

@celery_app.task
def run_vina_docking(prot_path: str, lig_path: str, work_dir: str, ligand_code: str, ligand_smiles: str = None):
    wd = Path(work_dir)
    try:
        # Define exact output targets in the working directory
        protein_pdbqt = wd / "protein.pdbqt"
        ligand_pdbqt = wd / "ligand.pdbqt"

        # 1. Prepare Ligand PDB or SMILES -> ligand.pdbqt
        if not lig_path and ligand_smiles:
            temp_lig_pdb = wd / f"{ligand_code}.pdb"
            subprocess.run(
                ["obabel", f"-:{ligand_smiles}", "-O", str(temp_lig_pdb), "--gen3d", "-h"],
                cwd=wd, capture_output=True, check=True
            )
            lig_path = str(temp_lig_pdb)

        subprocess.run(
            ["obabel", lig_path, "-O", str(ligand_pdbqt)],
            cwd=wd, capture_output=True, check=True
        )

        # ==========================================
        # FIX: Clean metadata tags from Ligand PDBQT
        # ==========================================
        with open(ligand_pdbqt, "r") as f:
            pdbqt_lines = f.readlines()
            
        with open(ligand_pdbqt, "w") as f:
            for line in pdbqt_lines:
                if not line.startswith(("COMPND", "AUTHOR", "HEADER")):
                    f.write(line)

        # 2. Prepare Protein PDB -> protein.pdbqt
        subprocess.run(
            ["obabel", prot_path, "-O", str(protein_pdbqt), "-p", "7.4", "-xr"],
            cwd=wd, capture_output=True, check=True
        )

        # 3. Create Vina Config
        with open(wd / "config.txt", "w") as f:
            f.write(
                "receptor = protein.pdbqt\n"
                "ligand = ligand.pdbqt\n"
                "center_x = 130\ncenter_y = 110\ncenter_z = 120\n"
                "size_x = 100\nsize_y = 100\nsize_z = 100\n"
                "exhaustiveness = 4\n"
                "cpu = 1\n"
            )

        # 4. Run AutoDock Vina
        subprocess.run(
            ["vina", "--config", "config.txt", "--out", "output.pdbqt"],
            cwd=wd, capture_output=True, check=True
        )

        # 5. Extract Best Pose and name it {ligand_code}AA.pdb
        best_pose_name = f"{ligand_code}AA.pdb"
        best_pose_path = wd / best_pose_name
        subprocess.run(
            ["obabel", "output.pdbqt", "-O", str(best_pose_path), "-f", "1", "-l", "1"],
            cwd=wd, capture_output=True, check=True
        )
        
        # ==========================================
        # 6. COMBINE PROTEIN AND LIGAND POSE
        # ==========================================
        complex_temp_path = wd / "complex.pdb"
        with open(prot_path, "r") as f_prot, open(best_pose_path, "r") as f_lig, open(complex_temp_path, "w") as f_out:
            for line in f_prot:
                if not line.startswith("END"):
                    f_out.write(line)
            f_out.write("TER\n")
            for line in f_lig:
                if not line.startswith("END"):
                    f_out.write(line)
            f_out.write("END\n")

        # ==========================================
        # 7. RENAME THE COMPLEX LIGAND
        # ==========================================
        complex_final_name = f"{ligand_code}_complex.pdb"
        complex_final_path = wd / complex_final_name
        
        lig_code_padded = ligand_code.ljust(4)[:4]
        
        with open(complex_temp_path, "r", encoding="utf-8") as f_in, open(complex_final_path, "w", encoding="utf-8") as f_out:
            for line in f_in:
                if line.startswith("ATOM") and (ligand_code in line[17:21] or "UNL" in line[17:21]):
                    m_line = "HETATM" + line[6:17] + lig_code_padded + line[21:]
                    f_out.write(m_line)
                else:
                    f_out.write(line)

        # 8. Parse output.pdbqt to generate Results Table
        scores = []
        if (wd / "output.pdbqt").exists():
            with open(wd / "output.pdbqt", "r") as f:
                mode = 1
                for line in f:
                    if line.startswith("REMARK VINA RESULT:"):
                        parts = line.split()
                        scores.append({
                            "Mode": mode,
                            "Affinity (kcal/mol)": float(parts[3]),
                            "Dist from best (rmsd l.b.)": float(parts[4]),
                            "Dist from best (rmsd u.b.)": float(parts[5])
                        })
                        mode += 1

        return {
            "status": "success", 
            "message": "Docking complete! Best pose saved.",
            "result_file": str(best_pose_path),
            "complex_file": str(complex_final_path), 
            "vina_scores": scores
        }
    except subprocess.CalledProcessError as e:
        error_output = e.stderr.decode() if isinstance(e.stderr, bytes) else str(e.stderr)
        return {"status": "failed", "message": f"Error: {error_output}"}

@celery_app.task
def generate_mol2_task(work_dir: str, input_pdb: str, ligand_name: str):
    wd = Path(work_dir)
    try:
        fixed_pdb = wd / f"{ligand_name}.pdb"
        subprocess.run(["obabel", input_pdb, "-O", str(fixed_pdb), "-p", "7.4", "-h"], cwd=wd, check=True)
        
        mol2_file = wd / f"{ligand_name}.mol2"
        subprocess.run(["obabel", str(fixed_pdb), "-O", str(mol2_file)], cwd=wd, check=True)
        
        with open(mol2_file, 'r') as f:
            lines = f.readlines()
            
        for i, line in enumerate(lines):
            if line.strip() == "@MOLECULE":
                lines[i+1] = f"{ligand_name}\n"
                break
                
        with open(mol2_file, 'w') as f:
            f.writelines(lines)
            
        return {"status": "success", "result_file": str(mol2_file)}
    except subprocess.CalledProcessError as e:
        return {"status": "failed", "message": f"Error: {e.stderr}"}

@celery_app.task
def run_aa_prep_task(work_dir: str, mol2_path: str, str_path: str, ligand_name: str = "LIG"):
    """Runs CGenFF to generate All-Atomic topologies."""
    wd = Path(work_dir)
    try:
        subprocess.run(
            [
                "python", "/app/tools/cgenff_charmm2gmx.py", 
                ligand_name, mol2_path, str_path, "/app/tools/charmm36-feb2026_cgenff-5.0.ff"
            ],
            cwd=wd, capture_output=True, text=True, check=True
        )
        
        # ======================================================================
        # FIX: AUTOMATIC LONE PAIR (LP) REMOVAL & CHARGE BALANCING
        # ======================================================================
        itp_file = wd / f"{ligand_name.lower()}.itp"
        if itp_file.exists():
            with open(itp_file, 'r') as f:
                lines = f.readlines()
                
            lp_indices = set()
            index_mapping = {}
            current_idx = 1
            last_heavy_idx = None
            charge_shifts = {}
            
            # Pass 1: Identify LPs, map indices, and calculate charge shifts
            in_atoms = False
            for line in lines:
                if line.strip().startswith('[ atoms ]'):
                    in_atoms = True
                    continue
                elif in_atoms and line.strip().startswith('['):
                    in_atoms = False
                elif in_atoms and line.strip() and not line.strip().startswith(';'):
                    parts = line.split()
                    if len(parts) >= 7:
                        old_idx = parts[0]
                        atom_name = parts[4]
                        charge = float(parts[6])
                        
                        if atom_name.startswith('LP'):
                            lp_indices.add(old_idx)
                            if last_heavy_idx:
                                # Apply the LP charge to the nearest preceding heavy atom
                                charge_shifts[last_heavy_idx] = charge_shifts.get(last_heavy_idx, 0.0) + charge
                        else:
                            if not atom_name.startswith('H'):
                                last_heavy_idx = old_idx
                            index_mapping[old_idx] = str(current_idx)
                            current_idx += 1

            # Pass 2: Rewrite file without LPs, maintaining valid indices and charge
            out_lines = []
            in_section = None
            for line in lines:
                stripped = line.strip()
                if stripped.startswith('['):
                    in_section = stripped.strip('[] ').lower()
                    out_lines.append(line)
                    continue
                if not stripped or stripped.startswith(';'):
                    out_lines.append(line)
                    continue
                    
                parts = line.split()
                
                if in_section == 'atoms':
                    old_idx = parts[0]
                    if old_idx in lp_indices:
                        continue # Drop the lone pair
                        
                    parts[0] = index_mapping[old_idx]
                    parts[6] = f"{float(parts[6]) + charge_shifts.get(old_idx, 0.0):.3f}"
                    out_lines.append("\t".join(parts) + "\n")
                    
                # ADDED 'exclusions' to the array to prevent Gromacs Invalid Atomnr errors
                elif in_section in ['bonds', 'pairs', 'angles', 'dihedrals', 'impropers', 'virtual_sites2', 'virtual_sites3', 'virtual_sites4', 'exclusions']:
                    # Drop any bond/angle/virtual_site/exclusion defining the lone pair
                    if any(p in lp_indices for p in parts if p.isdigit()):
                        continue
                    # Safely renumber remaining atoms so Gromacs doesn't throw index bounds errors
                    for i in range(len(parts)):
                        if parts[i].isdigit() and parts[i] in index_mapping:
                            parts[i] = index_mapping[parts[i]]
                    out_lines.append("\t".join(parts) + "\n")
                else:
                    out_lines.append(line)
                    
            with open(itp_file, 'w') as f:
                f.writelines(out_lines)
        # ======================================================================

        zip_path = wd / f"{ligand_name}_AA_topologies.zip"
        with zipfile.ZipFile(zip_path, 'w') as zipf:
            for ext in [".itp", ".prm"]:
                target_file = wd / f"{ligand_name.lower()}{ext}"
                if target_file.exists():
                    zipf.write(target_file, arcname=target_file.name)
                    
        return {"status": "success", "result_file": str(zip_path)}
    except subprocess.CalledProcessError as e:
        return {"status": "failed", "message": f"Error: {e.stderr}"}

@celery_app.task
def run_cg_prep_task(work_dir: str, smiles: str, ligand_name: str):
    wd = Path(work_dir)
    try:
        cg_gro = f"{ligand_name}_CG.gro"
        cg_itp = f"{ligand_name}_CG.itp"
        
        subprocess.run(
            ["python", "/app/tools/cg_param_m3_fixed.py", smiles, cg_gro, cg_itp, "1"],
            cwd=wd, capture_output=True, check=True
        )
        
        # 3.2 FIX GRO: Safely replace 'MOL' with the 3-letter code in the fixed-width GRO columns
        with open(wd / cg_gro, 'r') as f:
            gro_lines = f.readlines()
        with open(wd / cg_gro, 'w') as f:
            for line in gro_lines:
                if len(line) > 20 and "MOL" in line[5:10]:
                    line = line[:5] + line[5:10].replace("MOL", ligand_name.ljust(3)) + line[10:]
                f.write(line)
                
        # 3.2 FIX ITP: Replace 'MOL' with the 3-letter code globally
        with open(wd / cg_itp, 'r') as f:
            itp_content = f.read()
        itp_content = itp_content.replace('MOL', ligand_name)
        with open(wd / cg_itp, 'w') as f:
            f.write(itp_content)
        
        zip_path = wd / f"{ligand_name}_CG_topologies.zip"
        with zipfile.ZipFile(zip_path, 'w') as zipf:
            for f in [cg_gro, cg_itp]:
                target = wd / f
                if target.exists():
                    zipf.write(target, arcname=target.name)
                    
        return {"status": "success", "result_file": str(zip_path)}
    except subprocess.CalledProcessError as e:
        return {"status": "failed", "message": f"Error: {e.stderr}"}

# ==========================================
# CGMD PREPARATION TASKS
# ==========================================

def run_martinize(job_dir: str, protein_pdb: str):
    output_top = "system.top"
    output_cg_pdb = "CG.pdb"
    
    cmd = [
        "martinize2", "-f", protein_pdb, "-o", output_top, "-x", output_cg_pdb, 
        "-p", "backbone", "-ff", "martini3001", "-elastic", 
        "-ef", "500", "-el", "0.8", "-eu", "0.9", "-ss", "C", "-maxwarn", "10" # <-- ADDED -maxwarn 10
    ]
    
    subprocess.run(cmd, cwd=job_dir, capture_output=True, text=True, check=True)
    return {"top": output_top, "cg_pdb": output_cg_pdb}

def build_simulation_box(job_dir: str, protein_cg_pdb: str, ligand_gro: str, ligand_code: str):
    cg_gro = "CG.gro"
    
    # 1. Convert Protein CG PDB to GRO
    subprocess.run(["gmx", "editconf", "-f", protein_cg_pdb, "-o", cg_gro], cwd=job_dir, check=True)
    
    # 2. Insert Ligand into the Box
    pos_file = os.path.join(job_dir, "positions.dat")
    with open(pos_file, "w") as f:
        f.write("0 0 0\n")
        
    complex_gro = f"CG_{ligand_code}_CG.gro"
    subprocess.run([
        "gmx", "insert-molecules", "-f", cg_gro, "-ci", ligand_gro, 
        "-ip", "positions.dat", "-dr", "0", "0", "0", "-rot", "none", 
        "-nmol", "1", "-radius", "0", "-scale", "0", "-box", "5", "5", "5", "-o", complex_gro
    ], cwd=job_dir, check=True)
    
    # 3. Solvate the System with insane.py -> output as system.gro
    system_gro = "system.gro"
    insane_cmd = [
        "python", "/app/tools/insane.py", 
        "-f", complex_gro, "-o", system_gro, 
        "-pbc", "cubic", "-box", "10,10,10", 
        "-salt", "0.15", "-charge", "auto", "-sol", "W"
    ]
    
    subprocess.run(insane_cmd, cwd=job_dir, capture_output=True, text=True, check=True)
    
    # 4. Post-process system.gro to rename NA -> NA+ and CL -> CL- while keeping column width intact
    system_gro_path = os.path.join(job_dir, system_gro)
    with open(system_gro_path, 'r') as f:
        gro_lines = f.readlines()
        
    processed_gro_lines = []
    solvent_counts = {'W': 0, 'NA+': 0, 'CL-': 0, 'K+': 0}
    
    for i, line in enumerate(gro_lines):
        # Skip header (lines 0, 1) and box line (last line)
        if i >= 2 and i < len(gro_lines) - 1 and len(line) > 20:
            res_name = line[5:10].strip()
            if res_name == 'NA':
                line = line[:5] + '  NA+' + line[10:]
                solvent_counts['NA+'] += 1
            elif res_name == 'CL':
                line = line[:5] + '  CL-' + line[10:]
                solvent_counts['CL-'] += 1
            elif res_name == 'W':
                solvent_counts['W'] += 1
            elif res_name == 'K':
                line = line[:5] + '   K+' + line[10:]
                solvent_counts['K+'] += 1
        processed_gro_lines.append(line)
        
    with open(system_gro_path, 'w') as f:
        f.writelines(processed_gro_lines)
        
    solvents = [f"{k:<14} {v}" for k, v in solvent_counts.items() if v > 0]
            
    return system_gro, solvents

def update_topology(job_dir: str, top_file: str, ligand_itp: str, ligand_code: str, solvents: list):
    top_path = os.path.join(job_dir, top_file)
    
    with open(top_path, 'r') as file:
        content = file.read()
        
    # Martinize2 already writes '#include "molecule_0.itp"'.
    # We just replace '#include "martini.itp"' with the updated FF includes AND the ligand include.
    header_replacements = (
        '#include "martini_v3.0.0.itp"\n'
        '#include "martini_v3.0.0_ions_v1.itp"\n'
        '#include "martini_v3.0.0_solvents_v1.itp"\n'
        f'#include "{ligand_itp}"'
    )
    
    if '#include "martini.itp"' in content:
        content = content.replace('#include "martini.itp"', header_replacements)
    else:
        content = header_replacements + "\n\n" + content

    # Build exact molecules section layout requested (no blank line between molecule_0 and ligand)
    molecules_section = "[ molecules ]\nmolecule_0    1"
    
    new_molecules_block = f"[ molecules ]\nmolecule_0    1\n{ligand_code:<14} 1"
    for solvent in solvents:
        new_molecules_block += f"\n{solvent}"
        
    if molecules_section in content:
        content = content.replace(molecules_section, new_molecules_block)
    else:
        content += f"\n\n[ molecules ]\nmolecule_0    1\n{ligand_code:<14} 1"
        for solvent in solvents:
            content += f"\n{solvent}"
            
    with open(top_path, 'w') as file:
        file.write(content)

@celery_app.task
def run_cgmd_prep_workflow(work_dir: str, protein_pdb: str, ligand_itp: str, ligand_gro: str, ligand_code: str):
    wd = Path(work_dir)
    try:
        martinize_res = run_martinize(work_dir, protein_pdb)
        system_gro, solvents = build_simulation_box(work_dir, martinize_res["cg_pdb"], ligand_gro, ligand_code)
        update_topology(work_dir, martinize_res["top"], ligand_itp, ligand_code, solvents)
        
        # Pack the requested files into a ZIP
        zip_path = wd / "cgmd_prep_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            files_to_zip = ["system.top", "system.gro", "molecule_0.itp", "CG.pdb"]
            for f in files_to_zip:
                target_file = wd / f
                if target_file.exists():
                    zipf.write(target_file, arcname=f)
        
        return {
            "status": "success",
            "message": "CGMD Preparation Complete",
            "result_file": str(zip_path), # <-- ADDED FOR DOWNLOAD
            "files": {
                "topology": martinize_res["top"],
                "cg_protein": martinize_res["cg_pdb"],
                "system": system_gro
            }
        }
    except subprocess.CalledProcessError as e:
        error_output = e.stderr if e.stderr else e.stdout
        return {"status": "failed", "message": f"Command Failed: {error_output}"}
    except Exception as e:
        return {"status": "failed", "message": str(e)}

# ==========================================
# BACKMAPPING TASKS
# ==========================================

INITRAM_SCRIPT = r"""#!/bin/bash

PROGRAM=initram.sh
VERSION=0.7-purepython
AUTHOR="Tsjerk A. Wassenaar, PhD (Updated for GMX 2019+)"

DEPENDENCIES=(backward.py gmx python3)
SDIR=$( [[ \(0 !=\){0%/*} ]] && cd ${0%/*}; pwd )

INP=
TOP=
OUT=backmapped.gro
OTP=backmapped.top
NDX=backmapped.ndx
RAW=projected.gro
BW=0-backward.gro
CG=martini
AA=gromos54a7
NP=0
KICK=0.5
KEEP=false
TRJ=false
POSRE=true

# Read steps from environment variables (Python memory feature)
EMSTEPS=${EM_STEPS:-3000}
NBSTEPS=${NB_STEPS:-1000}
MDSTEPS=${MD_STEPS:-500}
DT=0.0002,0.0005,0.001,0.002
MDP=()

while [ -n "$1" ]; do
  case $1 in
        -f)    INP=$2   ; shift 2; continue ;;
        -p)    TOP=$2   ; shift 2; continue ;;
       -po)    OTP=$2   ; shift 2; continue ;;
        -o)    OUT=$2   ; shift 2; continue ;;
     -from)    CG=$2    ; shift 2; continue ;;
       -to)    AA=$2    ; shift 2; continue ;;
     -kick)    KICK=$2  ; shift 2; continue ;;
     -keep)    KEEP=true; shift  ; continue ;;
      -nopr)    POSRE=false; shift; continue ;;
          *)    shift ;;
  esac
done

[[ -z $INP ]] && echo "FATAL ERROR: MISSING INPUT STRUCTURE FILE" && exit 1
[[ -z $TOP ]] && echo "FATAL ERROR: MISSING INPUT TOPOLOGY FILE" && exit 1

$POSRE && MDPDEF=-DPOSRES || MDPDEF=
GARBAGE=()
trash() { for item in \(@; do GARBAGE[\){#GARBAGE[@]}]=$item; done; }

echo "=========================================================="
echo " Running Initram v0.7 with Pure-Python Spatial De-Clasher"
echo "=========================================================="

GRO=$BW
B="\(SDIR/backward.py -f\)INP -raw \(RAW -o\)GRO -kick \(KICK -sol -p\)TOP -po \(OTP -n\)NDX -from \(CG -to\)AA"
echo \(B;\)B || exit 1
trash \(RAW\)GRO

# --- ZERO-DEPENDENCY PURE PYTHON SPATIAL DE-CLASHER ---
python3 - << 'EOF'
import math
import random

gro_file = '0-backward.gro'
with open(gro_file, 'r') as f:
    lines = f.readlines()

header = lines[:2]
atom_lines = lines[2:-1]
box_line = lines[-1]

coords = []
for line in atom_lines:
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    coords.append([x, y, z])

min_dist = 0.105
min_dist_sq = min_dist * min_dist

for iteration in range(120):
    grid = {}
    clashes_found = 0
    for idx, (x, y, z) in enumerate(coords):
        cell = (int(x / min_dist), int(y / min_dist), int(z / min_dist))
        if cell not in grid:
            grid[cell] = []
        grid[cell].append(idx)

    for cell, indices in grid.items():
        cx, cy, cz = cell
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor = (cx + dx, cy + dy, cz + dz)
                    if neighbor in grid:
                        for i in indices:
                            for j in grid[neighbor]:
                                if i < j:
                                    x1, y1, z1 = coords[i]
                                    x2, y2, z2 = coords[j]
                                    rx, ry, rz = x2 - x1, y2 - y1, z2 - z1
                                    d2 = rx*rx + ry*ry + rz*rz
                                    if d2 < min_dist_sq:
                                        clashes_found += 1
                                        dist = math.sqrt(d2)
                                        if dist < 1e-5:
                                            rx = random.uniform(-1, 1)
                                            ry = random.uniform(-1, 1)
                                            rz = random.uniform(-1, 1)
                                            dist = math.sqrt(rx*rx + ry*ry + rz*rz)
                                        ux, uy, uz = rx/dist, ry/dist, rz/dist
                                        overlap = (min_dist - dist) / 2.0
                                        coords[i][0] -= ux * overlap
                                        coords[i][1] -= uy * overlap
                                        coords[i][2] -= uz * overlap
                                        coords[j][0] += ux * overlap
                                        coords[j][1] += uy * overlap
                                        coords[j][2] += uz * overlap

    if clashes_found == 0:
        print(f"-> All spatial clashes resolved in {iteration} iterations!")
        break

new_atom_lines = []
for idx, line in enumerate(atom_lines):
    x, y, z = coords[idx]
    new_line = line[:20] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[44:]
    new_atom_lines.append(new_line)

with open(gro_file, 'w') as f:
    f.writelines(header + new_atom_lines + [box_line])
EOF

i=1

###############################################
## STEP 2: ENERGY MINIMIZATION PASS 1        ##
###############################################
if [ "$EMSTEPS" -gt 0 ]; then
    BASE=$i-EM
    mdp=$BASE.mdp

cat << __MDP__ > $mdp
define                    = -DFLEXIBLE
integrator                = steep
nsteps                    = $EMSTEPS
emstep                    = 0.001
emtol                     = 100.0
pbc                       = xyz
cutoff-scheme             = Verlet
nstlist                   = 10
rlist                     = 1.0
rcoulomb                  = 1.0
rvdw                      = 1.0
coulombtype               = Cut-off
vdwtype                   = Cut-off
constraints               = none
nstxout-compressed        = 0
__MDP__

    G="gmx grompp -f \(mdp -c\)GRO -n \(NDX -p\)OTP -o $BASE -maxwarn 2"
    echo \(G;\)G || exit 1
    M="gmx mdrun -deffnm \(BASE -v -nt\)NP"
    echo \(M;\)M || { echo "FATAL: EM1 failed."; exit 1; }

    trash $BASE.*
    GRO=$((i++))-EM.gro
fi

############################################
## STEP 3: ENERGY MINIMIZATION PASS 2     ##
############################################
if [ "$NBSTEPS" -gt 0 ]; then
    BASE=$i-EM
cat << __MDP__ > $BASE.mdp
define                    = -DFLEXIBLE
integrator                = steep
nsteps                    = $NBSTEPS
emstep                    = 0.005
emtol                     = 10.0
pbc                       = xyz
cutoff-scheme             = Verlet
nstlist                   = 10
rlist                     = 1.0
rcoulomb                  = 1.0
rvdw                      = 1.0
coulombtype               = Cut-off
vdwtype                   = Cut-off
constraints               = none
nstxout-compressed        = 0
__MDP__
    mdp=$BASE.mdp

    G="gmx grompp -f \(mdp -c\)GRO -n \(NDX -p\)OTP -o $BASE -maxwarn 2"
    echo \(G;\)G || exit 1
    M="gmx mdrun -deffnm \(BASE -v -nt\)NP"
    echo \(M;\)M || { echo "FATAL: EM2 failed."; exit 1; }

    trash $BASE.*
    GRO=$i-EM.gro
fi

##########################################
## STEP 4: MD EQUILIBRATION CYCLES      ##
##########################################
if [ "$MDSTEPS" -gt 0 ]; then
    ifs=$IFS
    IFS=,
    DT=($DT)
    IFS=$ifs

    for DELTA_T in ${DT[@]}; do
        BASE=\(((++i))-mdpr-\)DELTA_T
        mdp=$BASE.mdp

cat << __MDP__ > $mdp
define                    = $MDPDEF
integrator                = md
nsteps                    = $MDSTEPS
dt                        = $DELTA_T
pbc                       = xyz
cutoff-scheme             = Verlet
nstlist                   = 20
rcoulomb                  = 1.0
rlist                     = 1.0
rvdw                      = 1.0
tcoupl                    = v-rescale
ref_t                     = 300
tau_t                     = 0.1
tc_grps                   = System
gen_vel                   = yes
gen_temp                  = 300
constraints               = h-bonds
nstxout-compressed        = 0
__MDP__

      G="gmx grompp -f \(mdp -c\)GRO -r \(BW -p\)OTP -o $BASE -maxwarn 2"
      echo \(G;\)G || exit 1
      M="gmx mdrun -deffnm \(BASE -v -nt\)NP"
      echo \(M;\)M || { echo "FATAL: Step $BASE (MD) failed."; exit 1; }
      trash $BASE.*
      GRO=$BASE.gro
    done
fi

cp \(GRO\)OUT
rm -f ${GARBAGE[@]}
echo "Backmapping successfully completed!"
"""

@celery_app.task
def run_backmapping_workflow(work_dir: str, protein_pdb: str, ligand_map: str, 
                             xtc_file: str, tpr_file: str, mdp_file: str, ndx_file: str, 
                             lig_itp_file: str, lig_prm_file: str, lig_cg_itp_file: str, system_gro: str,
                             ligand_code: str, loop_from: int, loop_until: int):
    wd = Path(work_dir)
    try:
        # 1. Prepare environment tools
        shutil.copytree("/app/tools/charmm36-feb2026_cgenff-5.0.ff", wd / "charmm36-feb2026_cgenff-5.0.ff")
        shutil.copy("/app/tools/backmapping/backward.py", wd / "backward.py")
        
        # Copy the ENTIRE Mapping folder so backward.py has all standard amino acid maps
        shutil.copytree("/app/tools/backmapping/Mapping", wd / "Mapping")
        
        # Move the uploaded ligand map into the Mapping folder
        shutil.move(wd / ligand_map, wd / "Mapping" / ligand_map)
        
        # 2. Convert protein.pdb to Gromacs aa.gro and generate topol.top
        subprocess.run([
            "gmx", "pdb2gmx", "-f", protein_pdb, "-o", "aa.gro", 
            "-water", "tip3p", "-ignh", "-ff", "charmm36-feb2026_cgenff-5.0"
        ], cwd=wd, input="1\n", text=True, check=True)
        
        # 3. Manually add ligand info into topol.top
        topol_path = wd / "topol.top"
        with open(topol_path, "r") as f:
            topol_data = f.read()
            
        # Add forcefield params near the top using regex to catch local path variations (e.g. ./)
        topol_data = re.sub(
            r'(#include ".*?forcefield\.itp")',
            f'\\1\n; additional params for {ligand_code}\n#include "{lig_prm_file}"\n',
            topol_data
        )
        
        # Add ligand topology before the system directive
        itp_insertion = f'; Include {ligand_code} topology\n#include "{lig_itp_file}"\n\n[ system ]'
        topol_data = topol_data.replace('[ system ]', itp_insertion)
        
        # Append ligand to the molecules list
        topol_data += f"{ligand_code:<14} 1\n"
        
        with open(topol_path, "w") as f:
            f.write(topol_data)
            
        # 4. Write initram-v5.sh and set permissions
        initram_path = wd / "initram-v5.sh"
        with open(initram_path, "w") as f:
            f.write(INITRAM_SCRIPT)
        
        subprocess.run(["chmod", "+x", "initram-v5.sh", "backward.py"], cwd=wd, check=True)
        
        # 5. --- PARSE CORRECT LIGAND ATOM NAMES FROM CG ITP ---
        itp_atom_names = []
        with open(wd / lig_cg_itp_file, "r") as f:
            in_atoms = False
            for line in f:
                strip_line = line.strip()
                if strip_line.startswith('[') and strip_line.endswith(']'):
                    in_atoms = ('atoms' in strip_line.lower())
                    continue
                if in_atoms and strip_line and not strip_line.startswith(';'):
                    parts = strip_line.split()
                    if len(parts) >= 5 and parts[0].isdigit():
                        itp_atom_names.append(parts[4])
        
        # ---> DEFINE DIRECTORIES HERE <---
        out_xtc_dir = wd / "all_xtc"
        out_gro_dir = wd / "gro"
        out_xtc_dir.mkdir(exist_ok=True)
        out_gro_dir.mkdir(exist_ok=True)
        
        # 6. Execute looping logic natively in Python
        for i in range(loop_from, loop_until + 1):
            # Extract frame using system_gro
            subprocess.run([
                "gmx", "trjconv", "-f", xtc_file, "-s", system_gro, "-n", ndx_file, 
                "-o", f"trj_{i}ns.gro", "-dump", f"{i}000"
            ], cwd=wd, input="16\n", text=True, check=True)
            
            # --- PATCH THE EXTRACTED FRAME WITH CORRECT ATOM NAMES ---
            with open(wd / f"trj_{i}ns.gro", "r") as f:
                gro_lines = f.readlines()
                
            lig_idx = 0
            with open(wd / f"trj_{i}ns.gro", "w") as f:
                for line in gro_lines:
                    # Identify ligand lines in the .gro file
                    if len(line) > 20 and ligand_code in line[5:10]:
                        if lig_idx < len(itp_atom_names):
                            correct_name = itp_atom_names[lig_idx].rjust(5)
                            line = line[:10] + correct_name + line[15:]
                            lig_idx += 1
                    f.write(line)
            # ---------------------------------------------------------
            
            # EXPEDITE SCRIPT USING MEMORY
            env = os.environ.copy()
            if i == loop_from:
                # First frame: Run the full EM and MD equilibration loop
                env["EM_STEPS"] = "3000"
                env["NB_STEPS"] = "1000"
                env["MD_STEPS"] = "500"
            else:
                # Subsequent frames: Fast EM pass
                env["EM_STEPS"] = "2000"
                env["NB_STEPS"] = "0"
                env["MD_STEPS"] = "0"

            # 1st Pass: Run Backmapping via initram-v5
            subprocess.run([
                "./initram-v5.sh", "-f", f"trj_{i}ns.gro", "-o", f"aa_{i}ns.gro", 
                "-to", "charmm36", "-p", "topol.top"
            ], cwd=wd, env=env, check=True)
            
            # ========================================================
            # ADDED: STRUCTURAL SCREENING FOR COLLAPSED LIGAND
            # ========================================================
            if i != loop_from:
                collapsed = False
                lig_coords = []
                # Read the newly generated backmapped coordinates
                with open(wd / f"aa_{i}ns.gro", "r") as f:
                    for line in f.readlines()[2:-1]: # skip header and box
                        if len(line) > 20 and ligand_code in line[5:10]:
                            try:
                                x, y, z = float(line[20:28]), float(line[28:36]), float(line[36:44])
                                lig_coords.append((x, y, z))
                            except ValueError:
                                pass
                
                # Check for distances < 0.08 nm (0.0064 nm^2 squared distance). 
                # If atoms are this close, they collapsed into the CG bead center.
                for c1 in range(len(lig_coords)):
                    for c2 in range(c1 + 1, len(lig_coords)):
                        dx = lig_coords[c1][0] - lig_coords[c2][0]
                        dy = lig_coords[c1][1] - lig_coords[c2][1]
                        dz = lig_coords[c1][2] - lig_coords[c2][2]
                        if (dx*dx + dy*dy + dz*dz) < 0.0064:
                            collapsed = True
                            break
                    if collapsed:
                        break
                        
                if collapsed:
                    # Fallback: Rerun this specific frame with full simulated annealing
                    env["EM_STEPS"] = "3000"
                    env["NB_STEPS"] = "1000"
                    env["MD_STEPS"] = "500"
                    subprocess.run([
                        "./initram-v5.sh", "-f", f"trj_{i}ns.gro", "-o", f"aa_{i}ns.gro", 
                        "-to", "charmm36", "-p", "topol.top"
                    ], cwd=wd, env=env, check=True)
            
            # Grep logic replacement (Update Header)
            with open(wd / f"trj_{i}ns.gro", "r") as f:
                header_line = f.readline().strip()
                sys_match = re.search(r"system\s+(.*)", header_line)
                sys_val = sys_match.group(1) if sys_match else header_line
                
            with open(wd / f"aa_{i}ns.gro", "r") as f:
                aa_lines = f.readlines()
            aa_lines[0] = f"Protein - PGE {sys_val}\n"
            with open(wd / f"aa_{i}ns.gro", "w") as f:
                f.writelines(aa_lines)
                
            # Convert final GRO to XTC
            subprocess.run([
                "gmx", "trjconv", "-f", f"aa_{i}ns.gro", "-o", f"aa_{i}ns.xtc"
            ], cwd=wd, check=True)
            
            # Store in output directories
            shutil.copy(wd / f"aa_{i}ns.xtc", out_xtc_dir / f"aa_{i}ns.xtc")
            shutil.copy(wd / f"aa_{i}ns.gro", out_gro_dir / f"aa_{i}ns.gro")

        # 7. Concatenate the xtc trajectories
        xtc_files = [f.name for f in out_xtc_dir.glob("*.xtc")]
        if xtc_files:
            # Change output to combine.xtc instead of all.xtc
            trjcat_cmd = ["gmx", "trjcat", "-f"] + xtc_files + ["-o", "combine.xtc"]
            subprocess.run(trjcat_cmd, cwd=out_xtc_dir, check=True)

        # 8. Copy the earliest loop and generate the all-atomic .tpr
        earliest_gro = f"aa_{loop_from}ns.gro"
        source_gro = out_gro_dir / earliest_gro
        target_gro = wd / earliest_gro
        
        if source_gro.exists():
            shutil.copy2(source_gro, target_gro)
            
        grompp_cmd = [
            "gmx", "grompp", 
            "-f", mdp_file, 
            "-c", earliest_gro, 
            "-p", "topol.top", 
            "-o", "dynamic.tpr", 
            "-maxwarn", "5"
        ]
        subprocess.run(grompp_cmd, cwd=wd, check=True)

        # 9. Copy and paste the generated .tpr into all_xtc folder
        source_tpr = wd / "dynamic.tpr"
        target_tpr = out_xtc_dir / "dynamic.tpr"
        
        if source_tpr.exists():
            shutil.copy2(source_tpr, target_tpr)

        # 9.5 Fix Periodic Boundary Condition (PBC) Artifacts
        if (out_xtc_dir / "combine.xtc").exists() and target_tpr.exists():
            trjconv_pbc_cmd = [
                "gmx", "trjconv", 
                "-s", "dynamic.tpr", 
                "-f", "combine.xtc", 
                "-o", "all.xtc", 
                "-pbc", "mol", 
                "-center"
            ]
            # Input "1\n0\n" selects "Protein" for centering and "System" for output
            subprocess.run(trjconv_pbc_cmd, cwd=out_xtc_dir, input="1\n0\n", text=True, check=True)
            
            # (Optional) Remove combine.xtc after to save storage space
            (out_xtc_dir / "combine.xtc").unlink()

        # 10. Generate index file (all atomic)
        make_ndx_cmd = ["gmx", "make_ndx", "-f", "dynamic.tpr", "-o", "index.ndx"]
        subprocess.run(make_ndx_cmd, cwd=out_xtc_dir, input="q\n", text=True, check=True)
            
        # --- ADDED: CLEANUP REDUNDANT FILES FROM ROOT DIRECTORY ---
        for f in wd.glob("aa_*ns.gro"):
            f.unlink()
        for f in wd.glob("aa_*ns.xtc"):
            f.unlink()
        for f in wd.glob("trj_*ns.gro"):
            f.unlink()
            
        # 11. Zip the results
        zip_path = wd / "backmapped_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(out_xtc_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    zipf.write(file_path, arcname=f"all_xtc/{file}")
            for root, dirs, files in os.walk(out_gro_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    zipf.write(file_path, arcname=f"gro/{file}")
                    
        return {
            "status": "success",
            "message": "Backmapping completed successfully.",
            "result_file": str(zip_path)
        }
    except subprocess.CalledProcessError as e:
        error_output = e.stderr if e.stderr else e.stdout
        return {"status": "failed", "message": f"Command Failed: {error_output}"}
    except Exception as e:
        return {"status": "failed", "message": str(e)}

@celery_app.task
def run_plip_analysis_task(job_id: str):
    """Run the complete PLIP analysis pipeline for Tab 5."""
    work_dir = f"/app/shared_data/{job_id}/analysis/plip"
    os.makedirs(work_dir, exist_ok=True)
    
    backmap_dir = f"/app/shared_data/{job_id}/all_xtc"
    
    tpr_file = f"{backmap_dir}/dynamic.tpr"
    ndx_file = f"{backmap_dir}/index.ndx"
    xtc_file = f"{backmap_dir}/all.xtc"
    
    try:
        # Step 1: Convert trajectories to concatenated PDB
        trjconv_cmd = f"echo '0' | gmx trjconv -s {tpr_file} -n {ndx_file} -f {xtc_file} -dt 1000 -o concatenated.pdb"
        subprocess.run(trjconv_cmd, shell=True, cwd=work_dir, check=True)
        
        # ==================================================================
        # FIX: Bypass Git GnuTLS issues by downloading the repo natively
        # ==================================================================
        import requests # (Ensure this is imported at the top of worker.py)
        
        repo_url = "https://github.com/tasyriqomar/ColabMD-Edu_Protein-Ligand/archive/refs/heads/main.zip"
        zip_out_path = os.path.join(work_dir, "repo.zip")
        
        # Download the zip
        r = requests.get(repo_url)
        with open(zip_out_path, "wb") as f:
            f.write(r.content)
            
        # Extract the zip
        with zipfile.ZipFile(zip_out_path, 'r') as zip_ref:
            zip_ref.extractall(work_dir)
            
        # Move the Analysis folder and clean up
        extracted_folder = os.path.join(work_dir, "ColabMD-Edu_Protein-Ligand-main")
        analysis_dir = os.path.join(work_dir, "Analysis")
        
        shutil.copytree(os.path.join(extracted_folder, "Analysis"), analysis_dir, dirs_exist_ok=True)
        shutil.rmtree(extracted_folder)
        os.remove(zip_out_path)
        # ==================================================================
        
        subprocess.run(f"mv concatenated.pdb {analysis_dir}/", shell=True, cwd=work_dir, check=True)
        
        # ------------------------------------------------------------------
        # ON-THE-FLY SCRIPT PATCHING
        # Fix hardcoded Google Drive paths, remove PyMOL dependencies, and fix plots
        # ------------------------------------------------------------------
        colab_base_path = "/content/drive/MyDrive/ColabMD-Edu_Protein-Ligand"
        scripts_to_patch = [
            "plip.sh", 
            "convert_xml_to_json.py", 
            "process_json.py", 
            "percentage_interactions.py", 
            "split_pdb.py",
            "type_interactions_color.py",
            "residue_interactions_color_csv_a.py",
            "timeline_interaction_color.py"
        ]
        
        for script_name in scripts_to_patch:
            script_path = os.path.join(analysis_dir, script_name)
            if os.path.exists(script_path):
                with open(script_path, "r") as f:
                    content = f.read()
                
                # Replace Colab paths with current Docker path
                content = content.replace(f"{colab_base_path}/Analysis", analysis_dir)
                content = content.replace(colab_base_path, os.path.dirname(analysis_dir))
                
                # Remove PyMOL flags from the bash script to prevent crashes
                if script_name == "plip.sh":
                    content = content.replace("--pics", "").replace("-p", "").replace("--pymol", "").replace("-y", "")
                
                # NEW: Tilt the X-axis labels vertically (90 degrees) and reduce font size to prevent overlapping
                if script_name == "timeline_interaction_color.py":
                    # Inject a Matplotlib rotation rule right before the layout is tightened/saved
                    content = content.replace("plt.tight_layout()", "plt.xticks(rotation=90, fontsize=7)\n    plt.tight_layout()")
                    
                with open(script_path, "w") as f:
                    f.write(content)
        # ------------------------------------------------------------------
        
        # Step 3: Split the concatenated.pdb
        subprocess.run(["python", "split_pdb.py"], cwd=analysis_dir, check=True)
        
        # Step 4: Run the PLIP bash script
        subprocess.run(["bash", "plip.sh"], cwd=analysis_dir, check=True)
        
        # Step 5: Convert XML to JSON
        subprocess.run(["python", "convert_xml_to_json.py"], cwd=analysis_dir, check=True)
        
        # Copy the extracted JSON up to the main Analysis directory so BOTH scripts find it
        subprocess.run("cp split_pdbs/extracted_data.json .", shell=True, cwd=analysis_dir, check=True)
        
        # Step 6: Process JSON to generate CSVs
        subprocess.run(["python", "process_json.py"], cwd=analysis_dir, check=True)
        
        # Copy the generated CSV folder into split_pdbs where the timeline script expects it
        subprocess.run("cp -r bond_csv_files split_pdbs/", shell=True, cwd=analysis_dir, check=False)
        
        # Step 7: Generate tables, figures, and plots
        subprocess.run(["python", "percentage_interactions.py"], cwd=analysis_dir, check=True)
        subprocess.run(["python", "type_interactions_color.py"], cwd=analysis_dir, check=True)
        subprocess.run(["python", "residue_interactions_color_csv_a.py"], cwd=analysis_dir, check=True)
        subprocess.run(["python", "timeline_interaction_color.py"], cwd=analysis_dir, check=True)
        
        # Step 8: Zip the results
        zip_path = f"{work_dir}/plip_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(analysis_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, analysis_dir)
                    zipf.write(file_path, arcname=arcname)
                    
        return {
            "status": "success", 
            "message": "PLIP Pipeline completed successfully.",
            "result_file": str(zip_path)
        }
        
    except subprocess.CalledProcessError as e:
        return {"status": "error", "message": f"PLIP Pipeline failed at command: {e.cmd}\nError: {e.stderr}"}

# Add to colabmd_service/tasks/worker.py
@celery_app.task
def run_mmpbsa_task(job_id: str, ligand_name: str, protein_idx: int = 1, ligand_idx: int = 13, start_frame: int = 1, end_frame: int = 0):
    job_dir = Path(f"/app/shared_data/{job_id}")
    analysis_dir = job_dir / "analysis"
    mmgbsa_dir = analysis_dir / "MMGBSA"
    mmgbsa_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        import requests
        import re 
        
        # 1. Copy required files from Backmapping root & all_xtc
        xtc_dir = job_dir / "all_xtc"
        files_to_copy = [
            (xtc_dir / "all.xtc", mmgbsa_dir / "all.xtc"),
            (xtc_dir / "dynamic.tpr", mmgbsa_dir / "dynamic.tpr"),
            (xtc_dir / "index.ndx", mmgbsa_dir / "index.ndx"),
            (job_dir / "aa.gro", mmgbsa_dir / "aa.gro"),
            (job_dir / "topol.top", mmgbsa_dir / "topol.top")
        ]
        
        for src, dst in files_to_copy:
            if src.exists():
                shutil.copy2(src, dst)
                
        # Copy ALL .itp and .prm files regardless of upper/lowercase
        for f_path in job_dir.iterdir():
            if f_path.suffix.lower() in [".itp", ".prm"]:
                shutil.copy2(f_path, mmgbsa_dir / f_path.name)
        
        # Copy forcefield and mapping directories
        for d_name in ["charmm36-feb2026_cgenff-5.0.ff", "Mapping"]:
            src_d = job_dir / d_name
            dest_d = mmgbsa_dir / d_name
            if src_d.exists() and not dest_d.exists():
                shutil.copytree(src_d, dest_d)
                
        # 2. Download required python and input scripts
        scripts = {
            "mmpbsa.in": "https://raw.githubusercontent.com/tasyriqomar/ColabMD-Edu_Protein-Ligand/refs/heads/main/mmpbsa.in",
            "per-residue.py": "https://raw.githubusercontent.com/tasyriqomar/ColabMD-Edu_Protein-Ligand/refs/heads/main/per-residue.py",
            "heatmap.py": "https://raw.githubusercontent.com/tasyriqomar/ColabMD-Edu_Protein-Ligand/refs/heads/main/heatmap.py"
        }
        
        for s_name, s_url in scripts.items():
            r = requests.get(s_url)
            with open(mmgbsa_dir / s_name, "wb") as f:
                f.write(r.content)
                
        # ====================================================================
        # 3. DYNAMICALLY UPDATE mmpbsa.in & PYTHON SCRIPTS
        # ====================================================================
        # Calculate End Frame
        xtc_files = list(xtc_dir.glob("aa_*ns.xtc"))
        detected_end = len(xtc_files) if xtc_files else 1
        
        # --- ADDED: USE REQUESTED UI FRAMES (Fallback to detected if 0) ---
        final_end_frame = end_frame if end_frame > 0 else detected_end
        
        # Parse aa.gro for Protein Residue Count and Start/End Indices
        protein_res_nums = set()
        ligand_res_nums = set()
        solvent_res_names = {'SOL', 'W', 'NA', 'CL', 'K', 'NA+', 'CL-', 'K+', 'ION'}
        
        if (mmgbsa_dir / "aa.gro").exists():
            with open(mmgbsa_dir / "aa.gro", "r") as f:
                lines = f.readlines()[2:-1] # Skip header and box vectors
                for line in lines:
                    try:
                        res_num_str = line[0:5].strip()
                        res_name = line[5:10].strip()
                        if res_num_str:
                            res_num = int(res_num_str)
                            if ligand_name in res_name:
                                ligand_res_nums.add(res_num)
                            elif res_name not in solvent_res_names:
                                protein_res_nums.add(res_num)
                    except ValueError:
                        continue
                        
        if protein_res_nums:
            protein_start_res = min(protein_res_nums)
            protein_end_res = max(protein_res_nums)
        else:
            protein_start_res = 1
            protein_end_res = 496
            
        if ligand_res_nums:
            ligand_res_num = min(ligand_res_nums)
        else:
            ligand_res_num = protein_end_res + 1 
            
        total_residues = (ligand_res_num - protein_start_res) + 1
        
        # Rewrite mmpbsa.in
        mmpbsa_file = mmgbsa_dir / "mmpbsa.in"
        with open(mmpbsa_file, "r") as f:
            mmpbsa_text = f.read()
            
        # --- UPDATED: INJECT FINAL FRAME VARIABLES ---
        mmpbsa_text = re.sub(r"startframe\s*=\s*\d+", f"startframe={start_frame}", mmpbsa_text)
        mmpbsa_text = re.sub(r"endframe\s*=\s*\d+", f"endframe={final_end_frame}", mmpbsa_text)
        mmpbsa_text = re.sub(r'print_res\s*=\s*".*?"', f'print_res="A/{protein_start_res}-{protein_end_res} B/{ligand_res_num}"', mmpbsa_text)
        
        with open(mmpbsa_file, "w") as f:
            f.write(mmpbsa_text)

        # PATCH: Fix hardcoded Google Drive paths in per-residue.py
        per_res_script = mmgbsa_dir / "per-residue.py"
        if per_res_script.exists():
            with open(per_res_script, "r") as f:
                pr_text = f.read()
            pr_text = re.sub(r'base_dir\s*=\s*".*?"', 'base_dir = "./"', pr_text)
            with open(per_res_script, "w") as f:
                f.write(pr_text)

        # PATCH: Fix paths and hardcoded limits in heatmap.py
        heatmap_script = mmgbsa_dir / "heatmap.py"
        if heatmap_script.exists():
            with open(heatmap_script, "r") as f:
                hm_text = f.read()
            hm_text = re.sub(r'base_dir\s*=\s*".*?"', 'base_dir = "./"', hm_text)
            hm_text = re.sub(r'num_frames\s*=\s*\d+', f'num_frames={end_frame}', hm_text)
            hm_text = re.sub(r'num_residues\s*=\s*\d+', f'num_residues={total_residues}', hm_text)
            with open(heatmap_script, "w") as f:
                f.write(hm_text)
        # ====================================================================
        
        # 4. Execute gmx_MMPBSA
        mmpbsa_cmd = [
            "gmx_MMPBSA", "-O", "-i", "mmpbsa.in", "-cs", "dynamic.tpr",
            "-ct", "all.xtc", "-ci", "index.ndx", "-cg", str(protein_idx), str(ligand_idx),
            "-cp", "topol.top", "-o", "FINAL_RESULTS_MMPBSA.dat",
            "-eo", "FINAL_RESULTS_MMPBSA.csv", "-do", "FINAL_DECOMP_MMPBSA.dat",
            "-deo", "FINAL_DECOMP_MMPBSA.csv", "-nogui"
        ]
        
        try:
            # First attempt to run the full pipeline
            subprocess.run(mmpbsa_cmd, cwd=mmgbsa_dir, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            error_output = e.stderr if e.stderr else e.stdout
            # If it crashed because of the AmberTools formatting bug
            if "ValueError: could not convert string to float: '*********'" in error_output:
                print("Detected AmberTools '*********' formatting bug. Patching output files and resuming...")
                
                # FIX: Find ALL intermediate files (decomp files are .out, not .mdout)
                for bad_file in mmgbsa_dir.glob("_GMXMMPBSA_*"):
                    if bad_file.is_file():
                        try:
                            with open(bad_file, "r") as f:
                                mdout_text = f.read()
                            
                            # If the file contains the broken asterisks, patch it
                            if "*********" in mdout_text:
                                # Replace exactly 9 asterisks with a safe 9-character float
                                patched_text = mdout_text.replace("*********", " 9999.999")
                                with open(bad_file, "w") as f:
                                    f.write(patched_text)
                        except UnicodeDecodeError:
                            # Safely skip any binary files (like temporary trajectories)
                            pass
                
                # Re-run gmx_MMPBSA with the --rewrite-output flag. 
                rewrite_cmd = [
                    "gmx_MMPBSA", "-O", "-i", "mmpbsa.in", "-cs", "dynamic.tpr",
                    "-ct", "all.xtc", "-ci", "index.ndx", "-cg", str(protein_idx), str(ligand_idx),
                    "-cp", "topol.top", "-o", "FINAL_RESULTS_MMPBSA.dat",
                    "-eo", "FINAL_RESULTS_MMPBSA.csv", "-do", "FINAL_DECOMP_MMPBSA.dat",
                    "-deo", "FINAL_DECOMP_MMPBSA.csv", "-nogui", "--rewrite-output"
                ]
                subprocess.run(rewrite_cmd, cwd=mmgbsa_dir, check=True, capture_output=True, text=True)
            else:
                # If it crashed for a different reason, bubble it up to the frontend
                raise e
        
        # 5. Generate Figures
        subprocess.run(["python", "per-residue.py"], cwd=mmgbsa_dir, check=True, capture_output=True, text=True)
        subprocess.run(["python", "heatmap.py"], cwd=mmgbsa_dir, check=True, capture_output=True, text=True)
        
        # Extract the Delta results from the .dat file
        summary_data = ""
        res_dat = mmgbsa_dir / "FINAL_RESULTS_MMPBSA.dat"
        if res_dat.exists():
            with open(res_dat, "r") as f:
                lines = f.readlines()
                summary_data = "".join(lines[-20:])
                
        # 6. Zip results
        zip_path = analysis_dir / "mmgbsa_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(mmgbsa_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, mmgbsa_dir)
                    zipf.write(file_path, arcname=arcname)
                    
        return {
            "status": "success",
            "message": "MM-GBSA Pipeline completed successfully.",
            "result_file": str(zip_path),
            "summary_data": summary_data
        }
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr if e.stderr else e.stdout
        return {"status": "error", "message": f"Command Failed: {' '.join(e.cmd)}\nError Details:\n{error_msg}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

@celery_app.task
def run_trajectory_analysis_task(job_id: str, rmsd_group: int, rmsf_group: int, rg_group: int):
    """Run RMSD, RMSF, and Rg calculations on the backmapped trajectories."""
    import matplotlib.pyplot as plt
    
    job_dir = Path(f"/app/shared_data/{job_id}")
    xtc_dir = job_dir / "all_xtc"
    analysis_dir = job_dir / "analysis" / "trajectory"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    
    # Required input files
    tpr_file = xtc_dir / "dynamic.tpr"
    xtc_file = xtc_dir / "all.xtc"
    ndx_file = xtc_dir / "index.ndx"
    
    if not (tpr_file.exists() and xtc_file.exists() and ndx_file.exists()):
        return {"status": "error", "message": "Missing required trajectory files (all.xtc, dynamic.tpr, or index.ndx). Please ensure Backmapping completed successfully."}
        
    try:
        # 1. RMSD Calculation (Native support for -tu ns)
        rmsd_cmd = f"echo '{rmsd_group} {rmsd_group}' | gmx rms -s {tpr_file} -f {xtc_file} -n {ndx_file} -o {analysis_dir}/rmsd.xvg -tu ns"
        subprocess.run(rmsd_cmd, shell=True, check=True, capture_output=True, text=True)
        
        # 2. RMSF Calculation (-res calculates average per residue)
        rmsf_cmd = f"echo '{rmsf_group}' | gmx rmsf -s {tpr_file} -f {xtc_file} -n {ndx_file} -o {analysis_dir}/rmsf.xvg -res"
        subprocess.run(rmsf_cmd, shell=True, check=True, capture_output=True, text=True)
        
        # 3. Radius of Gyration Calculation (Outputs in ps, NO -tu flag)
        rg_cmd = f"echo '{rg_group}' | gmx gyrate -s {tpr_file} -f {xtc_file} -n {ndx_file} -o {analysis_dir}/gyration.xvg"
        subprocess.run(rg_cmd, shell=True, check=True, capture_output=True, text=True)
        
        # --- Helper function to plot XVG files ---
        # Added x_scale to safely convert ps to ns where needed
        def plot_xvg(xvg_path, title, xlabel, ylabel, out_png, x_scale=1.0):
            x, y = [], []
            if Path(xvg_path).exists():
                with open(xvg_path, 'r') as f:
                    for line in f:
                        if not line.startswith(('@', '#')):
                            parts = line.split()
                            if len(parts) >= 2:
                                # Multiply x by the scale (e.g., 0.001 for ps to ns)
                                x.append(float(parts[0]) * x_scale)
                                y.append(float(parts[1]))
                plt.figure(figsize=(8, 5))
                plt.plot(x, y, color='blue', linewidth=1.5)
                plt.title(title, fontsize=14, fontweight='bold')
                plt.xlabel(xlabel, fontsize=12)
                plt.ylabel(ylabel, fontsize=12)
                plt.grid(True, linestyle='--', alpha=0.7)
                plt.tight_layout()
                plt.savefig(out_png, dpi=300)
                plt.close()

        # 4. Generate Plots
        # RMSD is already in ns, scale is 1.0
        plot_xvg(analysis_dir / "rmsd.xvg", "RMSD Over Time", "Time (ns)", "RMSD (nm)", analysis_dir / "RMSD.png", x_scale=1.0)
        
        # RMSF uses Residue Numbers on the X-axis, scale is 1.0
        plot_xvg(analysis_dir / "rmsf.xvg", "Per-Residue RMSF", "Residue Number", "RMSF (nm)", analysis_dir / "RMSF.png", x_scale=1.0)
        
        # Gyration is in ps, we apply x_scale=0.001 to convert to ns!
        plot_xvg(analysis_dir / "gyration.xvg", "Radius of Gyration", "Time (ns)", "Rg (nm)", analysis_dir / "Gyration.png", x_scale=0.001)
        
        # 5. Zip results
        zip_path = analysis_dir / "trajectory_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file_path in analysis_dir.glob("*"):
                if file_path.name != "trajectory_results.zip":
                    zipf.write(file_path, arcname=file_path.name)
                    
        return {
            "status": "success",
            "message": "Trajectory Analysis completed successfully.",
            "result_file": str(zip_path)
        }
        
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr if e.stderr else e.stdout
        return {"status": "error", "message": f"GROMACS Command Failed: {e.cmd}\nDetails:\n{error_msg}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

@celery_app.task
def run_pca_fel_task(job_id: str, backbone_group: int):
    import requests
    job_dir = Path(f"/app/shared_data/{job_id}")
    xtc_dir = job_dir / "all_xtc"
    pca_dir = job_dir / "analysis" / "PCA-FEL"
    pca_dir.mkdir(parents=True, exist_ok=True)
    
    tpr_file = xtc_dir / "dynamic.tpr"
    xtc_file = xtc_dir / "all.xtc"
    ndx_file = xtc_dir / "index.ndx"
    
    if not (tpr_file.exists() and xtc_file.exists() and ndx_file.exists()):
        return {"status": "error", "message": "Missing required trajectory files. Please ensure Backmapping completed."}
        
    try:
        # Define group selection inputs (fit and analysis are the same group)
        group_input = f"{backbone_group} {backbone_group}\n"
        
        # 1. Covariance Analysis
        covar_cmd = f"gmx covar -f {xtc_file} -s {tpr_file} -n {ndx_file} -o {pca_dir}/eigenval.xvg -v {pca_dir}/eigenvec.trr -av {pca_dir}/average.pdb -l {pca_dir}/covar.log -b 0 -tu ns"
        subprocess.run(covar_cmd, input=group_input, shell=True, check=True, text=True, cwd=pca_dir)
        
        # 2. Eigenvector Components (1-2, 1-3, 2-3)
        projections = [
            ("1", "2"), ("1", "3"), ("2", "3")
        ]
        
        for first, last in projections:
            anaeig_cmd = f"gmx anaeig -v {pca_dir}/eigenvec.trr -f {xtc_file} -s {tpr_file} -n {ndx_file} -comp {pca_dir}/eigcomp{first}-{last}.xvg -rmsf {pca_dir}/eigrmsf{first}-{last}.xvg -2d {pca_dir}/2dproj{first}-{last}.xvg -b 0 -tu ns -first {first} -last {last}"
            subprocess.run(anaeig_cmd, input=group_input, shell=True, check=True, text=True, cwd=pca_dir)
            
        # 3. FEL Sham Execution
        sham_cmd = f"gmx sham -f {pca_dir}/2dproj1-2.xvg -notime -ls {pca_dir}/gibb1-2.xpm"
        subprocess.run(sham_cmd, shell=True, check=True, text=True, cwd=pca_dir)
        
        # 4. Download Scripts
        scripts = {
            "xpm2txt.py": "https://raw.githubusercontent.com/tasyriqomar/ColabMD-Edu_Protein-Ligand/refs/heads/main/xpm2txt.py",
            "PCA.py": "https://raw.githubusercontent.com/tasyriqomar/ColabMD-Edu_Protein-Ligand/refs/heads/main/PCA.py",
            "FEL.py": "https://raw.githubusercontent.com/tasyriqomar/ColabMD-Edu_Protein-Ligand/refs/heads/main/FEL.py"
        }
        
        for name, url in scripts.items():
            r = requests.get(url)
            with open(pca_dir / name, "wb") as f:
                f.write(r.content)
                
        # 5. Patch Python Scripts (Remove hardcoded Colab Paths)
        colab_path = "/content/drive/MyDrive/ColabMD-Edu_Protein-Ligand/Analysis/PCA-FEL/"
        for py_script in ["PCA.py", "FEL.py"]:
            script_path = pca_dir / py_script
            if script_path.exists():
                with open(script_path, "r") as f:
                    content = f.read()
                # Point file lookups to the current directory
                content = content.replace(colab_path, "./")
                with open(script_path, "w") as f:
                    f.write(content)
                    
        # 6. Patch xpm2txt.py to be Python 3 compatible natively
        xpm_script = pca_dir / "xpm2txt.py"
        if xpm_script.exists():
            with open(xpm_script, "r") as f:
                lines = f.readlines()
            
            patched_lines = []
            for line in lines:
                # Fix Python 2 print statements (e.g., "print USAGE" -> "print(USAGE)")
                stripped = line.strip()
                if stripped.startswith("print ") and not stripped.startswith("print("):
                    indent = line[:len(line) - len(line.lstrip())]
                    expr = stripped[6:].strip()
                    line = f"{indent}print({expr})\n"
                
                # Fix Python 3 map() subscriptability issue
                line = line.replace(
                    "map(float, line.split()[2:-2])", 
                    "list(map(float, line.split()[2:-2]))"
                )
                patched_lines.append(line)
            
            with open(xpm_script, "w") as f:
                f.writelines(patched_lines)
        
        # 7. Execute the scripts natively in Python 3
        subprocess.run(["python3", "xpm2txt.py", "-f", "gibb1-2.xpm", "-o", "FEL.dat"], cwd=pca_dir, check=True)
        subprocess.run(["python3", "PCA.py"], cwd=pca_dir, check=True)
        subprocess.run(["python3", "FEL.py"], cwd=pca_dir, check=True)
        
        # 8. Zip Results
        zip_path = pca_dir / "pca_fel_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(pca_dir):
                for file in files:
                    if file != "pca_fel_results.zip":
                        file_path = os.path.join(root, file)
                        zipf.write(file_path, arcname=file)
                        
        return {
            "status": "success",
            "message": "PCA & FEL Analysis completed successfully.",
            "result_file": str(zip_path)
        }
        
    except subprocess.CalledProcessError as e:
        return {"status": "error", "message": f"Command Failed: {e.cmd}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

@celery_app.task
def run_ttclust_task(job_id: str, n_clusters: int):
    import os
    import zipfile
    import subprocess
    from pathlib import Path

    job_dir = Path(f"/app/shared_data/{job_id}")
    xtc_dir = job_dir / "all_xtc"
    analysis_dir = job_dir / "analysis" / "ttclust"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    
    xtc_file = xtc_dir / "all.xtc"
    gro_dir = job_dir / "gro"
    gro_files = list(gro_dir.glob("aa_*ns.gro"))
    
    if not xtc_file.exists() or not gro_files:
        return {"status": "error", "message": "Missing all.xtc or full-system .gro trajectory files. Please ensure Backmapping completed successfully."}
        
    gro_file = gro_files[0]
    
    try:
        # 1. Install TTClust on the fly
        subprocess.run(["pip", "install", "ttclust"], check=True)
        
        # 2. Run TTClust 
        ttclust_cmd = [
            "ttclust", 
            "-f", str(xtc_file), 
            "-t", str(gro_file),
            "-n", str(n_clusters)
        ]
        
        # FIX: Force Matplotlib to use a headless backend to prevent display crashes
        env = os.environ.copy()
        env["MPLBACKEND"] = "Agg"
        
        # Execute without check=True to bypass falsely triggered CalledProcessError
        process = subprocess.run(
            ttclust_cmd, 
            cwd=analysis_dir, 
            env=env,
            text=True, 
            capture_output=True
        )
        
        # 3. Verify if TTClust actually produced the expected outputs (plots/data)
        output_files = list(analysis_dir.glob("**/*.png")) + list(analysis_dir.glob("**/*.csv"))
        if not output_files:
            error_msg = process.stderr if process.stderr else process.stdout
            return {
                "status": "error", 
                "message": f"TTClust failed to generate output files.\nLog:\n{error_msg}"
            }
        
        # 4. Zip results
        zip_path = analysis_dir / "ttclust_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(analysis_dir):
                for file in files:
                    if file != "ttclust_results.zip":
                        file_path = os.path.join(root, file)
                        arcname = os.path.relpath(file_path, analysis_dir)
                        zipf.write(file_path, arcname=arcname)
                        
        return {
            "status": "success",
            "message": "TTClust Analysis completed successfully.",
            "result_file": str(zip_path)
        }
        
    except Exception as e:
        return {"status": "error", "message": str(e)}
