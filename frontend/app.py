import streamlit as st
import requests
import time
import zipfile
import io
import re
import os  # <-- ensure os is imported for path handling
import networkx as nx
import numpy as np
import math
import pandas as pd  # <-- ADD THIS IMPORT
from stmol import showmol
import py3Dmol

API_URL = "http://127.0.0.1:8000/api/v1"

# 1. Update Title of Web
st.set_page_config(page_title="ColabMD Web Service", page_icon="🧬")
st.title("🧬 ColabMD Web Service")

# Initialize session state for persistent data across UI interactions
if "docking_task_id" not in st.session_state:
    st.session_state.docking_task_id = None
if "best_pose_bytes" not in st.session_state:
    st.session_state.best_pose_bytes = None
if "mol2_bytes" not in st.session_state:
    st.session_state.mol2_bytes = None
if "ligand_code" not in st.session_state: 
    st.session_state.ligand_code = "A20"
if "vina_scores" not in st.session_state: 
    st.session_state.vina_scores = None
if "aa_zip_bytes" not in st.session_state: 
    st.session_state.aa_zip_bytes = None
if "cg_zip_bytes" not in st.session_state: 
    st.session_state.cg_zip_bytes = None
if "complex_bytes" not in st.session_state: 
    st.session_state.complex_bytes = None
if "mapped_cg_itp" not in st.session_state: 
    st.session_state.mapped_cg_itp = None
if "mapped_cg_preview" not in st.session_state:
    st.session_state.mapped_cg_preview = None
if "cg_gro_text" not in st.session_state:
    st.session_state.cg_gro_text = None
if "map_file_text" not in st.session_state:
    st.session_state.map_file_text = None
if "job_id" not in st.session_state:
    st.session_state.job_id = ""
if "plip_zip_bytes" not in st.session_state:
    st.session_state.plip_zip_bytes = None

# Initialize session state for CGMD and MM-GBSA
if "cgmd_zip_bytes" not in st.session_state:
    st.session_state.cgmd_zip_bytes = None
if "mmgbsa_zip_bytes" not in st.session_state:
    st.session_state.mmgbsa_zip_bytes = None
if "mmgbsa_summary" not in st.session_state:
    st.session_state.mmgbsa_summary = ""

# Initialize session state for Trajectory Analysis
if "traj_zip_bytes" not in st.session_state:
    st.session_state.traj_zip_bytes = None

# Initialize session state for PCA and FEL
if "pca_fel_zip_bytes" not in st.session_state:
    st.session_state.pca_fel_zip_bytes = None

# Initialize session state for TTClust
if "ttclust_zip_bytes" not in st.session_state:
    st.session_state.ttclust_zip_bytes = None

# Create the Tabs (Updated to rearrange Tabs 5 through 9)
tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs([
    "1. Docking", 
    "2. Ligand Preparation", 
    "3. CGMD Preparation", 
    "4. Backmapping", 
    "5. GROMACS", 
    "6. TTClust",
    "7. PLIP",
    "8. MM-GBSA",
    "9. PCA & FEL"
])

# ==========================================
# TAB 1: DOCKING
# ==========================================
with tab1:
    st.header("AutoDock Vina Docking")
    
    st.session_state.ligand_code = st.text_input(
        "Ligand 3-Letter Code (e.g., A20, LIG)", 
        value=st.session_state.ligand_code, 
        max_chars=3
    ).upper()

    # --- ADDED ALPHAFOLD LINK HERE ---
    st.markdown("Don't have a protein structure? Generate one using [AlphaFold2.ipynb](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=G4yBrceuFbf3)")
    
    prot_file = st.file_uploader("Upload Protein (PDB)", type=["pdb"], key="prot_dock")
    
    lig_input_method = st.radio("Ligand Input Method", ["Upload Ligand (PDB)", "Ligand SMILES string"])
    
    lig_file = None
    lig_smiles = None
    
    if lig_input_method == "Upload Ligand (PDB)":
        lig_file = st.file_uploader("Upload Ligand (PDB)", type=["pdb"], key="lig_dock")
    else:
        lig_smiles = st.text_input("Enter Ligand SMILES string")

    can_run = prot_file and (lig_file or lig_smiles)
    if st.button("Run Docking") and can_run:
        st.session_state.docking_task_id = None
        st.session_state.best_pose_bytes = None
        st.session_state.vina_scores = None
        
        with st.spinner("Preparing payload and queueing docking task..."):
            files = {"protein": (prot_file.name, prot_file.getvalue(), "chemical/x-pdb")}
            data = {"ligand_code": st.session_state.ligand_code}
            
            if lig_file:
                files["ligand"] = (lig_file.name, lig_file.getvalue(), "chemical/x-pdb")
            if lig_smiles:
                data["ligand_smiles"] = lig_smiles
                
            res = requests.post(f"{API_URL}/run-docking", files=files, data=data)
            
        if res.status_code == 200:
            st.session_state.docking_task_id = res.json()["task_id"]
            st.success(f"Docking queued! ID: `{st.session_state.docking_task_id}`")
            status_placeholder = st.empty()
            
            while True:
                status_res = requests.get(f"{API_URL}/status/{st.session_state.docking_task_id}").json()
                status = status_res["status"]
                
                if status == "SUCCESS":
                    task_result = status_res.get("result", {})
                    if task_result.get("status") == "success":
                        status_placeholder.success("✅ Docking complete!")
                        st.balloons()
                        
                        st.session_state.vina_scores = task_result.get("vina_scores", [])
                        
                        # Extract the actual job folder name from the complex_file path
                        complex_file_path = task_result.get("complex_file", "")
                        if complex_file_path:
                            # e.g. "/app/shared_data/job_1785741076/A20_complex.pdb" -> "job_1785741076"
                            st.session_state.job_id = complex_file_path.replace('\\', '/').split('/')[-2]
                        
                        download_res = requests.get(f"{API_URL}/download/{st.session_state.docking_task_id}")
                        if download_res.status_code == 200:
                            st.session_state.best_pose_bytes = download_res.content
                            
                        complex_res = requests.get(f"{API_URL}/download-complex/{st.session_state.docking_task_id}")
                        if complex_res.status_code == 200:
                            st.session_state.complex_bytes = complex_res.content
                    else:
                        status_placeholder.error(f"❌ Script Failed: {task_result.get('message')}")
                    break
                elif status == "FAILED":
                    status_placeholder.error("❌ Task Failed.")
                    break
                else:
                    status_placeholder.info(f"⏳ Status: {status} (Vina is running...)")
                    time.sleep(3)

    if st.session_state.vina_scores:
        st.subheader("Docking Results")
        st.dataframe(st.session_state.vina_scores, width="stretch", hide_index=True)
        
    if st.session_state.best_pose_bytes:
        best_pose_name = f"{st.session_state.ligand_code}AA.pdb"
        st.download_button(
            label=f"⬇️ Download Best Pose ({best_pose_name})",
            data=st.session_state.best_pose_bytes,
            file_name=best_pose_name,
            mime="chemical/x-pdb",
            key="best_pose_tab1"
        )
        
    if st.session_state.complex_bytes:
        complex_name = f"{st.session_state.ligand_code}_complex.pdb"
        st.download_button(
            label=f"⬇️ Download Complex ({complex_name})",
            data=st.session_state.complex_bytes,
            file_name=complex_name,
            mime="chemical/x-pdb",
            key="complex_download"
        )

        st.markdown("---")
        st.subheader("Interactive 3D Viewer")
        pdb_string = st.session_state.complex_bytes.decode("utf-8")
        view = py3Dmol.view(width=700, height=500)
        view.addModel(pdb_string, "pdb")
        view.setStyle({'hetflag': False}, {'cartoon': {'color': 'lightblue'}})
        view.setStyle({'resn': st.session_state.ligand_code}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})
        view.zoomTo()
        showmol(view, height=500, width=700)

# ==========================================
# TAB 2: LIGAND PREPARATION
# ==========================================
with tab2:
    st.header("CGenFF & Martini 3 Parameterization")
    st.info(f"Active Ligand Code: **{st.session_state.ligand_code}** (Change this in Tab 1)")
    st.markdown("---")
    
    st.subheader("All-atomic (AA)")
    st.info(f"Step 1: Generate {st.session_state.ligand_code}.mol2 from the best docked pose.")
    
    if st.button("Generate .mol2 file"):
        if not st.session_state.docking_task_id:
            st.warning("Please run docking first!")
        else:
            with st.spinner("Running OpenBabel conversion..."):
                res = requests.post(f"{API_URL}/generate-mol2", json={"task_id": st.session_state.docking_task_id, "ligand_name": st.session_state.ligand_code})
                if res.status_code == 200:
                    task_id = res.json()["task_id"]
                    while True:
                        status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                        if status_res["status"] == "SUCCESS":
                            task_result = status_res.get("result", {})
                            if task_result.get("status") == "success":
                                dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                st.session_state.mol2_bytes = dl.content
                                st.success(f"✅ {st.session_state.ligand_code}.mol2 generated successfully!")
                            else:
                                st.error(f"❌ Script Failed: {task_result.get('message')}")
                            break
                        elif status_res["status"] == "FAILED":
                            st.error("❌ Celery Task Failed")
                            break
                        time.sleep(2)

    if st.session_state.mol2_bytes:
        st.download_button(
            label=f"⬇️ Download {st.session_state.ligand_code}.mol2", 
            data=st.session_state.mol2_bytes, 
            file_name=f"{st.session_state.ligand_code}.mol2",
            key="mol2_download"
        )
    
    st.markdown("Ensure you have processed your `.mol2` file through the [CGenFF web server](https://cgenff.com/) first.")
    
    str_file = st.file_uploader("Upload CGenFF Output (.str)", type=["str"], key="str_up")
    up_mol2 = st.file_uploader(f"Upload {st.session_state.ligand_code}.mol2", type=["mol2"], key="mol2_up")
    
    if st.button("Generate Topologies AA") and str_file and up_mol2:
        with st.spinner("Generating All-Atom topologies..."):
            files = {
                "mol2_file": (up_mol2.name, up_mol2.getvalue(), "chemical/x-mol2"),
                "str_file": (str_file.name, str_file.getvalue(), "text/plain")
            }
            data = {"ligand_name": st.session_state.ligand_code}
            res = requests.post(f"{API_URL}/run-aa-prep", files=files, data=data)
            
            if res.status_code == 200:
                task_id = res.json()["task_id"]
                while True:
                    status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                    if status_res["status"] == "SUCCESS":
                        task_result = status_res.get("result", {})
                        if task_result.get("status") == "success":
                            dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                            st.session_state.aa_zip_bytes = dl.content
                            st.success("✅ AA Topologies Generated!")
                        else:
                            st.error(f"❌ Script Failed: {task_result.get('message')}")
                        break
                    elif status_res["status"] == "FAILED":
                        st.error("❌ Celery Task Failed")
                        break
                    time.sleep(2)
                
    if st.session_state.aa_zip_bytes:
        st.download_button(
            label="📦 Download AA Topologies", 
            data=st.session_state.aa_zip_bytes, 
            file_name=f"{st.session_state.ligand_code}_AA_topologies.zip", 
            mime="application/zip",
            key="aa_topologies_download"
        )

    st.markdown("---")
    st.subheader("Coarse-grained (CG)")
    
    lig_smiles = st.text_input("Ligand SMILES String", key="cg_smiles")
        
    if st.button("Generate Topologies CG") and lig_smiles:
        with st.spinner("Generating Coarse-Grained topologies..."):
            res = requests.post(f"{API_URL}/run-cg-prep", data={"ligand_name": st.session_state.ligand_code, "smiles": lig_smiles})
            if res.status_code == 200:
                task_id = res.json()["task_id"]
                while True:
                    status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                    if status_res["status"] == "SUCCESS":
                        task_result = status_res.get("result", {})
                        if task_result.get("status") == "success":
                            dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                            st.session_state.cg_zip_bytes = dl.content
                            st.success("✅ CG Topologies Generated!")
                        else:
                            st.error(f"❌ Script Failed: {task_result.get('message')}")
                        break
                    elif status_res["status"] == "FAILED":
                        st.error("❌ Celery Task Failed")
                        break
                    time.sleep(2)
                
    if st.session_state.cg_zip_bytes:
        st.download_button(
            label="📦 Download CG Topologies", 
            data=st.session_state.cg_zip_bytes, 
            file_name=f"{st.session_state.ligand_code}_CG_topologies.zip", 
            mime="application/zip",
            key="cg_topologies_download"
        )

    st.markdown("---")
    st.subheader("Topology Mapping & Preview")
    
    if st.button("Generate MAP Preview"):
        if not st.session_state.aa_zip_bytes or not st.session_state.cg_zip_bytes:
            st.warning("⚠️ Please generate both AA and CG topologies first!")
        else:
            with st.spinner("Executing Robust BFS Graph Mapping..."):
                try:
                    with zipfile.ZipFile(io.BytesIO(st.session_state.aa_zip_bytes)) as aa_zip:
                        aa_itp_name = [n for n in aa_zip.namelist() if n.endswith('.itp')][0]
                        aa_lines = aa_zip.read(aa_itp_name).decode('utf-8').splitlines(True)

                    with zipfile.ZipFile(io.BytesIO(st.session_state.cg_zip_bytes)) as cg_zip:
                        cg_itp_name = [n for n in cg_zip.namelist() if n.endswith('.itp')][0]
                        cg_lines = cg_zip.read(cg_itp_name).decode('utf-8').splitlines(True)
                    
                    atoms_dict = {}
                    mol_graph = nx.Graph()
                    current_section = None

                    for line in aa_lines:
                        stripped = line.strip()
                        if not stripped or stripped.startswith(';'): continue
                        if stripped.startswith('[') and stripped.endswith(']'):
                            current_section = stripped.replace(" ", "").lower()
                            continue

                        if current_section == '[atoms]':
                            tokens = stripped.split()
                            if len(tokens) >= 5:
                                idx = int(tokens[0])
                                name = tokens[4]
                                atoms_dict[idx] = name
                                mol_graph.add_node(idx)

                        elif current_section == '[bonds]':
                            tokens = stripped.split()
                            if len(tokens) >= 2:
                                mol_graph.add_edge(int(tokens[0]), int(tokens[1]))

                    heavy_atoms = [idx for idx, name in atoms_dict.items() if not name.upper().startswith('H')]
                    G_heavy = mol_graph.subgraph(heavy_atoms).copy()
                    
                    cg_beads = []
                    cg_graph = nx.Graph()
                    in_atoms_block = False
                    current_section = None

                    for line in cg_lines:
                        stripped = line.strip()
                        if not stripped: continue
                        if stripped.startswith('[') and stripped.endswith(']'):
                            current_section = stripped.replace(" ", "").lower()
                            in_atoms_block = (current_section == '[atoms]')
                            continue

                        if in_atoms_block and not stripped.startswith(';'):
                            parts = line.split(';', 1)
                            core = parts[0].split()
                            comment = parts[1].strip() if len(parts) > 1 else ""
                            if core and core[0].isdigit():
                                b_num = int(core[0])
                                cg_beads.append({'num': b_num, 'core': parts[0].rstrip(), 'comment': comment})
                                cg_graph.add_node(b_num)

                        elif current_section in ('[bonds]', '[constraints]'):
                            tokens = stripped.split()
                            if tokens and tokens[0].isdigit() and tokens[1].isdigit():
                                cg_graph.add_edge(int(tokens[0]), int(tokens[1]))

                    # --- ROBUST GRAPH MAPPING ALGORITHM ---
                    final_bead_assignments = {b['num']: [] for b in cg_beads}
                    mapped_aa = set()

                    # 1. Parse target sizes for each bead based on comments
                    bead_targets = {}
                    for b in cg_beads:
                        clean = re.sub(r'[Hh0-9()=\-#\[\]\s;]', '', b['comment'])
                        elements = [c for c in clean if c.isupper()]
                        bead_targets[b['num']] = len(elements) if elements else 1

                    # 2. Map Rings (Anchors)
                    aa_rings = [set(c) for c in nx.cycle_basis(G_heavy) if len(c) >= 5]
                    cg_rings = [set(c) for c in nx.cycle_basis(cg_graph) if len(c) == 3]

                    for cgr in cg_rings:
                        for aar in aa_rings:
                            if not aar.issubset(mapped_aa):
                                sub_g = G_heavy.subgraph(aar)
                                try:
                                    cycle_nodes = list(nx.find_cycle(sub_g))
                                    ordered_aar = [u for u, v in cycle_nodes]
                                except nx.NetworkXNoCycle:
                                    ordered_aar = list(aar)
                                
                                if not ordered_aar:
                                    ordered_aar = list(aar)

                                cgr_list = list(cgr)
                                chunk_size = max(1, len(ordered_aar) // len(cgr_list))
                                
                                for i, b in enumerate(cgr_list):
                                    if i == len(cgr_list) - 1:
                                        atoms = ordered_aar[i*chunk_size:]
                                    else:
                                        atoms = ordered_aar[i*chunk_size : (i+1)*chunk_size]
                                    final_bead_assignments[b].extend(atoms)
                                    mapped_aa.update(atoms)
                                break

                    # 3. BFS Expansion for Chains/Unmapped Beads
                    unmapped_beads = [b['num'] for b in cg_beads if not final_bead_assignments[b['num']]]
                    progress = True
                    while unmapped_beads and progress:
                        progress = False
                        for bead in list(unmapped_beads):
                            # Find mapped neighbors
                            mapped_neighbors = [n for n in cg_graph.neighbors(bead) if final_bead_assignments[n]]
                            if mapped_neighbors:
                                anchor_bead = mapped_neighbors[0]
                                anchor_atoms = final_bead_assignments[anchor_bead]
                                
                                available_aa = set()
                                for aa in anchor_atoms:
                                    for neighbor in G_heavy.neighbors(aa):
                                        if neighbor not in mapped_aa:
                                            available_aa.add(neighbor)
                                
                                if available_aa:
                                    target = bead_targets.get(bead, 3)
                                    collected = []
                                    queue = list(available_aa)
                                    
                                    while queue and len(collected) < target:
                                        curr = queue.pop(0)
                                        if curr not in mapped_aa and curr not in collected:
                                            collected.append(curr)
                                            for n in G_heavy.neighbors(curr):
                                                if n not in mapped_aa and n not in collected:
                                                    queue.append(n)
                                                    
                                    if collected:
                                        final_bead_assignments[bead].extend(collected)
                                        mapped_aa.update(collected)
                                        unmapped_beads.remove(bead)
                                        progress = True

                    # 4. Fallback for completely disjoint beads
                    if unmapped_beads:
                        remaining_aa = list(set(G_heavy.nodes()) - mapped_aa)
                        for bead in unmapped_beads:
                            target = bead_targets.get(bead, 3)
                            chunk = remaining_aa[:target]
                            final_bead_assignments[bead].extend(chunk)
                            mapped_aa.update(chunk)
                            remaining_aa = remaining_aa[target:]

                    # 5. Universal Orphan Rescue
                    unmapped_heavy = set(G_heavy.nodes()) - mapped_aa
                    for orphan in unmapped_heavy:
                        closest_bead = None
                        min_dist = float('inf')
                        for b, atoms in final_bead_assignments.items():
                            for a in atoms:
                                try:
                                    d = nx.shortest_path_length(G_heavy, orphan, a)
                                    if d < min_dist:
                                        min_dist = d
                                        closest_bead = b
                                except nx.NetworkXNoPath:
                                    pass
                        if closest_bead:
                            final_bead_assignments[closest_bead].append(orphan)
                        else:
                            first_bead = list(final_bead_assignments.keys())[0]
                            final_bead_assignments[first_bead].append(orphan)

                    # Append hydrogens back and rebuild .itp lines
                    updated_cg_lines = []
                    in_atoms_block = False

                    for line in cg_lines:
                        stripped = line.strip()
                        if stripped.startswith('[') and stripped.endswith(']'):
                            current_section = stripped.replace(" ", "").lower()
                            in_atoms_block = (current_section == '[atoms]')
                            updated_cg_lines.append(line)
                            continue

                        if in_atoms_block and stripped and not stripped.startswith(';'):
                            parts = line.split(';', 1)
                            core_line = parts[0].rstrip()
                            tokens = core_line.split()

                            if tokens and tokens[0].isdigit():
                                bead_num = int(tokens[0])
                                if bead_num in final_bead_assignments:
                                    assigned_heavy = final_bead_assignments[bead_num]
                                    full_set = set(assigned_heavy)

                                    for h_node in assigned_heavy:
                                        for neighbor in mol_graph.neighbors(h_node):
                                            if atoms_dict[neighbor].upper().startswith('H'):
                                                full_set.add(neighbor)

                                    sorted_names = [atoms_dict[i] for i in sorted(list(full_set))]
                                    mapped_string = " ".join(sorted_names)

                                    new_line = f"{core_line:<50} ; {mapped_string}\n"
                                    updated_cg_lines.append(new_line)
                                    continue

                        updated_cg_lines.append(line)
                    
                    st.session_state.mapped_cg_itp = "".join(updated_cg_lines)
                    
                    preview_lines = []
                    preview_started = False
                    for line in updated_cg_lines:
                        if '[atoms]' in line.replace(" ", "").lower():
                            preview_started = True
                        elif preview_started and line.startswith('['):
                            break
                        if preview_started:
                            preview_lines.append(line)
                            
                    st.session_state.mapped_cg_preview = "".join(preview_lines)
                    st.success("✅ Mapping generated successfully! (Using BFS Topo-Mapper)")
                    
                except Exception as e:
                    st.error(f"❌ Could not process mapping: {str(e)}")

    if st.session_state.mapped_cg_itp:
        st.download_button(
            label="⬇️ Download Mapped CG Topology (.itp)",
            data=st.session_state.mapped_cg_itp,
            file_name=f"{st.session_state.ligand_code}_CG_mapped.itp",
            mime="text/plain",
            key="mapped_cg_download"
        )
        st.text_area("MAP PREVIEW ([atoms] section)", value=st.session_state.mapped_cg_preview, height=400)


    # --- 3.4 GRO Generation Section ---
    st.markdown("---")
    st.subheader("Generate `.gro` with PDB Coordinates")
    
    if st.button("Generate .gro File"):
        if not st.session_state.mapped_cg_itp or not st.session_state.complex_bytes or not st.session_state.aa_zip_bytes:
            st.warning(f"⚠️ Please ensure you have generated the MAP Preview, AA Topologies, and have a docked {st.session_state.ligand_code}_complex.pdb!")
        else:
            try:
                with st.spinner(f"Aligning PDB heavy atoms to exact docked coordinates from {st.session_state.ligand_code}_complex.pdb..."):
                    
                    with zipfile.ZipFile(io.BytesIO(st.session_state.aa_zip_bytes)) as aa_zip:
                        aa_itp_name = [n for n in aa_zip.namelist() if n.endswith('.itp')][0]
                        aa_lines = aa_zip.read(aa_itp_name).decode('utf-8').splitlines()

                    aa_heavy_names = []
                    in_atoms_aa = False
                    for line in aa_lines:
                        line_clean = line.strip()
                        if not line_clean or line_clean.startswith(';'):
                            continue
                        if line_clean.startswith('[') and line_clean.endswith(']'):
                            in_atoms_aa = ('atoms' in line_clean.lower())
                            continue
                        if in_atoms_aa:
                            parts = line_clean.split()
                            if len(parts) >= 5:
                                atom_name = parts[4]
                                if not atom_name.upper().startswith('H'):
                                    aa_heavy_names.append(atom_name)

                    pdb_string = st.session_state.complex_bytes.decode("utf-8")
                    pdb_heavy_coords = []
                    
                    for line in pdb_string.splitlines():
                        if (line.startswith("HETATM") or line.startswith("ATOM")) and line[17:21].strip() == st.session_state.ligand_code:
                            atom_name_pdb = line[12:16].strip()
                            if not atom_name_pdb.upper().startswith('H'):
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
                                pdb_heavy_coords.append(np.array([x, y, z]))

                    if len(aa_heavy_names) != len(pdb_heavy_coords):
                        raise ValueError(f"Heavy atom mismatch! The AA Topologies have {len(aa_heavy_names)} heavy atoms, but the docked {st.session_state.ligand_code}_complex.pdb has {len(pdb_heavy_coords)}. Ensure the topologies match the ligand.")

                    name_to_coord = dict(zip(aa_heavy_names, pdb_heavy_coords))
                    
                    gro_template = []
                    bead_mapping = {}
                    molecule_name = st.session_state.ligand_code
                    in_atoms_section = False
                    
                    for line in st.session_state.mapped_cg_itp.splitlines():
                        line_clean = line.strip()
                        if not line_clean:
                            continue
                            
                        if line_clean.startswith("["):
                            section_name = line_clean.replace("[", "").replace("]", "").strip()
                            in_atoms_section = (section_name == "atoms")
                            continue
                            
                        if in_atoms_section and not line_clean.startswith(";"):
                            parts = line_clean.split(";")
                            core_data = parts[0].split()
                            
                            if len(core_data) < 5:
                                continue
                                
                            bead_id = int(core_data[0])
                            res_name = core_data[3]
                            bead_name = core_data[4]
                            
                            molecule_name = res_name
                            gro_template.append((res_name, bead_name, bead_id))
                            
                            if len(parts) > 1:
                                comment_text = parts[1]
                                tokens = re.findall(r'\w+', comment_text)
                                
                                mapped_coords = []
                                for token in tokens:
                                    if token in name_to_coord:
                                        mapped_coords.append(name_to_coord[token])
                                        
                                bead_mapping[bead_id] = mapped_coords
                                
                    computed_coords = {}
                    for bead_id, coord_list in bead_mapping.items():
                        if not coord_list:
                            raise ValueError(f"Bead {bead_id} has no valid heavy atoms mapped.")
                            
                        computed_coords[bead_id] = np.mean(coord_list, axis=0) / 10.0
                        
                    gro_lines = []
                    gro_lines.append(molecule_name)
                    gro_lines.append(f"{len(gro_template):>5}")
                    for resname, atomname, atomid in gro_template:
                        coord = computed_coords[atomid]
                        line = f"{1:>5}{resname:<5}{atomname:>5}{atomid:>5}{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                        gro_lines.append(line)
                    gro_lines.append(" 5.00000 5.00000 5.00000")
                    
                    st.session_state.cg_gro_text = "\n".join(gro_lines)
                    st.success(f"✅ {st.session_state.ligand_code}_CG.gro file generated successfully using coordinates from {st.session_state.ligand_code}_complex.pdb!")
                    
            except Exception as e:
                st.error(f"❌ Could not generate GRO: {str(e)}")
                
    if st.session_state.cg_gro_text:
        st.download_button(
            label="⬇️ Download Mapped CG Coordinate (.gro)",
            data=st.session_state.cg_gro_text,
            file_name=f"{st.session_state.ligand_code}_CG.gro",
            mime="text/plain",
            key="cg_gro_download"
        )
        st.text_area("GRO Preview", value=st.session_state.cg_gro_text, height=300)

    # --- 3.5 .map File Generation Section ---
    st.markdown("---")
    st.subheader("Generate `.map` File")
    
    if st.button("Generate .map File"):
        if not st.session_state.mapped_cg_itp or not st.session_state.aa_zip_bytes:
            st.warning("⚠️ Please ensure you have generated both the AA Topologies and the MAP Preview!")
        else:
            with st.spinner("Generating .map file..."):
                try:
                    # 1. Parse rules dynamically from your mapped preview text
                    mapping_rules = {}
                    martini_beads = []
                    in_atoms_cg = False
                    
                    for line in st.session_state.mapped_cg_itp.splitlines():
                        stripped = line.strip()
                        if not stripped: continue
                        
                        if stripped.startswith('[') and stripped.endswith(']'):
                            current_section = stripped.replace('[', '').replace(']', '').strip().lower()
                            in_atoms_cg = (current_section == 'atoms')
                            continue
                            
                        if in_atoms_cg and ';' in line:
                            main_part, comment_part = line.split(';', 1)
                            main_parts = main_part.strip().split()
                            
                            if len(main_parts) >= 5:
                                cg_bead = main_parts[4]
                                if cg_bead not in martini_beads:
                                    martini_beads.append(cg_bead)
                                    
                                aa_atoms = comment_part.strip().split()
                                for aa_atom in aa_atoms:
                                    mapping_rules[aa_atom] = cg_bead

                    # 2. Parse sequential atom list from your All-Atom file inside the zip
                    aa_atoms_list = []
                    in_atoms_aa = False
                    with zipfile.ZipFile(io.BytesIO(st.session_state.aa_zip_bytes)) as aa_zip:
                        aa_itp_name = [n for n in aa_zip.namelist() if n.endswith('.itp')][0]
                        aa_lines = aa_zip.read(aa_itp_name).decode('utf-8').splitlines()

                    for line in aa_lines:
                        clean_line = line.split(';')[0].split('#')[0].strip()
                        if not clean_line: continue
                        
                        if clean_line.startswith('[') and clean_line.endswith(']'):
                            current_section = clean_line.replace('[', '').replace(']', '').strip().lower()
                            in_atoms_aa = (current_section == 'atoms')
                            continue
                            
                        if in_atoms_aa:
                            parts = clean_line.split()
                            if len(parts) >= 5:
                                atom_idx = parts[0]
                                atom_name = parts[4]
                                aa_atoms_list.append((atom_idx, atom_name))
                                
                    # 3. Write out the beautifully formatted .map file string
                    map_lines = []
                    map_lines.append("[ molecule ]")
                    map_lines.append(f"{st.session_state.ligand_code}\n")
                    
                    map_lines.append("[ martini ]")
                    map_lines.append(" ".join(martini_beads) + "\n")
                    
                    map_lines.append("[ mapping ]")
                    map_lines.append("charmm36\n")
                    
                    map_lines.append("[ atoms ]")
                    unmapped_count = 0
                    for idx, atom_name in aa_atoms_list:
                        target_bead = mapping_rules.get(atom_name, "UNMAPPED")
                        if target_bead == "UNMAPPED":
                            unmapped_count += 1
                        map_lines.append(f"{idx:<5} {atom_name:<6} {target_bead}")
                        
                    st.session_state.map_file_text = "\n".join(map_lines)
                    
                    if unmapped_count > 0:
                        st.warning(f"⚠️ Warning: {unmapped_count} atoms in AA.itp were not found in your preview comments. Check the UNMAPPED tags in the file.")
                    else:
                        st.success(f"✅ {st.session_state.ligand_code}.map file generated successfully!")

                except Exception as e:
                    st.error(f"❌ Could not generate .map file: {str(e)}")

    if st.session_state.map_file_text:
        st.download_button(
            label="⬇️ Download .map File",
            data=st.session_state.map_file_text,
            file_name=f"{st.session_state.ligand_code}.map",
            mime="text/plain",
            key="map_file_download"
        )
        st.text_area(".map File Preview", value=st.session_state.map_file_text, height=300)

# ==========================================
# TAB 3: CGMD PREPARATION (Continued execution)
# ==========================================
with tab3:
    st.header("Coarse-Grained MD (CGMD) Preparation")
    st.markdown("Upload your Protein PDB, Ligand ITP, and Ligand GRO to prepare the CG system.")
    
    prot_file_cg = st.file_uploader("1. Upload Protein PDB (Original)", type=["pdb"], key="prot_cg")
    itp_file_cg = st.file_uploader("2. Upload Ligand CG Topology (.itp)", type=["itp"], key="itp_cg_up")
    gro_file_cg = st.file_uploader("3. Upload Ligand CG Coordinate (.gro)", type=["gro"], key="gro_cg_up")
    
    ligand_code = st.text_input("Ligand Code", value=st.session_state.ligand_code, key="lig_code_cg")
    
    if st.button("Prepare CGMD System"):
        st.session_state.cgmd_zip_bytes = None # Reset state on new run
        
        if not prot_file_cg or not itp_file_cg or not gro_file_cg:
            st.warning("⚠️ Please upload the Protein PDB, Ligand ITP, and Ligand GRO files.")
        else:
            with st.spinner("Queueing Martinize2 and GROMACS tools..."):
                files = {
                    "protein_pdb": (prot_file_cg.name, prot_file_cg.getvalue(), "chemical/x-pdb"),
                    "ligand_itp": (itp_file_cg.name, itp_file_cg.getvalue(), "text/plain"),
                    "ligand_gro": (gro_file_cg.name, gro_file_cg.getvalue(), "text/plain")
                }
                data = {"ligand_code": ligand_code}
                
                # Call the FastAPI backend
                response = requests.post(f"{API_URL}/run-cgmd-prep", files=files, data=data)
                
                if response.status_code == 200:
                    task_id = response.json()["task_id"]
                    status_placeholder = st.empty()
                    
                    # Poll Celery status
                    while True:
                        status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                        status = status_res["status"]
                        
                        if status == "SUCCESS":
                            task_result = status_res.get("result", {})
                            if task_result.get("status") == "success":
                                status_placeholder.success("✅ CGMD Preparation Complete!")
                                
                                # Download the result ZIP
                                dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                if dl.status_code == 200:
                                    st.session_state.cgmd_zip_bytes = dl.content
                                    time.sleep(1)
                                    st.rerun() # Refresh to show the download button safely
                                else:
                                    st.error("❌ Failed to download CGMD files.")
                            else:
                                status_placeholder.error(f"❌ Failed: {task_result.get('message')}")
                            break
                        elif status == "FAILED":
                            status_placeholder.error("❌ Celery Task Failed")
                            break
                        else:
                            status_placeholder.info(f"⏳ Status: {status} (Processing Martinize2/GROMACS...)")
                            time.sleep(3)
                else:
                    st.error(f"Error {response.status_code}: {response.text}")

    # Render the download button if the bytes exist
    if st.session_state.cgmd_zip_bytes:
        st.markdown("---")
        st.download_button(
            label="⬇️ Download CGMD Preparation Files (.zip)",
            data=st.session_state.cgmd_zip_bytes,
            file_name=f"{st.session_state.ligand_code}_cgmd_prep_files.zip",
            mime="application/zip",
            key="cgmd_prep_download"
        )
        st.info("📦 **Included files:** `system.top`, `system.gro`, `molecule_0.itp`, and `CG.pdb`")

    # --- ADDED ColabMD-Edu:Protein-Ligand Dynamics LINK HERE ---
    st.markdown("Proceed with Coarse-Grained MD Simulation in GPU Environment? Generate one using [ColabMD-Edu_Protein-Ligand_Run-CGMD.ipynb](https://github.com/tasyriqomar/ColabMD-Edu/blob/main/ColabMD_Edu_Protein_Ligand_Run_CGMD.ipynb)")

# ==========================================
# TAB 4: BACKMAPPING (CG to AA)
# ==========================================
with tab4:
    st.header("Backmapping (CG to AA)")
    st.markdown("Upload your CGMD simulation files and topologies to convert back to all-atom resolution.")
    
    col1, col2 = st.columns(2)
    with col1:
        protein_pdb_bm = st.file_uploader("1. Upload protein.pdb (AA)", type=["pdb"], key="prot_bm")
        ligand_map = st.file_uploader("2. Upload ligand.map", type=["map"], key="map_bm")
        xtc_file = st.file_uploader("3. Upload dynamic_fitted.xtc", type=["xtc"], key="xtc_bm")
        tpr_file = st.file_uploader("4. Upload dynamic.tpr", type=["tpr"], key="tpr_bm")
        mdp_file = st.file_uploader("5. Upload dynamic.mdp", type=["mdp"], key="mdp_bm")
    with col2:
        ndx_file = st.file_uploader("6. Upload index.ndx", type=["ndx", "txt"], key="ndx_bm")
        lig_itp_file = st.file_uploader("7. Upload ligand.itp (CGenFF AA)", type=["itp"], key="itp_bm")
        lig_prm_file = st.file_uploader("8. Upload ligand.prm (CGenFF AA)", type=["prm"], key="prm_bm")
        lig_cg_itp_file = st.file_uploader("9. Upload ligand CG .itp (from Tab 3)", type=["itp"], key="cg_itp_bm")
        system_gro = st.file_uploader("10. Upload system.gro (CGMD)", type=["gro"], key="sys_gro_bm")
        
        st.markdown("#### Frame Loop Settings")
        ligand_code = st.text_input("Ligand Code", value=st.session_state.ligand_code, key="lig_code_bm")
        loop_from = st.number_input("Loop from (ns):", min_value=1, value=1)
        loop_until = st.number_input("Loop until (ns):", min_value=10, value=10)

    all_files = [protein_pdb_bm, ligand_map, xtc_file, tpr_file, mdp_file, ndx_file, lig_itp_file, lig_prm_file, lig_cg_itp_file, system_gro]
    
    if st.button("Run Backmapping"):
        if not all(all_files):
            st.warning("⚠️ Please upload all 10 required files.")
        else:
            with st.spinner(f"Running Backmapping loop from {loop_from} to {loop_until}..."):
                files = {
                    "protein_pdb": (protein_pdb_bm.name, protein_pdb_bm.getvalue(), "chemical/x-pdb"),
                    "ligand_map": (ligand_map.name, ligand_map.getvalue(), "text/plain"),
                    "xtc_file": (xtc_file.name, xtc_file.getvalue(), "application/octet-stream"),
                    "tpr_file": (tpr_file.name, tpr_file.getvalue(), "application/octet-stream"),
                    "mdp_file": (mdp_file.name, mdp_file.getvalue(), "text/plain"),
                    "ndx_file": (ndx_file.name, ndx_file.getvalue(), "text/plain"),
                    "lig_itp_file": (lig_itp_file.name, lig_itp_file.getvalue(), "text/plain"),
                    "lig_prm_file": (lig_prm_file.name, lig_prm_file.getvalue(), "text/plain"),
                    "lig_cg_itp_file": (lig_cg_itp_file.name, lig_cg_itp_file.getvalue(), "text/plain"),
                    "system_gro": (system_gro.name, system_gro.getvalue(), "text/plain")
                }
                data = {
                    "ligand_code": ligand_code,
                    "loop_from": loop_from,
                    "loop_until": loop_until
                }
                
                response = requests.post(f"{API_URL}/run-backmapping", files=files, data=data)
                
                if response.status_code == 200:
                    task_id = response.json()["task_id"]
                    status_placeholder = st.empty()
                    
                    while True:
                        status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                        status = status_res["status"]
                        
                        if status == "SUCCESS":
                            task_result = status_res.get("result", {})
                            if task_result.get("status") == "success":
                                status_placeholder.success(f"✅ Backmapping Complete for frames {loop_from} to {loop_until}!")
                                
                                # ==========================================
                                # FIX: EXTRACT AND SAVE JOB ID FOR TAB 5 & 6
                                # ==========================================
                                res_file_path = task_result.get("result_file", "")
                                if res_file_path:
                                    # Extract Job ID and force UI fields to update
                                    extracted_job_id = res_file_path.replace('\\', '/').split('/')[-2]
                                    st.session_state.job_id = extracted_job_id
                                    st.session_state.plip_job_id = extracted_job_id
                                    st.session_state.mmgbsa_job_id = extracted_job_id
                                
                                dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                
                                st.download_button(
                                    label="⬇️ Download Backmapped Trajectories (.zip)",
                                    data=dl.content,
                                    file_name=f"{ligand_code}_backmapped_results.zip",
                                    mime="application/zip",
                                    key="backmap_download"
                                )
                                st.info(
                                    "**Containing:**\n"
                                    "1. `gro` folder\n"
                                    "2. `all_xtc` folder\n"
                                    "3. `all.xtc` file\n"
                                    "4. `dynamic.tpr` file\n"
                                    "5. `index.ndx` file"
                                )
                                
                                # Force the app to refresh immediately so Tabs 5 & 6 show the ID
                                time.sleep(1)
                                st.rerun()
                            else:
                                status_placeholder.error(f"❌ Failed: {task_result.get('message')}")
                            break
                        elif status == "FAILED":
                            status_placeholder.error("❌ Celery Task Failed")
                            break
                        else:
                            status_placeholder.info(f"⏳ Status: {status} (Processing Trajectories...)")
                            time.sleep(3)
                else:
                    st.error(f"Error {response.status_code}: {response.text}")

# ==========================================
# TAB 5: GROMACS (Trajectory Analysis)
# ==========================================
with tab5:
    st.header("Step 5: GROMACS Trajectory Analysis (RMSD, RMSF, Rg)")
    
    job_id_input_traj = st.text_input("Job ID", value=st.session_state.get("job_id", ""), key="traj_job_id")
    
    if not job_id_input_traj:
        st.info("⏳ Waiting for Backmapping (Tab 4) to finish to detect Job ID, or paste your Job ID above to proceed.")
    else:
        st.session_state.job_id = job_id_input_traj
        st.success(f"**Linked Job ID:** `{job_id_input_traj}`")
        
        st.markdown("### GROMACS Group Selection")
        col1, col2, col3 = st.columns(3)
        with col1:
            rmsd_group = st.number_input("RMSD Group (e.g., 4=Backbone, 1=Protein)", value=4, step=1)
        with col2:
            rmsf_group = st.number_input("RMSF Group (e.g., 1=Protein, 4=Backbone)", value=1, step=1)
        with col3:
            rg_group = st.number_input("Rg Group (e.g., 1=Protein)", value=1, step=1)
            
        if st.button("Run GROMACS Analysis"):
            st.session_state.traj_zip_bytes = None
            
            with st.spinner("Calculating RMSD, RMSF, and Rg..."):
                payload = {
                    "job_id": job_id_input_traj,
                    "rmsd_group": rmsd_group,
                    "rmsf_group": rmsf_group,
                    "rg_group": rg_group
                }
                
                try:
                    response = requests.post(f"{API_URL}/analysis/trajectory", json=payload)
                    if response.status_code == 200:
                        task_id = response.json()["task_id"]
                        status_placeholder = st.empty()
                        
                        while True:
                            status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                            status = status_res["status"]
                            
                            if status == "SUCCESS":
                                task_result = status_res.get("result", {})
                                if task_result.get("status") == "success":
                                    status_placeholder.success("✅ GROMACS Analysis Complete!")
                                    
                                    dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                    if dl.status_code == 200:
                                        st.session_state.traj_zip_bytes = dl.content
                                        time.sleep(1)
                                        st.rerun()
                                else:
                                    status_placeholder.error(f"❌ Failed: {task_result.get('message')}")
                                break
                            elif status == "FAILED":
                                status_placeholder.error("❌ Celery Task Failed")
                                break
                            else:
                                status_placeholder.info(f"⏳ Status: {status} (Processing trajectories...)")
                                time.sleep(3)
                    else:
                        st.error("Failed to connect to API.")
                except Exception as e:
                    st.error(f"Connection Error: {str(e)}")

    if st.session_state.get("traj_zip_bytes"):
        st.markdown("---")
        st.download_button(
            label="⬇️ Download GROMACS Results (.zip)",
            data=st.session_state.traj_zip_bytes,
            file_name=f"{job_id_input_traj}_gromacs_analysis.zip",
            mime="application/zip",
            key="traj_download"
        )
        
        st.subheader("Analysis Plots")
        try:
            with zipfile.ZipFile(io.BytesIO(st.session_state.traj_zip_bytes)) as z:
                png_files = [f for f in z.namelist() if f.endswith('.png')]
                if not png_files:
                    st.info("No plots (.png) found in the results zip.")
                else:
                    cols = st.columns(3)
                    for i, png in enumerate(sorted(png_files)):
                        with cols[i % 3]:
                            st.image(z.read(png), caption=os.path.basename(png), use_container_width=True)
        except Exception as e:
            st.warning(f"Could not load images from zip: {e}")

# ==========================================
# TAB 6: TTClust
# ==========================================
with tab6:
    st.header("Step 6: Trajectory Clustering (TTClust)")
    
    job_id_input_ttclust = st.text_input("Job ID", value=st.session_state.get("job_id", ""), key="ttclust_job_id")
    
    if not job_id_input_ttclust:
        st.info("⏳ Waiting for Backmapping (Tab 4) to finish to detect Job ID, or paste your Job ID above to proceed.")
    else:
        st.session_state.job_id = job_id_input_ttclust
        st.success(f"**Linked Job ID:** `{job_id_input_ttclust}`")
        
        st.info("Since TTClust's auto-detection feature has a known package bug, please explicitly specify the number of clusters you want to generate.")
        
        n_clusters = st.number_input("Number of Clusters to Generate", min_value=2, max_value=50, value=5, step=1)
        
        if st.button("Run TTClust Analysis"):
            st.session_state.ttclust_zip_bytes = None
            
            with st.spinner(f"Executing TTClust Trajectory Clustering into {n_clusters} clusters..."):
                payload = {
                    "job_id": job_id_input_ttclust,
                    "n_clusters": n_clusters
                }
                
                try:
                    response = requests.post(f"{API_URL}/analysis/ttclust", json=payload)
                    if response.status_code == 200:
                        task_id = response.json()["task_id"]
                        status_placeholder = st.empty()
                        
                        while True:
                            status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                            status = status_res["status"]
                            
                            if status == "SUCCESS":
                                task_result = status_res.get("result", {})
                                if task_result.get("status") == "success":
                                    status_placeholder.success("✅ TTClust Analysis Complete!")
                                    
                                    dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                    if dl.status_code == 200:
                                        st.session_state.ttclust_zip_bytes = dl.content
                                        time.sleep(1)
                                        st.rerun()
                                else:
                                    status_placeholder.error(f"❌ Failed: {task_result.get('message')}")
                                break
                            elif status == "FAILED":
                                status_placeholder.error("❌ Celery Task Failed")
                                break
                            else:
                                status_placeholder.info(f"⏳ Status: {status} (Processing TTClust pipelines...)")
                                time.sleep(3)
                    else:
                        st.error("Failed to connect to API.")
                except Exception as e:
                    st.error(f"Connection Error: {str(e)}")

    if st.session_state.get("ttclust_zip_bytes"):
        st.markdown("---")
        st.download_button(
            label="⬇️ Download TTClust Results (.zip)",
            data=st.session_state.ttclust_zip_bytes,
            file_name=f"{job_id_input_ttclust}_ttclust_results.zip",
            mime="application/zip",
            key="ttclust_download"
        )
        
        st.subheader("Clustering Plots")
        try:
            with zipfile.ZipFile(io.BytesIO(st.session_state.ttclust_zip_bytes)) as z:
                png_files = [f for f in z.namelist() if f.endswith('.png')]
                if not png_files:
                    st.info("No plots (.png) found in the results zip.")
                else:
                    cols = st.columns(2)
                    for i, png in enumerate(sorted(png_files)):
                        with cols[i % 2]:
                            st.image(z.read(png), caption=os.path.basename(png), use_container_width=True)
        except Exception as e:
            st.warning(f"Could not load images from zip: {e}")

# ==========================================
# TAB 7: PLIP
# ==========================================
with tab7:
    st.header("Step 7: PLIP (Protein-Ligand Interaction Profiler)")
    
    job_id_input = st.text_input("Job ID", value=st.session_state.get("job_id", ""), key="plip_job_id")
    
    if not job_id_input:
        st.info("⏳ Waiting for Backmapping (Tab 4) to finish to detect Job ID, or paste your Job ID above to proceed.")
    else:
        st.session_state.job_id = job_id_input 
        st.success(f"**Linked Job ID:** `{job_id_input}`")
        
        if st.button("Run PLIP Pipeline"):
            st.session_state.plip_zip_bytes = None
            
            with st.spinner("Executing full PLIP pipeline... This may take several minutes."):
                try:
                    response = requests.post(f"{API_URL}/analysis/plip", json={"job_id": job_id_input})
                    
                    if response.status_code == 200:
                        task_id = response.json().get("task_id")
                        if task_id:
                            status_placeholder = st.empty()
                            while True:
                                status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                                status = status_res["status"]
                                
                                if status == "SUCCESS":
                                    task_result = status_res.get("result", {})
                                    if task_result.get("status") == "success":
                                        status_placeholder.success("✅ PLIP Pipeline Complete!")
                                        dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                        if dl.status_code == 200:
                                            st.session_state.plip_zip_bytes = dl.content
                                            time.sleep(1)
                                            st.rerun()
                                        else:
                                            st.error(f"❌ Failed to fetch results from backend. Code: {dl.status_code}")
                                    else:
                                        status_placeholder.error(f"❌ Failed: {task_result.get('message')}")
                                    break
                                elif status == "FAILED":
                                    status_placeholder.error("❌ Celery Task Failed")
                                    break
                                else:
                                    status_placeholder.info(f"⏳ Status: {status} (Processing PLIP & Generating Plots...)")
                                    time.sleep(3)
                    else:
                        st.error(f"Error {response.status_code}: {response.text}")
                except Exception as e:
                    st.error(f"Failed to connect to backend: {str(e)}")

    if st.session_state.get("plip_zip_bytes"):
        st.markdown("---")
        st.success("✅ PLIP Analysis Results are ready!")
        
        st.download_button(
            label="⬇️ Download All PLIP Analysis Results (.zip)",
            data=st.session_state.plip_zip_bytes,
            file_name=f"{job_id_input}_plip_results.zip",
            mime="application/zip",
            key="plip_download"
        )
        
        st.markdown("---")
        st.subheader("📊 Results Preview")
        try:
            with zipfile.ZipFile(io.BytesIO(st.session_state.plip_zip_bytes)) as z:
                png_files = [f for f in z.namelist() if f.endswith('.png') or f.endswith('.jpg')]
                csv_files = [f for f in z.namelist() if f.endswith('.csv')]
                
                if png_files:
                    st.markdown("#### 🖼️ Interaction Plots")
                    cols = st.columns(2) 
                    for i, png in enumerate(png_files):
                        with cols[i % 2]:
                            st.image(z.read(png), caption=os.path.basename(png), use_container_width=True)
                            
                st.markdown("---")
                if csv_files:
                    st.markdown("#### 📑 Data Tables")
                    for csv_file in csv_files:
                        with st.expander(f"📄 View {os.path.basename(csv_file)}"):
                            try:
                                import pandas as pd
                                df = pd.read_csv(io.BytesIO(z.read(csv_file)))
                                st.dataframe(df, use_container_width=True)
                            except Exception:
                                st.write("Could not preview this CSV.")
        except Exception as e:
            st.error(f"Error reading zip file: {e}")

# ==========================================
# TAB 8: MM-GBSA
# ==========================================
with tab8:
    st.header("Step 8: MM-GBSA & Decomposition Analysis")
    
    job_id_input_mmgbsa = st.text_input("Job ID", value=st.session_state.get("job_id", ""), key="mmgbsa_job_id")
    
    if not job_id_input_mmgbsa:
        st.info("⏳ Waiting for Backmapping (Tab 4) to finish to detect Job ID, or paste your Job ID above to proceed.")
    else:
        st.session_state.job_id = job_id_input_mmgbsa 
        st.success(f"**Linked Job ID:** `{job_id_input_mmgbsa}`")
        
        col1, col2 = st.columns(2)
        with col1:
            protein_idx = st.number_input("Protein Group Index", value=1, step=1)
            start_frame = st.number_input("Start Frame", value=1, step=1, min_value=1) 
        with col2:
            ligand_idx = st.number_input("Ligand Group Index", value=13, step=1)
            end_frame = st.number_input("End Frame (0 for auto-detect)", value=0, step=1, min_value=0)
            
        if st.button("Run MM-GBSA Analysis"):
            st.session_state.mmgbsa_zip_bytes = None
            st.session_state.mmgbsa_summary = ""
            
            with st.spinner("Executing MM-GBSA Analysis... This requires significant time."):
                payload = {
                    "job_id": job_id_input_mmgbsa,
                    "ligand_name": st.session_state.get("ligand_code", "A20"),
                    "protein_idx": protein_idx,
                    "ligand_idx": ligand_idx,
                    "start_frame": start_frame, 
                    "end_frame": end_frame
                }
                
                try:
                    response = requests.post(f"{API_URL}/analysis/mmpbsa", json=payload)
                    if response.status_code == 200:
                        task_id = response.json()["task_id"]
                        status_placeholder = st.empty()
                        
                        while True:
                            status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                            status = status_res["status"]
                            
                            if status == "SUCCESS":
                                task_result = status_res.get("result", {})
                                if task_result.get("status") == "success":
                                    status_placeholder.success("✅ MM-GBSA Complete!")
                                    st.session_state.mmgbsa_summary = task_result.get("summary_data", "")
                                    
                                    dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                    if dl.status_code == 200:
                                        st.session_state.mmgbsa_zip_bytes = dl.content
                                        time.sleep(1)
                                        st.rerun()
                                else:
                                    status_placeholder.error(f"❌ Failed: {task_result.get('message')}")
                                break
                            elif status == "FAILED":
                                status_placeholder.error("❌ Celery Task Failed")
                                break
                            else:
                                status_placeholder.info(f"⏳ Status: {status} (Processing trajectories...)")
                                time.sleep(3)
                    else:
                        st.error("Failed to connect to API.")
                except Exception as e:
                    st.error(f"Connection Error: {str(e)}")

    if st.session_state.get("mmgbsa_zip_bytes"):
        st.markdown("---")
        st.download_button(
            label="⬇️ Download MM-GBSA Results (.zip)",
            data=st.session_state.mmgbsa_zip_bytes,
            file_name=f"{job_id_input_mmgbsa}_mmgbsa_results.zip",
            mime="application/zip",
            key="mmgbsa_download"
        )
        
        if st.session_state.mmgbsa_summary:
            st.subheader("Delta (Complex - Receptor - Ligand)")
            st.code(st.session_state.mmgbsa_summary, language="text")
        
        st.subheader("Decomposition Plots")
        try:
            with zipfile.ZipFile(io.BytesIO(st.session_state.mmgbsa_zip_bytes)) as z:
                png_files = [f for f in z.namelist() if f.endswith('.png') or f.endswith('.jpg')]
                if not png_files:
                    st.info("No plots (.png) found in the results zip. Check if the python scripts generated them.")
                else:
                    cols = st.columns(2)
                    for i, png in enumerate(png_files):
                        with cols[i % 2]:
                            st.image(z.read(png), caption=os.path.basename(png), use_container_width=True)
        except Exception as e:
            st.warning(f"Could not load images from zip: {e}")

# ==========================================
# TAB 9: PCA and FEL Analysis
# ==========================================
with tab9:
    st.header("Step 9: Principal Component Analysis (PCA) & FEL")
    
    job_id_input_pca = st.text_input("Job ID", value=st.session_state.get("job_id", ""), key="pca_job_id")
    
    if not job_id_input_pca:
        st.info("⏳ Waiting for Backmapping (Tab 4) to finish to detect Job ID, or paste your Job ID above to proceed.")
    else:
        st.session_state.job_id = job_id_input_pca
        st.success(f"**Linked Job ID:** `{job_id_input_pca}`")
        
        st.markdown("### GROMACS Group Selection")
        backbone_group = st.number_input("Backbone Group Index (Default: 4)", value=4, step=1)
            
        if st.button("Run PCA & FEL Analysis"):
            st.session_state.pca_fel_zip_bytes = None
            
            with st.spinner("Executing Covariance, Eigenvectors, and Free Energy Landscape... This may take a while."):
                payload = {
                    "job_id": job_id_input_pca,
                    "backbone_group": backbone_group
                }
                
                try:
                    response = requests.post(f"{API_URL}/analysis/pca-fel", json=payload)
                    if response.status_code == 200:
                        task_id = response.json()["task_id"]
                        status_placeholder = st.empty()
                        
                        while True:
                            status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                            status = status_res["status"]
                            
                            if status == "SUCCESS":
                                task_result = status_res.get("result", {})
                                if task_result.get("status") == "success":
                                    status_placeholder.success("✅ PCA and FEL Analysis Complete!")
                                    
                                    dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                    if dl.status_code == 200:
                                        st.session_state.pca_fel_zip_bytes = dl.content
                                        time.sleep(1)
                                        st.rerun()
                                else:
                                    status_placeholder.error(f"❌ Failed: {task_result.get('message')}")
                                break
                            elif status == "FAILED":
                                status_placeholder.error("❌ Celery Task Failed")
                                break
                            else:
                                status_placeholder.info(f"⏳ Status: {status} (Processing PCA & FEL pipelines...)")
                                time.sleep(3)
                    else:
                        st.error("Failed to connect to API.")
                except Exception as e:
                    st.error(f"Connection Error: {str(e)}")

    if st.session_state.get("pca_fel_zip_bytes"):
        st.markdown("---")
        st.download_button(
            label="⬇️ Download PCA & FEL Results (.zip)",
            data=st.session_state.pca_fel_zip_bytes,
            file_name=f"{job_id_input_pca}_pca_fel_results.zip",
            mime="application/zip",
            key="pca_fel_download"
        )
        
        st.subheader("Generated Plots")
        try:
            with zipfile.ZipFile(io.BytesIO(st.session_state.pca_fel_zip_bytes)) as z:
                png_files = [f for f in z.namelist() if f.endswith('.png')]
                if not png_files:
                    st.info("No plots (.png) found in the results zip.")
                else:
                    cols = st.columns(2)
                    for i, png in enumerate(sorted(png_files)):
                        with cols[i % 2]:
                            st.image(z.read(png), caption=os.path.basename(png), use_container_width=True)
        except Exception as e:
            st.warning(f"Could not load images from zip: {e}")
