from fastapi import APIRouter, UploadFile, File, Form, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel
import time
import os
from pathlib import Path
from tasks.worker import generate_mol2_task, run_aa_prep_task, run_cg_prep_task, celery_app
from pydantic import BaseModel
from tasks.worker import run_cgmd_prep_workflow
from tasks.worker import run_backmapping_workflow

router = APIRouter()
SHARED_DIR = Path("/app/shared_data")

# 1. Add ligand_name to the Mol2 Request
class Mol2Request(BaseModel):
    task_id: str
    ligand_name: str

@router.post("/generate-mol2")
def generate_mol2(req: Mol2Request):
    docking_task = celery_app.AsyncResult(req.task_id)
    if not docking_task.ready() or docking_task.result.get("status") != "success":
        raise HTTPException(status_code=400, detail="Docking incomplete.")
        
    ligaa_pdb_path = docking_task.result.get("result_file")
    
    job_id = f"mol2_{int(time.time())}"
    wd = SHARED_DIR / job_id
    wd.mkdir(parents=True, exist_ok=True)
    
    # Pass the dynamic ligand_name here
    task = generate_mol2_task.delay(str(wd), ligaa_pdb_path, req.ligand_name)
    return {"message": "MOL2 queued.", "task_id": task.id}

@router.post("/run-aa-prep")
async def run_aa_prep(
    ligand_name: str = Form(...),  # 2. Accept ligand_name as a Form field
    mol2_file: UploadFile = File(...),
    str_file: UploadFile = File(...)
):
    job_id = f"aa_{int(time.time())}"
    wd = SHARED_DIR / job_id
    wd.mkdir(parents=True, exist_ok=True)
    
    mol2_path = wd / mol2_file.filename
    str_path = wd / str_file.filename
    
    with open(mol2_path, "wb") as f: f.write(await mol2_file.read())
    with open(str_path, "wb") as f: f.write(await str_file.read())
        
    # Pass the dynamic ligand_name here
    task = run_aa_prep_task.delay(str(wd), str(mol2_path), str(str_path), ligand_name)
    return {"message": "AA prep queued.", "task_id": task.id}

@router.post("/run-cg-prep")
async def run_cg_prep(
    ligand_name: str = Form(...),
    smiles: str = Form(...)
):
    job_id = f"cg_{int(time.time())}"
    wd = SHARED_DIR / job_id
    wd.mkdir(parents=True, exist_ok=True)
    
    task = run_cg_prep_task.delay(str(wd), smiles, ligand_name)
    return {"message": "CG prep queued.", "task_id": task.id}

@router.get("/download-prep-result/{task_id}")
def download_prep_result(task_id: str):
    task_result = celery_app.AsyncResult(task_id)
    if task_result.ready() and task_result.result.get("status") == "success":
        file_path = task_result.result.get("result_file")
        if file_path and os.path.exists(file_path):
            mt = "application/zip" if file_path.endswith(".zip") else "chemical/x-mol2"
            return FileResponse(path=file_path, filename=Path(file_path).name, media_type=mt)
    return {"error": "File not ready or task failed."}

class CGMDPrepRequest(BaseModel):
    job_id: str
    protein_pdb: str
    ligand_itp: str

@router.post("/run-cgmd-prep")
async def start_cgmd_prep(
    protein_pdb: UploadFile = File(...),
    ligand_itp: UploadFile = File(...),
    ligand_gro: UploadFile = File(...),
    ligand_code: str = Form(...)
):
    # 1. Create a fresh, unique folder for this CGMD job
    timestamp = str(int(time.time()))
    work_dir = f"/app/shared_data/cgmd_{timestamp}"
    os.makedirs(work_dir, exist_ok=True)
    
    # 2. Save the uploaded files
    for upload in [protein_pdb, ligand_itp, ligand_gro]:
        file_path = os.path.join(work_dir, upload.filename)
        with open(file_path, "wb") as f:
            f.write(await upload.read())
        
    # 3. Dispatch to the Celery worker
    task = run_cgmd_prep_workflow.delay(
        work_dir, 
        protein_pdb.filename, 
        ligand_itp.filename, 
        ligand_gro.filename, 
        ligand_code
    )
    
    return {"message": "CGMD Prep queued", "task_id": task.id}

@router.post("/run-backmapping")
async def start_backmapping(
    protein_pdb: UploadFile = File(...),
    ligand_map: UploadFile = File(...),
    xtc_file: UploadFile = File(...),
    tpr_file: UploadFile = File(...),
    mdp_file: UploadFile = File(...),
    ndx_file: UploadFile = File(...),
    lig_itp_file: UploadFile = File(...),
    lig_prm_file: UploadFile = File(...),
    lig_cg_itp_file: UploadFile = File(...),
    system_gro: UploadFile = File(...),
    ligand_code: str = Form(...),
    loop_from: int = Form(...),
    loop_until: int = Form(...)
):
    timestamp = str(int(time.time()))
    work_dir = f"/app/shared_data/backmap_{timestamp}"
    os.makedirs(work_dir, exist_ok=True)
    
    # Save all files to the working directory
    file_list = [protein_pdb, ligand_map, xtc_file, tpr_file, mdp_file, ndx_file, lig_itp_file, lig_prm_file, lig_cg_itp_file, system_gro]
    for upload in file_list:
        file_path = os.path.join(work_dir, upload.filename)
        with open(file_path, "wb") as f:
            f.write(await upload.read())
            
    task = run_backmapping_workflow.delay(
        work_dir, protein_pdb.filename, ligand_map.filename, 
        xtc_file.filename, tpr_file.filename, mdp_file.filename, ndx_file.filename, 
        lig_itp_file.filename, lig_prm_file.filename, lig_cg_itp_file.filename, system_gro.filename,
        ligand_code, loop_from, loop_until
    )
    
    return {"message": "Backmapping queued", "task_id": task.id}