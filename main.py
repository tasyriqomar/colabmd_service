from fastapi.responses import FileResponse
from fastapi import FastAPI, UploadFile, File, Form
from typing import Optional
from tasks.worker import run_vina_docking, celery_app
import os
import time

# 1. Import the routers from your routers folder
from routers import preparation
from routers import analysis  # <-- ADDED THIS IMPORT

app = FastAPI(title="ColabMD-Edu Web Service")

SHARED_DIR = "/app/shared_data"
os.makedirs(SHARED_DIR, exist_ok=True)

# 2. Register the routers right here!
app.include_router(preparation.router, prefix="/api/v1")
app.include_router(analysis.router, prefix="/api/v1")  # <-- ADDED THIS TO REGISTER TAB 5

@app.post("/api/v1/run-docking")
async def run_docking(
    ligand_code: str = Form(...),
    protein: UploadFile = File(...), 
    ligand: Optional[UploadFile] = File(None),
    ligand_smiles: Optional[str] = Form(None)
):
    # Create a unique job directory
    job_id = f"job_{int(time.time())}"
    job_dir = os.path.join(SHARED_DIR, job_id)
    os.makedirs(job_dir, exist_ok=True)
    
    # Save the protein
    prot_path = os.path.join(job_dir, protein.filename)
    with open(prot_path, "wb") as f: f.write(await protein.read())
    
    # Save ligand if uploaded
    lig_path = None
    if ligand:
        lig_path = os.path.join(job_dir, ligand.filename)
        with open(lig_path, "wb") as f: f.write(await ligand.read())
        
    # Dispatch task with the new parameters
    task = run_vina_docking.delay(prot_path, lig_path, job_dir, ligand_code, ligand_smiles)
    
    return {"message": "Docking queued.", "task_id": task.id}

@app.get("/api/v1/status/{task_id}")
def get_status(task_id: str):
    task_result = celery_app.AsyncResult(task_id)
    return {
        "status": task_result.status,
        "result": task_result.result if task_result.ready() else None
    }

@app.get("/api/v1/download/{task_id}")
def download_result(task_id: str):
    task_result = celery_app.AsyncResult(task_id)
    
    # Check if the task was successful and has a result_file
    if task_result.ready() and task_result.status == "SUCCESS":
        file_path = task_result.result.get("result_file")
        
        # Verify the file actually exists on the shared drive
        if file_path and os.path.exists(file_path):
            return FileResponse(
                path=file_path, 
                filename="ligand_pose1.pdb", 
                media_type="chemical/x-pdb"
            )
            
    return {"error": "File not ready, failed, or missing."}

@app.get("/api/v1/download-complex/{task_id}")
def download_complex(task_id: str):
    task_result = celery_app.AsyncResult(task_id)
    
    if task_result.ready() and task_result.status == "SUCCESS":
        result_data = task_result.result
        if isinstance(result_data, dict) and result_data.get("status") == "success":
            file_path = result_data.get("complex_file")
            
            if file_path and os.path.exists(file_path):
                filename = os.path.basename(file_path)
                return FileResponse(
                    path=file_path, 
                    filename=filename, 
                    media_type="chemical/x-pdb"
                )
            
    return {"error": "File not ready, failed, or missing."}

@app.get("/api/v1/download-prep-result/{task_id}")
def download_prep_result(task_id: str):
    task_result = celery_app.AsyncResult(task_id)
    if task_result.ready() and task_result.status == "SUCCESS":
        result_data = task_result.result
        if isinstance(result_data, dict) and "result_file" in result_data:
            file_path = result_data.get("result_file")
            if file_path and os.path.exists(file_path):
                filename = os.path.basename(file_path)
                return FileResponse(path=file_path, filename=filename, media_type="application/zip")
    return {"error": "File not ready, failed, or missing."}