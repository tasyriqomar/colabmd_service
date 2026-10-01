from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from tasks.worker import run_plip_analysis_task, run_mmpbsa_task, run_trajectory_analysis_task, run_pca_fel_task, run_ttclust_task # Added import

router = APIRouter(prefix="/analysis", tags=["Analysis"])

class AnalysisRequest(BaseModel):
    job_id: str

@router.post("/plip")
async def start_plip(req: AnalysisRequest):
    try:
        # Trigger the task in Celery using .delay()
        task = run_plip_analysis_task.delay(req.job_id)
        
        # Return the task_id so Streamlit can track the progress
        return {
            "message": "PLIP interaction profiling queued.", 
            "task_id": task.id
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# Add to colabmd_service/routers/analysis.py
from tasks.worker import run_mmpbsa_task # Ensure this is imported

# --- UPDATE THE PYDANTIC MODEL ---
class MMPBSARequest(BaseModel):
    job_id: str
    ligand_name: str
    protein_idx: int = 1
    ligand_idx: int = 13
    start_frame: int = 1       # <-- ADDED
    end_frame: int = 0         # <-- ADDED

@router.post("/mmpbsa")
async def start_mmpbsa(req: MMPBSARequest):
    try:
        # --- PASS THE NEW VARIABLES TO THE WORKER ---
        task = run_mmpbsa_task.delay(
            req.job_id, 
            req.ligand_name, 
            req.protein_idx, 
            req.ligand_idx,
            req.start_frame,   # <-- ADDED
            req.end_frame      # <-- ADDED
        )
        return {"task_id": task.id, "message": "MM-GBSA task submitted successfully."}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# NEW: Trajectory Analysis Request
class TrajectoryRequest(BaseModel):
    job_id: str
    rmsd_group: int = 4
    rmsf_group: int = 1
    rg_group: int = 1

@router.post("/trajectory")
async def start_trajectory_analysis(req: TrajectoryRequest):
    try:
        task = run_trajectory_analysis_task.delay(
            req.job_id, 
            req.rmsd_group, 
            req.rmsf_group, 
            req.rg_group
        )
        return {"task_id": task.id, "message": "Trajectory analysis task submitted successfully."}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# NEW: PCA & FEL Request Model
class PCAFELRequest(BaseModel):
    job_id: str
    backbone_group: int = 4

@router.post("/pca-fel")
async def start_pca_fel(req: PCAFELRequest):
    try:
        task = run_pca_fel_task.delay(
            req.job_id, 
            req.backbone_group
        )
        return {"task_id": task.id, "message": "PCA & FEL task submitted successfully."}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# NEW: TTClust Request
class TTClustRequest(BaseModel):
    job_id: str
    n_clusters: int = 5  # <-- ADDED FIELD

@router.post("/ttclust")
async def start_ttclust(req: TTClustRequest):
    try:
        task = run_ttclust_task.delay(req.job_id, req.n_clusters)  # <-- PASS TO WORKER
        return {"task_id": task.id, "message": "TTClust task submitted successfully."}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))