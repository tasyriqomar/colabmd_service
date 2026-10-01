from fastapi import APIRouter, BackgroundTasks, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel
import os
from tasks.worker import run_martinize, update_topology, build_simulation_box

# Adjust import based on your actual project structure
from colabmd_service.tasks.worker import run_backmapping_task

router = APIRouter()

class MartinizeRequest(BaseModel):
    job_id: str
    protein_pdb: str
    ligand_itp: str

@router.post("/run-cgmd-prep")
async def prep_cgmd(req: MartinizeRequest):
    job_dir = os.path.join("shared_data", req.job_id)
    
    try:
        # Step 3.1: Run Martinize2[cite: 1]
        martinize_res = run_martinize(job_dir, req.protein_pdb)
        
        # Step 3.2: Update Topology[cite: 1]
        update_topology(job_dir, martinize_res["top"], req.ligand_itp)
        
        # Step 3.3: Generate Box with editconf[cite: 1]
        # In the notebook, this targets the initial structure[cite: 1]
        box_res = build_simulation_box(job_dir, martinize_res["cg_pdb"], "system_box.gro")
        
        return {
            "message": "CGMD Preparation Complete",
            "files": {
                "topology": martinize_res["top"],
                "cg_protein": martinize_res["cg_pdb"],
                "system_box": box_res["gro"]
            }
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/api/simulation/backmap")
async def trigger_backmapping(base_dir: str, background_tasks: BackgroundTasks):
    """
    Endpoint to start the backmapping process.
    """
    if not os.path.exists(base_dir):
        raise HTTPException(status_code=404, detail="Base directory not found.")
    
    # Run the backmapping command (can be moved to background task if it takes too long)
    try:
        zip_path = run_backmapping_task(base_dir)
        return {"status": "success", "message": "Backmapping completed.", "zip_path": zip_path}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Backmapping failed: {str(e)}")

@router.get("/api/simulation/download_backmap")
async def download_backmapped_results(base_dir: str):
    """
    Endpoint to download the generated ZIP file.
    """
    zip_path = os.path.join(base_dir, "backmapped_results.zip")
    
    if not os.path.exists(zip_path):
        raise HTTPException(status_code=404, detail="Zip file not found. Please run backmapping first.")
        
    return FileResponse(
        path=zip_path,
        media_type='application/zip',
        filename="backmapped_results.zip"
    )