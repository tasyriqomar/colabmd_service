FROM python:3.12-slim

WORKDIR /app

# 1. Install system dependencies (removed gromacs from here)
RUN apt-get update && apt-get install -y \
    openbabel \
    wget \
    git \
    dssp \
    && rm -rf /var/lib/apt/lists/*

# 1.5 Install Miniconda and Dependencies (Auto-detects Architecture for Apple/Intel)
RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "x86_64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"; \
    elif [ "$ARCH" = "aarch64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"; \
    else \
        echo "Unsupported architecture: $ARCH" && exit 1; \
    fi && \
    wget $MINICONDA_URL -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

ENV PATH="/opt/conda/bin:$PATH"

# Install mamba for faster resolution, then install GROMACS, AmberTools, and gmx_MMPBSA
# Accept Anaconda Terms of Service, install mamba, then install specific GROMACS version
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main || true && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r || true && \
    conda install -c conda-forge mamba -y && \
    mamba install -c conda-forge "gromacs=2019.6" ambertools mpi4py gmx_mmpbsa -y

# 2. Download and install AutoDock Vina (from your notebook)
RUN wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64 -O /usr/local/bin/vina && \
    chmod +x /usr/local/bin/vina

# 3. Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 4. Fetch Ligand Prep Scripts and Forcefields
WORKDIR /app/tools

# Download the CGenFF to GROMACS conversion script
RUN wget https://raw.githubusercontent.com/Lemkul-Lab/cgenff_charmm2gmx/main/cgenff_charmm2gmx.py

# Download the Python 3 & Martini 3 compatible insane3.py and save it as /app/tools/insane.py
RUN wget https://raw.githubusercontent.com/pstansfeld/MemProtMD/main/insane3.py -O /app/tools/insane.py

# Clone the ColabMD-Edu repo to extract the CHARMM36 forcefield, then clean up
RUN git clone https://github.com/tasyriqomar/ColabMD-Edu_Protein-Ligand.git && \
    mv ColabMD-Edu_Protein-Ligand/Ligand-preparation/CG/charmm36-feb2026_cgenff-5.0.ff . && \
    rm -rf ColabMD-Edu_Protein-Ligand

# Download Martini 3 parameterization script and its required dat file
RUN wget https://raw.githubusercontent.com/tasyriqomar/ColabMD-Edu_Protein-Ligand/refs/heads/main/Ligand-preparation/CG/cg_param_m3_fixed.py && \
    wget https://raw.githubusercontent.com/tasyriqomar/ColabMD-Edu_Protein-Ligand/refs/heads/main/Ligand-preparation/CG/fragments-exp.dat && \
    sed -i 's/int(shared)/int(np.ravel(shared)[0])/g' cg_param_m3_fixed.py

# Download the entire Backmapping folder (including all map files) from your repository
RUN git clone https://github.com/tasyriqomar/ColabMD-Edu_Protein-Ligand.git temp_repo && \
    mv temp_repo/Backmapping /app/tools/backmapping && \
    chmod +x /app/tools/backmapping/backward.py && \
    rm -rf temp_repo

# Reset the working directory for the application code
WORKDIR /app

# 5. Copy the app
COPY . .

CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0", "--port", "8000"]
