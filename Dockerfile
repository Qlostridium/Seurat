FROM qrouchon/r-base-plus
RUN apt-get update

RUN apt-get install -y libgit2-dev libpng-dev python3 python3-pip && rm -rf /var/lib/apt/lists/*

# Python dependencies
RUN pip3 install --upgrade setuptools
RUN pip3 install numpy scipy scikit-learn numba umap-learn

RUN R -e 'BiocManager::install(c("SingleCellExperiment","multtest"),update = TRUE, ask = FALSE);install.packages(c("Seurat"));'
