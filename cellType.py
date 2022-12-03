import numpy as np


class CellType:
    def __init__(self, cluster : np.ndarray):
        self.cluster = cluster
        self.cellType = ["Fibroblasts", "Epithelium", "Bcells", "Monocytes",
                        "Macrophages", "IL17", "T cells", "GamadeltaT"] 
        self.cellTypeDict = {}

    def getCellType(self, cluster : int):
        return self.cellTypeDict(cluster)
