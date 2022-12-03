import numpy as np
import pandas as pd


class CellType:
    def __init__(self, cluster : pd.DataFrame, marker : list):
        self.cluster = cluster
        self.cellType = ["Fibroblasts", "Epithelium", "Bcells", "Monocytes",
                        "Macrophages", "IL17", "T cells", "GamadeltaT"] 
        self.cellTypeDict = {}

    def getCellType(self, cluster : int):
        return self.cellTypeDict(cluster)

    def setCellType(self):
        pass
        
if __name__ == "__main__":
    print(len(CellType(None, None).cellType))