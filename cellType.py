import numpy as np
import pandas as pd
import pickle
from sklearn_som.som import SOM
import pandas as pd

class CellType:
    def __init__(self, som : SOM, df : pd.DataFrame, col : list):
        self.som = som
        self.df = df
        self.df["predictions"] = self.som.predict(self.df.values)
        self.df2 = df.drop('predictions', axis=1)[col].iloc[0:0]
        self.col = col
        self.importantcol = ["PanKeratin_Mean", "Ecad_Mean","aSMA_Mean",
                            "CD45_Mean", "CD3_Mean","CD68_Mean","CD20_Mean",
                            "CD4_Mean", "IL17a_Mean"]
        self.labels = ['epithelium', 'fibroblasts', 'T-Cells', 'macrophage',
                        'B-Cells', 'IL17']
        self.plotMe = None
        self.plotlabel = None

    def proportion(self):
        for i in range(self.som.n * self.som.m):
            mask = self.df['predictions'] == i
            myList = []
            for column in self.df.drop('predictions', axis=1)[mask][self.col]:
                myList.append(self.df[mask][column].mean())
            self.df2.loc[len(self.df2)] = myList
        self.plotMe = self.df2[self.importantcol]
    
    def fixeClusterLabel(self):     
        if self.plotMe is None:
            self.proportion()   
        for col in self.df2[self.importantcol]:
            string = col
            string += "-"
            count = 0
            l = []
            for x in self.df2[col].tolist():
                if x - self.df2[col].mean() > self.df2[col].std()*0.75:
                    string += str(count)
                    string += "-"
                    l.append(True)
                else:
                    l.append(False)
                count = count + 1
            col += '_'    
            self.df2[col] = l

        epithelium_mask = (self.df2['PanKeratin_Mean_'] == True) & (self.df2['Ecad_Mean_'] == True)
        fibroblast_mask = self.df2['aSMA_Mean_'] == True
        tcells_mask = (self.df2['CD45_Mean_'] == True) & (self.df2['CD3_Mean_'] == True)
        macrophage_mask = (self.df2['CD45_Mean_'] == True) & (self.df2['CD68_Mean_'] == True)
        bcells_mask = (self.df2['CD45_Mean_'] == True) & (self.df2['CD20_Mean_'] == True)
        il17_mask1 = (self.df2['CD4_Mean_'] == True) & (self.df2['IL17a_Mean_'] == True)
        il17_mask = tcells_mask&il17_mask1
        masks = [epithelium_mask,fibroblast_mask,tcells_mask,macrophage_mask,bcells_mask,il17_mask]
        for mask, label in zip(masks,self.labels):
            self.df2[label] = 0
            self.df2.loc[mask, label] = 1    
        self.plotlabel = self.df2[self.labels]

    def findLabel(self):
        if self.plotlabel is None:
            self.fixeClusterLabel()
        
        # prediction = self.som.predict(data)

        self.celltype = [[] for i in range(len(self.plotlabel.columns))]
        for i in range(len(self.plotlabel.columns)):
            tr = self.plotlabel[self.plotlabel.columns[i]] == 1
            index = self.plotlabel[tr].index
            self.celltype[i].extend(index)
        
        self.diClusterToLabel = dict(zip(self.labels, self.celltype))

                    
    