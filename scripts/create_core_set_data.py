import pandas as pd
import math

tmp = open('../data/exp_data/CoreSet.dat', 'r')
lines = tmp.readlines()

PDB_code_list = []
Affinity_Data_Type_list = []
Affinity_Data_Value_list = []
pKi_pKd_IC50_list = []
Delta_G_list = []
K_reciprocal_list = []

for line in lines:
    test = str(line).split()
    # print(test)
    if test[0] == "#code":
        pass
    else:
        PDB_code_list.append(test[0])
        Affinity_Data_Type, Affinity_Data_Value = test[4].split(sep="=")
        Affinity_Data_Type_list.append(Affinity_Data_Type)

        if Affinity_Data_Value[-2:] == "nM":
            Affinity_Data_Value = Affinity_Data_Value[0:-2]
        elif Affinity_Data_Value[-2:] == "uM":
            Affinity_Data_Value = Affinity_Data_Value[0:-2]
            Affinity_Data_Value = float(Affinity_Data_Value)*1000
            Affinity_Data_Value = int(Affinity_Data_Value)
        elif Affinity_Data_Value[-2:] == "pM":
            Affinity_Data_Value = Affinity_Data_Value[0:-2]
            Affinity_Data_Value = float(Affinity_Data_Value)/1000
            Affinity_Data_Value = round(Affinity_Data_Value,5)
        elif Affinity_Data_Value[-2:] == "mM":
            Affinity_Data_Value = Affinity_Data_Value[0:-2]
            Affinity_Data_Value = float(Affinity_Data_Value)*1000000
            Affinity_Data_Value = int(Affinity_Data_Value)
        
        Affinity_Data_Value_list.append(Affinity_Data_Value)
        
        pKi_pKd_IC50_list.append(test[3])

        Delta_G = round(0.001987*298*math.log(float(Affinity_Data_Value)/1000000000),2)
        Delta_G_list.append(Delta_G)

        K_reciprocal = round(1/float(Affinity_Data_Value),5)
        K_reciprocal_list.append(K_reciprocal)

# print(len(PDB_code_list),len(Affinity_Data_Type_list),len(Affinity_Data_Value_list),len(pKi_pKd_IC50_list),len(Delta_G_list),len(K_reciprocal_list))
data = {'PDB code':PDB_code_list,'Affinity Data Type':Affinity_Data_Type_list,'Affinity Data Value':Affinity_Data_Value_list,'pKd pKi pIC50':pKi_pKd_IC50_list,'delta G':Delta_G_list,'1/K':K_reciprocal_list}
df = pd.DataFrame(data)
df.to_csv("../data/exp_data/PDBbind_core_set_all.csv")


df_kd = df[df['Affinity Data Type'] == "Kd"]
df_ki = df[df['Affinity Data Type'] == "Ki"]

df_kd.to_csv("../data/exp_data/PDBbind_core_set_kd.csv")
df_ki.to_csv("../data/exp_data/PDBbind_core_set_ki.csv")

print('Done')