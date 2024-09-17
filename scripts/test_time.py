import pandas as pd
import phantomdragon.functions as ph
import time

modeltype_list = []
scoretype_list = []
datatype_list = []
mae_list = []
mse_list = []
sd_list = []
pear_list = []
coef_list = []
confidence_interval_list =[]
add_info_list = []
timelist =[]

datatypes = ["all","ki","kd"]
modeltypes = ["linearRegression","Ridge","Lasso","ElasticNet","SVR","DecisionTree","RandomForest","XGBoost"]
scoretypes = ["pKd pKi pIC50"]
descriptor = ["GRADE","X-GRADE"]
test_val = ["core","general"]

for modeltype in modeltypes:
    for score in scoretypes:
        for k in datatypes:
            for des in descriptor:
                for t in test_val:

                    tmp_modeltype_list = []
                    tmp_scoretype_list = []
                    tmp_datatype_list = []
                    tmp_mae_list = []
                    tmp_mse_list = []
                    tmp_sd_list = []
                    tmp_pear_list = []
                    tmp_coef_list = []
                    tmp_confidence_interval_list =[]
                    tmp_add_info_list = []
                    tmp_timelist =[]
                    blub = ph.parameterCollector(
                                        add_information=f"{des}_{t}",
                                        modeltype=modeltype,
                                        scoretype=score,
                                    )
                    blub = ph.parameterCollector(add_information=f"{des}",modeltype=modeltype,scoretype=score)
                    x_train,x_test,y_train,y_test = ph.prepare_data(score,
                    f"../data/Descriptors/PDBbind_refined_set_{des}.csv",
                    f"../data/Descriptors/PDBbind_{t}_set_{des}.csv",
                    f"../data/exp_data/PDBbind_refined_set_{k}.csv",
                    f"../data/exp_data/PDBbind_{t}_set_all.csv",
                    f"{des}")
                    blub.set_trainingdata(x_train,y_train)
                    blub.set_testingdata(x_test,y_test)
                    blub.set_datatype(k)
                    blub.train_and_save_model(savepath="../models/")
                    for i in range(10):
                        start = time.time()
                        blub.phantomtest(loadpath="../models/")
                        end = time.time()
                        # blub.plot_phantomtest("../plots/")
                        modeltype,scoret,datatype,mae,mse,sd,pearsonr,confidence_interval,r_2,add_info = blub.get_stats()
                        tmp_modeltype_list.append(modeltype)
                        tmp_scoretype_list.append(scoret)
                        tmp_datatype_list.append(datatype)
                        tmp_mae_list.append(mae)
                        tmp_mse_list.append(mse)
                        tmp_sd_list.append(sd)
                        tmp_pear_list.append(pearsonr)
                        tmp_confidence_interval_list.append(confidence_interval)
                        tmp_coef_list.append(r_2)
                        tmp_add_info_list.append(add_info)
                        tmp_timelist.append(end-start)
                    modeltype_list.append(tmp_modeltype_list[0])
                    scoretype_list.append(tmp_scoretype_list[0])
                    datatype_list.append(tmp_datatype_list[0])
                    mae_list.append(sum(tmp_mae_list)/10)
                    mse_list.append(sum(tmp_mse_list)/10)
                    coef_list.append(sum(tmp_coef_list)/10)
                    sd_list.append(sum(tmp_sd_list)/10)
                    pear_list.append(sum(tmp_pear_list)/10)
                    confidence_interval_list.append(tmp_confidence_interval_list[0])
                    add_info_list.append(tmp_add_info_list[0])
                    timelist.append(sum(tmp_timelist)/10)


                print(k,modeltype,score,"done")
                print("Training Set size",len(x_train))
                print("Testing Set size",len(x_test))

#print(len(modeltype_list),len(featuretype_list),len(datatype_list),len(mse_list),len(pear_list),len(coef_list),len(add_info_list))
data = {'Modeltype':modeltype_list,'Scoretype':scoretype_list,'Datatype':datatype_list,'Mean absolute error (mae)':mae_list,'Mean squared error (mse)':mse_list,'Standard Diviation (SD)':sd_list,'Pearson correlation coefficient (r)':pear_list,'90% Confidence interval':confidence_interval_list,'Coefficient of determination (r²)':coef_list,'add. information':add_info_list,'Time':timelist}
df = pd.DataFrame(data)
df.to_csv("../results/time_GRADE.csv")