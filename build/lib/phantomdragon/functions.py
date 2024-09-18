import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import linear_model
from sklearn import svm
from sklearn import tree
from sklearn import ensemble
from sklearn import preprocessing
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.utils import shuffle
from scipy import stats
import joblib
import xgboost as xgb
import warnings

def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return r_value**2        

def filter_and_sort_features(exp_data, features,identifier="PDB code"):
    """
    Parameters
    ----------
    exp_data : TYPE PANDAS Dataframe
        A pandas Dataframe of the avavible experimental scores.
    features : TYPE PANDAS Dataframe
        A pandas Dataframe of the avaviable features.
        It is assumed that exp_data is corresponding to a subset of features.

    Returns
    -------
    df : TYPE PANDAS Dataframe
        Gives back a sorted subset of score_features.
        This subset corresponds to the subset of the experimental data.
        It is sorted and filtered according to the "PDB code" bracket.

    """
    df = features[features[identifier].isin(exp_data[identifier])]
    df = df.sort_values(identifier)
    df = df.reset_index(drop=True)
    return df


def combine(list1, list2, list3):
    """
    Parameters
    ----------
    list1 : TYPE List
    list2 : TYPE List
    list3 : TYPE List

    Returns
    -------
    combinations: TYPE List
        List that contains a list of all possible combinations between the three input lists
    """
    combinations = []
    for lst1 in list1:
        for lst2 in list2:
            for lst3 in list3:
                combinations.append([lst1, lst2, lst3])
    return combinations


def prepare_data(
    scoretype,
    featurepath_train,
    featurepath_test,
    experimentpath_train,
    experimentpath_test,
    add_information,
    sh=True,
    identifier="PDB code",
    experiment_identifier="Affinity Data Type",
    polynomial=False,
):
    
    if " " in featurepath_train or " " in featurepath_test or " " in experimentpath_train or " " in experimentpath_test:
        featurepath_train = featurepath_train.replace(" ","_")
        featurepath_test = featurepath_test.replace(" ","_")
        experimentpath_train = experimentpath_train.replace(" ","_")
        experimentpath_test = experimentpath_test.replace(" ","_")
    #Ignoring warnings of converters and dtype. Converters override dtype. This is desired behavior
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=pd.errors.ParserWarning)
        features_train = pd.read_csv(featurepath_train, dtype='float64',converters={identifier: str,})
        features_test = pd.read_csv(featurepath_test, dtype='float64',converters={identifier: str})
        experiment_train = pd.read_csv(experimentpath_train, dtype='float64',converters={identifier: str, experiment_identifier: str})
        experiment_test = pd.read_csv(experimentpath_test, dtype='float64',converters={identifier: str, experiment_identifier: str})
    experiment_train = experiment_train.sort_values(identifier)
    experiment_train = experiment_train.reset_index()
    experiment_train = experiment_train.drop("index", axis=1)
    experiment_test = experiment_test.sort_values(identifier)
    experiment_test = experiment_test.reset_index()
    experiment_test = experiment_test.drop("index", axis=1)
    features_train = filter_and_sort_features(experiment_train, features_train)
    features_test = filter_and_sort_features(experiment_test, features_test)
    experiment_train = filter_and_sort_features(features_train, experiment_train)
    experiment_test = filter_and_sort_features(features_test,experiment_test)

    if features_test[identifier].tolist() == experiment_test[identifier].tolist():
        PDB_codes = features_test[identifier].tolist()
    else:
        raise ValueError("Testfeatures and Testlables have different Values")

    if all(x in PDB_codes for x in features_train[identifier]):
        pass
    else:
        bool_list = []
        for i in range(len(PDB_codes)):
            if PDB_codes[i] in list(features_train[identifier]):
                bool_list.append(True)
            else:
                bool_list.append(False)
        PDB_codes = np.array(PDB_codes)
        PDB_codes = list(PDB_codes[bool_list])

    features_train = features_train.set_index(features_train[identifier])
    features_train = features_train.transpose()

    for PDB_code in PDB_codes:
        features_train = features_train.drop(PDB_code, axis=1)

    features_train = features_train.transpose()
    features_train = features_train.reset_index(drop=True)

    experiment_train = experiment_train.set_index(experiment_train[identifier])
    experiment_train = experiment_train.transpose()

    experiment_train = experiment_train.drop(PDB_codes, axis=1)
    experiment_train = experiment_train.transpose()

    scores_test = np.array(experiment_test[scoretype])
    scores_train = np.array(experiment_train[scoretype])

    if add_information == "basic":
        drop = [identifier, " HW-HW_SUM", " HW-HW_MAX", " ES", " VDW_ATT", " VDW_REP"]
    elif add_information == "-w":
        drop = [identifier, " H-H_SUM", " H-H_MAX", " ES", " VDW_ATT", " VDW_REP"]
    elif add_information == "basic-el":
        drop = [identifier, " HW-HW_SUM", " HW-HW_MAX", " VDW_ATT", " VDW_REP"]
    elif add_information == "-w-el":
        drop = [identifier, " H-H_SUM", " H-H_MAX", " VDW_ATT", " VDW_REP"]
    elif add_information == "basic-el-vdw":
        drop = [identifier, " HW-HW_SUM", " HW-HW_MAX"]
    elif add_information == "-w-el-vdw":
        drop = [identifier, " H-H_SUM", " H-H_MAX"]
    elif add_information == "basic-vdw":
        drop = [identifier, " HW-HW_SUM", " HW-HW_MAX", " ES", " VDW_ATT", " VDW_REP"]
    elif add_information == "-w-vdw":
        drop = [identifier, " HW-HW_SUM", " HW-HW_MAX", " ES", " VDW_ATT", " VDW_REP"]
    else:
        drop = [identifier]

    features_train = features_train.drop(drop, axis=1)
    features_test = features_test.drop(drop, axis=1)
    features_train = features_train.to_numpy()
    features_test = features_test.to_numpy()

    if polynomial == True:
        features_train = np.hstack((features_train, features_train**2))
        features_test = np.hstack((features_test, features_test**2))

    scaler = preprocessing.StandardScaler().fit(features_train)

    if sh == True:
        features_train, scores_train = shuffle(features_train, scores_train)
        features_test, scores_test = shuffle(features_test, scores_test)


    return features_train, features_test, scores_train, scores_test


class parameterCollector:
    def __init__(self, modeltype, add_information, scoretype):
        self.modeltype = modeltype
        self.add_information = add_information
        self.scoretype = scoretype
        self.mse = None
        self.r = None
        self.r_2 = None

    def get_data(self):
        x_train, y_train = self.get_trainingdata()
        x_test, y_test = self.get_testingdata()
        return x_train, x_test, y_train, y_test

    def get_trainingdata(self):
        return self.features_train, self.scores_train

    def get_testingdata(self):
        return self.features_test, self.scores_test
    
    def get_predicted_values(self):
        return self.scores_pre

    def set_trainingdata(self, features_train, scores_train):
        self.features_train = features_train
        self.scores_train = scores_train

    def set_testingdata(self, features_test, scores_test):
        self.features_test = features_test
        self.scores_test = scores_test

    def set_datatype(self, datatype):
        self.datatype = datatype

    def set_trainingscoretype(self, trainingscores_type):
        self.trainingscores_type = trainingscores_type

    def load_data(
        self, featurepath, experimentpath, split=0.2, sh=True, identifier="PDB code"
    ):
        if " " in featurepath or " " in experimentpath:
            featurepath = featurepath.replace(" ","_")
            experimentpath = experimentpath.replace(" ","_")
        
        features = pd.read_csv(featurepath, dtype={identifier: str})
        experiment = pd.read_csv(experimentpath, dtype={identifier: str})
        experiment = experiment.sort_values(identifier)
        experiment = experiment.reset_index(drop=True)
        scores = np.array(experiment[self.scoretype])

        if self.add_information == "basic":
            drop = [
                identifier,
                " HW-HW_SUM",
                " HW-HW_MAX",
                " ES",
                " VDW_ATT",
                " VDW_REP",
            ]
        elif self.add_information == "-w":
            drop = [identifier, " H-H_SUM", " H-H_MAX", " ES", " VDW_ATT", " VDW_REP"]
        elif self.add_information == "basic-el":
            drop = [identifier, " HW-HW_SUM", " HW-HW_MAX", " VDW_ATT", " VDW_REP"]
        elif self.add_information == "-w-el":
            drop = [identifier, " H-H_SUM", " H-H_MAX", " VDW_ATT", " VDW_REP"]
        elif self.add_information == "basic-el-vdw":
            drop = [identifier, " HW-HW_SUM", " HW-HW_MAX"]
        elif self.add_information == "-w-el-vdw":
            drop = [identifier, " H-H_SUM", " H-H_MAX"]
        elif self.add_information == "basic-vdw":
            drop = [
                identifier,
                " HW-HW_SUM",
                " HW-HW_MAX",
                " ES",
                " VDW_ATT",
                " VDW_REP",
            ]
        elif self.add_information == "-w-vdw":
            drop = [
                identifier,
                " HW-HW_SUM",
                " HW-HW_MAX",
                " ES",
                " VDW_ATT",
                " VDW_REP",
            ]
        else:
            drop = [identifier]
        
        features = filter_and_sort_features(experiment, features)
        features = features.drop(drop, axis=1)
        features = features.to_numpy()

        if sh == True:
            features, scores = shuffle(features, scores)

        scaler = preprocessing.StandardScaler().fit(features)
        features = scaler.transform(features)

        x_train, x_test, y_train, y_test = train_test_split(
            features, scores, test_size=split
        )

        self.features_train = x_train
        self.features_test = x_test
        self.scores_train = y_train
        self.scores_test = y_test

    def train_and_save_model(self, savepath="", additional_marker="",hyperparametersearch=False):
        if self.modeltype == "linearRegression":
            reg = linear_model.LinearRegression()
        elif self.modeltype == "Ridge":
            reg = linear_model.RidgeCV()
        elif self.modeltype == "Lasso":
            reg = linear_model.LassoCV()
        elif self.modeltype == "ElasticNet":
            reg = linear_model.ElasticNetCV(cv=10, tol=0.00001, n_alphas=1000, eps=0.0001, max_iter=3000)
        elif self.modeltype == "SVR":
            reg = svm.SVR()
        elif self.modeltype == "DecisionTree":
            reg = tree.DecisionTreeRegressor()
        elif self.modeltype == "RandomForest":
            reg = ensemble.RandomForestRegressor()
        elif self.modeltype == "XGBoost":
            reg = xgb.XGBRegressor(max_depth=3,n_estimators=1000)
        else:
            raise ValueError("Unexpected string in modeltype")

        if hyperparametersearch ==True:
            if self.modeltype == "SVR":
                params = {
                    "kernel": ["linear", "poly", "rbf", "sigmoid"],
                    "C": [1, 2.5, 5],
                    "gamma": ["scale", "auto"],
                    "degree": [2, 3, 4],
                    "epsilon": [0.001, 0.01, 0.1, 1, 10, 100],
                }
                gs_reg = GridSearchCV(reg, params)
                gs_reg.fit(self.features_train, self.scores_train)
                reg.set_params(**gs_reg.best_params_)
            elif self.modeltype == "DecisionTree":
                params = {
                    "criterion": [
                        "squared_error",
                        "friedman_mse",
                        "absolute_error",
                        "poisson",
                    ],
                    "splitter": ["best", "random"],
                    "max_features": ["auto", "sqrt", "log2"],
                }
                gs_reg = GridSearchCV(reg, params)
                gs_reg.fit(self.features_train, self.scores_train)
                reg.set_params(**gs_reg.best_params_)
            elif self.modeltype == "RandomForest":
                params = {
                    "n_estimators": [100, 200, 300],
                    "criterion": [
                        "squared_error",
                        "friedman_mse",
                        "absolute_error",
                        "poisson",
                    ],
                    "max_features": ["auto", "sqrt", "log2"],
                }
                gs_reg = GridSearchCV(reg, params)
                gs_reg.fit(self.features_train, self.scores_train)
                reg.set_params(**gs_reg.best_params_)

        pipe = Pipeline([('scaler',preprocessing.StandardScaler()),(self.modeltype,reg)])

        pipe.fit(self.features_train, self.scores_train)

        if "/" in self.scoretype:
            self.scoretype = self.scoretype.replace("/", "div")

        joblib.dump(
            pipe,
            f"{savepath}{self.modeltype}_{self.scoretype}_{self.datatype}_{self.add_information}{additional_marker}.sav",
        )

        if "div" in self.scoretype:
            self.scoretype = self.scoretype.replace("div", "/")

    def phantomtest(
        self,
        testing_features=0,
        testing_scores=0,
        loadpath="",
        confidence_level=0.9,
        additional_marker="",
        return_values=False,
    ):
        
        if testing_features != 0:
            self.features_test = testing_features
        if testing_scores != 0:
            self.scores_test = testing_scores

        if "/" in self.scoretype:
            self.scoretype = self.scoretype.replace("/", "div")

        reg = joblib.load(
            f"{loadpath}{self.modeltype}_{self.scoretype}_{self.datatype}_{self.add_information}{additional_marker}.sav"
        )

        if "div" in self.scoretype:
            self.scoretype = self.scoretype.replace("div", "/")
       
        scores_pre = reg.predict(self.features_test)

        self.scores_pre = scores_pre
        self.mae = mean_absolute_error(self.scores_test, self.scores_pre)
        self.mse = mean_squared_error(self.scores_test, self.scores_pre)
        self.sd = np.std(self.scores_pre)
        self.r = round(stats.pearsonr(self.scores_test, self.scores_pre).statistic, 6)
        conf_int = stats.pearsonr(
            self.scores_test, self.scores_pre
        ).confidence_interval(confidence_level=confidence_level)
        conf_int_low = round(conf_int.low, 3)
        conf_int_high = round(conf_int.high, 3)
        self.conf_int = f"[{conf_int_low} ~ {conf_int_high}]"
        self.spearman_r = round(stats.spearmanr(self.scores_test,self.scores_pre).correlation,6)
        #self.r_2 = round(r2_score(self.scores_test, self.scores_pre),6)
        self.r_2 = round(rsquared(self.scores_test, self.scores_pre),6)

        if return_values == True:
            return self.scores_pre

    def plot_phantomtest(self, savepath,name=None):
        res = stats.linregress(self.scores_test, self.scores_pre)
        plt.plot(self.scores_test, self.scores_pre, 'o',alpha=0.2)

        plt.plot(self.scores_test, res.intercept + res.slope*self.scores_test, 'k', label=f"r\u00b2 = {self.r_2}")

        # plt.title(f"{self.scoretype}, {self.modeltype}, {self.add_information}, {self.datatype}")  
        plt.legend()
        plt.xlabel("y_real")
        plt.ylabel("y_predicted")

        if (
            " " in self.scoretype
            or " " in self.modeltype
            or " " in self.add_information
            or "/" in self.scoretype
        ):
            self.scoretype = self.scoretype.replace(" ", "_")
            self.scoretype = self.scoretype.replace("/", "div")
            self.modeltype = self.modeltype.replace(" ", "_")
            self.add_information = self.add_information.replace(" ", "_")

        if name != None:
            plt.savefig(
                f"{savepath}{name}.png"
            )
        else:
            plt.savefig(
            f"{savepath}{self.datatype}_{self.scoretype}_{self.modeltype}_{self.add_information}.png"
            )

        if (
            "_" in self.scoretype
            or "_" in self.modeltype
            or "_" in self.add_information
            or "div" in self.scoretype
        ):
            self.scoretype = self.scoretype.replace("_", " ")
            self.scoretype = self.scoretype.replace("div", "/")
            self.modeltype = self.modeltype.replace("_", " ")
            self.add_information = self.add_information.replace("_", " ")

        plt.close()

    def phantomscore(self, features_test, loadpath, identifier="PDB code",):
        if isinstance(features_test, str):
            self.features_test = pd.read_csv(features_test, dtype={identifier: str})
        
        PDB_codes = self.features_test[identifier]

        if self.add_information == "basic":
            drop = [
                identifier,
                " HW-HW_SUM",
                " HW-HW_MAX",
                " ES",
                " VDW_ATT",
                " VDW_REP",
            ]
        elif self.add_information == "-w":
            drop = [identifier, " H-H_SUM", " H-H_MAX", " ES", " VDW_ATT", " VDW_REP"]
        elif self.add_information == "basic-el":
            drop = [identifier, " HW-HW_SUM", " HW-HW_MAX", " VDW_ATT", " VDW_REP"]
        elif self.add_information == "-w-el":
            drop = [identifier, " H-H_SUM", " H-H_MAX", " VDW_ATT", " VDW_REP"]
        elif self.add_information == "basic-el-vdw":
            drop = [identifier, " HW-HW_SUM", " HW-HW_MAX"]
        elif self.add_information == "-w-el-vdw":
            drop = [identifier, " H-H_SUM", " H-H_MAX"]
        elif self.add_information == "basic-vdw":
            drop = [
                identifier,
                " HW-HW_SUM",
                " HW-HW_MAX",
                " ES",
                " VDW_ATT",
                " VDW_REP",
            ]
        elif self.add_information == "-w-vdw":
            drop = [
                identifier,
                " HW-HW_SUM",
                " HW-HW_MAX",
                " ES",
                " VDW_ATT",
                " VDW_REP",
            ]
        else:
            drop = [identifier]
            
        self.features_test = self.features_test.drop(drop, axis=1)
        self.features_test = self.features_test.to_numpy()
        
        if "/" in self.scoretype:
            self.scoretype = self.scoretype.replace("/", "div")

        reg = joblib.load(
            f"{loadpath}{self.modeltype}_{self.scoretype}_{self.datatype}_{self.add_information}.sav"
        )

        if "div" in self.scoretype:
            self.scoretype = self.scoretype.replace("div", "/")
            
        scores_pre = reg.predict(self.features_test)
        self.scores_pre = scores_pre

        (PDBs, scores) = PDB_codes, self.scores_pre
        return (PDBs, scores)

    def get_stats(self,spearman=False):
        if spearman == False:
            return (
                self.modeltype,
                self.scoretype,
                self.datatype,
                self.mae,
                self.mse,
                self.sd,
                self.r,
                self.conf_int,
                self.r_2,
                self.add_information,
            )
        else:
            return (
                self.modeltype,
                self.scoretype,
                self.datatype,
                self.mae,
                self.mse,
                self.sd,
                self.r,
                self.conf_int,
                self.r_2,
                self.spearman_r,
                self.add_information,
            )
