import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Printer:
    # Paths
    path_to_solutions = os.path.join("data", "solutions")
    path_to_calibration = os.path.join("data", "calibration", "parameters")
    path_to_trial = os.path.join("data", "comparisons", "trial")
    path_to_size = os.path.join("data", "comparisons", "size")
    path_to_polygons = os.path.join("data", "polygons")
    path_to_iterations = os.path.join("data", "iterations")
    path_to_results = os.path.join("data", "comparisons", "results")
    plots_folder = "plots"

    # Colors
    genetic_color = "#4f6cf1"
    genetic_accent = "#B0BCF2"
    tabu_color = "#F24F4F"
    tabu_accent = "#F2B0B0"
    cplex_color = "#3BB273"
    cplex_accent = "#82B298"

    @staticmethod
    def __split_solutions(data, num_trials, mean=True):
        objs, times = zip(*[(float(v.split(":")[0]), float(v.split(":")[1]))
                            for s in [d.split(" ") for d in data] for v in s if v != ""])
        if mean:
            objs = np.mean(np.reshape(objs, (-1, num_trials)), axis=1)
            times = np.mean(np.reshape(times, (-1, num_trials)), axis=1)
        return objs, times

    def print_times_and_obj(self):
        # Set the list of solutions to be plotted
        # solutions = ["random_polygons.csv"]
        solutions = os.listdir(self.path_to_solutions)

        for i, file_name in enumerate(solutions):
            base_name = str(file_name.split(".")[0])
            solver_name, sizes, target_time, num_trials, problem_name = base_name.split("__")
            print("Printing plots for problem instance of type {p} for target time {t}...\n".format(p=problem_name,
                                                                                                    t=target_time))
            data = pd.read_csv(os.path.join(self.path_to_solutions, file_name))
            objs, times = self.__split_solutions(data["solutions"], int(num_trials.split("_")[0]))
            title = "{n} trials on {p} for target time {t}".format(n=num_trials.split("_")[0],
                                                                   p=problem_name,
                                                                   t=target_time)

            fig, ax = plt.subplots()
            ax.plot(data["size"], times)
            ax.set(xlabel='Size', ylabel='Time (s)',
                   title='Time over increasing size\n' + title)
            os.makedirs(self.plots_folder, exist_ok=True)
            fig.savefig(os.path.join(self.plots_folder, "solutions", base_name + "__time.png"))
            plt.clf()

            fig, ax = plt.subplots()
            ax.plot(data["size"], objs)
            ax.set(xlabel='Size', ylabel='Obj value',
                   title='Obj values over increasing size\n' + title)
            os.makedirs(self.plots_folder, exist_ok=True)
            fig.savefig(os.path.join(self.plots_folder, "solutions", base_name + "__objs.png"))
            plt.clf()

    def print_polygons(self):
        solutions = os.listdir(self.path_to_polygons)
        for file_name in solutions:
            data = pd.read_csv(os.path.join(self.path_to_polygons, file_name))

            print("Plotting polygons for {} holes...\n".format(data["num_vertices"].sum()))

            vertices = [i for ls in [vs.split(" ") for vs in data["vertices"]] for i in ls if i != ""]
            x, y = zip(*[(float(v.split(":")[0]), float(v.split(":")[1])) for v in vertices])

            board_size = 200

            fig, ax = plt.subplots()
            ax.plot(x, y, 'ro')
            ax.plot(range(0, board_size + 1), [board_size] * (board_size + 1))
            ax.plot([board_size] * (board_size + 1), range(0, board_size + 1))
            ax.set(title='Grid of random polygons')
            ax.grid()

            os.makedirs(self.plots_folder, exist_ok=True)
            fig.savefig(os.path.join(self.plots_folder, "polygons", str(file_name.split(".")[0]) + ".png"))

    def print_iterations(self):
        iterations = os.listdir(self.path_to_iterations)
        for file_name in iterations:
            data = pd.read_csv(os.path.join(self.path_to_iterations, file_name))
            base_name = str(file_name.split(".")[0])
            print("Plotting iterations for {}...\n".format(base_name))
            y = [float(v) for v in data["obj_value"]]

            fig, ax = plt.subplots()
            ax.plot(range(0, len(y)), y)
            ax.set(title=str(len(data)) + ' iterations for ' + base_name)
            ax.grid()

            os.makedirs(self.plots_folder, exist_ok=True)
            fig.savefig(os.path.join(self.plots_folder, "iterations", base_name + ".png"))

    def print_calibration(self):
        calibrations = os.listdir(self.path_to_calibration)
        exact_solution = 629.8
        for dirname in calibrations:
            df = pd.DataFrame(columns=["value", "mte", "mov", "me", "score"])
            print("Printing calibration of parameter: " + dirname)
            for filename in sorted(os.listdir(os.path.join(self.path_to_calibration, dirname))):
                print(str(filename))
                data = pd.read_csv(os.path.join(self.path_to_calibration, dirname, filename))
                data = [d.split(":") for d in data.loc[0, "solutions"].split(" ") if d != ""]
                ov, te = zip(*[(float(d[0]), float(d[1])) for d in data])
                mov = np.array(ov).mean()
                mte = np.array(te).mean()
                me = (mov - exact_solution) * 100 / exact_solution
                df = df.append({
                    "value": int(str(filename).split("_")[1].split(".")[0]),
                    "mte": mte,
                    "mov": mov,
                    "me": me,
                    "score": 1 / (mte * me + 1)
                }, ignore_index=True)

            df.to_csv(os.path.join(self.path_to_calibration, dirname, dirname + ".csv"))

    def __print_size_comparison(self, data, base_name):
        base_size = int(base_name.split("__")[1].split("_")[0])
        target_size = int(base_name.split("__")[1].split("_")[1])
        instance = str(base_name.split("__")[2])
        num_trials = int(base_name.split("__")[0].split("_")[0])

        gm_objs, gm_times = self.__split_solutions(data["genetic"], num_trials)
        tm_objs, tm_times = self.__split_solutions(data["tabu"], num_trials)
        cm_objs, cm_times = self.__split_solutions(data["cplex"], num_trials)
        cmm_objs = [cm_objs.mean()] * len(data["size"])
        cmm_times = [cm_times.mean()] * len(data["size"])
        tmm_objs = [tm_objs.mean()] * len(data["size"])
        tmm_times = [tm_times.mean()] * len(data["size"])
        gmm_objs = [gm_objs.mean()] * len(data["size"])
        gmm_times = [gm_times.mean()] * len(data["size"])

        print("Plotting comparison for {n} comparisons {d}, size {b} to {t}...\n".format(n=num_trials,
                                                                                         d=instance,
                                                                                         b=base_size,
                                                                                         t=target_size))
        # OBJ FUNCTION

        # Tabu
        fig, ax = plt.subplots()
        ax.plot(data["size"], cm_objs, self.cplex_accent, linestyle="dashed", label="CPLEX")
        ax.plot(data["size"], cmm_objs, self.cplex_color, label="Mean of CPLEX")
        ax.plot(data["size"], tm_objs, self.tabu_accent, linestyle="dashed", label="Tabu Search")
        ax.plot(data["size"], tmm_objs, self.tabu_color, label="Mean of Tabu")
        ax.set(xlabel='Size', ylabel='Objective function',
               title='Tabu Search - Obj value \n {n} comparisons on {d} - Size {b} to {t}'.format(n=num_trials,
                                                                                                  d=instance,
                                                                                                  b=base_size,
                                                                                                  t=target_size))
        ax.legend()
        os.makedirs(self.plots_folder, exist_ok=True)
        fig.savefig(os.path.join(self.plots_folder, "comparisons", "size", "tabu_" + base_name + "_size_objs.png"))
        plt.clf()

        # Genetic
        fig, ax = plt.subplots()
        ax.plot(data["size"], cm_objs, self.cplex_accent, linestyle="dashed", label="CPLEX")
        ax.plot(data["size"], cmm_objs, self.cplex_color, label="Mean of CPLEX")
        ax.plot(data["size"], gm_objs, self.genetic_accent, linestyle="dashed", label="Genetic Algorithm")
        ax.plot(data["size"], gmm_objs, self.genetic_color, label="Mean of Genetic")
        ax.set(xlabel='Size', ylabel='Objective function',
               title='Genetic Algorithm - Obj value \n {n} comparisons on {d} - Size {b} to {t}'.format(
                   n=num_trials,
                   d=instance,
                   b=base_size,
                   t=target_size))
        ax.legend()
        os.makedirs(self.plots_folder, exist_ok=True)
        fig.savefig(os.path.join(self.plots_folder, "comparisons", "size", "genetic_" + base_name + "_size_objs.png"))
        plt.clf()

        # Comparison
        fig, ax = plt.subplots()
        ax.plot(data["size"], cm_objs, self.cplex_accent, linestyle="dashed", label="CPLEX")
        ax.plot(data["size"], cmm_objs, self.cplex_color, label="Mean of CPLEX")
        ax.plot(data["size"], tm_objs, self.tabu_accent, linestyle="dashed", label="Tabu Search")
        ax.plot(data["size"], tmm_objs, self.tabu_color, label="Mean of Tabu")
        ax.plot(data["size"], gm_objs, self.genetic_accent, linestyle="dashed", label="Genetic Algorithm")
        ax.plot(data["size"], gmm_objs, self.genetic_color, label="Mean of Genetic")
        ax.set(xlabel='Size', ylabel='Objective function',
               title='Comparison - Obj value \n {n} comparisons on {d} - Size {b} to {t}'.format(n=num_trials,
                                                                                                 d=instance,
                                                                                                 b=base_size,
                                                                                                 t=target_size))
        ax.legend()
        os.makedirs(self.plots_folder, exist_ok=True)
        fig.savefig(os.path.join(self.plots_folder, "comparisons", "size", base_name + "_size_objs.png"))
        plt.clf()

        # TIMES
        fig, ax = plt.subplots()
        ax.plot(data["size"], cm_times, self.cplex_accent, linestyle='dashed', label="CPLEX")
        ax.plot(data["size"], cmm_times, self.cplex_color, label="Mean of CPLEX")
        ax.plot(data["size"], tm_times, self.tabu_accent, linestyle='dashed', label="Tabu Search")
        ax.plot(data["size"], tmm_times, self.tabu_color, label="Mean of Tabu")
        ax.plot(data["size"], gm_times, self.genetic_accent, linestyle='dashed', label="Genetic Algorithm")
        ax.plot(data["size"], gmm_times, self.genetic_color, label="Mean of Genetic")
        ax.set(xlabel='Size', ylabel='Time of execution',
               title='Time of execution \n {n} comparisons on {d} - Size {b} to {t}'.format(n=num_trials,
                                                                                            d=instance,
                                                                                            b=base_size,
                                                                                            t=target_size))
        ax.legend()
        ax.grid()
        os.makedirs(self.plots_folder, exist_ok=True)
        fig.savefig(os.path.join(self.plots_folder, "comparisons", "size", base_name + "_size_times.png"))

    def __analyze_size_comparison(self):
        comparisons = os.listdir(self.path_to_size)
        for file_name in comparisons:
            data = pd.read_csv(os.path.join(self.path_to_size, file_name))
            base_name = str(file_name.split(".")[0])
            num_trials = int(base_name.split("__")[0].split("_")[0])

            c_objs, c_times = self.__split_solutions(data["cplex"], num_trials, mean=False)
            c_analysis = pd.DataFrame({
                "size": data["size"],
                "mean_obj": np.mean(np.reshape(c_objs, (-1, num_trials)), axis=1),
                "std_obj": np.std(np.reshape(c_objs, (-1, num_trials)), axis=1),
                "min_obj": np.min(np.reshape(c_objs, (-1, num_trials)), axis=1),
                "max_obj": np.max(np.reshape(c_objs, (-1, num_trials)), axis=1),
                "mean_time": np.mean(np.reshape(c_times, (-1, num_trials)), axis=1),
                "std_time": np.std(np.reshape(c_times, (-1, num_trials)), axis=1),
                "min_time": np.min(np.reshape(c_times, (-1, num_trials)), axis=1),
                "max_time": np.max(np.reshape(c_times, (-1, num_trials)), axis=1)
            })
            c_analysis.to_csv(os.path.join(self.path_to_results,
                                           "size_cplex_on_{}_analysis.csv".format(str(base_name.split("__")[2]))),
                              index=None)

            t_objs, t_times = self.__split_solutions(data["tabu"], num_trials, mean=False)
            t_analysis = pd.DataFrame({
                "size": data["size"],
                "mean_obj": np.mean(np.reshape(t_objs, (-1, num_trials)), axis=1),
                "std_obj": np.std(np.reshape(t_objs, (-1, num_trials)), axis=1),
                "min_obj": np.min(np.reshape(t_objs, (-1, num_trials)), axis=1),
                "max_obj": np.max(np.reshape(t_objs, (-1, num_trials)), axis=1),
                "mean_time": np.mean(np.reshape(t_times, (-1, num_trials)), axis=1),
                "std_time": np.std(np.reshape(t_times, (-1, num_trials)), axis=1),
                "min_time": np.min(np.reshape(t_times, (-1, num_trials)), axis=1),
                "max_time": np.max(np.reshape(t_times, (-1, num_trials)), axis=1)
            })
            t_analysis.to_csv(os.path.join(self.path_to_results,
                                           "size_tabu_on_{}_analysis.csv".format(str(base_name.split("__")[2]))),
                              index=None)

            g_objs, g_times = self.__split_solutions(data["genetic"], num_trials, mean=False)
            g_analysis = pd.DataFrame({
                "size": data["size"],
                "mean_obj": np.mean(np.reshape(g_objs, (-1, num_trials)), axis=1),
                "std_obj": np.std(np.reshape(g_objs, (-1, num_trials)), axis=1),
                "min_obj": np.min(np.reshape(g_objs, (-1, num_trials)), axis=1),
                "max_obj": np.max(np.reshape(g_objs, (-1, num_trials)), axis=1),
                "mean_time": np.mean(np.reshape(g_times, (-1, num_trials)), axis=1),
                "std_time": np.std(np.reshape(g_times, (-1, num_trials)), axis=1),
                "min_time": np.min(np.reshape(g_times, (-1, num_trials)), axis=1),
                "max_time": np.max(np.reshape(g_times, (-1, num_trials)), axis=1)
            })
            g_analysis.to_csv(os.path.join(self.path_to_results,
                                           "size_genetic_on_{}_analysis.csv".format(str(base_name.split("__")[2]))),
                              index=None)
            self.__print_size_comparison(data, base_name)

    def __print_trial_comparison(self, data, base_name):
        num_trials = int(base_name.split("__")[0].split("_")[0])
        instance = str(base_name.split("__")[1])
        r = range(0, num_trials)

        gm_obj = [data["genetic_obj"].mean()] * num_trials
        tm_obj = [data["tabu_obj"].mean()] * num_trials
        gm_time = [data["genetic_time"].mean()] * num_trials
        tm_time = [data["tabu_time"].mean()] * num_trials
        cm_time = [data["cplex_time"].mean()] * num_trials

        print("Plotting comparison for {n} trial on {d}...\n".format(n=num_trials, d=instance))

        # OBJ FUNCTION

        # Tabu
        fig, ax = plt.subplots()
        ax.plot(r, data["cplex_obj"], self.cplex_color, label="CPLEX")
        ax.plot(r, data["tabu_obj"], self.tabu_accent, linestyle="dashed", label="Tabu Search")
        ax.plot(r, tm_obj, self.tabu_color, label="Mean of Tabu")
        ax.set(xlabel='Trial', ylabel='Objective function',
               title='Tabu Search - Obj value - {n} trial on {d}'.format(n=num_trials, d=instance))
        ax.legend()
        os.makedirs(self.plots_folder, exist_ok=True)
        fig.savefig(os.path.join(self.plots_folder, "comparisons", "trial", "tabu_" + base_name + "_trials_objs.png"))
        plt.clf()

        # Genetic
        fig, ax = plt.subplots()
        ax.plot(r, data["cplex_obj"], self.cplex_color, label="CPLEX")
        ax.plot(r, data["genetic_obj"], self.genetic_accent, linestyle="dashed", label="Genetic Algorithm")
        ax.plot(r, gm_obj, self.genetic_color, label="Mean of Genetic")
        ax.set(xlabel='Trial', ylabel='Objective function',
               title='Genetic Algorithm - Obj value - {n} trial on {d}'.format(n=num_trials, d=instance))
        ax.legend()
        os.makedirs(self.plots_folder, exist_ok=True)
        fig.savefig(
            os.path.join(self.plots_folder, "comparisons", "trial", "genetic_" + base_name + "_trials_objs.png"))
        plt.clf()

        # Comparison
        fig, ax = plt.subplots()
        ax.plot(r, data["cplex_obj"], self.cplex_color, label="CPLEX")
        ax.plot(r, data["tabu_obj"], self.tabu_accent, linestyle="dashed", label="Tabu Search")
        ax.plot(r, tm_obj, self.tabu_color, label="Mean of Tabu")
        ax.plot(r, data["genetic_obj"], self.genetic_accent, linestyle="dashed", label="Genetic Algorithm")
        ax.plot(r, gm_obj, self.genetic_color, label="Mean of Genetic")
        ax.set(xlabel='Trial', ylabel='Objective function',
               title='Comparison on obj value - {n} trial on {d}'.format(n=num_trials, d=instance))
        ax.legend()
        os.makedirs(self.plots_folder, exist_ok=True)
        fig.savefig(os.path.join(self.plots_folder, "comparisons", "trial", base_name + "_trial_objs.png"))
        plt.clf()

        # TIMES
        fig, ax = plt.subplots()
        ax.plot(r, data["cplex_time"], self.cplex_accent, linestyle='dashed', label="CPLEX")
        ax.plot(r, cm_time, self.cplex_color, label="Mean of CPLEX")
        ax.plot(r, data["tabu_time"], self.tabu_accent, linestyle='dashed', label="Tabu Search")
        ax.plot(r, tm_time, self.tabu_color, label="Mean of Tabu")
        ax.plot(r, data["genetic_time"], self.genetic_accent, linestyle='dashed', label="Genetic Algorithm")
        ax.plot(r, gm_time, self.genetic_color, label="Mean of Genetic")
        ax.set(xlabel='Trial', ylabel='Time of execution',
               title='Comparison on time of execution - {n} trial on {d}'.format(n=num_trials, d=instance))
        ax.legend()
        ax.grid()
        os.makedirs(self.plots_folder, exist_ok=True)
        fig.savefig(os.path.join(self.plots_folder, "comparisons", "trial", base_name + "_trial_times.png"))

    def __analyze_trial_comparison(self):
        comparisons = os.listdir(self.path_to_trial)
        for file_name in comparisons:
            data = pd.read_csv(os.path.join(self.path_to_trial, file_name), usecols=["tabu_obj",
                                                                                     "tabu_time",
                                                                                     "genetic_obj",
                                                                                     "genetic_time",
                                                                                     "cplex_obj",
                                                                                     "cplex_time"])
            base_name = str(file_name.split(".")[0])
            analysis = pd.DataFrame({
                "method_metrics": ["cplex_obj", "cplex_time",
                                   "tabu_obj", "tabu_time",
                                   "genetic_obj", "genetic_time"],
                "mean": [data["cplex_obj"].mean(), data["cplex_time"].mean(),
                         data["tabu_obj"].mean(), data["tabu_time"].mean(),
                         data["genetic_obj"].mean(), data["genetic_time"].mean()],
                "std": [data["cplex_obj"].std(), data["cplex_time"].std(),
                        data["tabu_obj"].std(), data["tabu_time"].std(),
                        data["genetic_obj"].std(), data["genetic_time"].std()],
                "min": [data["cplex_obj"].min(), data["cplex_time"].min(),
                        data["tabu_obj"].min(), data["tabu_time"].min(),
                        data["genetic_obj"].min(), data["genetic_time"].min()],
                "max": [data["cplex_obj"].max(), data["cplex_time"].max(),
                        data["tabu_obj"].max(), data["tabu_time"].max(),
                        data["genetic_obj"].max(), data["genetic_time"].max()]
            })
            analysis.to_csv(os.path.join(self.path_to_results,
                                         "trial_on_{}_analysis.csv".format(str(base_name.split("__")[1]))),
                            index=None)
            self.__print_trial_comparison(data, base_name)

    def print_comparison(self):
        # self.__analyze_trial_comparison()
        self.__analyze_size_comparison()
