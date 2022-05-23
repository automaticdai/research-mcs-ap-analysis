########################################################################################################################
# Analysis of Mixed-Criticality System Approximated Computing (MCS-AP)
# Dr. Xiaotian Dai
# University of York
# 2021-2022
########################################################################################################################

import math
import random
import pprint
import logging
import sys
from util import UUniFast

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)  # logging.debug() / logging.info()
random.seed(10)  # fix the random seed for reproducibility

########################################################################################################################
# Parameters
########################################################################################################################
number_of_tasks = 10     # total number of tasks
number_of_hi_tasks = 5  # control proportion of hi-tasks to be ~50%
periods = [100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000]  # hyper-period = 10000

C_low_Mi = 0.8
C_high_multiplier = 3

number_of_trials = 50

# Mi = 0, 1, 2, 3, 4, 5, ..., 23
QoC_candidates = [1.0, 0.95, 0.9, 0.85, 0.8]
C_low_candidates = [1.0, 0.9, 0.8, 0.7, 0.6]

########################################################################################################################
# Variables
########################################################################################################################
# task = type {0:LO, 1:HI}, Ci(LO), Ci(MI), Ci(HI), Ti
# with Di = Ti
# the taskset is indexed and ordered by priority

# taskset = []
# taskidx = []
# taskidx_high = []
# taskidx_low = []


########################################################################################################################
# Functions
########################################################################################################################
# check low RTA
def rta_low(i, taskset):
    Ci = taskset[i]["C_low"]
    Ti = taskset[i]["T_i"]
    Di = Ti
    typei = taskset[i]["Type"]

    Ri = -100000
    Ri_new = Ci

    while abs(Ri_new - Ri) > 1:
        Ri = Ri_new
        Ij = 0
        for j in range(len(taskset)):
            typej = taskset[j]["Type"]
            if j < i:
                Cj = taskset[j]["C_low"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri / Tj) * Cj
        Ri_new = Ci + Ij

        if Ri_new > Di:
            Ri = Ri_new
            break

    return Ri


# check low -> mid
def rta_low_to_mid(i, taskset):
    Ri_LO = rta_low(i, taskset)

    Ci = taskset[i]["C_low"]
    Ci_mid = taskset[i]["C_mid"]
    Ci_high = taskset[i]["C_high"]

    Ti = taskset[i]["T_i"]
    Di = Ti
    typei = taskset[i]["Type"]

    Ri = -100000
    Ri_new = Ci_mid

    while (abs(Ri_new - Ri) > 1):
        Ri = Ri_new
        Ij = 0
        for j in range(len(taskset)):
            typej = taskset[j]["Type"]
            if j < i and typej == 0:  # j = hpL(i)
                Cj_mid = taskset[j]["C_mid"]
                Cj_low = taskset[j]["C_low"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri_LO / Tj) * Cj_low + (math.ceil(Ri / Tj) - math.ceil(Ri_LO / Tj)) * Cj_mid
            elif j < i and typej == 1:  # j = hpH(i)
                Cj_mid = taskset[j]["C_mid"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri / Tj) * Cj_mid
        Ri_new = Ci_mid + Ij

        if Ri_new > Di:
            Ri = Ri_new
            break

    return Ri


# check mid RTA
def rta_mid(i, taskset):
    Ci = taskset[i]["C_low"]
    Ci_mid = taskset[i]["C_mid"]
    Ti = taskset[i]["T_i"]
    Di = Ti
    typei = taskset[i]["Type"]

    Ri = -100000
    Ri_new = Ci_mid

    while (abs(Ri_new - Ri) > 1):
        Ri = Ri_new
        Ij = 0
        for j in range(len(taskset)):
            typej = taskset[j]["Type"]
            if j < i:
                if typej == 0:  # j = hpL(i)
                    Cj_mid = taskset[j]["C_mid"]
                    Tj = taskset[j]["T_i"]
                    Ij += math.ceil(Ri / Tj) * Cj_mid
                elif typej == 1:  # j = hpH(i)
                    Cj_mid = taskset[j]["C_mid"]
                    Tj = taskset[j]["T_i"]
                    Ij += math.ceil(Ri / Tj) * Cj_mid
        Ri_new = Ci_mid + Ij

        if Ri_new > Di:
            Ri = Ri_new
            break

    return Ri


# check mid -> high
def rta_mid_to_high(i, taskset):
    Ri_MI = rta_mid(i, taskset)
    Ri_LO = rta_low(i, taskset)

    Ci = taskset[i]["C_low"]
    Ci_mid = taskset[i]["C_mid"]
    Ci_high = taskset[i]["C_high"]

    Ti = taskset[i]["T_i"]
    Di = Ti
    typei = taskset[i]["Type"]

    Ri = -100000
    Ri_new = Ci_high

    while abs(Ri_new - Ri) > 1:
        Ri = Ri_new
        Ij = 0
        for j in range(len(taskset)):
            typej = taskset[j]["Type"]
            if j < i:
                if typej == 0:  # j = hpL(i)
                    Cj_low = taskset[j]["C_low"]
                    Cj_mid = taskset[j]["C_mid"]
                    Tj = taskset[j]["T_i"]
                    Ij += math.ceil(Ri_LO / Tj) * Cj_low + (math.ceil(Ri_MI / Tj) - math.ceil(Ri_LO / Tj)) * Cj_mid
                elif typej == 1:  # j = hpH(i)
                    Cj_high = taskset[j]["C_high"]
                    Tj = taskset[j]["T_i"]
                    Ij += math.ceil(Ri / Tj) * Cj_high
        Ri_new = Ci_high + Ij

        if Ri_new > Di:
            Ri = Ri_new
            break

    return Ri


# check low -> high (for dual-criticality only)
def rta_low_to_high(i, taskset):
    Ri_LO = rta_low(i, taskset)

    Ci = taskset[i]["C_low"]
    Ci_high = taskset[i]["C_high"]

    Ti = taskset[i]["T_i"]
    Di = Ti
    typei = taskset[i]["Type"]

    Ri = -100000
    Ri_new = Ci_high

    while abs(Ri_new - Ri) > 1:
        Ri = Ri_new
        Ij = 0
        for j in range(len(taskset)):
            if j < i:
                typej = taskset[j]["Type"]
                if typej == 0:  # j = hpL(i)
                    Cj_low = taskset[j]["C_low"]
                    Tj = taskset[j]["T_i"]
                    Ij += math.ceil(Ri_LO / Tj) * Cj_low
                elif typej == 1:  # j = hpH(i)
                    Cj_high = taskset[j]["C_high"]
                    Tj = taskset[j]["T_i"]
                    Ij += math.ceil(Ri / Tj) * Cj_high
        Ri_new = Ci_high + Ij

        if Ri_new > Di:
            Ri = Ri_new
            break

    return Ri


# check high RTA
def rta_high(i, taskset):
    Ci = taskset[i]["C_low"]
    Ci_high = taskset[i]["C_high"]
    Ti = taskset[i]["T_i"]
    Di = Ti
    typei = taskset[i]["Type"]

    Ri = -100000
    Ri_new = Ci_high

    while abs(Ri_new - Ri) > 1:
        Ri = Ri_new
        Ij = 0
        for j in range(len(taskset)):
            typej = taskset[j]["Type"]
            if j < i and typej == 1:  # j = hpH(i)
                Cj_high = taskset[j]["C_high"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri / Tj) * Cj_high
        Ri_new = Ci_high + Ij

        if Ri_new > Di:
            Ri = Ri_new
            break

    return Ri


########################################################################################################################
# Trials
########################################################################################################################
def task_gen(target_util):
    taskset = []
    taskidx = []
    taskidx_high = []
    taskidx_low = []

    # generate the taskset
    U = UUniFast(number_of_tasks, target_util)
    # print(U)

    sum_U = 0

    for U_i in U:
        T_i = random.choice(periods)
        C_i = math.floor(U_i * T_i)
        if C_i < 1:
            C_i = 1
        sum_U += C_i / T_i
        task_i = {"Type": 0, "C_low": C_i, "C_mid": 0, "C_high": 0, "T_i": T_i}

        # append through deadline monotonic
        found = False
        for j in range(len(taskset)):
            if T_i <= taskset[j]["T_i"]:
                taskset.insert(j, task_i)
                found = True
                break
        if not taskset or not found:
            taskset.append(task_i)

    # print("Total util: ", sum_U)

    # select the HI tasks
    for i in range(number_of_tasks):
        taskidx.append(i)
    idx_high_ = random.sample(taskidx, number_of_hi_tasks)

    for i in range(number_of_tasks):
        if i in idx_high_:
            taskset[i]["Type"] = 1
            taskset[i]["C_high"] = math.ceil(taskset[i]["C_low"] * C_high_multiplier)
            taskidx_high.append(i)
        else:
            taskset[i]["C_mid"] = math.ceil(taskset[i]["C_low"] * C_low_Mi)
            taskidx_low.append(i)

    # check schedulability of low mode; if not then discard (discard skipped)
    for i in range(number_of_tasks):
        Di = taskset[i]["T_i"]
        Ri = rta_low(i, taskset)
        if (Ri > Di):
            # print("Error #1.")
            pass

    return taskset, taskidx, taskidx_high, taskidx_low


def trial_searching(taskset_, taskidx, taskidx_high, taskidx_low, search_gamma=True, gamma_in=-1):
    taskset = taskset_.copy()

    # --------------------------------------------------------------------------------------------
    # searching
    # --------------------------------------------------------------------------------------------
    if search_gamma:
        # generate the gamma set
        # gamma_set = []
        #
        # step = 0.01
        # gamma_this = 0
        # while gamma_this < 1:
        #     gamma_set.append(gamma_this)
        #     gamma_this += step
        gamma_set = [x * 0.01 for x in range(1, 100)]

        # calculate C_i^MD based on gamma
        L = 0
        R = len(gamma_set) - 1

        while L <= (R - 2):
            M = math.ceil((R + L) / 2)
            gamma = gamma_set[M]

            # calculate C_i^MID
            for i in taskidx_high:
                taskset[i]["C_mid"] = taskset[i]["C_low"] + math.ceil((taskset[i]["C_high"] - taskset[i]["C_low"]) * gamma)

            # test scheduablity
            scheduable = True

            # check sched(MID) for all tasks
            for i in taskidx:
                Di = taskset[i]["T_i"]
                typei = taskset[i]["Type"]
                if typei == 0:      # for low task, only stable mode sched(MID) needs to be guaranteed
                    Ri = rta_mid(i, taskset)
                else:               # for high task, that are two transiation modes
                    Ri_l2m = rta_low_to_mid(i, taskset)
                    Ri_m2h = rta_mid_to_high(i, taskset)
                    Ri = max(Ri_l2m, Ri_m2h)
                if (Ri >= Di):
                    scheduable = False
                elif (Ri <= Di and typei == 0):
                    pass

            if scheduable:
                L = M
            else:
                R = M
    else:
        gamma = gamma_in

        # calculate C_i^MID
        for i in taskidx_high:
            taskset[i]["C_mid"] = taskset[i]["C_low"] + math.ceil((taskset[i]["C_high"] - taskset[i]["C_low"]) * gamma)
            # check sched(MID) for all tasks
            for i in taskidx:
                Di = taskset[i]["T_i"]
                typei = taskset[i]["Type"]
                if typei == 0:      # for low task, only stable mode sched(MID) needs to be guaranteed
                    Ri = rta_mid(i, taskset)
                else:               # for high task, that are two transiation modes
                    Ri_l2m = rta_low_to_mid(i, taskset)
                    Ri_m2h = rta_mid_to_high(i, taskset)
                    Ri = max(Ri_l2m, Ri_m2h)
                if (Ri >= Di):
                    return 0, 0, 0, 0, gamma  #  scheduable = False
                elif (Ri <= Di and typei == 0):
                    pass

    # --------------------------------------------------------------------------------------------
    # pprint.pprint(taskset)

    # select and overrun one of the HI tasks
    overrun_task_idx = random.choice(taskidx_high)
    # print("Overrun idx: ", overrun_task_idx)

    # check schedulability of HI or AP
    overrun_Ci_low = taskset[overrun_task_idx]["C_low"]
    overrun_Ci_mid = taskset[overrun_task_idx]["C_mid"]
    overrun_Ci_high = taskset[overrun_task_idx]["C_high"]

    overrun_Ci = int(random.randrange(overrun_Ci_low, overrun_Ci_high))
    # print("overrun_Ci: ", overrun_Ci)

    # METHOD 1
    method1_feasible = True
    # check sched(HI) for low-to-high
    for i in taskidx_high:
        Di = taskset[i]["T_i"]
        Ri = rta_low_to_high(i, taskset)
        if (Ri > Di):
            method1_feasible = False
            # print("Error #2.")

    # METHOD 2
    method2_feasible = True
    method2_feasible_high = True
    low_all = len(taskidx_low)
    low_count = 0
    if overrun_Ci > overrun_Ci_mid:
        # check sched(HI) for high task only
        for i in taskidx_high:
            Di = taskset[i]["T_i"]
            Ri = rta_mid_to_high(i, taskset)
            if (Ri > Di):
                method2_feasible = False
                # print("Error #3.")
    else:  # overrun_Ci > taskset[overrun_task_idx]["C_low"]:
        # check sched(MID) for all tasks
        for i in taskidx:
            Di = taskset[i]["T_i"]
            typei = taskset[i]["Type"]

            if typei == 0:      # for low task, only stable mode sched(MID) needs to be guaranteed
                Ri = rta_mid(i, taskset)
            else:
                Ri = rta_low_to_mid(i, taskset)
            if (Ri > Di):
                method2_feasible = False
                if (typei == 1):
                    method2_feasible_high = False
                # print("Error #4.")
            elif (Ri <= Di and typei == 0):
                low_count += 1

    # if any high task is infeasible, then the result is invalid
    if not method2_feasible_high:
        low_count = 0

    return method1_feasible, method2_feasible, low_count, low_all, gamma


def taskset_gen(target_util):
    pass


########################################################################################################################
# Main starts here
########################################################################################################################
if __name__ == "__main__":
    for util_this in [0.10, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]:
        target_util = util_this
        count1 = 0
        count2 = 0
        count_low = 0
        count_low_1 = 0
        count_low_2 = 0
        count_low_3 = 0
        count_low_4 = 0
        count_low_all = 0

        for k in range(number_of_trials):
            # generate the taskset
            taskset, taskidx, taskidx_high, taskidx_low = task_gen(util_this)

            # Optimal
            method1_feasible, method2_feasible, low_count, low_all, gamma_this = trial_searching(taskset, taskidx, taskidx_high, taskidx_low)
            #print(round(gamma_this, 2))

            if method1_feasible:
                count1 += 1
            if method2_feasible:
                count2 += 1
            count_low += low_count
            count_low_all += low_all

            # compare fixed gammas
            method1_feasible, method2_feasible, low_count, low_all, _ = trial_searching(taskset, taskidx, taskidx_high, taskidx_low, search_gamma=False, gamma_in=0.1)
            if method2_feasible:
                count_low_1 += low_count

            method1_feasible, method2_feasible, low_count, low_all, _ = trial_searching(taskset, taskidx, taskidx_high, taskidx_low, search_gamma=False, gamma_in=0.3)
            if method2_feasible:
                count_low_2 += low_count

            method1_feasible, method2_feasible, low_count, low_all, _ = trial_searching(taskset, taskidx, taskidx_high, taskidx_low, search_gamma=False, gamma_in=0.5)
            if method2_feasible:
                count_low_3 += low_count

            method1_feasible, method2_feasible, low_count, low_all, _ = trial_searching(taskset, taskidx, taskidx_high, taskidx_low, search_gamma=False, gamma_in=0.7)
            if method2_feasible:
                count_low_4 += low_count

        print("{0},{1},{2},{3},{4},{5},{6},{7},{8}".format(util_this, count1, count2, count_low_1, count_low_2, count_low_3, count_low_4, count_low, count_low_all))
