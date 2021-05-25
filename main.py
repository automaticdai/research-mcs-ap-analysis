########################################################################################################################
# Analysis of Mixed-Criticality Approximated Computing
# Dr. Xiaotian Dai
# University of York
########################################################################################################################

import math
import random
import pprint
import logging, sys
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
#logging.debug()
#logging.info()

########################################################################################################################
# Parameters
########################################################################################################################
number_of_tasks = 5
number_of_hi_tasks = math.ceil(number_of_tasks / 2)
periods = [250, 500, 1000, 2000, 5000, 10000]   # hyperperiod = 10000

C_low_Mi = 0.7
C_high_multiplier = 3

number_of_trials = 10000


########################################################################################################################
# Variables
########################################################################################################################
# task = type {0:LO, 1:HI}, Ci(LO), Ci(MI), Ci(HI), Ti
# Di = Ti
# the taskset is ordered by priority
taskset = []
taskidx = []
taskidx_high = []
taskidx_low = []


########################################################################################################################
# Functions
########################################################################################################################
# UUniFast
def UUniFast(n, U_avg):
    USet = []
    sumU = U_avg
    for i in range(n-1):
        nextSumU = sumU*math.pow(random.random(), 1/(n-i))
        USet.append(sumU-nextSumU)
        sumU = nextSumU
    USet.append(sumU)
    return USet



# check low RTA
def rta_low(i, taskset):
    Ci = taskset[i]["C_low"]
    Ti = taskset[i]["T_i"]
    Di = Ti
    typei = taskset[i]["Type"]

    Ri = -100000
    Ri_new = Ci

    while (abs(Ri_new - Ri) > 1):
        Ri = Ri_new
        Ij = 0
        for j in taskidx:
            typej = taskset[j]["Type"]
            if j < i:
                Cj = taskset[j]["C_low"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri/Tj) * Cj
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
    Ci_high = taskset[i]("C_high")

    Ti = taskset[i]["T_i"]
    Di = Ti
    typei = taskset[i]["Type"]

    Ri = -100000
    Ri_new = Ci_mid

    while (abs(Ri_new - Ri) > 1):
        Ri = Ri_new
        Ij = 0
        for j in taskidx:
            typej = taskset[j]["Type"]
            if j < i and typej == 0:  # j = hpL(i)
                Cj_mid = taskset[j]["C_mid"]
                Cj_low = taskset[j]["C_low"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri_LO / Tj) * Cj_low + math.ceil((Ri - Ri_LO) / Tj) * Cj_mid
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
        for j in taskidx:
            typej = taskset[j]["Type"]
            if j < i and typej == 0:    # j = hpL(i)
                Cj_mid = taskset[j]["C_mid"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri/Tj) * Cj_mid
            elif j < i and typej == 1:  # j = hpH(i)
                Cj_mid = taskset[j]["C_mid"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri/Tj) * Cj_mid
        Ri_new = Ci_mid + Ij

        if Ri_new > Di:
            Ri = Ri_new
            break

    return Ri



# check mid -> high
def rta_mid_to_high(i, taskset):
    Ri_MI = rta_low(i, taskset)

    Ci = taskset[i]["C_low"]
    Ci_mid = taskset[i]["C_mid"]
    Ci_high = taskset[i]("C_high")

    Ti = taskset[i]["T_i"]
    Di = Ti
    typei = taskset[i]["Type"]

    Ri = -100000
    Ri_new = Ci_mid

    while (abs(Ri_new - Ri) > 1):
        Ri = Ri_new
        Ij = 0
        for j in taskidx:
            typej = taskset[j]["Type"]
            if j < i and typej == 0:  # j = hpL(i)
                Cj_mid = taskset[j]["C_mid"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri_MI / Tj) * Cj_mid
            elif j < i and typej == 1:  # j = hpH(i)
                Cj_high = taskset[j]["C_high"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri / Tj) * Cj_high
        Ri_new = Ci_mid + Ij

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

    while (abs(Ri_new - Ri) > 1):
        Ri = Ri_new
        Ij = 0
        for j in taskidx:
            typej = taskset[j]["Type"]
            if j < i and typej == 1:  # j = hpH(i)
                Cj = taskset[j]["C_low"]
                Cj_mid = taskset[j]["C_mid"]
                Cj_high = taskset[j]["C_high"]
                Tj = taskset[j]["T_i"]
                Ij += math.ceil(Ri/Tj) * Cj_high
        Ri_new = Ci_high + Ij

        if Ri_new > Di:
            Ri = Ri_new
            break

    return Ri

########################################################################################################################
# Trial
########################################################################################################################

def trial(target_util, gamma):
    taskset = []
    taskidx = []
    taskidx_high = []
    taskidx_low = []

    # generate the taskset
    U = UUniFast(number_of_tasks, target_util)
    #print(U)

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

    #print("Total util: ", sum_U)

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
            #print("Error #1.")
            pass


    # calculate C_i^MD based on gamma
    for i in taskidx_high:
        taskset[i]["C_mid"] = taskset[i]["C_low"] + math.ceil((taskset[i]["C_high"] - taskset[i]["C_low"]) * gamma)

    #pprint.pprint(taskset)


    # select and overrun one of the HI tasks
    overrun_task_idx = random.choice(taskidx_high)
    #print("Overrun idx: ", overrun_task_idx)


    # check schedulability of HI or AP
    overrun_Ci_low = taskset[overrun_task_idx]["C_low"]
    overrun_Ci_mid = taskset[overrun_task_idx]["C_mid"]
    overrun_Ci_high = taskset[overrun_task_idx]["C_high"]

    overrun_Ci = int(random.randrange(overrun_Ci_low, overrun_Ci_high))
    #print("overrun_Ci: ", overrun_Ci)

    # METHOD 1
    method1_feasible = True
    # check sched(HI) for low task only
    for i in taskidx_high:
        Di = taskset[i]["T_i"]
        Ri = rta_high(i, taskset)
        if (Ri > Di):
            method1_feasible = False
            #print("Error #2.")

    # METHOD 2
    method2_feasible = True
    method2_feasible_high = True
    low_all = len(taskidx_low)
    low_count = 0
    if overrun_Ci > overrun_Ci_mid:
        # check sched(HI) for high task only
        for i in taskidx_high:
            Di = taskset[i]["T_i"]
            Ri = rta_high(i, taskset)
            if (Ri > Di):
                method2_feasible = False
                #print("Error #3.")
    else:  # overrun_Ci > taskset[overrun_task_idx]["C_low"]:
        # check sched(MID) for all tasks
        for i in taskidx:
            Di = taskset[i]["T_i"]
            typei = taskset[i]["Type"]
            Ri = rta_mid(i, taskset)
            if (Ri > Di):
                method2_feasible = False
                if (typei == 1):
                    method2_feasible_high = False
                #print("Error #4.")
            elif (Ri <= Di and typei == 0):
                low_count += 1
    # if high task is not feasible, then the result is invalid
    if not method2_feasible_high:
        low_count = 0

    #print(method1_feasible, method2_feasible)

    return (not method1_feasible and method2_feasible, low_count, low_all)



def trial_searching(target_util):
    taskset = []
    taskidx = []
    taskidx_high = []
    taskidx_low = []

    # generate the taskset
    U = UUniFast(number_of_tasks, target_util)
    #print(U)

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

    #print("Total util: ", sum_U)

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
            #print("Error #1.")
            pass


    # --------------------------------------------------------------------------------------------
    # searching
    # --------------------------------------------------------------------------------------------
    # generate the gamma set
    gamma_set = []

    step = 0.01
    gamma_this = 0

    while (gamma_this < 1):
        gamma_set.append(gamma_this)
        gamma_this += step

    # calculate C_i^MD based on gamma
    L = 0
    R = len(gamma_set) - 1

    while (L <= R-2):
        M = math.ceil((R + L) / 2)
        gamma = gamma_set[M]

        # calculate C_i^MID
        for i in taskidx_high:
            taskset[i]["C_mid"] = taskset[i]["C_low"] + math.ceil((taskset[i]["C_high"] - taskset[i]["C_low"]) * gamma)

        method2_feasible = True
        method2_feasible_high = True

        # check sched(MID) for all tasks
        for i in taskidx:
            Di = taskset[i]["T_i"]
            typei = taskset[i]["Type"]
            Ri = rta_mid(i, taskset)
            if (Ri >= Di):
                method2_feasible = False
                if (typei == 1):
                    method2_feasible_high = False
            elif (Ri <= Di and typei == 0):
                pass

        if method2_feasible_high and method2_feasible:
            L = M + 1
        else:
            R = M


    # --------------------------------------------------------------------------------------------

    #pprint.pprint(taskset)


    # select and overrun one of the HI tasks
    overrun_task_idx = random.choice(taskidx_high)
    #print("Overrun idx: ", overrun_task_idx)


    # check schedulability of HI or AP
    overrun_Ci_low = taskset[overrun_task_idx]["C_low"]
    overrun_Ci_mid = taskset[overrun_task_idx]["C_mid"]
    overrun_Ci_high = taskset[overrun_task_idx]["C_high"]

    overrun_Ci = int(random.randrange(overrun_Ci_low, overrun_Ci_high))
    #print("overrun_Ci: ", overrun_Ci)

    # METHOD 1
    method1_feasible = True
    # check sched(HI) for low task only
    for i in taskidx_high:
        Di = taskset[i]["T_i"]
        Ri = rta_high(i, taskset)
        if (Ri > Di):
            method1_feasible = False
            #print("Error #2.")

    # METHOD 2
    method2_feasible = True
    method2_feasible_high = True
    low_all = len(taskidx_low)
    low_count = 0
    if overrun_Ci > overrun_Ci_mid:
        # check sched(HI) for high task only
        for i in taskidx_high:
            Di = taskset[i]["T_i"]
            Ri = rta_high(i, taskset)
            if (Ri > Di):
                method2_feasible = False
                #print("Error #3.")
    else:  # overrun_Ci > taskset[overrun_task_idx]["C_low"]:
        # check sched(MID) for all tasks
        for i in taskidx:
            Di = taskset[i]["T_i"]
            typei = taskset[i]["Type"]
            Ri = rta_mid(i, taskset)
            if (Ri > Di):
                method2_feasible = False
                if (typei == 1):
                    method2_feasible_high = False
                #print("Error #4.")
            elif (Ri <= Di and typei == 0):
                low_count += 1
    # if high task is not feasible, then the result is invalid
    if not method2_feasible_high:
        low_count = 0

    #print(method1_feasible, method2_feasible)


    return (method1_feasible, method2_feasible, low_count, low_all, gamma)



########################################################################################################################
# Main starts here
########################################################################################################################
if __name__ == "__main__":
    # Exp 1.
    # for util_this in [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]:
    #     for gamma_this in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    #         target_util = util_this
    #         count = 0
    #         count_low = 0
    #         count_low_all = 0
    #
    #         for k in range(number_of_trials):
    #             ret, low_count, low_all = trial(util_this, gamma_this)
    #             if ret:
    #                 count += 1
    #             count_low += low_count
    #             count_low_all += low_all
    #         print(util_this, gamma_this, count, count_low, count_low_all)

    # Exp 2.
    for util_this in [0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950]:
        target_util = util_this
        count1 = 0
        count2 = 0
        count_low = 0
        count_low_all = 0

        for k in range(number_of_trials):
            ret_method1, ret_method2, low_count, low_all, gamma_this = trial_searching(util_this)
            #print(round(gamma_this, 2))
            if ret_method1:
                count1 += 1
            if ret_method2:
                count2 += 1
            count_low += low_count
            count_low_all += low_all
        print(util_this, count1, count2, count_low, count_low_all)
