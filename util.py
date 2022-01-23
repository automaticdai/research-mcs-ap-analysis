########################################################################################################################
# Analysis of Mixed-Criticality System Approximated Computing (MCS-AP)
# Dr. Xiaotian Dai
# University of York
# 2021-2022
########################################################################################################################

import math
import random

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