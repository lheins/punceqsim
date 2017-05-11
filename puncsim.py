#!/usr/bin/env python

# Simulation of evolution of a multi-allele phenotype
# in a population due to genetic drift and selection
# Allows to test how different parameterizations affect
# the recovery of a punctuated equilibrium pattern

from math import sin
from random import random, gauss, choice, uniform

#Variables

# Determines how much fitness matters for sampling into the next generation
selectionStrength = 0.5

# Starting number of individuals in population
popsize = 10000

# Mutation rate for all alleles in population
mutRate = 0.01

# Standard deviation of the mutation value
mutStrength = 0.3

# Number of generations to run the simulation
numGens = 1000

# Determines whether population size periodically changes
changePop = True

# How large a population to expand to
bigPop = 1000

# How small a population to shrink to
smallPop = 30

# Number of generations at larger size
bigGens = 1000

# Number of generations at smaller size
smallGens = 100

# TODO: For better optimization, NumPy random functions would be more appropriate

# Core Functions
def computeFitness(aVal, bVal):
    fitness = sin(7*aVal*bVal) + 1
    return fitness

def computeMeanFit(popVec):
    fits = [computeFitness(popVec[2*i], popVec[2*i+1]) for i in range(popsize)]
    return sum(fits)/popsize

def calcMeanA(popVec):
    AVals = [popVec[2*i] for i in range(popsize)]
    return sum(AVals)/popsize

def calcMeanB(popVec):
    BVals = [popVec[2*i+1] for i in range(popsize)]
    return sum(BVals)/popsize

def mutate(popVec):
    newVec = popVec
    for i in range(popsize):
        AVal = popVec[2*i]
        BVal = popVec[2*i + 1]
        if (random() < mutRate):
            newVec[2*i] = gauss(calcMeanA(popVec), mutStrength)
        else:
            newVec[2*i] = popVec[2*i]
        if (random() < mutRate):
            newVec[2*i+1] = gauss(calcMeanB(popVec), mutStrength)
        else:
            newVec[2*i+1] = popVec[2*i+1]

    return newVec

def selectPop(popVec):
    meanfit = computeMeanFit(popVec)
    nextPop = list(range(popsize*2))
    counter = 0
    while(counter < popsize):
        ind = choice(range(popsize))
        indfit = computeFitness(popVec[2*ind], popVec[2*ind + 1])
        if(uniform(0.0, (indfit+meanfit)**selectionStrength) < indfit**selectionStrength):
            nextPop[2*counter] = popVec[2*ind]
            nextPop[2*counter+1] = popVec[2*ind + 1]
            counter += 1

    return nextPop

def resizePop(popVec):
    newVec = list(range(popsize*2))
    inds = list(range(popsize))
    for i in range(popsize):
        inds[i] = choice(range(len(popVec)/2))
    for i in range(popsize):
        newVec[2*i] = popVec[2*inds[i]]
        newVec[2*i+1] = popVec[2*inds[i]+1]

    return newVec

# Simulation

print("Starting simulation for population of size", popsize)
testPop = list(range(popsize*2))
fitVec = list(range(numGens))
AMeans = list(range(numGens))
BMeans = list(range(numGens))
for i in range(popsize*2):
    testPop[i] = 0
switcher = 0
curr_counter = 0
if(changePop):
    popsize = bigPop
    for i in range(numGens):
        if(not i % 50):
            print("Currently on generation", i)
        if (not switcher):
            if curr_counter == bigGens:
                popsize = smallPop
                switcher = 1
                curr_counter = 0
                testPop = resizePop(testPop)
        elif curr_counter == smallGens:
            popsize = bigPop
            switcher = 0
            curr_counter = 0
            testPop = resizePop(testPop)
        testPop = mutate(testPop)
        testPop = selectPop(testPop)
        fitVec[i] = computeMeanFit(testPop)
        AMeans[i] = calcMeanA(testPop)
        BMeans[i] = calcMeanB(testPop)

else:
    for j in range(numGens):
         testPop = mutate(testPop)
         testPop = selectPop(testPop)
         fitVec[j] = computeMeanFit(testPop)
         AMeans[j] = calcMeanA(testPop)
         BMeans[j] = calcMeanB(testPop)

print("A Allele Vector: ")
print(AMeans)
print("B Allele Vector: ")
print(BMeans)
