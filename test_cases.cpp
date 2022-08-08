/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "test_cases.h"
#include "simulation_parameters.h"
#include "particles.h"
#include "force_calculation_methods.h"

void compareForceCalculationMethodsAcceleration(float dt, float simTime,
        int numItersStored, int numThreads, const char* folderName,
        const char* infoFileName, const char* fieldsName, 
        const char* trianglesName){
    
    ParticleSystem particleSystem;
    particleSystem.params.initialize(dt, simTime, numItersStored, numThreads,
            folderName);
    particleSystem.params.readInputInfoFile(infoFileName);
    particleSystem.params.createOutputInfoFile();
    particleSystem.params.initializeDataStorageFiles();
    
    particleSystem.params.timerInfo = TimerOptions::minimized;
    
    particleSystem.initialize(Initialization::singleEmitter,
            InjectionTimingOptions::constantRate, folderName);
    
    particleSystem.initializeBackgroundEField(BackgroundEfield::singleEmitter, 
            folderName, fieldsName, trianglesName);
    
    DirectForceCalculation df(&particleSystem);
    
    float tempTheta[] = {0.1, 0.5, 1.0};
    int numThetaVals = sizeof(tempTheta) / sizeof(tempTheta[0]);
    
    int tempExpansions[] = {7, 5, 3};
    int numExpansionVals = sizeof(tempExpansions) / sizeof(tempExpansions[0]);
//    int numExpansionVals = 0;
    Array<ForceCalculationMethod*> fca(numThetaVals + 2 * numExpansionVals);
    
    MultipoleMethodConstants c;
    MultipoleMethodVariables v;
    FmmType t = FmmType::fmm;
    
    BarnesHutForceCalculation bh_100(&particleSystem);
    bh_100.setTheta(tempTheta[0]);
    bh_100.setDomainToggle(DomainToggle::rectangularPrismCells);
    BarnesHutForceCalculation bh_500(&particleSystem);
    bh_500.setTheta(tempTheta[1]);
    bh_500.setDomainToggle(DomainToggle::rectangularPrismCells);
    BarnesHutForceCalculation bh_1000(&particleSystem);
    bh_1000.setTheta(tempTheta[2]);
    bh_1000.setDomainToggle(DomainToggle::rectangularPrismCells);
    
    fca[0] = &bh_100;
    fca[1] = &bh_500;
    fca[2] = &bh_1000;
    
    MultipoleMethodForceCalculation fmm_7_1(&particleSystem, &c, &v, t);
    fmm_7_1.setNumExpansions(tempExpansions[0]);
    fmm_7_1.setMaxLevel(-1);
    MultipoleMethodForceCalculation fmm_7_2(&particleSystem, &c, &v, t);
    fmm_7_2.setNumExpansions(tempExpansions[0]);
    fmm_7_2.setMaxLevel(0);
    fmm_7_1.olderSibling = &fmm_7_2;
    
    MultipoleMethodForceCalculation fmm_5_1(&particleSystem, &c, &v, t);
    fmm_5_1.setNumExpansions(tempExpansions[1]);
    fmm_5_1.setMaxLevel(-1);
    MultipoleMethodForceCalculation fmm_5_2(&particleSystem, &c, &v, t);
    fmm_5_2.setNumExpansions(tempExpansions[1]);
    fmm_5_2.setMaxLevel(0);
    fmm_5_1.olderSibling = &fmm_5_2;
    
    MultipoleMethodForceCalculation fmm_3_1(&particleSystem, &c, &v, t);
    fmm_3_1.setNumExpansions(tempExpansions[2]);
    fmm_3_1.setMaxLevel(-1);
    MultipoleMethodForceCalculation fmm_3_2(&particleSystem, &c, &v, t);
    fmm_3_2.setNumExpansions(tempExpansions[2]);
    fmm_3_2.setMaxLevel(0);
    fmm_3_1.olderSibling = &fmm_3_2;
    
    fca[numThetaVals] = &fmm_7_1;
    fca[numThetaVals + 1] = &fmm_7_2;
    fca[numThetaVals + 2] = &fmm_5_1;
    fca[numThetaVals + 3] = &fmm_5_2;
    fca[numThetaVals + 4] = &fmm_3_1;
    fca[numThetaVals + 5] = &fmm_3_2;
    
    LeapfrogDKD leap(particleSystem, df);
    leap.setComparisons(fca);
    
    leap.runSimulation();
}

void compareForceCalculationMethodsEvolution(float dt, float simTime,
        int numItersStored, int numThreads, const char* folderName,
        const char* infoFileName, const char* fieldsName, 
        const char* trianglesName){
    
    ParticleSystem particleSystem;
    particleSystem.params.initialize(dt, simTime, numItersStored, numThreads,
            folderName);
    particleSystem.params.readInputInfoFile(infoFileName);
    particleSystem.params.createOutputInfoFile();
    particleSystem.params.initializeDataStorageFiles();
    
    particleSystem.params.timerInfo = TimerOptions::minimized;
    
    particleSystem.initialize(Initialization::singleEmitter,
            InjectionTimingOptions::constantRate, folderName);
    
    particleSystem.initializeBackgroundEField(BackgroundEfield::singleEmitter, 
            folderName, fieldsName, trianglesName);
    
    DirectForceCalculation df(&particleSystem);
    
    LeapfrogDKD leapDF(particleSystem, df);
    leapDF.runSimulation();
    
    float tempTheta = 0.5;
    BarnesHutForceCalculation bh_500(&particleSystem);
    bh_500.setTheta(tempTheta);
    bh_500.setDomainToggle(DomainToggle::rectangularPrismCells);
    
    LeapfrogDKD leapBH(particleSystem, bh_500);
    leapBH.runSimulation();
    
    int tempExpansions = 5;
    
    MultipoleMethodConstants c;
    MultipoleMethodVariables v;
    FmmType t = FmmType::fmm;
    Array<ForceCalculationMethod*> fca(1);
    
    MultipoleMethodForceCalculation fmm_5_1(&particleSystem, &c, &v, t);
    fmm_5_1.setNumExpansions(tempExpansions);
    fmm_5_1.setMaxLevel(-1);
    MultipoleMethodForceCalculation fmm_5_2(&particleSystem, &c, &v, t);
    fmm_5_2.setNumExpansions(tempExpansions);
    fmm_5_2.setMaxLevel(0);
    fmm_5_1.olderSibling = &fmm_5_2;
    
    fca[0] = &fmm_5_2;
    
    LeapfrogDKD leapFMM(particleSystem, fmm_5_1);
    leapFMM.setComparisons(fca);
    leapFMM.runSimulation();
}

void compareForceCalclationMethodsEndSim(float dt, float simTime,
        int numItersStored, int numThreads, const char* folderName,
        const char* infoFileName, const char* fieldsName, 
        const char* trianglesName){
    int i, j;
    
    ParticleSystem particleSystem;
    particleSystem.params.initialize(dt, simTime, numItersStored, numThreads,
            folderName);
    particleSystem.params.readInputInfoFile(infoFileName);
    particleSystem.params.createOutputInfoFile();
    particleSystem.params.initializeDataStorageFiles();
    
    particleSystem.params.timerInfo = TimerOptions::minimized;
    
    particleSystem.initialize(Initialization::singleEmitter,
            InjectionTimingOptions::constantRate, folderName);
    
    particleSystem.initializeBackgroundEField(BackgroundEfield::singleEmitter, 
            folderName, fieldsName, trianglesName);
    
    // Set timestep and dN to zero
    // Set numIters to number of trials to perform
    particleSystem.params.timeStep = 0;
    particleSystem.params.dN = 0;
    particleSystem.params.numIters = 10;
    
    char* initConditionsName = concat(folderName, "IC_endSim.txt");
    FILE* ICp = fopen(initConditionsName,"r");
    for (i = 0; i < particleSystem.chargedParticles.maxNumParticles; i++){
        fscanf(ICp,"%g, %g, %g",&particleSystem.chargedParticles.position[i].x,
                &particleSystem.chargedParticles.position[i].y,
                &particleSystem.chargedParticles.position[i].z);
    }
    fclose(ICp);
    
    float theta_array[] = {1.0, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01, 0};
    int p_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int maxLevel_array[] = {-1, 0, 1, 2, 3, 4, 5, 6};
    
    int sizeThetaArray = sizeof(theta_array) / sizeof(theta_array[0]);
    int size_pArray = sizeof(p_array) / sizeof(p_array[0]);
    int sizeMaxLevelArray = sizeof(maxLevel_array) / sizeof(maxLevel_array[0]);
    
    MultipoleMethodConstants c;
    MultipoleMethodVariables v;
    FmmType t = FmmType::fmm;
    Array<ForceCalculationMethod*> fca(1);
    DirectForceCalculation df(&particleSystem);
    
    for (i = 0; i < sizeThetaArray; i++){
        printf("\nStarting BH with theta = %f\n", theta_array[i]);
        BarnesHutForceCalculation bh(&particleSystem);
        fca[0] = &bh;
        LeapfrogDKD leap(particleSystem, df);
        leap.setComparisons(fca);
        leap.runSimulation();
    }
    
    for (i = 0; i < size_pArray; i++){
        for (j = 0; j < sizeMaxLevelArray; j++){
            printf("\nStarting FMM with p = %d, maxLevel = %d\n",p_array[i],maxLevel_array[j]);
            MultipoleMethodForceCalculation fmm(&particleSystem, &c, &v, t);
            fmm.setNumExpansions(p_array[i]);
            fmm.setMaxLevel(maxLevel_array[j]);
            fca[0] = &fmm;
            LeapfrogDKD leap(particleSystem, df);
            leap.setComparisons(fca);
            leap.runSimulation();
        }
    }
}