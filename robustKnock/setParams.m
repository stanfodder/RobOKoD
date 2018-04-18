function newProb=setParams(prob)
    newProb=prob;
    newProb.Alg = 2; % Depth First, then Breadth (Default Depth First)
    newProb.optParam.MaxIter = 100000; % Must increase iterations from default 500
    newProb.optParam.IterPrint = 0;
    newProb.PriLev = 1;    
    newProb.MIP.cpxControl.EPINT=10e-9;
    newProb.MIP.cpxControl.EPOPT=1e-8;
    newProb.MIP.cpxControl.EPRHS=1e-8;
    newProb.MIP.cpxControl.EPGAP=1e-6;
    newProb.MIP.cpxControl.EPAGAP=1e-8;
    newProb.MIP.cpxControl.NUMERICALEMPHASIS =1;
    newProb.MIP.cpxControl.SCAIND=1;
    %newProb.MIP.cpxControl.TILIM=120;
    
