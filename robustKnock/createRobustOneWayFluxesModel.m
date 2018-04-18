function [model,yInd, notYInd,  m,n,notY, coupled, coupledYInd, coupledNotYInd]=createRobustOneWayFluxesModel(originalModel, chemicalInd, coupledFlag)

tmpModel=originalModel;
[~, rxn_num] = size(tmpModel.S);
yInd=[];
notYInd=[];
coupled=[];
coupledYInd=[];
coupledNotYInd=[];
coupledInd=[];
%vEqual=[];
%seperateFlux=[];

notGeneRelatedNum=0;
notKnockablesNum=0;
noInfluenceNum=0;

%create variables>=0
%matrices to change: S,c,lb,ub,rxns,rxnGeneMat, rxnNames,  in tmpModel
tmp_rxn_num=rxn_num;
if (coupledFlag)
    for i=1:rxn_num
        lbound=tmpModel.lb(i);
        ubound=tmpModel.ub(i);
        if(lbound < 0)
            sLeft=tmpModel.S(:,1:i-1);
            sRight=tmpModel.S(:,i+1:tmp_rxn_num);
            sVar=tmpModel.S(:,i);

            rxnGeneMatLeft=tmpModel.rxnGeneMat(1:i-1,:);
            rxnGeneMatRight=tmpModel.rxnGeneMat(i+1:tmp_rxn_num,:);
            rxnGeneMatVar=tmpModel.rxnGeneMat(i,:);

            cLeft=tmpModel.c(1:i-1);
            cRight=tmpModel.c(i+1:tmp_rxn_num);
            cVar=tmpModel.c(i);
            
            organismObjectiveLeft=tmpModel.organismObjective(1:i-1);
            organismObjectiveRight=tmpModel.organismObjective(i+1:tmp_rxn_num);
            organismObjectiveVar=tmpModel.organismObjective(i);
            
            C_chemicalLeft=tmpModel.C_chemical(1:i-1);
            C_chemicalRight=tmpModel.C_chemical(i+1:tmp_rxn_num);
            C_chemicalVar=tmpModel.C_chemical(i);
       
            lbLeft=tmpModel.lb(1:i-1);
            lbRight=tmpModel.lb(i+1:tmp_rxn_num);

            ubLeft=tmpModel.ub(1:i-1);
            ubRight=tmpModel.ub(i+1:tmp_rxn_num);

            rxnsLeft=tmpModel.rxns(1:i-1);
            rxnsRight=tmpModel.rxns(i+1:tmp_rxn_num);
            rxnsVar=tmpModel.rxns(i);

            rxnNamesLeft=tmpModel.rxnNames(1:i-1);
            rxnNamesRight=tmpModel.rxnNames(i+1:tmp_rxn_num);
            rxnNamesVar=tmpModel.rxnNames(i);

            if(ubound>0)    %need to create 2 reactions
                if (i == chemicalInd)
                    tmpModel.lb(i)=0;  % only secretion from this chemical!
                else
                    s1=sVar;
                    s2=-sVar;
                    c1=cVar;
                    c2=-cVar;
                    lb1=0;
                    lb2=0;
                    ub1=ubound;
                    ub2=-lbound;
                    organismObjective1=organismObjectiveVar;
                    organismObjective2=-organismObjectiveVar;
                    C_chemical1=C_chemicalVar;
                    C_chemical2=-C_chemicalVar;                    

                    Snew=[sLeft,s1,sRight,s2];
                    cnew=[cLeft;c1;cRight;c2];
                    lbnew=[lbLeft;lb1;lbRight;lb2];
                    ubnew=[ubLeft;ub1;ubRight;ub2];
                    rxnsNew=[rxnsLeft;rxnsVar;rxnsRight;strcat(rxnsVar,'2')];
                    rxnNamesNew=[rxnNamesLeft;rxnNamesVar;rxnNamesRight;strcat(rxnNamesVar,'2')];
                    rxnGeneMatNew=[rxnGeneMatLeft;rxnGeneMatVar;rxnGeneMatRight;rxnGeneMatVar];
                    organismObjectiveNew=[organismObjectiveLeft;organismObjective1;organismObjectiveRight;organismObjective2];
                    C_chemicalNew=[C_chemicalLeft;C_chemical1;C_chemicalRight;C_chemical2];

                    tmpModel.S=Snew;
                    tmpModel.c=cnew;
                    tmpModel.lb=lbnew;
                    tmpModel.ub=ubnew;
                    tmpModel.rxns=rxnsNew;
                    tmpModel.rxnNames=rxnNamesNew;
                    tmpModel.rxnGeneMat=rxnGeneMatNew;
                    tmpModel.organismObjective=organismObjectiveNew;
                    tmpModel.C_chemical=C_chemicalNew;

                    tmp_rxn_num=tmp_rxn_num+1;
                    coupled(end+1, 1:2)=[i,tmp_rxn_num];
                end
            end
            if(ubound<=0)    %need to change directions
                Snew=[sLeft,-sVar,sRight];
                cnew=[cLeft;-cVar;cRight];
                lbnew=[lbLeft;-ubound;lbRight];
                ubnew=[ubLeft;-lbound;ubRight];
                organismObjectiveNew=[organismObjectiveLeft;-organismObjectiveVar;organismObjectiveRight];
                C_chemicalNew=[C_chemicalLeft;C_chemicalVar;C_chemicalRight];
                
                tmpModel.S=Snew;
                tmpModel.c=cnew;
                tmpModel.lb=lbnew;
                tmpModel.ub=ubnew;
                tmpModel.organismObjective=organismObjectiveNew;
                tmpModel.C_chemical=C_chemicalNew;
            end
        end
    end
end
% construct constraints
model = tmpModel;
[n m] = size(model.S);

for i=1:rxn_num
    if (find(notYInd == i) )
        %go to the next index
    else
        [xCoupled,y]=find(coupled==i);
        if(xCoupled)
            coupledInd=coupled(xCoupled, y+1);
        end
        if (max(tmpModel.rxnGeneMat(i,:)) == 0)     %find if the flux is gene related
            notYInd(end+1)=i;
            notGeneRelatedNum=notGeneRelatedNum+1;
            if(xCoupled)
                coupledNotYInd(end+1)=coupledInd;
            end
        else             %find if the gene is knockable
            c=zeros(m,1);
            c(i)=1;
            ub=model.ub;
            ub(i)=0;
            lb=model.lb;
            lb(i)=0;

            if (xCoupled)
                c(coupledInd)=-1;
                ub(coupledInd)=0;
                lb(coupledInd)=0;
            end

            Prob = lpAssign( -c, model.S, model.row_lb, model.row_ub, lb, ub);
            Result = tomRun('cplex', Prob,0);

            if (Result.ExitFlag ~=0)
                notYInd(end+1)=i;
                if(xCoupled)
                    coupledNotYInd(end+1)=coupledInd;
                end
                notKnockablesNum=notKnockablesNum+1;
            else           %check if it makes a difference
                ub(i)=model.ub(i);
                lb(i)=model.lb(i);
                if (coupledInd)
                    ub(coupledInd)=model.ub(coupledInd);
                    lb(coupledInd)=model.lb(coupledInd);
                end

                Prob2 = lpAssign( -c, model.S, model.row_lb, model.row_ub, lb, ub);
                Result2 = tomRun('cplex', Prob2,0);

                Prob3 = lpAssign( c, model.S, model.row_lb, model.row_ub, lb, ub);
                Result3 = tomRun('cplex', Prob3,0);

                if (Result3.f_k == Result2.f_k)
                    notYInd(end+1)=i;
                    if(xCoupled)
                        coupledNotYInd(end+1)=coupledInd;
                    end
                    noInfluenceNum=noInfluenceNum+1;
                else
                    yInd(end+1)=i;
                    if(xCoupled)
                        coupledYInd(end+1)=coupledInd;
                    end
                end
            end
        end
    end
end

notY.notGeneRelatedNum=notGeneRelatedNum;
notY.notKnockablesNum=notKnockablesNum;
notY.noInfluenceNum=noInfluenceNum;

notYInd=sort(notYInd)';
yInd=sort(yInd)';
coupledNotYInd=sort(coupledNotYInd)';
coupledYInd=sort(coupledYInd)';

