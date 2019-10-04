function simulation
% ===========================================================================
% This is a main function for simulating X individuals repeating the five words and
% recording their performance

% ===========================================================================            
clc; clear 

% Adding path directory
addpath d:/PhD/Code/spm/toolbox/DEM/
addpath d:/PhD/Code/spm/
addpath d:/PhD/Code/degeneracy/
addpath D:\PhD\Code\code

n = 10; % number of trials
X = 50; % number of subjects
E = 10; % ep responses

Nf = 4;
Ng = 3;


% FULL SPECIFICATION WITH REDUNDANT LEVEL
% ================================
policy_lesion = 1;  
t = 4;
content_lesion = 1; 
agency_lesion = 1;
bet = 0;
alpha =16;
pre = 2;
c = 1/4;
red = 1;

[REDmdp, REDa, REDMDP, REDqa, REDMDP_WS_O, REDMDP_WS_S, REDMDP_WS_H] = WORD_SIMULATION(n, X, policy_lesion, content_lesion, agency_lesion, bet, alpha, pre, pre, c, t,red);
[RMno_o, RMno_av_o, RMno_s_o] = calculate_average(REDMDP_WS_O, X, n, 4);

% Component of Free Energy:
for i = 1:X
    for j = 1:n
        [REDMDPentropy{i}{j},  REDMDPenergy{i}{j}, REDMDPcost{i}{j}, REDMDPaccuracy{i}{j}, REDMDPredundancy{i}{j}]  = free_energy_decomp(REDMDP{i}(j));
    end
end

avREDentropy = average_quantity(REDMDPentropy, Nf, t,n,X);
avREDenergy = average_quantity(REDMDPenergy, Nf, t,n,X);
avREDcost = average_quantity(REDMDPcost, Nf, t,n,X);
avREDaccuracy = average_quantity(REDMDPaccuracy, Nf, t,n,X);
avREDredundancy = average_quantity(REDMDPredundancy, Nf, t,n,X);
total(1,:) = [avREDredundancy - avREDaccuracy, avREDredundancy, avREDaccuracy, avREDentropy, avREDenergy, avREDcost];

tw_Emdp = GEN_MDP_SPEC(policy_lesion, content_lesion, agency_lesion, bet, alpha, pre, pre, c, t, red);   
tw_EMDP  = spm_MDP_VB_X(tw_Emdp);  
spm_figure('GetWin',['Figure ' num2str(1)]); clf
for i = 1:Nf
spm_MDP_VB_LFP(tw_EMDP(1),[],i);
end 



% FULL SPECIFICATION WITHOUT REDUNDANT LEVEL
% ================================
red = 0;

[mdp, a, MDP, qa, MDP_WS_O, MDP_WS_S, RMDP_WS_H] = WORD_SIMULATION(n, X, policy_lesion, content_lesion, agency_lesion, bet, alpha, pre, pre, c, t,red);
[Mno_o, Mno_av_o, Mno_s_o] = calculate_average(MDP_WS_O, X, n, 4);

% Component of Free Energy:
for i = 1:X
    for j = 1:n
        [MDPentropy{i}{j},  MDPenergy{i}{j}, MDPcost{i}{j}, MDPaccuracy{i}{j}, MDPredundancy{i}{j}]  = free_energy_decomp(MDP{i}(j));
    end
end

aventropy = average_quantity(MDPentropy, Nf, t,n,X);
avenergy = average_quantity(MDPenergy, Nf, t,n,X);
avcost = average_quantity(MDPcost, Nf, t,n,X);
avaccuracy = average_quantity(MDPaccuracy, Nf, t,n,X);
avredundancy = average_quantity(MDPredundancy, Nf, t,n,X);
total(2,:) = [avredundancy - avaccuracy, avredundancy, avaccuracy, aventropy, avenergy, avcost];

tw_ERmdp = GEN_MDP_SPEC(policy_lesion, content_lesion, agency_lesion, bet, alpha, pre, pre, c, t, red);   
tw_ERMDP  = spm_MDP_VB_X(tw_ERmdp);  
spm_figure('GetWin',['Figure ' num2str(1)]); clf
for i = 1:Nf
spm_MDP_VB_LFP(tw_ERMDP(1),[],i);
end


% Lesions: 
% =========================================
% B WITH REDUNDANT LEVEL
% ================================
content_lesion = 2; 
agency_lesion = 1;
alpha =4;
pre = 0.8;
red = 1;

[bREDmdp, bREDa, bREDMDP, bREDqa, bREDMDP_WS_O, bREDMDP_WS_S, bREDMDP_WS_H] = WORD_SIMULATION(n, X, policy_lesion, content_lesion,  agency_lesion, bet, alpha, pre,pre, c, t,red);
[bRMno_o, bRMno_av_o, bRMno_s_o] = calculate_average(bREDMDP_WS_O, X, n, 4);

% Component of Free Energy:
for i = 1:X
    for j = 1:n
        [bREDMDPentropy{i}{j}, bREDMDPenergy{i}{j}, bREDMDPcost{i}{j}, bREDMDPaccuracy{i}{j}, bREDMDPredundancy{i}{j}]  = free_energy_decomp(bREDMDP{i}(j));
    end
end

avbREDentropy = average_quantity(bREDMDPentropy, Nf, t,n, X);
avbREDenergy = average_quantity(bREDMDPenergy, Nf, t,n, X);
avbREDcost = average_quantity(bREDMDPcost, Nf, t,n, X);
avbREDaccuracy = average_quantity(bREDMDPaccuracy, Nf, t,n, X);
avbREDredundancy = average_quantity(bREDMDPredundancy, Nf, t,n, X);
total(3,:) = [avbREDredundancy - avbREDaccuracy avbREDredundancy, avbREDaccuracy, avbREDentropy, avbREDenergy, avbREDcost];
    
% =========================================
% B WITHOUT REDUNDANT LEVEL
% ================================
red = 0;

[bmdp, ba, bMDP, bqa, bMDP_WS_O, bMDP_WS_S, bMDP_WS_H] = WORD_SIMULATION(n, X, policy_lesion, content_lesion,  agency_lesion, bet, alpha, pre, pre, c, t,red);
[bMno_o, bMno_av_o, bMno_s_o] = calculate_average(bMDP_WS_O, X, n, 4);

% Component of Free Energy:
for i = 1:X
    for j = 1:n
        [bMDPentropy{i}{j}, bMDPenergy{i}{j}, bMDPcost{i}{j}, bMDPaccuracy{i}{j}, bMDPredundancy{i}{j}]  = free_energy_decomp(bMDP{i}(j));
    end
end

avbentropy = average_quantity(bMDPentropy, Nf, t,n,X);
avbenergy = average_quantity(bMDPenergy, Nf, t,n,X);
avbcost = average_quantity(bMDPcost, Nf, t,n,X);
avbaccuracy = average_quantity(bMDPaccuracy, Nf, t,n,X);
avbredundancy = average_quantity(bMDPredundancy, Nf, t,n,X);
total(4,:) = [avbredundancy - avbaccuracy, avbredundancy, avbaccuracy, avbentropy, avbenergy, avbcost];
        
% =========================================
% A WITH REDUNDANT LEVEL
% ================================

content_lesion = 1; 
agency_lesion = 2;
alpha =16;
pre2 = 0.9;
red = 1;

[aREDmdp, aREDa, aREDMDP, aREDqa, aREDMDP_WS_O, aREDMDP_WS_S, aREDMDP_WS_H] = WORD_SIMULATION(n, X, policy_lesion, content_lesion,  agency_lesion, bet, alpha, pre, pre2, c, t,red);
[aRMno_o, aRMno_av_o, aRMno_s_o] = calculate_average(aREDMDP_WS_O, X, n, 4);

% Component of Free Energy:
for i = 1:X
    for j = 1:n
        [aREDMDPentropy{i}{j}, aREDMDPenergy{i}{j}, aREDMDPcost{i}{j}, aREDMDPaccuracy{i}{j}, aREDMDPredundancy{i}{j}]  = free_energy_decomp(aREDMDP{i}(j));
    end
end

avaREDentropy = average_quantity(aREDMDPentropy, Nf, t,n,X);
avaREDenergy = average_quantity(aREDMDPenergy, Nf, t,n,X);
avaREDcost = average_quantity(aREDMDPcost, Nf, t,n,X);
avaREDaccuracy = average_quantity(aREDMDPaccuracy, Nf, t,n,X);
avaREDredundancy = average_quantity(aREDMDPredundancy, Nf, t,n,X);
total(5,:) = [avaREDredundancy - avaREDaccuracy, avaREDredundancy, avaREDaccuracy, avaREDentropy, avaREDenergy, avaREDcost];

% =========================================
% A WITHOUT REDUNDANT LEVEL
% ================================
red = 0;

[amdp, aa, aMDP, aqa, aMDP_WS_O, aMDP_WS_S, aMDP_WS_H] = WORD_SIMULATION(n, X, policy_lesion, content_lesion,  agency_lesion, bet, alpha, pre,  pre2, c, t,red);
[aMno_o, aMno_av_o, aMno_s_o] = calculate_average(aMDP_WS_O, X, n, 4);

% Component of Free Energy:
for i = 1:X
    for j = 1:n
        [aMDPentropy{i}{j}, aMDPenergy{i}{j}, aMDPcost{i}{j}, aMDPaccuracy{i}{j}, aMDPredundancy{i}{j}]  = free_energy_decomp(aMDP{i}(j));
    end
end

avaentropy = average_quantity(aMDPentropy, Nf, t,n,X);
avaenergy = average_quantity(aMDPenergy, Nf, t,n,X);
avacost = average_quantity(aMDPcost, Nf, t,n,X);
avaaccuracy = average_quantity(aMDPaccuracy, Nf, t,n,X);
avaredundancy = average_quantity(aMDPredundancy, Nf, t,n,X);
total(6,:) = [avaredundancy - avaaccuracy, avaredundancy, avaaccuracy, avaentropy, avaenergy, avacost];
    
% =========================================
% A & B WITH REDUNDANT LEVEL
% ================================

content_lesion = 2; 
agency_lesion = 2;
alpha =4;
red = 1;

[abREDmdp, abREDa, abREDMDP, abREDqa, abREDMDP_WS_O, abREDMDP_WS_S, abREDMDP_WS_H] = WORD_SIMULATION(n, X, policy_lesion, content_lesion,  agency_lesion, bet, alpha, pre, pre2, c, t,red);
[abRMno_o, abRMno_av_o, abRMno_s_o] = calculate_average(abREDMDP_WS_O, X, n, 4);

% Component of Free Energy:
for i = 1:X
    for j = 1:n
        [abREDMDPentropy{i}{j}, abREDMDPenergy{i}{j}, abREDMDPcost{i}{j}, abREDMDPaccuracy{i}{j}, abREDMDPredundancy{i}{j}]  = free_energy_decomp(abREDMDP{i}(j));
    end
end

avabREDentropy = average_quantity(abREDMDPentropy, Nf, t,n,X);
avabREDenergy = average_quantity(abREDMDPenergy, Nf, t,n,X);
avabREDcost = average_quantity(abREDMDPcost, Nf, t,n, X);
avabREDaccuracy = average_quantity(abREDMDPaccuracy, Nf, t,n,X);
avabREDredundancy = average_quantity(abREDMDPredundancy, Nf, t,n,X);

total(7,:) = [avabREDredundancy - avabREDaccuracy, avabREDredundancy, avabREDaccuracy, avabREDentropy, avabREDenergy, avabREDcost];

% =========================================
% A & B WITHOUT REDUNDANT LEVEL
% ================================
red = 0;

[abmdp, aba, abMDP, abqa, abMDP_WS_O, abMDP_WS_S, abMDP_WS_H] = WORD_SIMULATION(n, X, policy_lesion, content_lesion,  agency_lesion, bet, alpha, pre, pre2, c, t,red);
[abMno_o, abMno_av_o, abMno_s_o] = calculate_average(abMDP_WS_O, X, n, 4);
ab_mdp = WORD_REACTION(n,X,policy_lesion, content_lesion, agency_lesion, bet, alpha, pre, pre, c, t, red);

% Component of Free Energy:
for i = 1:X
    for j = 1:n
        [abMDPentropy{i}{j}, abMDPenergy{i}{j}, abMDPcost{i}{j}, abMDPaccuracy{i}{j}, abMDPredundancy{i}{j}]  = free_energy_decomp(abMDP{i}(j));
    end
end

avabentropy = average_quantity(abMDPentropy, Nf, t,n,X);
avabenergy = average_quantity(abMDPenergy, Nf, t,n,X);
avabcost = average_quantity(abMDPcost, Nf, t,n,X);
avabaccuracy = average_quantity(abMDPaccuracy, Nf, t,n,X);
avabredundancy = average_quantity(abMDPredundancy, Nf, t,n,X);
total(8,:) = [avabredundancy - avabaccuracy, avabredundancy, avabaccuracy, avabentropy, avabenergy, avabcost];


% Output:  for behavioural graphics: 

final_o = spm_unvec([RMno_av_o Mno_av_o bRMno_av_o bMno_av_o aRMno_av_o aMno_av_o abRMno_av_o abMno_av_o], zeros(X*8, 1));
t = array2table(final_o); 
writetable(t, '~\red_deg_output2.csv');

return



