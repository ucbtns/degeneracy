function [allmdp, alla, allMDP, allqa, mdp_outcomes, mdp_states, mdp_H] = WORD_SIMULATION(n, X,policy, content, agency, bet, alpha, pre,pre2, c, t, redundancy, violation)

% ===========================================================================
% This functionaims to simulate X individuals repeating the five words and
% recording their performance:
%  - Simulate both multiple agents individuals via different seeds                                     
% ===========================================================================            

% default parameters
if ~exist('n','var'), n=10; end
if ~exist('X','var'), X=50; end
if ~exist('policy','var'), policy=1; end
if ~exist('content','var'), content=1; end
if ~exist('agency','var'), agency=1; end
if ~exist('bet','var'), bet=0; end
if ~exist('alpha','var'), alpha=16; end
if ~exist('pre','var'), pre=0.8; end
if ~exist('pre2','var'), pre2=0.8; end
if ~exist('c','var'), c=1/2; end
if ~exist('t','var'), t= 3; end
if ~exist('redundancy','var'), redundancy= 1; end
if ~exist('violation', 'var'), violation =0;end

% Initialise:
mdp_outcomes = zeros(n*X, (3*3 + 1));
mdp_states = zeros(n*X, (3*3 + 1));
mdp_H = zeros(n*X, (t+1));

num = 1;

allmdp = {};
allMDP = {};
alla = {};
allqa = {};

for x = 1:X
    
    % setting the seed:
    rng(x);
    clear mdp MDP
    
    % specifying the set-up:
    mdp = GEN_MDP_SPEC(policy, content, agency, bet, alpha, pre, pre2, c, t, redundancy, violation);

    a = mdp.a;
    [mdp(1:n)] = deal(mdp);      
    
    MDP  = spm_MDP_VB_X(mdp);
    
    % storing results:
    mdp_o = zeros(n, (3*3 + 1));
    mdp_s = zeros(n, (3*3 + 1));
    mdp_h = zeros(n, size(MDP(1).H,2) + 1);
    
    for i = 1:n 
        mdp_o(i,:) = [x reshape(transpose(MDP(i).o), 3*3,1)'];
        mdp_s(i,:) = [x reshape(transpose(MDP(i).s), 3*3,1)'];
        mdp_h(i,:) = [x MDP(i).H];
    end
   
    mdp_outcomes(num:num+(n-1), :) = mdp_o;
    mdp_states(num:num+(n-1), :) = mdp_s;
    mdp_H(num:num+(n-1), :) = mdp_h;
    
    allmdp{x} = mdp;
    allMDP{x} = MDP;
    alla{x} = a;
    allqa{x} = MDP.a;
    
    num = num + n;    
end

return 