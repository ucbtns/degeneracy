function mdp = GEN_MDP_SPEC(policy, content, agency, bet, alpha, z, w, c, t, redundancy, violation)

% ===========================================================================
% This function specifies the generative model (and process) for a simple
% speech repetition task. There are few parametres that can be tweaked to
% enable few things to be tested.

% Parametres:
% --------------------------------------------------------------------------------   
            % For redundancy: 
            % redudancy ==1 is redundant model specification; 0 reduced generative model form
            
            % For degeneracy: 
            % policy == 1 is deep policy; Otherwise one-shot policy
            % content == 1 is normal transition for content; Otherwise uncertain transition
            % agency == 1 is normal, other uncertain likelihood
            % bet = beta cofficient
            % z = level associated with lesion 
            
            % Additional: 
            % alpha: precision over action selection
            % c: (lower bound) how much uncertainty to resolve
            % t: the total duration of the trial
            
% ===========================================================================   
% Default HP:
%clear;clc
if ~exist('policy','var'), policy=1; end
if ~exist('content','var'), content=1; end
if ~exist('agency','var'), agency=1; end
if ~exist('bet','var'), bet=0; end
if ~exist('alpha','var'), alpha=4; end
if ~exist('z','var'), z=0.9; end
if ~exist('w','var'), w=0.8; end
if ~exist('c','var'), c=1/4; end
if ~exist('t','var'), t=4; end
if ~exist('redundancy','var'), redundancy=0; end
if ~exist('violation', 'var'), violation =0;end

fprintf('Policy Lesion: %i\n', policy)
fprintf('Target Lesion: %i\n', content)
fprintf('Agency Lesion: %i\n', agency)
fprintf('Redundancy Parameter: %i\n', redundancy)

% MODEL:
%--------------------------------------------------------------------------
% Prior beliefs about initial states (count): D
% This shows that we think the agent could be in any of the starting
% stages: Categorical distribution
%--------------------------------------------------------------------------

D{1} = [10 10 10 10 10]';          % Content: what I said?:  {'red' 'red' 'square' 'triangle', 'blue'} 
D{2} = [10 10 10 10]';                   % Target: what is the true context?: {'square' 'triangle', 'blue' 'red'} 
D{3} = [10 0 0]';                      % Time: {1, 2, 3}
D{4} = [10 0]';                         % Agency: {Listening, Speaking}


if violation == 1
    D{1} = [10 10 10 10 0]';
end
d = D; 


%--------------------------------------------------------------------------
% Probabilistic mapping from hidden states to outcomes: A Matrix, P(o|s)
% It is the likelihood where the elements are the probability of an outcome 
% under every combination of hidden states. We specify this as both the
% Categorical (A) & Dirchlet (a; prior)
%--------------------------------------------------------------------------

% Number of factors:
Nf = numel(D); 

% Number of factor in each state
for f = 1:Nf
    Ns(f) = numel(D{f}); 
end

No    = [2 3 5]; % Outcome modalities; 
                 % Proprioceptive: {'Other', 'Me'}; 
                 % Critic: {'Right', 'Wrong', 'Nothing'} 
                 % Auditory: {'square' 'triangle', 'blue' 'red' 'none'} 

% Total outcome levels: 3
Ng    = numel(No); 

% Matrix with all possible combinations (outcomes;states)'
e = 0.1e-4;
for g = 1:Ng
    A{g} = ones([No(g),Ns])*e; 
end

for f1 = 1:Ns(1) 
    for f2 = 1:Ns(2)  
        for f3 = 1:Ns(3)  
            for f4 = 1:Ns(4)
             
                       if f4 == 1 
                            A{1}(f4,:,:,:,f4) = 1;
                        elseif f4 == 2 
                            A{1}(f4,:,:,:,f4) = 1;
                       end 
                        
                         if f3 == 1 || f3 == 2
                            A{2}(3,:,:,f3,:) = 1;
                            
                        elseif (f2+1) == f1 
                            A{2}(1,f1,f2,3,2) = 1;
                        elseif (f2 ~= f1 && f1 <=4)
                            A{2}(2,f1,f2,3,2) = 1;                       
                        end   
                        A{2}(3, 1:5, :,3,1)= 1;
                        
                        A{2}(1,1,1,3,2) = 1;
                        A{2}(2,5,1:3,3,2) = 1;
                        
                        if f2 == f3 && f2 > 1
                            A{2}(2,f2,f3,3,2) = 1;
                        end

                        A{2}(2,4,4,3,2) = 1;
                        
                        if f3 == 1 && f1 <= 5 && (f4 == 1 || f4 ==2)
                            A{3}(:,f1, :,f3, f4) = [ones(4,4)*e+eye(4,4); ones(1,4)*e];
                        
                        elseif f3 == 2 && (f4 == 1|| f4 ==2)
                            A{3}(5,:,:,f3,f4) = 1;
                        
                        elseif (f3 == 3) && (f4 == 1 || f4 ==2)
                                A{3}(:,2:end,f2,f3,f4) =[ones(4,4)*e+eye(4,4); ones(1,4)*e];
                                A{3}(1,1,f2,f3,f4) =1;
                        end                  
            end                       
        end
    end       
end

for g = 1:Ng
    a{g} = A{g}*10;
end

if redundancy == 0
    for g = 1:Ng
        a{g}(:,2,:,:,:) = 10;
    end    
end


if agency ~=  1
    for f1 = 1:Ns(1) 
        for f2 = 1:Ns(2)  
            for f3 = 1:Ns(3)  
                for f4 = 1:Ns(4)                         
                          a{2}(:,:,f2,f3,f4) =  spm_softmax(z*log(A{2}(:,:,f2,f3,f4)+0.1))*10;                        
                end
            end
        end
    end
end
   

%--------------------------------------------------------------------------
% Transitions from t to t+1: P(S_t| S_t-1, pi)
% The B(u) matrices encode action-specific transitions
%--------------------------------------------------------------------------
for f = 1:Nf
        B{f} = e*ones(Ns(f));
end

 for i = 1:5
    B{1}(i,:,i) = 1;
 end 
    
B{1}(:,:,2:end) = B{1}(:,:,2:end) + e;

% State 2: Same context in play #u=1
B{2} = e*ones(4,4)+eye(4);

% State 3: Moving everytime except for last #u=1
B{3} = circshift( e*ones(3,3)+eye(3),1);
B{3}(:,3) = circshift(B{3}(:,3),2); 

% State 4: Moving in both directions #u=2
B{4}(:,:,1) = circshift( e*ones(2,2)+eye(2),1); % listen to speak; speak to listen
B{4}(:,:,2) =  (e*ones(2,2)+eye(2));                 % listen to listen; speak to speak

for f = 1:Nf
    b{f} = 10*B{f};
end

% Content lesion
if content ~= 1
    b{2} = spm_softmax(w*log(eye(4)+0.1));
end 


%--------------------------------------------------------------------------
% Priors: (utility) C; ln(P(o))
% Here, the agent prefers being correct when receiving feedback and only 
% likes speaking post the first two timepoints. 
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g}  = zeros(No(g),t);
end

C{2}(1,:) =  c; % does want to be right
C{2}(2,:) = -c; % does not want to be wrong

C{1}(1,1:2) = c; % does not want to speak before time point 3
C{1}(2,3:end) = c; 
     

%------------------------------------------------------------------------------
% Allowable policies (of depth T).  These are just sequences of actions
% Time * Poilicy * Factor
%------------------------------------------------------------------------------

% Deep Policies:

V(:,:,1) =     [1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5; 
                    1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5;
                    1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5;
                ];
V(:,:,2) = 1;
V(:,:,3) = 1;
V(:,:,4) =    [1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2;
                    2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2;
                    2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
                    ];
     
% One-shot Policies:
U(:,:,1) = [1 2 3 4 5 1 2 3 4 5];
U(:,:,2) = 1;
U(:,:,3) = 1;
U(:,:,4) = [1 1 1 1 1 2 2 2 2 2];

 
%--------------------------------------------------------------------------
% Specify the generative model
%--------------------------------------------------------------------------
% Hyper-parameters:
mdp.beta = exp(bet);            % Precision over precision (Gamma hyperprior)   
mdp.alpha = alpha;            % Precision over action selection: default 512

% Parameters:
mdp.T = t;                      % Time step

% = Likelihood
mdp.a = a;                      % Internal model
mdp.A = A;                      % Observation 

% = Transition
mdp.b = b;                      % Internal model
mdp.B = B;                      % Observation 

mdp.C = C;                      % Preferred outcomes
mdp.d = d; 
mdp.D = D;                      % Prior over initial states

% Specifying policy based on lesion type:
if policy == 1
    mdp.V = V;    
else
    mdp.U = U; 
end

% Labels:
mdp.Bname = {'Content', 'Context','Time','Agency'};
mdp.Aname = {'Proprioceptive', 'Feedback', 'Auditory'};
mdp.label.modality = {'Proprioceptive', 'Feedback', 'Auditory'};
mdp.label.outcome{1} = {'Other', 'Me'};
mdp.label.outcome{2} = {'Right', 'Wrong', 'None'};
mdp.label.outcome{3} = {'Square', 'Triange', 'Blue', 'Red','none'};

% Checking model:
mdp  = spm_MDP_check(mdp);
return 



