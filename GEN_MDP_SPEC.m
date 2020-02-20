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
if ~exist('policy','var'), policy=1; end
if ~exist('content','var'), content=1; end
if ~exist('agency','var'), agency=1; end
if ~exist('bet','var'), bet=0; end
if ~exist('alpha','var'), alpha=8; end
if ~exist('z','var'),  z=0.85; end
if ~exist('w','var'), w=0.8; end
if ~exist('c','var'), c=1/2; end
if ~exist('t','var'), t=3; end
if ~exist('redundancy','var'), redundancy=1; end
if ~exist('violation', 'var'), violation =0;end

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
D{2} = [10 10 10 10]';                   % Target: what is the true context?:   { 'red' 'square' 'triangle', 'blue'} 
D{3} = [128 0 0]';                      % Time: {1, 2, 3}


if violation == 1
    D{1} = [10 10 10 10 0]';
end

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

No    = [2 3 4]; % Outcome modalities; 
                 % Proprioceptive: {'Other', 'Me'}; 
                 % Critic: {'Right', 'Wrong', 'Nothing'} 
                 % Auditory:   { 'red' 'square' 'triangle', 'blue'} 

% Total outcome levels: 3
Ng    = numel(No); 

% Matrix with all possible combinations (outcomes;states)'
e = 0;
for g = 1:Ng
    A{g} = ones([No(g),Ns])*e; 
end

for f1 = 1:Ns(1) 
    for f2 = 1:Ns(2)  
        for f3 = 1:Ns(3)  
             
                       if (f3 == 1 )
                            A{1}(1,:,:,f3) = 1;
                        elseif (f3 == 2 || f3 == 3)
                            A{1}(2,:,:,f3) = 1;
                       end 
                        
                         if f3 == 1 || f3 == 2 
                            A{2}(3,:,:,f3) = 1;
                            
                        elseif (f2+1) == f1 
                            A{2}(1,f1,f2,3) = 1;
                        elseif (f2 ~= f1)
                            A{2}(2,f1,f2,3) = 1;                       
                        end   
                        
                        A{2}(1,1,1,3) = 1;
                        
                        if f2 == f3 && f2 > 1
                            A{2}(2,f2,f3,3) = 1;
                        end

                        A{2}(2,4,4,3) = 1;
                        
                        if f3 == 1 
                            A{3}(:,f1, :,f3) = [ones(4,4)*e+eye(4,4)];                   
                        else
                                A{3}(:,2:end,f2,f3) =[ones(4,4)*e+eye(4,4)];
                                A{3}(1,1,f2,f3) =1;
                        end                                        
        end
    end       
end

if (content ~= 1 || agency ~=1 )
    % Switching off learning:
    for g = 1:Ng
        a{g} = A{g}*10;
    end
else
    % Learning:
    for g = 1:Ng
        a{g} = (A{g})+1/2;
     end
end


if redundancy == 0
    % Increasing the connection strength; 
    % disconnection from the model
    for g = 1:Ng
        a{g}(:,2,:,:,:) = 10;
    end    
end


if agency ~=  1
    for f1 = 1:Ns(1) 
        for f2 = 1:Ns(2)  
            for f3 = 1:Ns(3)                         
                          a{2}(:,:,f2,f3) =  spm_softmax(z*log(A{2}(:,:,f2,f3)+0.5));                        
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
    
B{1}(:,:,2:end) = B{1}(:,:,2:end) ;

% State 2: Same context in play #u=1
B{2} = eye(4);

% State 3: Moving everytime except for last #u=1
B{3} = circshift(eye(3),1);
B{3}(:,3) = circshift(B{3}(:,3),2); 

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

%------------------------------------------------------------------------------
% Allowable policies (of depth T).  These are just sequences of actions
% Time * Poilicy * Factor
%------------------------------------------------------------------------------

% Deep Policies:
V(:,:,1) =     [1 2 3 4 5; 
                    1 2 3 4 5;
                ];
V(:,:,2) = 1;
V(:,:,3) = 1;
 
%--------------------------------------------------------------------------
% Specify the generative model
%--------------------------------------------------------------------------
% Hyper-parameters:
mdp.beta = exp(bet);            % Precision over precision (Gamma hyperprior)   
mdp.alpha = alpha;            % Precision over action selection: default 512

if (content ~= 1 ||  agency ~= 1)
    mdp.eta = 0;
end

% Parameters:
mdp.T = t;                      % Time step

% = Likelihood
mdp.a = a;                      % Internal model
mdp.A = A;                      % Observation model

% = Transition
mdp.b = b;                      % Internal model
mdp.B = B;                      % Observation model

mdp.C = C;                      % Preferred outcomes
mdp.D = D;                      % Prior over initial states

% Specifying policy based on lesion type:

 mdp.V = V;    

% Labels:
mdp.Bname = {'Content', 'Context','Time'};
mdp.Aname = {'Proprioceptive', 'Feedback', 'Auditory'};
mdp.label.modality = {'Proprioceptive', 'Feedback', 'Auditory'};
mdp.label.name{1} = {'red' 'red' 'square' 'triangle', 'blue'} ;
mdp.label.name{2} = { 'red' 'square' 'triangle', 'blue'} ;
mdp.label.name{3} = { '1' '2' '3',} ;
mdp.label.outcome{1} = {'Other', 'Me'};
mdp.label.outcome{2} = {'Right', 'Wrong', 'None'};
mdp.label.outcome{3} =   { 'red' 'square' 'triangle', 'blue'} ;


% Checking model:
mdp  = spm_MDP_check(mdp);

return 


