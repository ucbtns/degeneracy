function [entropy,  energy, cost, accuracy, redundancy] = free_energy_decomp_pol(MDP)

% ========================================
% MDP: output from spm_MDP_X_VB
% Calculation of the Free Energy Decompositions:
%               - Degeneracy / Entropy:  E_Q(s)[-ln Q(s)] 
%               - Cost:  E_Q(s)[-ln P(s)] 
%               - Redundacy / Complexity: D_KL[Q(S)|P(S));
%               - Accuracy: E_Q(s)[lnP(O|S);  
%               - Energy: E_Q(s)[lnP(O,S)]

% Core: Q(s_t| s_t-1, pi), P(s_t| pi),  P(O|S_t, s_t-1, pi)
% Equations: F = Complexity (Red) - Accuracy
%                       = Energy + Entropy (Deg)
% Complexity   = Cost - Degeneracy

% ========================================
%--------------------------------------------------------------------------
% preclude numerical overflow
%--------------------------------------------------------------------------
spm_log = @(x)log(x + exp(-16));
% ========================================

% Specifying model components:
Nf = numel(MDP.D); 
for f = 1:Nf
    Ns(f) = numel(MDP.D{f}); 
end
t= MDP.T;

for f = 1:Nf
    Nu(f) = numel(MDP.b{f})/(Ns(f).^2); 
end

Ng = numel(MDP.a);    
for g = 1:Ng
    No(g) = size(MDP.a{g},1);
end

pg  = 1;                                                                     % prior precision
qg  = MDP.w;                                                           % posterior precision

% ========================================
% Degeneracy / Entropy = E_Q(s)[-ln Q(s)]: 
% ---------------------------------------------
% Accumulating hidden state (approx. beliefs) about each time point:

% Taking the beliefs about t, at all timepoints; Q(s) from bayesian
% averages: MDP.X :

for f = 1:Nf
   q_s{f} = (MDP.X{f}*diag(qg));
end


% Entropy for each factor, at each time point:
entropy = zeros(Nf,t);
for f = 1:Nf
    for i = 1:t
       sx = q_s{f}(:,i);
        entropy(f,i) = sx'*(-spm_log(sx));
    end
end

% Cost = E_Q(s)[ln P(s)]: 
% ---------------------------------------------
% Accumulating hidden state prob about each time point: P(s)

p_s = {};
for f = 1:Nf
        ps = zeros(Ns(f),t);
        for i = 1:t
            if i == 1 
                ps(:,i) = MDP.d{f}; 
            else        
                jz = 0;
                if Nu(f) > 1
                      j = MDP.s(f,i);
                     b = MDP.b{f}(:,:,j);                
                     jz = jz+ (b*q_s{f}(:,(i-1)));
                     jz = sum(jz,2);
                else
                    b = MDP.b{f};                
                    jz =jz + (b*q_s{f}(:,(i-1)));
                end             
              ps(:,i) = jz;
            end
        end
       p_s{f} = spm_norm(ps);
end

% Cost for each factor, at each time point:
cost = zeros(Nf,t);
for f = 1:Nf
    for i = 1:t
        px = p_s{f}(:,i);
        sx = q_s{f}(:,i);
        cost(f,i) = sx'*(-spm_log(px));
    end
end

% Redundacy / Complexity: D_KL[Q(S)|P(S)); 
% ---------------------------------------------
redundancy = zeros(Nf,t);
for f = 1:Nf
    for i = 1:t
        px = p_s{f}(:,i);
        sx = q_s{f}(:,i);
        redundancy(f,i) = sum(sx'*(spm_log(sx) - spm_log(px)));
    end
end


% Accuracy = E_Q(s)[lnP(O|S); 
% ---------------------------------------------
% Likelihood of hidden states for all the timepoints :
for i = 1:t
        L{i} = 1;
        for g = 1:Ng
            L{i} = L{i}.*spm_dot(spm_norm(MDP.a{g}),MDP.O{g,i});
      end
end       

% (current) Posterior over hidden factors
for i = 1:t
    for f= 1:Nf
        xq{f,i} = MDP.X{f}(:,i);
    end
end

% for timepoint; marginal likelihood over outcome factors
for i = 1:t
     for f= 1:Nf
         qL{f}(:,i) = spm_dot(L{i}, xq(:,i),f);
         lnqL{f}(:,i) = spm_log(qL{f}(:,i));
    end
end

% (negative) Accuracy for each factor, at each time point :
accuracy = zeros(Nf,t);

for f = 1:Nf
    for i = 1:t
        ql = lnqL{f}(:,i);
        sx =q_s{f}(:,i);
        accuracy(f,i) = sx'*ql;
    end
end

% Energy = E_Q(s)[lnP(O,S) = E_Q(s)[ln(P(o|s)P(s))] 
%                                             = E_Q(s)[lnP(o|s) + ln(P(s)]
% ---------------------------------------------
for i = 1:t
     for f = 1:Nf
         os{f}(:,i) = lnqL{f}(:,i) + spm_log(p_s{f}(:,i));        
    end
end

energy = zeros(Nf,t);

for f = 1:Nf
    for i = 1:t
        pos =  os{f}(:,i) ;
        sx =q_s{f}(:,i);
        energy(f,i) = sx'*(-pos);
    end
end

return 
