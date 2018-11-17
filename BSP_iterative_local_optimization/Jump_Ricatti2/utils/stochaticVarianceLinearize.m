function [Fi, Gi, ei] = stochaticVarianceLinearize(b_nom, u_nom)
    % Linearizing W
    global N nState mControl delx delu
    % Fi is a set of matrices, each matrix represents linearization of
    % column i of W
%     Fi = zeros(N, nState+nState*(1+nState)/2, nState+nState*(1+nState)/2, nState);
%     Gi = zeros(N, nState+nState*(1+nState)/2, mControl, nState);
%     ei = zeros(N, nState+nState*(1+nState)/2, nState);

    Fi = cell(N, nState);
    Gi = cell(N, nState);
    ei = cell(N, nState);
    
    for t = 1:N %for all time t
        x = b_nom(1:nState, t);
        sigma = b_nom(nState+1:end, t);
        b = b_nom(:, t);
        u = u_nom(:,t);
        for i = 1:nState+nState*(1+nState)/2 % over all belief terms
            db = zeros(nState+nState*(1+nState)/2,1);
            db(i) = delx;
            dWdb = (stochasticVarianceTermW(b+db,u) - stochasticVarianceTermW(b,u))./delx;
            for j = 1:nState % Filling for each column individually
%                 Fi(t,:,i,j) = dWdb(:,j);
                Fi{t,j}(:,i) = dWdb(:,j);
            end
        end
        for i = 1:mControl % over all belief terms
            du = zeros(mControl, 1);
            du(i) = delu;
            dWdu = (stochasticVarianceTermW(b,u+du) - stochasticVarianceTermW(b,u))./delu;
            for j = 1:nState % Filling for each column individually
%                 Gi(t,:,i,j) = dWdu(:,j);
                Gi{t,j}(:,i) = dWdu(:,j);
            end
        end
        etemp = stochasticVarianceTermW(b,u);
        for j = 1:nState
%             ei(t,:,:) = stochasticVarianceTermW(b,u);
            ei{t,j} = etemp(:,j);
        end
    end
end