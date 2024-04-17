%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: MAP Adaptation Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SuperVector, pointsPerMixtureComponent, mixtureWeights, meanVectors, posteriors] = AdaptGMM(UBM_GMM, SeqData)

    n_mixtures = UBM_GMM.NumComponents;

    % Mixture weights
    Pi_k = UBM_GMM.ComponentProportion';

    % Mixture means
    Mu_k = UBM_GMM.mu;

    % Mixture covariance matrices
    sigma_k = UBM_GMM.Sigma;
    for i = 1:n_mixtures
        Sigma_k(:,:,i) = diag(sigma_k(:,:,i));
    end
   
    % no of datapoints x dim. of each datapoint
    feature_set = SeqData;
    
    n_points = size(feature_set,1);

    %Expectation-Maximization Algorithm------------------------------------
        %E-step------------------------------------------------------------
            %Calculating Responsibilities
            gamma = posterior(UBM_GMM, feature_set);

            %Estimating weights of mixture components
            N_k = sum(gamma,1);

        %M-step------------------------------------------------------------
            %New Weights
            New_Pi_k = N_k./n_points;
            New_Pi_k = New_Pi_k';

            %New Means
            for k = 1:n_mixtures
                New_Mu_k(k,:) = (gamma(:,k)'*feature_set)./N_k(k);
            end

    %MAP Adaptation -- mean only
        for k = 1:n_mixtures
            alpha_k = N_k(k)/(N_k(k)+16);
            Mu_k_final(k,:) = alpha_k*New_Mu_k(k,:) + (1-alpha_k)*Mu_k(k,:);
        end

    %Super_Mean_Vector
        super_mean_vector = [];
        for k = 1:n_mixtures
            super_mean_vector = [super_mean_vector, Mu_k_final(k,:)];
        end
    %----------------------------------------------------------------------

    pointsPerMixtureComponent = N_k;
    SuperVector = super_mean_vector;
    mixtureWeights = New_Pi_k';
    meanVectors = Mu_k_final;
    posteriors = gamma;
            
end

