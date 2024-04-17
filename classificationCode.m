%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: Classification Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [accuracy, clNames, sensitivity, specificity, precision, recall, f1_score, cMat] = classificationCode(featureVectorsTrain, classLabelsTrain, featureVectorsTest, classLabelsTest, classifiers, seed)
    % Classification using supervised learning on feature vectors

    % Initialize cell arrays to store confusion matrices for each classifier
    cMat = cell(length(classifiers), 1);

    % Get unique class labels
    uniqueLabels = unique(classLabelsTrain);
    n = length(unique(uniqueLabels));
    ord = 1:n; % Use ordinal values directly

    % Initialize clNames cell array to store classifier names
    clNames = cell(length(classifiers), 1);

    % Initialize accuracy cell array to store accuracy for each classifier
    accuracy = zeros(length(classifiers), 1); % Initialize as numeric array

    % Initialize weighted sensitivity, specificity, precision, and F1 score arrays
    sensitivity = zeros(length(classifiers), 1);
    specificity = zeros(length(classifiers), 1);
    precision = zeros(length(classifiers), 1);
    recall = zeros(length(classifiers), 1);
    f1_score = zeros(length(classifiers), 1);

    % Initialize weights for each class based on the number of true instances
    classWeights = histc(classLabelsTrain, uniqueLabels) / length(classLabelsTrain);

    % Iterate over classifiers
    for i = 1:length(classifiers)

        % Random seed
        rng(seed,'twister');
        
        % Classifier
        classifierName = classifiers{i, 1};

        % Train the classifier on the training data
        classifier = trainClassifier(featureVectorsTrain, classLabelsTrain, classifierName);

        % Test the classifier on the testing data
        predictedLabels = testClassifier(classifier, featureVectorsTest);

        % Compute confusion matrix
        cMat{i} = confusionmat(classLabelsTest, predictedLabels, 'Order', ord);

        % Compute accuracy
        cMatrix = cMat{i};
        accuracy(i) = round((trace(cMatrix)/length(classLabelsTest))*100,1);

        % Reset the metrics for the current classifier
        weightedSens = 0;
        weightedSpec = 0;
        weightedPrec = 0;
        weightedRecall = 0;
        weightedF1 = 0;

        % Compute weighted sensitivity, specificity, precision, recall, and F1 score for each class
        for j = 1:n
            TP = cMatrix(j, j);
            FN = sum(cMatrix(j, :)) - TP;
            FP = sum(cMatrix(:, j)) - TP;
            TN = sum(cMatrix(:)) - TP - FN - FP;

            classSens = TP / (TP + FN);
            classSpec = TN / (TN + FP);
            classPrec = TP / (TP + FP);
            classRecall = classSens;  % Recall is the same as sensitivity
            classF1 = 2 * TP / (2 * TP + FP + FN);

            % Update weighted metrics
            weightedSens = weightedSens + classSens * classWeights(j);
            weightedSpec = weightedSpec + classSpec * classWeights(j);
            weightedPrec = weightedPrec + classPrec * classWeights(j);
            weightedRecall = weightedRecall + classRecall * classWeights(j);
            weightedF1 = weightedF1 + classF1 * classWeights(j);
        end

        % Store weighted metrics
        sensitivity(i) = weightedSens;
        specificity(i) = weightedSpec;
        precision(i) = weightedPrec;
        recall(i) = weightedRecall;
        f1_score(i) = weightedF1;

        % Store classifier name
        clNames{i} = classifiers{i, 1};
    end

    % Display classification results (you can customize this part)
    displayClassificationResults(clNames, accuracy, precision, recall, specificity, f1_score);
end

function classifier = trainClassifier(featureVectorsTrain, classLabelsTrain, classifierName)
    % Train a classifier on the training data

    cn = unique(classLabelsTrain);
    n = length(cn);

    switch classifierName
        case 'LinearDiscriminant'
            classifier = fitcdiscr( ...
                featureVectorsTrain, ...
                classLabelsTrain, ...
                'DiscrimType', 'diaglinear', ...
                'Gamma', 0, ...
                'FillCoeffs', 'off', ...
                'ClassNames', cn);

        case 'LinearSVM'
            if(n==2)
                classifier = fitcsvm(...
                featureVectorsTrain, ...
                classLabelsTrain, ...
                'KernelFunction', 'linear', ...
                'PolynomialOrder', [], ...
                'KernelScale', 'auto', ...
                'BoxConstraint', 1, ...
                'Standardize', true, ...
                'ClassNames', cn);
            else
                template = templateSVM(...
                'KernelFunction', 'linear', ...
                'PolynomialOrder', [], ...
                'KernelScale', 'auto', ...
                'BoxConstraint', 1, ...
                'Standardize', true);

                classifier = fitcecoc(...
                featureVectorsTrain, ...
                classLabelsTrain, ...
                'Learners', template, ...
                'Coding', 'onevsone', ...
                'ClassNames', cn);
            end

        % case 'QuadraticSVM'
        %     if (n==2)
        %     classifier = fitcsvm( ...
        %         featureVectorsTrain, ...
        %         classLabelsTrain, ...
        %         'KernelFunction', 'polynomial', ...
        %         'PolynomialOrder', 2, ...
        %         'KernelScale', 'auto', ...
        %         'BoxConstraint', 1, ...
        %         'Standardize', true, ...
        %         'ClassNames', cn);
        %     else
        %         template = templateSVM(...
        %         'KernelFunction', 'polynomial', ...
        %         'PolynomialOrder', 2, ...
        %         'KernelScale', 'auto', ...
        %         'BoxConstraint', 1, ...
        %         'Standardize', true);
        %         classifier = fitcecoc(...
        %         featureVectorsTrain, ...
        %         classLabelsTrain, ...
        %         'Learners', template, ...
        %         'Coding', 'onevsone', ...
        %         'ClassNames', cn);
        %     end

        % case 'FineKNN'
        %     classifier = fitcknn( ...
        %         featureVectorsTrain, ...
        %         classLabelsTrain, ...
        %         'Distance', 'Euclidean', ...
        %         'Exponent', [], ...
        %         'NumNeighbors', 1, ...
        %         'DistanceWeight', 'Equal', ...
        %         'Standardize', true, ...
        %         'ClassNames', cn);

        otherwise
            error('Unsupported classifier.');
    end
end

function predictedLabels = testClassifier(classifier, featureVectorsTest)
    % Test the classifier on the testing data

    predictedLabels = predict(classifier, featureVectorsTest);
end

function displayClassificationResults(clNames, accuracy, precision, recall, specificity, f1_score)
    % Display header
    fprintf('%-20s%-12s%-12s%-12s%-12s%-12s\n', 'Classifier', 'Accuracy', 'Precision', 'Recall', 'Specificity', 'F1 Score');
    
    % Loop through each classifier and display its metrics
    for i = 1:length(clNames)
        fprintf('%-20s%.2f%%\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n', clNames{i}, accuracy(i), precision(i), recall(i), specificity(i), f1_score(i));
    end
end
