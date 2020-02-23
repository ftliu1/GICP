function outputMatrix = kdNeighbour(kd_tree, target, k, nodeNumber)
%
% This function is used to find k nearest points in a kd tree structure to
% a desired point
%
% INPUT:
% kd_tree:          a built KD tree structure
% target:             the 1 x d point searching for nearest point
% k:                   the number of points needed finding
% nodeNumber:  give the postion in the tree where the current node is
%
% OUTPUT:
% outputDataset:       the k nearest points founded in kd_tree
%
%
% Fangtong Liu
% University of Michigan
% ftliu@umich.edu
% 2020/02/12

global kdBuildTree
global distanceArray    % to store k closest distance
global dataMatrix        % to sotre k closest point

%% Initialization
if nargin == 3
    kdBuildTree = kd_tree;
    clear kd_tree;
    distanceArray = [];
    dataMatrix = [];
    nodeNumber = 1;
end

%% Find k nearest points
% if node input is the parent of root
if ~nodeNumber
    leftChild = kdBuildTree(1).left;
    rightChild = kdBuildTree(1).right;
    if kdBuildTree(leftChild).visit && kdBuildTree(rightChild).visit
        outputMatrix = dataMatrix;
        return;
    else
        maxDis = max(distanceArray);
        [dis_median, intersection] = checkIntersect(kdBuildTree, target, 1, maxDis);
        % if points stored is less than k or the target circle with maxDis
        % intersects the spliting hyperplane of root, explore the other side
        if (length(distanceArray) < k || intersection) && ~kdBuildTree(1).visit
            % update distanceArray and dataMatrix
            tempDis = getDistance(kdBuildTree(1).data, target);
            [dataMatrix, distanceArray] = updateNeighbour(tempDis, kdBuildTree(1).data, distanceArray, dataMatrix, k);
            % update visit
            kdBuildTree(1).visit = 1;        
            % if target is left to parent 
            if dis_median <= 0
                kdNeighbour(kdBuildTree, target, k, rightChild);
            % if target is right to parent
            elseif dis_median > 0
            kdNeighbour(kdBuildTree, target, k, leftChild);
            end
        % if not intersect
        else
            kdBuildTree(1).visit = 1;
            return;
        end
    end
end

% if current node is not visited or not a leaf, find next node
if ~kdBuildTree(nodeNumber).visit
    while ~kdBuildTree(nodeNumber).leaf && ~kdBuildTree(nodeNumber).visit
        leftChild = kdBuildTree(nodeNumber).left;
        rightChild = kdBuildTree(nodeNumber).right;
        if ~kdBuildTree(leftChild).visit && target(kdBuildTree(nodeNumber).split) <= kdBuildTree(nodeNumber).median
            nodeNumber = leftChild;
        elseif ~kdBuildTree(rightChild).visit
            nodeNumber = rightChild;
        else
            break;
        end
    end
    
    % update visit
    kdBuildTree(nodeNumber).visit = 1;
    
    % update distanceArray and dataMatrix
    pointFound = kdBuildTree(nodeNumber).data;
    tempDis = getDistance(pointFound, target);
    [dataMatrix, distanceArray] = updateNeighbour(tempDis, pointFound, distanceArray, dataMatrix, k);
    
    %  find next node if the root is not visited
    parentNumber = kdBuildTree(nodeNumber).parent;
    if parentNumber
        leftChild = kdBuildTree(parentNumber).left;
        rightChild = kdBuildTree(parentNumber).right;
        % if both children has been visited, explore grandparent
        if kdBuildTree(leftChild).visit && kdBuildTree(rightChild).visit
            if kdBuildTree(parentNumber).parent == 1
                maxDis = max(distanceArray);
                [~, intersection] =  checkIntersect(kdBuildTree, target, 1, maxDis);
                if ~intersection
                    outputMatrix = dataMatrix;
                    return;
                end
            end
                kdNeighbour(kdBuildTree, target, k, kdBuildTree(parentNumber).parent);
        else
            maxDis = max(distanceArray);
            [dis_median, intersection] = checkIntersect(kdBuildTree, target, parentNumber, maxDis);
            % if points stored is less than k or if the target cicrle with
            % max distance stored as radius intersects the spliting hyperplane
            % check the other side of parent node
            if length(distanceArray) < k || intersection
                % update distanceArray and dataMatrix
                if ~kdBuildTree(parentNumber).visit
                    pointFound = kdBuildTree(parentNumber).data;
                    tempDis = getDistance(pointFound, target);
                    [dataMatrix, distanceArray] = updateNeighbour(tempDis, pointFound, distanceArray, dataMatrix, k);
                end
        
                % update visit
                kdBuildTree(parentNumber).visit = 1;
                % if target is left to parent 
                if dis_median <= 0
                    kdNeighbour(kdBuildTree, target, k, rightChild);
                % if target is right to parent
                elseif dis_median > 0 
                    kdNeighbour(kdBuildTree, target, k, leftChild);
                end
            % if not intersect and more than k nodes stored
            else
                % update visit
                kdBuildTree(parentNumber).visit = 1;
                if kdBuildTree(parentNumber).parent == 1
                    maxDis = max(distanceArray); 
                    [~, intersection] =  checkIntersect(kdBuildTree, target, 1, maxDis);
                    if ~intersection
                        outputMatrix = dataMatrix;                  
                        return;
                    end
                end
            kdNeighbour(kdBuildTree, target, k, kdBuildTree(parentNumber).parent);
            end
        end      
    % if the parent node is root
    else
        kdNeighbour(kdBuildTree, target, k, 0);
    end
end

%% Output
if nargin == 3
    outputMatrix = dataMatrix;
    clear global kdBuildTree;
    clear global dataMatrix;
    clear global distanceArray;
end
end

%%  This fuction is used to calcuate Euclidean distance between two point
function distance = getDistance(source, target)
temp = (source - target).^2;
distance = sqrt(sum(temp));
end

%% This function is used to update distanceArray and dataMatrix
function [dataMatrix, distanceArray] = updateNeighbour(tempDis, pointFound, distanceArray, dataMatrix, k)
     % if number of points stored is less than k, intersect new data
     if length(distanceArray) < k
         dataMatrix = [dataMatrix; pointFound];
         distanceArray = [distanceArray tempDis];
     % if number of points stored is more than k and tempDis is smaller
     % than the max distance sotored, deleted the old and intersect the new
     else
         [maxDistance, maxIndex] = max(distanceArray);
         if tempDis < maxDistance
             distanceArray(maxIndex) = tempDis;
             dataMatrix(maxIndex, :) = pointFound;
         end
     end
end

%% This function is used to check whether intersection exists between:
%  the circle whose center is target and whose radius is minDis
%  the spliting hyperplane of currentNode
function [dis_median, intersection] = checkIntersect(kdBuildTree, target, nodeNumber, maxDis)
split = kdBuildTree(nodeNumber).split;
median = kdBuildTree(nodeNumber).median;
dis_median = target(split)-median;
if abs(dis_median) <= maxDis
    intersection = 1;
else
    intersection = 0;
end
end

