function outputPoint = kdSearch(kd_tree, target, minDis, nodeNumber)
%
% This function is used to find the nearest point in a kd tree structure to
% a desired point
%
% INPUT:
% kd_tree:          a built KD tree structure
% target:             the 1 x d point searching for nearest point
% minDis:           the minimal distance from nodes in KD tree to target
% nodeNumber:  give the postion in the tree where the current node is
%
% minDis, nodeNumber do not need to assign
%
% OUTPUT:
% outputPoint:       the nearest point founded in kd_tree
%
%
% Fangtong Liu
% University of Michigan
% ftliu@umich.edu
% 2020/02/07

global kdBuildTree
global nearestPoint

%% Initialization
if nargin == 2
    nodeNumber = 1; % assign the root number
    kdBuildTree = kd_tree;
    clear kd_tree
    minDis = inf;
end

%% Find nearest point
% check if the current node is the parent of root
if ~nodeNumber
    leftChild = kdBuildTree(1).left;
    rightChild = kdBuildTree(1).right;
    if kdBuildTree(leftChild).visit && kdBuildTree(rightChild).visit
        outputPoint = nearestPoint;
        return;
    else
        % check whether target has intersection with root median
        [dis_median, intersection] = checkIntersect(kdBuildTree, target, 1, minDis);
        if intersection && ~kdBuildTree(1).visit
            % check if root is closer to the target
            tempDis = getDistance(kdBuildTree(1).data, target);   
            [nearestPoint, minDis] = updateDis(tempDis, kdBuildTree(1).data, minDis, nearestPoint);
            % update the root as visited
            kdBuildTree(1).visit = 1;       
            % if the target is left to the root, check its right child
            if dis_median <= 0
                kdSearch(kdBuildTree, target, minDis, rightChild);
            % if the target is right to the root, check its left child
            elseif dis_median > 0
                kdSearch(kdBuildTree, target, minDis, leftChild);
            end
        % if not intersect
        else
            kdBuildTree(1).visit = 1;
            return;
        end
    end
end

% if current node is not visit or is not a leaf, find the closed point
% along the tree
if (~kdBuildTree(nodeNumber).visit) 
    while (~kdBuildTree(nodeNumber).leaf) && (~kdBuildTree(nodeNumber).visit)
        leftChild = kdBuildTree(nodeNumber).left;
        rightChild = kdBuildTree(nodeNumber).right;
        % find the closest node by deep first search
        if ~kdBuildTree(leftChild).visit && target(kdBuildTree(nodeNumber).split) <= kdBuildTree(nodeNumber).median
            nodeNumber = leftChild;
        elseif ~kdBuildTree(rightChild).visit
            nodeNumber = rightChild;
        else
            break;
        end
    end
    
    % update visit status
    kdBuildTree(nodeNumber).visit = 1;   
    
    % update the nearest point and the minimal distance
     pointFound = kdBuildTree(nodeNumber).data;
     tempDis = getDistance(pointFound, target);   
     [nearestPoint, minDis] = updateDis(tempDis, pointFound, minDis, nearestPoint);
     
     % find next node if the root is not visited    
     parentNumber = kdBuildTree(nodeNumber).parent;
     if parentNumber
         leftChild = kdBuildTree(parentNumber).left;
         rightChild = kdBuildTree(parentNumber).right;
         % if both children been visited, explore grandparent
         if kdBuildTree(leftChild).visit && kdBuildTree(rightChild).visit
                 if kdBuildTree(parentNumber).parent == 1
                     [~, intersection] = checkIntersect(kdBuildTree, target, 1, minDis);
                     if ~intersection
                     outputPoint = nearestPoint;
                     return;
                     end
                 end
                 kdSearch(kdBuildTree, target, minDis, kdBuildTree(parentNumber).parent);
         else           
            [dis_median, intersection] = checkIntersect(kdBuildTree, target, parentNumber, minDis);
            if intersection 
                % check if parent node is closer to the target
                tempDis = getDistance(kdBuildTree(parentNumber).data, target);   
                [nearestPoint, minDis] = updateDis(tempDis, kdBuildTree(parentNumber).data, minDis, nearestPoint);
             
                % update the parent as visited
                kdBuildTree(parentNumber).visit = 1;
                
                % if the target is left to the parent, check its right child
                if dis_median <= 0 
                    kdSearch(kdBuildTree, target, minDis, rightChild);
                % if the target is right to the parent, check its left child
                elseif dis_median > 0 
                    kdSearch(kdBuildTree, target, minDis, leftChild);
                end
            % if not intersect
            else
                kdBuildTree(parentNumber).visit = 1;
                if kdBuildTree(parentNumber).parent == 1
                    [~, intersection] = checkIntersect(kdBuildTree, target, 1, minDis);
                    if ~intersection
                        outputPoint = nearestPoint;
                        return;
                    end
                end
                kdSearch(kdBuildTree, target, minDis, kdBuildTree(parentNumber).parent);
            end
         end
        % if currentNode is root
     else
         kdSearch(kdBuildTree, target, minDis, 0);
     end
end

% Output
if nargin == 2
    outputPoint = nearestPoint;
    clear global kdBuildTree;
    clear global nearestPoint;
end
end
%%  This fuction is used to calcuate Euclidean distance between two point
function distance = getDistance(source, target)
temp = (source - target).^2;
distance = sqrt(sum(temp));
end

%% This function is used to update closest point and minimal distance
function [nearestPoint, minDis] = updateDis(tempDis, pointFound, minDis, nearestPoint)
     if tempDis < minDis
         minDis = tempDis;
         nearestPoint = pointFound;
     end
end

%% This function is used to check whether intersection exists between:
%  the circle whose center is target and whose radius is minDis
%  the spliting hyperplane of currentNode
function [dis_median, intersection] = checkIntersect(kdBuildTree, target, nodeNumber, minDis)
split = kdBuildTree(nodeNumber).split;
median = kdBuildTree(nodeNumber).median;
dis_median = target(split)-median;
if abs(dis_median) <= minDis
    intersection = 1;
else
    intersection = 0;
end
end


