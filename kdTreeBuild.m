function kdBuildTree = kdTreeBuild(dataMatrix, parentIndex)
%
% This function is used to build a kd-tree structure for d dimensional data
%
% INPUT:
% data:  n x d matrix, while n is the number of points of the given
%          pointcloud, d is the dimension of the points
% parentIndex: do not need be given, will asign during iteration, parent
%                    position in kd_tree
%
% OUTPUT:
% count:     the total number of nodes
% kd_tree:  An array of structs, which represents a node in the tree.
%
%               The structure contains the following information:
%
%               left:           the position of left child node
%               right:         the position of right child node
%               leftVector:  the nodes vector of left child tree
%               rightVector: the nodes vector of right child tree
%               parent:       the positoin of parent node
%               spilt:          the spliting dimension
%               median:     the point vector, median of the data along
%                                current spliting dimension
%               data:          the current node vector, the first point
%                                with median data along spliting dimension
%               leaf:           "1" if it is a leaf, "0" if it is a node 
%               visit:          "1" if it is visited, "0" if it is not
%               
%               
% Fangtong Liu
% University of Michigan
% ftliu@umich.edu
% 2020/02/05

%% Initialize parameters
global kd_tree;
global count;

if nargin == 1
    % dimension of the data
    [n, ~] = size(dataMatrix);
    parentIndex = 0;
    count = 1;
    % initialize node
    node.left = [];
    node.right = [];
%     node.leftVector = [];
%     node.rightVector = [];
    node.parent = 0;
    node.data = [];
    node.split = 0;
    node.median = 0;
    node.leaf = 0;
    node.visit = 0;
    % initialize kd_tree
    kd_temp(1:n) = node;
    kd_tree = kd_temp;
else
    count = count + 1;
    [n, ~] = size(dataMatrix);
%     kd_tree(count).root = 0;
end

numberCurrent = count;
kd_tree(numberCurrent).parent = parentIndex;
%% If there is no point in dataMatrix
if n==0
    kd_tree(numberCurrent).data = [];
    kd_tree(numberCurrent).leaf = 1;
    kd_tree(numberCurrent).left = [];
    kd_tree(numberCurrent).right = [];
    kd_tree(numberCurrent).visit = 1;
    kdBuildTree = numberCurrent;
    return;
end
%% If only one point in dataMatrix
if n == 1
    kd_tree(numberCurrent).data = dataMatrix(1, :);
    kd_tree(numberCurrent).leaf = 1;
    kd_tree(numberCurrent).left = [];
    kd_tree(numberCurrent).right = [];
    kd_tree(numberCurrent).visit = 0;
%     kd_tree(numberCurrent).leftVector = [];
%     kd_tree(numberCurrent).rightVector = [];
    kdBuildTree = numberCurrent;
    return;
end
kd_tree(numberCurrent).visit = 0;
%% If there are two or more points left in dataMatrix
% find the spliting dimesion 
 variance = var(dataMatrix);
 [~, splitValue] = max(variance);
 kd_tree(numberCurrent).split = splitValue;
 % find the median of the desired dimension
 % if there are two points, get the smaller one as node, the bigger one as
 % right child leaf, if there are more than two points, get the middle one
 if n == 2
     if dataMatrix(1, splitValue) < dataMatrix(2, splitValue)
         medianIndex = 1;
         otherIndex = 2;
     else
         medianIndex = 2;
         otherIndex = 1;
     end

     kd_tree(numberCurrent).data = dataMatrix(medianIndex, :);
     kd_tree(numberCurrent).median = dataMatrix(medianIndex, splitValue);
%      kd_tree(numberCurrent).rightVector = (dataMatrix(otherIndex, :));
     kd_tree(numberCurrent).left = kdTreeBuild([], numberCurrent);
     kd_tree(numberCurrent).right = kdTreeBuild(dataMatrix(otherIndex, :), numberCurrent);
     
 else
      medianPosition = round(n/2);
      medianArray = sort(dataMatrix(:, splitValue));
      medianValue = medianArray(medianPosition);
      kd_tree(numberCurrent).median = medianValue;
      % assign node.data as the first point in dataMatrix with median value
      % along spliting point
      index_temp = find(dataMatrix(:, splitValue) == medianValue);
      index_node = index_temp(1);
      kd_tree(numberCurrent).data = dataMatrix(index_node, :);
      % find the index of points smaller than medianValue on spliting dimension
      index_left = find(dataMatrix(:, splitValue) < medianValue);
      % include the remaining points with medianValue
      if length(index_temp) > 1
          index_left = [index_left; index_temp(2:end)];
      end
      clear index_temp;
      % find the index of points smaller than medianValue on spliting dimension
      index_right = find(dataMatrix(:, splitValue) > medianValue);
%       kd_tree(numberCurrent).leftVector = dataMatrix(index_left, :); 
%       kd_tree(numberCurrent).rightVector = dataMatrix(index_right, :);
      kd_tree(numberCurrent).left = kdTreeBuild(dataMatrix(index_left, :), numberCurrent);
      kd_tree(numberCurrent).right = kdTreeBuild(dataMatrix(index_right, :), numberCurrent);         
 end

 %% Output 
 if nargin == 1
     kdBuildTree = kd_tree;
     clear global kd_tree;
     clear global count;
 else     
     kdBuildTree = numberCurrent;
end
end
