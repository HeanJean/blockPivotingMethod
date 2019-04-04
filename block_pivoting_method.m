function [SolutionX] = block_pivoting_method(MatrixA, VectorQ)
% This is the matlab script for "Ximing Fang, Zhijun Qiao, The block principle pivoting
% algorithm for the linear complementarity problem with an M-matrix."
% x^T (Ax + q) =  0, x >= 0, Ax + q >=0                  (1)
% input varies:  MatrixA : M-matrix
%                 VectorQ
% output varies: SolutionX : solution of equation (1)


%% Pseudocode 
% Input：MatrixA, VectorQ
% output：SolutionX



%% Matlab Code
LenP = length(VectorQ)
%pause

[FlagAllNegative, FlagAllNotNegative, NegativeLocation, NotNegativeLocation] = isAllNegative(VectorQ);


if FlagAllNegative
    SolutionX = (MatrixA) \ (- VectorQ);
    return;
end

if FlagAllNotNegative
    SolutionX = zeros(length(VectorQ), 1);
    return;
end  



[MatrixAnn, MatrixApn, MatrixApp, MatrixAnp] = extractElements(MatrixA, NegativeLocation, NotNegativeLocation);
VectorQn = VectorQ(NegativeLocation);
VectorQp = VectorQ(NotNegativeLocation);

%MatrixB1 = - inv(MatrixAnn) * MatrixAnp;
MatrixB1 = - MatrixAnn\MatrixAnp;
%VectorQ1 = - inv(MatrixAnn) * VectorQn;
VectorQ1 = - MatrixAnn\VectorQn;

%MatrixA1 = MatrixApp - MatrixApn * inv(MatrixAnn) * MatrixAnp;
MatrixA1 = MatrixApp - MatrixApn * (MatrixAnn \ MatrixAnp);
%VectorQ2 = - MatrixApn * inv(MatrixAnn) * VectorQn + VectorQp;
VectorQ2 = - MatrixApn * (MatrixAnn\VectorQn) + VectorQp;

[~, FlagAllNotNegative, ~, ~] = isAllNegative(VectorQ2);
if FlagAllNotNegative
    SolutionX = assignSolution(VectorQ1, NegativeLocation, length(VectorQ));
    return;
end  


SolutionPartXp = block_pivoting_method(MatrixA1, VectorQ2);
SolutionPartXn = MatrixB1 * SolutionPartXp + VectorQ1;


SolutionX = assignSolution2(SolutionPartXn, SolutionPartXp, NegativeLocation, LenP); 
return;
 end



function output = assignSolution2(Vector1, Vector2, Location, Len)
output = zeros(Len,1);
for ForI = 1 : length(Location)
    output(Location(ForI)) = Vector1(ForI);
end

Location2 = setdiff(1 : Len, Location);
for ForI = 1 : length(Location2)
    output(Location2(ForI)) = Vector2(ForI);
end


end



function output = assignSolution(Vector, Location, Len)
output = zeros(Len, 1);
for ForI = 1 : length(Location)
    output(Location(ForI)) = Vector(ForI);
end
end


function [FlagAllNegative, FlagAllNotNegative, NegativeLocation, NotNegativeLocation] = isAllNegative(Vector)
FlagAllNegative = 0;
FlagAllNotNegative = 0;

NegativeLocation = find(Vector < 0);

if length(NegativeLocation) == length(Vector)
    FlagAllNegative = 1;
end

if isempty(NegativeLocation)
    FlagAllNotNegative = 1;
end

NotNegativeLocation = setdiff(1 : length(Vector), NegativeLocation);
end

function [MatrixA, MatrixB, MatrixC, MatrixD] =  extractElements(InputMatrix, LocationA, LocationB)

LenA = length(LocationA);
MatrixA = zeros(LenA, LenA);


for ForI = 1 : LenA
    for ForJ = 1 : LenA
        MatrixA(ForI, ForJ) = InputMatrix(LocationA(ForI), LocationA(ForJ));
    end
end

LenB = length(LocationB);
MatrixB = zeros(LenB, LenA);
for ForI = 1 : LenB
    for ForJ = 1 : LenA
        MatrixB(ForI, ForJ) = InputMatrix(LocationB(ForI), LocationA(ForJ));
    end
end

MatrixC = zeros(LenB, LenB);
for ForI = 1 : LenB
    for ForJ = 1 : LenB
        MatrixC(ForI, ForJ) = InputMatrix(LocationB(ForI), LocationB(ForJ));
    end
end


MatrixD = zeros(LenA, LenB);
for ForI = 1 : LenA
    for ForJ = 1 : LenB
        MatrixD(ForI, ForJ) = InputMatrix(LocationA(ForI), LocationB(ForJ));
    end
end

end