format long%rat;


%% example 1 and 2
MatrixA = [
    1 0 -1 0 -1
    -1 2 0 -1 -1
    0 -1 3 -1 0
    -1 0 -1 4 -1
    0 -1 -1 0 5];

q_line = [-1 2 -1 2 1]';
q_hat = [-1 1 -1 0 1]';


Solution_MatrixA_q_line = block_pivoting_method(MatrixA, q_line);  %solution of example 1
Error =  norm(min(MatrixA * Solution_MatrixA_q_line + q_line, Solution_MatrixA_q_line)); % example 1

Solution_MatrixA_q_hat = block_pivoting_method(MatrixA, q_hat);  %solution of example 2
Error =   norm(min(MatrixA * Solution_MatrixA_q_hat + q_hat, Solution_MatrixA_q_hat)); % example 2
%hand_solution = [16/7 17/14 29/28 25/28 1/4]';


%% example tridiagonal M-matrix
Dim = 30;

Tridiagonal_Matrix=BlockPivotMatrix(Dim)+eye(Dim^2);
inv(Tridiagonal_Matrix);

Dim=Dim^2;
q_line = ones(Dim, 1);
Location_n_1 = 1 : 2  : Dim;
q_line(Location_n_1) = -1;
q_line;

tic
Solution_Tridiagonal_q_line = block_pivoting_method(Tridiagonal_Matrix, q_line);
toc

Error = norm(min(Tridiagonal_Matrix * Solution_Tridiagonal_q_line + q_line, Solution_Tridiagonal_q_line))

q_hat = -1 * ones(Dim, 1);
Location_n_1 = 2 : 3 : Dim;
q_hat(Location_n_1) = 1;
Location_n_0 = 3 : 3 : Dim;
q_hat(Location_n_0) = 0;

q_hat;


tic
Solution_Tridiagonal_q_hat = block_pivoting_method(Tridiagonal_Matrix, q_hat);
toc

Error = norm(min(Tridiagonal_Matrix * Solution_Tridiagonal_q_hat + q_hat, Solution_Tridiagonal_q_hat))

%% example tridiagonal M-matrix ,example 3
Dim=1000;
Tridiagonal_Matrix = diag(ones(1, Dim) * 2) + diag(ones(1, Dim - 1) * -1 , -1 ) + diag(ones(1, Dim - 1)*-1, 1);
q_hat=randn(Dim,1);
tic
Solution_Tridiagonal_q_hat = block_pivoting_method(Tridiagonal_Matrix, q_hat);
toc
Error = norm(min(Tridiagonal_Matrix * Solution_Tridiagonal_q_hat + q_hat, Solution_Tridiagonal_q_hat))