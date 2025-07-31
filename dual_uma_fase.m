clc, clear all

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

syms M
% Variaveis basicas
% xb = [3 4];
xb = [3 4];
% Variaveis nao basicas
% xn = [1 2];
xn = [1 2];
% Variaveis artificiais
fic = [];
% Matriz de coeficientes A
% A = [-2 -3 1 0; 1 -2 0 1];
A = [-1 -1 1 0; 3 1 0 1];
% Coeficientes das basicas na funcao objetivo z
% cb = [0 0];
cb = [0 0];
% Coeficientes das nao basicas na funcao objetivo z
% cn = [2 6];
cn = [-2 3];
% Coeficientes das basicas nas restricoes
B = A(:,xb);
% Coeficientes das nao basicas nas restricoes
N = A(:,xn);
% Vetor de resultados b
% b = [-12; -4];
b = [-3; 6];
% Se precisar usar o metodo de restricao artificial, use RA = 1
RA = 0;
if RA == 1
    % Mais uma variavel basica
    xb = [xb size(A,2)+1];
    % Coficiente 1 para xN e 0 para xB na restricao artificial
    A = [A; zeros(1,size(A,2))];
    A(end,xn) = 1;
    % Adicionando a presenca de x0 nas restricoes
    A(:,end+1) = [zeros(size(A,1)-1,1); 1];
    % Adicionando M em b
    b = [b; M];
    % Atualizar cb, B e N
    cb = [cb 0];
    B = A(:,xb);
    N = A(:,xn);
end

% Operacoes para construcao da tabela
op{1} = cb * inv(B) * N - cn;
op{2} = inv(B) * N;
op{3} = inv(B) * b;
op{4} = cb * inv(B) * b;
for i = 1:length(op)
    % Elimina imprecisao
    op{i}(subs(abs(op{i}),M,9999) < 1e-6) = 0;
end
fprintf('cb * inv(B) * N - cn =\n');
disp(op{1});
fprintf('inv(B) * N =\n');
disp(op{2});
fprintf('inv(B) * b =\n');
disp(op{3});
fprintf('cb * inv(B) * b =\n');
disp(op{4});

% Gera a tabela inicial
tab = geratab1fase(op, A, xb, xn, b);
disp('-------------------------');
xb
xn
tab

while 1
    % Encontra o pivot
    [p,xn,xb] = pivotdual(tab,xn,xb,fic);
    if all(p == 0)
        % Retornou pivot 0, fim do metodo
        break
    end
    p
    disp('-------------------------');
    xb
    xn
    
    % Pivot vira 1
    tab(p(1),:) = tab(p(1),:)/tab(p(1),p(2));
    for i = 1:size(tab,1)
        if i ~= p(1)
            % Zerando demais valores na coluna do pivot
            tab(i,:) = tab(p(1),:)*(-tab(i,p(2))) + tab(i,:);
        end
    end
    tab
end

function tab = geratab1fase(op, A, xb, xn, b)
    % Gera a tabela inicial
    %    z  xB        xN            RHS
    % z  1  0   cB B^-1 N - cN   cB B^-1 b
    % xB 0  I       B^-1 N         B^-1 b
    
    % Inicializa com zero
    tab = sym(zeros(size(A,1)+1, size(A,2)+2));
    
    % Preenche a linha de z
    tab(1,1) = 1;
    tab(1,[xn+1]) = op{1};
    tab(1,end) = op{4};
    
    % Preenche as linhas de xB
    tab(2:size(A,1)+1,[xn+1]) = op{2};
    tab(2:size(A,1)+1,[xb+1]) = eye(length(xb));
    tab(2:end,end) = op{3};
end

function [p,xn,xb] = pivotdual(tab,xn,xb,fic)
    %    z  xB        xN            RHS
    % z  1  0   cB B^-1 N - cN   cB B^-1 b
    % xB 0  I       B^-1 N         B^-1 b
    
    syms M
    % Se existir algum M, substituir por um valor suficientemente grande
    tabsub = subs(tab,M,9999);
    if all(tabsub(2:end,end) >= 0)
        % O metodo acaba quando os valores da coluna de RHS sao >= 0
        fprintf('end of method\n');
        if any(tabsub(1,xn+1) > 1e-6)
            fprintf(['dual problem is infeasible\n'...
                'primal may be unbounded or infeasible\n']);
        end
        p = [0 0];
        xn = xn;
        xb = xb;
        return
    else
        % Indice do menor valor na coluna de RHS
        idx = find(tabsub(2:end,end) == min(tabsub(2:end,end)));
        if length(idx) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(idx)
                if any(ismember(fic,xb(idx(i))) == 1)
                    % Se houver uma variavel artificial candidata a sair da
                    % base, sera selecionada
                    idx = idx(i);
                    break
                end
            end
            if length(idx) > 1
                % Existe mais de um valor igual e nao existem variaveis
                % artificiais para sair da base, entao escolhe qualquer um
                idx = idx(1);
            end
        end
        p(1) = idx + 1;
        
        % Variaveis nao basicas com yrj < 0
        idx = xn(tabsub(p(1),xn+1) < 0);
        if isempty(idx)
            % Checa se nao existem valores negativos para yrj
            fprintf('Problem may be not limited\n');
            p = [0 0];
            xn = xn;
            xb = xb;
            return
        end
        % Encontra o minimo valor de (zj - cj) / yrj
        idx = idx(find((tabsub(1,idx+1)./tabsub(p(1),idx+1)) == ...
            min(tabsub(1,idx+1)./tabsub(p(1),idx+1))));
        if length(idx) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(idx)
                if all(ismember(fic,idx(i)) == 0)
                    % Checa se preferivelmente existe alguma variavel que
                    % nao seja artificial para entrar na base
                    idx = idx(i);
                    break
                end
            end
            if length(idx) > 1
                % Aparentemente somente variaveis artificiais estao
                % disponiveis para entrar na base
                idx = idx(1);
            end
        end
        p(2) = idx + 1;
        
        % Atualizando xB e xN
        iaux = xn(xn==p(2)-1);
        xn(xn==p(2)-1) = xb(p(1)-1);
        xb(p(1)-1) = iaux;
    end
end