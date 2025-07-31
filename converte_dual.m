clc, clear all

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

% Matriz de coeficientes A
% A = [3 1; 5 2];
% A = [1 -6 1; 5 7 -2];
A = [1 1 2 1 3; 2 -2 3 1 1];
% Coeficientes na funcao objetivo z
% c = [6 8];
% c = [8 3 -2];
c = [2 3 5 2 3];
% Vetor de resultados b
% b = [4; 7];
% b = [2; -4];
b = [4; 3];
% Operadores (1: =, 2: <=, 3: >=, 4: irrestrito)
% o = [3 3 3 3];
% o = [3 1 2 3 4];
o = [3 3 3 3 3 3 3];
% Minimizacao (1) ou Maximizacao (0)
mima = 1;
% Solucao do dual (deixar vazio se nao quiser converter para primal)
sd = [4/5 3/5];

% Matriz de coeficientes A DUAL
Ad = A.';
% Coeficientes na funcao objetivo z DUAL
cd = b.';
% Vetor de resultados b DUAL
bd = c.';
% Convertendo os operadores
conv_var_rest = mima .* [4 3 2 1] + ~mima .* [4 2 3 1];
conv_rest_var = mima .* [4 2 3 1] + ~mima .* [4 3 2 1];
od = zeros(size(o));
% Numero de restricoes b ou bd
od(1:length(bd)) = conv_var_rest(o(length(b)+1:length(b)+length(c)));
% Numero de variaveis c ou cd
od(length(bd)+1:length(bd)+length(cd)) = conv_rest_var(o(1:length(b)));

% Operadores
os = {'=','<=','>=','irrestrito'};
mm = {'Maximize','Minimize'};
% Variaveis
x = sym('x',[1 length(c)]);
w = sym('w',[1 length(cd)]);
% Funcoes objetivos
z = sum(c.*x);
v = sum(cd.*w);
% Restricoes
for i = 1:length(b)
    P{i} = sum(A(i,:).*x);
end
for i = 1:length(bd)
    D{i} = sum(Ad(i,:).*w);
end

% Mostra o primal
fprintf('%s z(x) = %s\ns.a.\n', mm{mima+1}, char(z));
for i = 1:length(b)
    fprintf('%s %s %d\n', P{i}, os{o(i)}, b(i));
end
for i = 1:length(c)
    if o(i+length(b)) < 4
        fprintf('x%d %s 0\n', i, os{o(i+length(b))});
    else
        fprintf('x%d %s\n', i, os{o(i+length(b))});
    end
end
% Mostra o dual
fprintf('\n%s v(w) = %s\ns.a.\n',  mm{~mima+1}, char(v));
for i = 1:length(bd)
    fprintf('%s %s %d\n', D{i}, os{od(i)}, bd(i));
end
for i = 1:length(cd)
    if od(i+length(bd)) < 4
        fprintf('w%d %s 0\n', i, os{od(i+length(bd))});
    else
        fprintf('w%d %s\n', i, os{od(i+length(bd))});
    end
end

if ~isempty(sd)
    % Substituindo solucao do dual no dual para encontrar as restricoes
    % ativas i.e. as restricoes que atingem a igualdade
    ativaD = find(abs(sum(Ad.*sd,2)-bd) < 1e-6);
    
    % Variaveis nao nulas no primal <==> restricoes ativas no dual
    sistema = A(:,ativaD);
    % Restricoes nao ativas no primal <==> variaveis nulas no dual
    sistema(sd==0,:) = [];
    % Resolve o sistema de equações
    sol = linprog(ones(size(x(ativaD))),-eye(length(x(ativaD))),zeros(length(x(ativaD)),1),sistema,b(sd~=0));
    %sol = fmincon(@(x) sum(x(ativaD)),zeros(size(x(ativaD))),-eye(length(x(ativaD))),zeros(length(x(ativaD)),1),sistema,b(sd~=0));
    
    fprintf('\nRestricoes ativas no primal:');
    find(sd~=0)
    fprintf('Restricoes ativas no dual:');
    ativaD.'
    fprintf('Sistema de equacoes para ser resolvido:');
    sum(sistema.*repmat(x(ativaD),size(sistema,1),1),2)==b(sd~=0)
    fprintf('Variaveis nao nulas: ');
    x(ativaD)==sol(:).'
    fprintf('Variaveis nulas: ');
    x(setdiff(find(x),ativaD))==zeros(1,length(x)-length(ativaD))
    fprintf('Funcao objetivo: z(x) =');
    subs(sum(c.*x), [x(ativaD) x(setdiff(find(x),ativaD))], [sol(:).' zeros(1,length(x)-length(ativaD))])
end