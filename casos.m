clc, clear all

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

% Exemplo 1
% Variaveis basicas
xb = [1 2];
% Variaveis nao basicas
xn = [3 4];
% Variaveis artificiais
fic = [];
% Matriz de coeficientes A
A = [1 1 1 0; -1 2 0 1];
% Coeficientes das basicas na funcao objetivo z
cb = [-1 3];
% Coeficientes das nao basicas na funcao objetivo z
cn = [0 0];
% Coeficientes das basicas nas restricoes
B = A(:,xb);
% Coeficientes das nao basicas nas restricoes
N = A(:,xn);
% Vetor de resultados b
b = [6; 8];
% Tabela
tab = [1 0 0 -5/3 -2/3 -46/3;...
    0 1 0 2/3 -1/3 4/3;...
    0 0 1 1/3 1/3 14/3];
% Solucao
s = [4/3 14/3 0 0];

% Exemplo 2
% % Variaveis basicas
% xb = [1 4];
% % Variaveis nao basicas
% xn = [2 3];
% % Variaveis artificiais
% fic = [];
% % Matriz de coeficientes A
% A = [1 1 1 0; -1 2 0 1];
% % Coeficientes das basicas na funcao objetivo z
% cb = [-1 0];
% % Coeficientes das nao basicas na funcao objetivo z
% cn = [2 0];
% % Coeficientes das basicas nas restricoes
% B = A(:,xb);
% % Coeficientes das nao basicas nas restricoes
% N = A(:,xn);
% % Vetor de resultados b
% b = [6; 8];
% % Tabela
% tab = [1 0 -3 -1 0 -6;...
%     0 1 1 1 0 6;...
%     0 0 3 1 1 14];
% % Solucao
% s = [6 0 0 14];

% Se o primal simplex falhar, use dual = 1 que usaremos o dual simplex
dual = 0;

% 1. Introducao de uma nova variavel no problema
% 2. Adicionar uma restricao de desigualdade
% 3. Adicionar uma restricao de igualdade
% 4. Variacao do coeficiente de custo de uma variavel nao basica
% 5. Variacao do coeficiente de custo de uma variavel basica
% 6. Variacao de um elemento do vetor b
% 7. Variacao do coeficiente aij da coluna de uma variavel nao basica
% 8. Variacao do coeficiente aij da coluna de uma variavel basica
caso = 1;
switch caso
    case 1
        % MODIFIQUE cnew e anew
        cnew = -1; % Coeficiente da nova variavel
        anew = [-1; 2]; % Coeficientes da nova variavel na matriz A
        
        % Expandindo a tabela
        tab(:,end+1) = tab(:,end);
        % Adicionando coluna da nova variavel
        tab(1,end-1) = tab(1,[xb(cb==0) xn(cn==0)]+1) * anew - cnew;
        tab(2:end,end-1) = tab(2:end,[xb(cb==0) xn(cn==0)]+1) * anew;
        % Adicionar nova variavel nao basica
        xn = [xn size(tab,2)-2];
    case 2
        % MODIFIQUE anew e bnew
        % x1 - x2 <= 4
        % x1 + 2*x2 <= 9
        asnew = [1 2 0 0]; % Coeficientes da nova restricao na matriz A
        bnew = 9; % Resultado da restricao
        
        % MODIFIQUE a operacao da restricao (>=, <=)
        if ~(sum(s .* asnew) <= bnew)
            % Expandindo a tabela
            tab(:,end+1) = tab(:,end);
            tab(end+1,:) = zeros(1,size(tab,2));
            % Adicionar nova variavel basica
            xb = [xb size(tab,2)-2];
            % Adicionando nova linha
            tab(end,2:length(asnew)+1) = asnew;
            % Adicionando coluna da variavel de folga
            tab(:,end-1) = [zeros(length(xb),1); 1];
            % Adicionando valor em RHS
            tab(end,end) = bnew;
            disp('-------------------------');
            disp('Quadro pre-preparado');
            tab
            % Formatando no formato factivel
            for j = xb(1:end-1)
                if tab(end,j+1) > 0
                    tab(end,:) = tab(end,:) + (tab(j+1,:) * -tab(end,j+1));
                end
            end
        else
            error('>> A restricao nao elimina a solucao otima! <<');
        end
    case 3
        % MODIFIQUE anew e bnew
        % 2*x1 - x2 = 4
        asnew = [2 -1 0 0]; % Coeficientes da nova restricao na matriz A
        bnew = 4; % Resultado da restricao
        
        if ~(sum(s .* asnew) == bnew)
            % Expandindo a tabela
            for i = size(tab,2):-1:2
                tab(:,i+1) = tab(:,i);
            end
            tab(:,end+1) = tab(:,end);
            for i = size(tab,1):-1:2
                tab(i+1,:) = tab(i,:);
            end
            tab(end+1,:) = tab(end,:);
            % Adicionando colunas de x0 e da variavel artificial
            tab(:,2) = zeros(size(tab,1),1);
            tab(:,end-1) = zeros(size(tab,1),1);
            % Adicionando linhas de x0 e da variavel artificial
            tab(2,:) = [0 1 zeros(1,size(tab,2)-4) -1 0];
            tab(end,3:length(asnew)+4) = [asnew 1 bnew];
            disp('-------------------------');
            disp('Quadro pre-preparado');
            tab
            % Adicionar nova variavel basica
            xb = [xb size(tab,2)-3];
            % Formatando no formato factivel
            for j = xb(1:end-1)
                if tab(end,j+2) > 0
                    tab(end,:) = tab(end,:) + (tab(j+2,:) * -tab(end,j+2));
                end
            end
            tab(2,:) = tab(2,:) + tab(end,:);
        else
            error('>> A restricao nao elimina a solucao otima! <<');
        end
    case 4
        % MODIFIQUE xvar, cold e cnew
        % c2 = 2 para c2' = 1
        % c2 = 2 para c2' = -2
        xvar = 2; % Indice do x nao basico
        cold = 2; % Antigo coeficiente na funcao objetivo
        cnew = 1; % Novo coeficiente na funcao objetivo
        
        if cb * A(:,xvar) - cnew < 0
            error('>> O quadro continua otimo! <<');
        else
            % Atualiza o custo de x
            tab(1,xvar+1) = cb * A(:,xvar) - cnew;
        end
    case 5
        % MODIFIQUE xvar, cold e cnew
        % c2 = -3 para c2' = -2
        % c2 = -3 para c2' = 2
        xvar = 2; % Indice do x basico
        cold = -3; % Antigo coeficiente na funcao objetivo
        cnew = 2; % Novo coeficiente na funcao objetivo
        
        diffc = cnew - cold;
        for i = xn
            % Calcula todos cj barra
            cbs(i) = tab(1,i+1) + diffc*tab(xvar+1,i+1);
        end
        cbs(cbs < 1e-6) = 0;
        cbs = cbs(xn)
        
        % Calculando o z barra
        cb(xb==xvar) = cnew;
        zb = cb * s(xb).'
        
        if all(cbs <= 0)
            error('>> O quadro continua otimo! <<');
        else
            % Atualiza a tabela
            tab(1,xn+1) = cbs;
            tab(1,end) = zb;
        end
    case 6
        % MODIFIQUE nr, bold e bnew
        % b2 = 8 para b2' = 6
        % b2 = 8 para b2' = 15
        nr = 2; % Numero da restricao
        bold = 8; % Antigo valor de bi
        bnew = 15; % Novo valor de bi
        
        % Atualizando b
        b(nr) = bnew;
        
        % Calcula o b barra
        bb = tab(2:end,xn+1) * b
        
        % Calculando o z barra
        zb = cb * bb
        
        if all(bb > 0)
            error('>> O quadro continua otimo! <<');
        else
            % Atualiza a tabela
            tab(2:end,end) = bb;
            tab(1,end) = zb;
        end
    case 7
        % MODIFIQUE a e lc
        % a22 = 2 para a22 = 4
        % a2 = [1; 2] para a2 = [-3; 5]
        a = [-3; 5]; % Novo valor do coeficiente
        lc = {1:2; 2}; % Linha e coluna do coeficiente na matriz A
        
        al = A(:,lc{2});
        al(lc{1}) = a;
        % Calcula z linha
        if ~isempty(cb * al - cb(xb==lc{2}))
            zl = cb * al - cb(xb==lc{2});
        else
            zl = cb * al - cn(xn==lc{2});
        end
        
        if zl < 0
            error('>> O quadro continua otimo! <<');
        else
            % Calcula y'
            yl = tab(2:end,xn+1) * al;
            
            % Atualiza a tabela
            tab(1,lc{2}+1) = zl;
            tab(2:end,lc{2}+1) = yl;
        end
    case 8
        % MODIFIQUE a e lc
        % a21 = -1 para a21 = 1
        a = 1; % Novo valor do coeficiente
        lc = {2; 1}; % Linha e coluna do coeficiente na matriz A
        
        al = A(:,lc{2});
        al(lc{1}) = a;
        % Calcula c barra
        if ~isempty(cb * al - cb(xb==lc{2}))
            cba = tab(1,xn+1) * al - cb(xb==lc{2});
        else
            cba = tab(1,xn+1) * al - cn(xn==lc{2});
        end
        
        if all(cba <= 0)
            error('>> O quadro continua otimo! <<');
        else
            % Calcula y'
            yl = tab(2:end,xn+1) * al;
            
           % Expandindo a tabela
            tab(:,end+1) = tab(:,end);
            % Adicionando coluna da nova variavel
            tab(:,end-1) = [cba; yl];
        end
    otherwise
        fprintf('Esse caso nao existe')
end

disp('-------------------------');
xb
xn
tab
if caso ~= 3
    while 1
        % Encontra o pivot
        if dual == 1
            [p,xn,xb] = pivotdual(tab,xn,xb,fic);
        else
            [p,xn,xb] = pivot(tab,xn,xb,fic);
        end
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
else
    % Fase 1
    while 1
        % Encontra o pivot
        [p,xn,xb] = pivotfase1(tab,xn,xb,fic);
        if all(p == 0)
            % Retornou pivot 0, fim da fase 1
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
    % Fase 2
    % Elimina linha e coluna de x0
    tab(2,:) = [];
    tab(:,2) = [];
    % Variaveis artificiais que nao estao na base
    ficaux = fic;
    ficaux(find(ismember(xb,fic))) = [];
    % Elimina demais variaveis artificiais da tabela
    tab(:,ficaux+1) = [];
    xn(find(ismember(xn, ficaux))) = [];
    disp('-------------------------');
    tab
    while 1
        % Encontra o pivot
        [p,xn,xb] = pivotfase2(tab,xn,xb,fic);
        if all(p == 0)
            % Retornou pivot 0, fim da fase 2
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
end

function [p,xn,xb] = pivot(tab,xn,xb,fic)
    %    z  xB        xN            RHS
    % z  1  0   cB B^-1 N - cN   cB B^-1 b
    % xB 0  I       B^-1 N         B^-1 b
    
    if all(tab(1,xn+1) <= 1e-6)
        % O metodo acaba quando os valores das variaveis nao basicas na
        % linha de z sao <= 0
        fprintf('end of method\n');
        if any(tab(2:end,end) < 0)
            fprintf('primal problem is infeasible\n');
        end
        p = [0 0];
        xn = xn;
        xb = xb;
        return
    else
        % Maximo valor positivo das nao basicas na linha de z
        p(2) = max(tab(1,xn+1));
        % Armazena o indice
        iaux = find(tab(1,2:end-1)==p(2));
        iaux = iaux+1;
        if length(iaux) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(iaux)
                if ~isempty(find(tab(2:end,iaux(i)) > 0)) && all(ismember(fic,iaux(i)) == 0)
                    % Checa se existe ao menos um valor positivo na coluna
                    % da variavel nao basica selecionada e se
                    % preferivelmente nao eh uma variavel artificial
                    iaux = iaux(i);
                    break
                end
            end
            if length(iaux) > 1
                % Verifica se somente variaveis artificiais estao
                % disponiveis para entrar na base
                for i = 1:length(iaux)
                    if ~isempty(find(tab(2:end,iaux(i)) > 0))
                        % Checa se existe ao menos um valor positivo na
                        % coluna da variavel nao basica selecionada
                        iaux = iaux(i);
                        break
                    end
                end
                if length(iaux) > 1
                    % Nao existem valores positivos na coluna da variavel
                    % nao basica selecionada, entao o problema pode ser
                    % ilimitado
                    fprintf('Problem may be not limited\n');
                    p = [0 0];
                    xn = xn;
                    xb = xb;
                    return
                end
            end
        end
        p(2) = iaux(1);
        % Indices dos valores > 0 na coluna selecionada
        index = find(tab(2:end,p(2)) > 0);
        if isempty(index)
            % Checa se nao existem valores positivos na coluna selecionada
            fprintf('Problem may be not limited\n');
            p = [0 0];
            xn = xn;
            xb = xb;
            return
        end
        index = index + 1;
        % Indices dos menores valores de RHS dividido pelos valores
        % positivos da coluna selecionada
        iaux = find((tab(index,end)./tab(index,p(2))) == ...
            min(tab(index,end)./tab(index,p(2))));
        if length(iaux) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(iaux)
                if any(ismember(fic,iaux(i)) == 1)
                    % Se houver uma variavel artificial candidata a sair da
                    % base, sera selecionada
                    iaux = iaux(i);
                    break
                end
            end
        end
        p(1) = index(iaux(1));
        
        % Atualizando xB e xN
        iaux = xn(xn==p(2)-1);
        xn(xn==p(2)-1) = xb(p(1)-1);
        xb(p(1)-1) = iaux;
    end
end

function [p,xn,xb] = pivotdual(tab,xn,xb,fic)
    %    z  xB        xN            RHS
    % z  1  0   cB B^-1 N - cN   cB B^-1 b
    % xB 0  I       B^-1 N         B^-1 b
    
    if all(tab(2:end,end) >= 0)
        % O metodo acaba quando os valores da coluna de RHS sao >= 0
        fprintf('end of method\n');
        if any(tab(1,xn+1) > 1e-6)
            fprintf(['dual problem is infeasible\n'...
                'primal may be unbounded or infeasible\n']);
        end
        p = [0 0];
        xn = xn;
        xb = xb;
        return
    else
        % Indice do menor valor na coluna de RHS
        idx = find(tab(2:end,end) == min(tab(2:end,end)));
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
        idx = xn(tab(p(1),xn+1) < 0);
        if isempty(idx)
            % Checa se nao existem valores negativos para yrj
            fprintf('Problem may be not limited\n');
            p = [0 0];
            xn = xn;
            xb = xb;
            return
        end
        % Encontra o minimo valor de (zj - cj) / yrj
        idx = idx(find((tab(1,idx+1)./tab(p(1),idx+1)) == ...
            min(tab(1,idx+1)./tab(p(1),idx+1))));
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

function [p,xn,xb] = pivotfase1(tab,xn,xb,fic)
    %    z  x0  xB        xN            RHS
    % z  1  0   0   cB B^-1 N - cN   cB B^-1 b
    % x0 0  1   0  cB' B^-1 N - cN'  cB' B^-1 b
    % xB 0  0   I       B^-1 N         B^-1 b
    
    if tab(2,end) == 0
        % Checa se RHS na linha de x0 atingiu valor 0
        if ~isempty(intersect(xb,fic)) && any(tab(find(ismember(xb,fic))+2,end) ~= 0)
            % Se houver variavel artificial na base com valor diferente de
            % zero, o problema pode nao ter solucao
            fprintf('problem may have no solution\n');
        else
            % Sem variaveis artificiais na base, fim da fase 1
            fprintf('end of phase 1\n');
        end
        p = [0 0];
        xn = xn;
        xb = xb;
        return
    else
        % Maximo valor positivo das nao basicas na linha de x0
        p(2) = max(tab(2,xn+2));
        % Armazena o indice
        iaux = find(tab(2,3:end-1)==p(2));
        iaux = iaux+2;
        if length(iaux) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(iaux)
                if ~isempty(find(tab(3:end,iaux(i)) > 0)) && all(ismember(fic,iaux(i)) == 0)
                    % Checa se existe ao menos um valor positivo na coluna
                    % da variavel nao basica selecionada e se
                    % preferivelmente nao eh uma variavel artificial
                    iaux = iaux(i);
                    break
                end
            end
            if length(iaux) > 1
                % Verifica se somente variaveis artificiais estao
                % disponiveis para entrar na base
                for i = 1:length(iaux)
                    if ~isempty(find(tab(3:end,iaux(i)) > 0))
                        % Checa se existe ao menos um valor positivo na
                        % coluna da variavel nao basica selecionada
                        iaux = iaux(i);
                        break
                    end
                end
                if length(iaux) > 1
                    % Nao existem valores positivos na coluna da variavel
                    % nao basica selecionada, entao o problema pode ser
                    % ilimitado
                    fprintf('Problem may be not limited\n');
                    p = [0 0];
                    xn = xn;
                    xb = xb;
                    return
                end
            end
        end
        if length(iaux) > 1 && ~isempty(intersect(iaux,fic))
            % Double checa se existe mais de um valor igual e se sao
            % variaveis artificiais
            iaux = intersect(iaux,fic);
        end
        p(2) = iaux(1);
        % Indices dos valores > 0 na coluna selecionada
        index = find(tab(3:end,p(2)) > 0);
        if isempty(index)
            % Checa se nao existem valores positivos na coluna selecionada
            fprintf('Problem may be not limited\n');
            p = [0 0];
            xn = xn;
            xb = xb;
            return
        end
        index = index + 2;
        % Indices dos menores valores de RHS dividido pelos valores
        % positivos da coluna selecionada
        iaux = find((tab(index,end)./tab(index,p(2))) == ...
            min(tab(index,end)./tab(index,p(2))));
        if length(iaux) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(iaux)
                if any(ismember(fic,iaux(i)) == 1)
                    % Se houver uma variavel artificial candidata a sair da
                    % base, sera selecionada
                    iaux = iaux(i);
                    break
                end
            end
        end
        p(1) = index(iaux(1));
        
        % Atualizando xB e xN
        iaux = xn(xn==p(2)-2);
        xn(xn==p(2)-2) = xb(p(1)-2);
        xb(p(1)-2) = iaux;
    end
end

function [p,xn,xb] = pivotfase2(tab,xn,xb,fic)
    %    z  xB        xN            RHS
    % z  1  0   cB B^-1 N - cN   cB B^-1 b
    % xB 0  I       B^-1 N         B^-1 b
    
    if all(tab(1,xn+1) <= 1e-6)
        % A fase 2 acaba quando os valores das variaveis nao basicas na
        % linha de z sao <= 0
        fprintf('end of phase 2\n');
        if any(tab(2:end,end) < 0)
            fprintf('primal problem is infeasible\n');
        end
        p = [0 0];
        xn = xn;
        xb = xb;
        return
    else
        % Maximo valor positivo das nao basicas na linha de z
        p(2) = max(tab(1,xn+1));
        % Armazena o indice
        iaux = find(tab(1,2:end-1)==p(2));
        iaux = iaux+1;
        if length(iaux) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(iaux)
                if ~isempty(find(tab(2:end,iaux(i)) > 0)) && all(ismember(fic,iaux(i)) == 0)
                    % Checa se existe ao menos um valor positivo na coluna
                    % da variavel nao basica selecionada e se
                    % preferivelmente nao eh uma variavel artificial
                    iaux = iaux(i);
                    break
                end
            end
            if length(iaux) > 1
                % Verifica se somente variaveis artificiais estao
                % disponiveis para entrar na base
                for i = 1:length(iaux)
                    if ~isempty(find(tab(2:end,iaux(i)) > 0))
                        % Checa se existe ao menos um valor positivo na
                        % coluna da variavel nao basica selecionada
                        iaux = iaux(i);
                        break
                    end
                end
                if length(iaux) > 1
                    % Nao existem valores positivos na coluna da variavel
                    % nao basica selecionada, entao o problema pode ser
                    % ilimitado
                    fprintf('Problem may be not limited\n');
                    p = [0 0];
                    xn = xn;
                    xb = xb;
                    return
                end
            end
        end
        p(2) = iaux(1);
        % Indices dos valores > 0 na coluna selecionada
        index = find(tab(2:end,p(2)) > 0);
        if isempty(index)
            % Checa se nao existem valores positivos na coluna selecionada
            fprintf('Problem may be not limited\n');
            p = [0 0];
            xn = xn;
            xb = xb;
            return
        end
        index = index + 1;
        % Indices dos menores valores de RHS dividido pelos valores
        % positivos da coluna selecionada
        iaux = find((tab(index,end)./tab(index,p(2))) == ...
            min(tab(index,end)./tab(index,p(2))));
        if length(iaux) > 1
            % Checa se existe mais de um valor igual
            for i = 1:length(iaux)
                if any(ismember(fic,iaux(i)) == 1)
                    % Se houver uma variavel artificial candidata a sair da
                    % base, sera selecionada
                    iaux = iaux(i);
                    break
                end
            end
        end
        p(1) = index(iaux(1));
        
        % Atualizando xB e xN
        iaux = xn(xn==p(2)-1);
        xn(xn==p(2)-1) = xb(p(1)-1);
        xb(p(1)-1) = iaux;
    end
end