# Metodos-Dual-Simplex-Matlab
Este repositório contém uma coleção de scripts em MATLAB que implementam variações do método Dual Simplex para resolver problemas de programação linear (PL), códigos auxiliares para converter o modelo primal para o dual e também para 8 casos de análise de sensibilidade, quando são realizadas pequenas alterações no problema.

---

## Implementações do Método Dual Simplex em MATLAB
Esta seção expande o repositório https://github.com/Lucca-Rodrigues-main/Metodos-Simplex-Matlab com ferramentas para análise de dualidade e sensibilidade em problemas de programação linear. Os códigos permitem a conversão de um problema primal para seu dual, a aplicação do algoritmo Dual Simplex e a realização de análises pós-otimalidade.

<img width="512" height="266" alt="image" src="https://github.com/user-attachments/assets/edf3a094-8b96-4301-86f8-bbb1fa36a2c3" />

Os seguintes scripts e métodos foram adicionados:
* **Conversão Primal-Dual**
* **Método Dual Simplex de Uma Fase**
* **Método Dual Simplex de Duas Fases**
* **Análise de Sensibilidade (Pós-otimalidade)**

---

## Conversão Primal-Dual
Em programação linear, todo problema (chamado de **primal**) tem um problema correspondente associado a ele, chamado de **dual**. A relação entre eles é fundamental, pois a solução de um revela informações sobre a solução do outro (Teorema da Dualidade Forte) e as variáveis de um se relacionam com as restrições do outro (Teorema das Folgas Complementares). A conversão segue um conjunto de regras bem definidas para transformar as variáveis, restrições e a função objetivo do primal no formato dual.

### Implementação: `converte_dual.m`
Este script é uma ferramenta poderosa que automatiza a conversão de um problema de programação linear primal para sua forma dual. O usuário fornece a matriz de coeficientes (`A`), os vetores de custo (`c`) e de recursos (`b`), e as naturezas das restrições e variáveis. O código então:
- Constrói e exibe a formulação matemática completa tanto do problema primal quanto do seu correspondente dual.
- Oferece uma funcionalidade avançada que, a partir de uma solução ótima do problema dual, utiliza as relações de folgas complementares para encontrar a solução ótima do problema primal original.

---

## O Método Dual Simplex de Uma Fase
O algoritmo Dual Simplex é uma variação do método Simplex que é particularmente útil quando o quadro (tableau) Simplex é **primal infactível** (pelo menos um valor no lado direito das restrições, o vetor `b`, é negativo) mas **dual factível** (os custos reduzidos na linha da função objetivo já satisfazem a condição de otimalidade).

Em vez de buscar uma variável com custo reduzido favorável para entrar na base, o Dual Simplex primeiro seleciona uma variável básica com valor negativo para sair da base. Em seguida, aplica um critério de razão sobre a linha da função objetivo para determinar qual variável não básica entrará na base, restaurando a factibilidade primal passo a passo, enquanto mantém a factibilidade dual.

### Implementação: `dual_uma_fase.m`
Este script implementa o algoritmo Dual Simplex de forma tabular. Ele inicia com uma solução que é primal infactível e itera através de pivoteamentos, conforme as regras do Dual Simplex, até que a factibilidade primal seja alcançada (todos os elementos do RHS sejam não negativos), resultando na solução ótima.

#### Exemplo Resolvido no Código
O script está configurado para resolver o seguinte problema de minimização, que possui uma solução básica inicial infactível ($x_3 = -3$).

$$
\begin{aligned}
\text{Minimizar } Z = -2x_1 &+ 3x_2 \\
\text{Sujeito a:} \\
-x_1 - x_2 + x_3 &= -3 \\
3x_1 + x_2 + x_4 &= 6 \\
x_i &\ge 0, \quad \forall i \in \{1, ..., 4\}
\end{aligned}
$$

---

## O Método Dual Simplex de Duas Fases
Assim como no Simplex primal, pode haver casos em que a tabela inicial não é nem primal factível, nem dual factível. Para esses cenários, uma abordagem de duas fases para o Dual Simplex pode ser aplicada. A Fase 1 tem como objetivo encontrar uma solução básica que seja dual factível, utilizando uma função objetivo auxiliar. Uma vez que a factibilidade dual é alcançada, a Fase 2 prossegue com o algoritmo Dual Simplex padrão para encontrar a solução ótima e factível do problema.

### Implementação: `dual_duas_fases.m`
Este script implementa a variação de duas fases do método Dual Simplex. Ele utiliza uma variável artificial e uma função objetivo para a Fase 1 para manobrar a tabela até um estado dual factível. Em seguida, na Fase 2, ele otimiza a função objetivo original usando as regras do Dual Simplex para resolver o problema.

#### Exemplo Resolvido no Código
O problema exemplo no código é um caso que requer a abordagem de duas fases, contendo tanto uma restrição de igualdade (que necessita de uma variável artificial, $x_5$) quanto uma restrição que leva a uma infactibilidade no RHS.

$$
\begin{aligned}
\text{Minimizar } Z = -2x_1 &+ 3x_2 \\
\text{Sujeito a:} \\
-x_1 - x_2 + x_3 &= -3 \\
3x_1 + x_2 + x_4 &= 6 \\
x_1 + x_2 &= 3 \\
x_i &\ge 0, \quad \forall i \in \{1, ..., 4\}
\end{aligned}
$$

---

## Análise de Sensibilidade (Pós-otimalidade)
A análise de sensibilidade é o estudo de como as mudanças nos parâmetros de um problema de programação linear afetam sua solução ótima. Uma vez que a solução ótima é encontrada, em vez de resolver o problema do zero para cada nova alteração, podemos usar a tabela ótima final para analisar o impacto dessas mudanças de forma muito mais eficiente.

### Implementação: `casos.m`
O script `casos.m` é uma ferramenta de análise de sensibilidade que parte de uma solução ótima de um PPL e permite ao usuário investigar o efeito de diferentes tipos de alterações. Após a modificação, o script avalia se a solução ainda é ótima ou factível e, caso não seja, aplica o método Simplex (Primal ou Dual) para encontrar a nova solução ótima.

O código contempla os seguintes 8 casos de análise de sensibilidade:
1.  **Introdução de uma nova variável no problema:** Analisa o impacto da adição de uma nova atividade ou produto ao modelo.
2.  **Adicionar uma restrição de desigualdade:** Verifica como uma nova restrição do tipo $\le$ ou $\ge$ afeta a solução ótima existente.
3.  **Adicionar uma restrição de igualdade:** Verifica o impacto de uma nova restrição de igualdade, que é mais restritiva.
4.  **Variação do coeficiente de custo de uma variável não básica:** Determina a faixa de valores que o custo de uma variável que está fora da solução pode assumir sem alterar a solução ótima atual.
5.  **Variação do coeficiente de custo de uma variável básica:** Analisa como a mudança no custo de uma variável que já está na solução afeta a otimalidade e o valor da função objetivo.
6.  **Variação de um elemento do vetor b:** Estuda o efeito da alteração na disponibilidade de um recurso (lado direito de uma restrição).
7.  **Variação do coeficiente $\alpha_{ij}$ da coluna de uma variável não básica:** Analisa o impacto da mudança em um coeficiente tecnológico de uma variável não básica.
8.  **Variação do coeficiente $\alpha_{ij}$ da coluna de uma variável básica:** O caso mais complexo, que analisa a mudança em um coeficiente tecnológico de uma variável que já está na solução, podendo afetar a matriz base e a solução como um todo.
