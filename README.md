# plotarBCSTwoBandModel

BCS plotagem utilizando do model descrito no artigo 'A minimal two-band model for the superconducting Fe-pnictides'

Utiliza-se o Hamiltoniano como descrito no artigo 'A minimal two-band model for the superconducting 
Fe-pnictides'. Para testar-lo, apenas necessita os pacotes necessarios instalados, e inserir a faixa de valores que quer trabalhar. Mas dependendo do que quer, é necessario que defina as intenções de plotagem.

Caso queira mudar a precisao em calcularDelta_MI scipy.optimize.fsolve tem um parametro que por default no caso desse programa
calcula no máximo 300(Grande contribuidor das descontinuidades no gráfico), se não chegar na precisao desejada, outro paramentro
que é por default 1.49012e-8, retorna onde parou. A figura é guardar por default no mesmo local do codigo, e já é automaticamente
numerada, utilizando os parametros utilizados para definir o nome do arquivo.

Para analisar tempo temos uma classe Tempo que funciona assim: t.inicio() guarda o tempo, quando terminar t.fim(nome) guardar
o novo tempo e associa esses dois a um nome. Quando quiser saber quanto tempo passou só usar t.obterTempo(nome). Os methods de
fimCalculo e fimPlot tem interessações diferentes em relação a guardar quantidades, sendo o ultimo usado para guardar nomes.

A previsao de tempo é normalmente errado em 3-8 vezes quando os numeros de calculos são muito grandes, mas em casos pequenos sao bem
precisos. O segundo tempo é so o primeiro multiplicado por um fator de correção. Ele funciona calculando uma vez um fila de valores,
por exemplo uma serie de N, para um dado J. E ai multiplica pelo numero de J, e tambem o numero de particoes da temperatura. Não é o 
jeito mais preciso(talvez uma melhoria seria um que atualiza apos cada calculo e ai talvez mostre 'tempo previsto ate termino')

Exemplos de codigo (certas opcoes teriam que adicionar na funcao diretamente):

    for T in Temperaturas:
        listaMi, listaDelta = calcularDelta_MI([1, 1.1, 1.2], J, s0, t1, t2, t3, t4, T, Divisoes) 

#2D para n = 1:
    
    plot2D(listaJ, listaMi, 'J', 'μ', 'Grafico relacionando J e μ em T = ' + str(T) + 'K', False, 0, [0,2], True) 
    # O ultimo true é pq alem de salvar, ira mostrar o plot na tela

#3D:

    plot3D('wireframe', N, J, Mi, 'n', 'J', 'μ', 'μ em T = ' + str(T) + 'K')
    plot3D('curvas', N, J, Delta, 'n', 'J', 'Δ', 'Δ em T = ' + str(T) + 'K', 50)
    plotDensityPlot(N, J, Delta, 'n', 'J', 'Δ')
    plotDensityPlot(N, J, Mi, 'n', 'J', 'μ')
