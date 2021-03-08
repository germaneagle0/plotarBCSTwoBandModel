"""
Criado em 21 de fevereiro de 2020

O programa primeiramente pede para inserir a faixa de valores, na função obterConstantesIniciais.
Utiliza-se o Hamiltoniano como descrito no artigo 'A minimal two-band model for the superconducting 
Fe-pnictides'. É necessario que defina as intenções de plotagem.

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

Exemplos de codigo (certas opcoes teriam que adicionar na funcao diretamente)

for T in Temperaturas:
    listaMi, listaDelta = calcularDelta_MI([1, 1.1, 1.2], J, s0, t1, t2, t3, t4, T, Divisoes) 

#2D para n = 1:
    
    plot2D(listaJ, listaMi, 'J', 'μ', 'Grafico relacionando J e μ em T = ' + str(T) + 'K', False, 0, [0,2], True) # O ultimo true é pq alem de salvar, ira mostrar o plot na tela

#3D
    plot3D('wireframe', N, J, Mi, 'n', 'J', 'μ', 'μ em T = ' + str(T) + 'K')
    plot3D('curvas', N, J, Delta, 'n', 'J', 'Δ', 'Δ em T = ' + str(T) + 'K', 50)
    plotDensityPlot(N, J, Delta, 'n', 'J', 'Δ')
    plotDensityPlot(N, J, Mi, 'n', 'J', 'μ')

@author: Thiago Oliveira Jucá Vasconcelos
"""


from time import time

# Classe para guardar tempo de modo organizado
class Tempo:

    def __init__(self, digitos = 3):
        self.inicial = [time()]
        self.termino = []
        self.frase = []
        self.plot = {}
        self.digitos = digitos
        self.calculo = 0


    def inicio(self):
        self.inicial.append(time())


    def fim(self, frase, boolean = True):
        self.termino.append(time())
        self.frase.append(frase)

        if (boolean):
            self.printTempo(frase)

    def fimPlot(self, boolean = False, name = ''):
        try:
            self.plot[name] += 1
        except:
            self.plot[name] = 0

        self.fim('Plot ' + str(self.plot[name]), boolean)

    def fimCalculo(self, boolean = False):
        self.calculo += 1
        self.fim(str(self.calculo) + 'o Calculo', boolean)

    def obterTempo(self, frase):
        try:
            index = self.frase.index(frase)
            return self.termino[index] - self.inicial[index]
        except:
            print(f'{frase} não encontrado em {self.frase}')


    def printTempo(self, frase, optional = ''):
        tempo = 0
        
        try:
            if (frase == 'Total'):

                for palavra in self.frase:
                    tempo += self.obterTempo(palavra)
                
            else:
                tempo += self.obterTempo(frase)
            
            self.printFormat(frase, tempo, optional)
        
        except:
            print(f'{frase} não encontrado em {self.frase} e não é \'Total\'')

    def printFormat(self, frase, tempo, optional):
        print(f'{frase}:', round(tempo, self.digitos), '\u0008s', optional)

# Obter tempo para ser representado no console
tempo = Tempo()

# Pacotes para resolução de integrais
import scipy.optimize
from scipy.optimize import fsolve
from scipy.integrate import tplquad, quad, dblquad

# Para obter constante e funcoes matematicas
import math
import cmath
from math import sqrt

# Para plotar os graficos
import matplotlib.pyplot as plt
import numpy

# Para garantir que o usuário escape quando precionar CTRL C no questionário
import sys

# Remover alertas enquanto roda o programa
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

# Lidar com certos tipo de erro
    #numpy.seterr(divide='ignore', invalid='ignore')

#from mpl_toolkits import mplot3d

def main():

    # Constantes que definem a forma/estrutura cristalina
    t1 = -1 
    t2 = 1.3
    t3 = -0.85
    t4 = -0.85

    # Definir estimativas iniciais para as raizes das respectivas equacoes (eq1 e eq2) 
    s0 = numpy.array([1.5, 0.05])

    # Para determinar precisao do arrendondamento no print no terminal do tempo
    tempo.digitos = 3

    # Guardar tempo que o programa está começando inicializando, mas não imprimir ainda
    tempo.fim('Inicialização', False)

    # Questionar o usuario sobre os dados iniciais
    Numero, EnergiaJ, Divisoes, Temperaturas = obterConstantesIniciaisDoGrafico()

    # Obter uma lista de todos valores de temperaturas a ser analisado
    Temperaturas = numpy.linspace(*Temperaturas)

    # Obter as listas de todos N e J que seram analisados
    N, J = numpy.meshgrid(numpy.linspace(*Numero), numpy.linspace(*EnergiaJ))

    # Imprimir frase
    tempo.printTempo('Inicialização')

    # Imprimir previsão de tempo de calculo
    preverTempoDeCalculo(Numero, EnergiaJ, s0, t1, t2, t3, t4, Temperaturas, Divisoes)

    # Lista contendo todos calculos feitos, caso queira calcular primeiro, depois mostrar, talvez queira usar sua transposta(com np.array(Deltas) pode usar Deltas.transpose())

    Deltas = []
    Mis = []

    ## Calculo e Plotagem##
    print('Calculando...\n')

    # Plot 3D/HeatMap
    for T in Temperaturas:
        Mi, Delta = calcularDelta_MI(N, J, s0, t1, t2, t3, t4, T, Divisoes)

        Mis.append(Mi)
        Deltas.append(Delta)

        plot3D('superficie', N, J, Delta, 'n', 'J', 'Δ', f'Δ (T = {str(round(T, tempo.digitos))}K)', 0, True, True)

        plotDensityPlot(N, J, Delta, 'n', 'J', f'Δ (T = {str(round(T, tempo.digitos))}K)')
        plotDensityPlot(N, J, Mi, 'n', 'J', f'μ (T = {str(round(T, tempo.digitos))}K)')

    ## Plotagem pós calculo caso queira

    # Calcula toda tempo em que o programa esta rodando, ou seja exclui qualquer tempo de interação com usuario
    tempo.printTempo('Total')


def plotDensityPlot(listaX, listaY, gridZ, nomeX, nomeY, nomeZ, boolean = False, show = False):
    tempo.inicio()

    plt.pcolormesh(listaY, listaX, gridZ, cmap=plt.get_cmap('coolwarm'), alpha = 1, shading='gouraud')
    
    plt.xlabel(nomeY)
    plt.ylabel(nomeX)

    clb = plt.colorbar()
    clb.set_label(nomeZ)

    name = f'density-{nomeZ[0]}x{nomeY}x{nomeX}'
    tempo.fimPlot(boolean, name)
    plt.savefig(name + '-' + str(tempo.plot[name]) + '.png')
    if show:
        plt.show()
    plt.close()


def plot3D(opcoes, listaX, listaY, gridZ, nomeX = 'x', nomeY = 'y', nomeZ = 'z', titulo = '', niveis = 0, boolean = False,  show = False):
    tempo.inicio()

    # Definir Eixos e Aparencia 3D
    ax = plt.axes(projection='3d')

    if opcoes == 'superficie':
        ax = plt.axes(projection='3d')
        ax.view_init(20,-140)
        ax.plot_surface(listaX, listaY, gridZ, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    elif opcoes == 'curvas':
        ax.contour3D(listaX, listaY, gridZ, niveis, cmap='binary')
        ax.view_init(20,-140)
    elif opcoes == 'wireframe':
        ax.plot_wireframe(listaX, listaY, gridZ, color='black')
        ax.view_init(20,-140)
    
    ax.set_title(titulo)
    ax.set_xlabel(nomeX)
    ax.set_ylabel(nomeY)
    ax.set_zlabel(nomeZ)

    name = f'{opcoes}-{nomeZ}x{nomeY}x{nomeX}'

    tempo.fimPlot(boolean, name)
    plt.savefig(name + '-' + str(tempo.plot[name]) + '.png')
    if show:
        plt.show()

def plot2D(listaX, listaY, nomeX = 'x', nomeY = 'y', titulo = '', boolean = False, limX = 0, limY = 0,  show = False):
    tempo.inicio()

    # Definir como será apresentado
    plt.plot( listaX,listaY, 'go') # green bolinha
    plt.plot( listaX,listaY, 'k:', color='blue') # linha pontilha orange, pode botar os dois ao mesmo tempo sem problema

    plt.xlabel(nomeX)
    plt.ylabel(nomeY)
    plt.title(titulo)

    # Definir o intervalo a ser inicializado
    if limX != 0:
        plt.xlim(*limX)
    if limY != 0:
        plt.ylim(*limY)

    name = f'2D-{nomeY}x{nomeX}'

    tempo.fimPlot(boolean, name)
    plt.savefig(name + '-' + str(tempo.plot[name]) + '.png')
    
    if show:
        plt.show()
    # Voce pode ficar montando uma figura com varias curvas, ai teria que tirar esse plt.close()(plt.show tb fecha a figura)
    plt.close()



# Definimos a funçao responsavel por resolver o problema, retorna o grid mi, delta, e os arrays de N e J
def calcularDelta_MI(N, J, s0, t1, t2, t3, t4, T, p_l, boolean = False):
    tempo.inicio()

    # Função que define e resolve, os sistemas além de preencher a lista
    @numpy.vectorize
    def solving(n,j):

        # Definindo o sistema de equações
        def sistema(variaveis):
                (u, delta) = variaveis
                equacoes = numgapT1s_nTs(t1,t2,t3,t4,u,delta,T,p_l)
                eq1 = numpy.sum(equacoes['gap']) - 1 / j
                eq2 = numpy.sum(equacoes['nts']) - n
                return [eq1, eq2]

        # Resolver sistema
        s = scipy.optimize.fsolve(sistema, s0) # Se quiser que desaparece aquelas aberações no gráfico, tem como aumentar o numero de interacoes, diminuir a precisao necessaria
        s = [abs(s[0]),abs(s[1])] # Talvez seja util modificar s0(bota ele fora de def main(), e bota 'global s0' dentro dessa def solving(), ai modificar ele)
        # Guardar resultado em listas
        lista_mi.append(s[0])

        # para mudar para mi so troca return s[1] para s[0]
        return s[1] 
    

    # Initializando lista mi, a lista delta sera automaticamente gerada na função solving
    lista_mi=[]
    grid_Delta = solving(N, J)

    # Solving function de tal modo que o primeiro valor é repetido duas vezes
    lista_mi.pop(0)

    # Transformar array 1D em um array de array (2D)
    grid_Mi = numpy.reshape(lista_mi, (len(N), -1))

    tempo.fimCalculo(boolean)

    return [numpy.array(grid_Mi), numpy.array(grid_Delta)]


def preverTempoDeCalculo(N, J, s0, t1, t2, t3, t4, T, p_l):

    erro = 1.5
    temp =  tempo.obterTempo(tempo.frase[-1]) * len(T)
    
    n,j = numpy.meshgrid(numpy.linspace(1, 1, 1), numpy.linspace(1, 1, 1))

    if len(N) > len(J):
        calcularDelta_MI(n, J, s0, t1, t2, t3, t4, T[0], p_l)
        print(f'\n=>ESTIMATIVA INFERIOR PARA O TEMPO DE CALCULO TOTAL: {round(temp * N[2], tempo.digitos)}\u0008s a {round(temp * N[2] * erro, tempo.digitos)}\u0008s\n')
    else:
        calcularDelta_MI(N, j, s0, t1, t2, t3, t4, T[0], p_l)
        print(f'\n=>ESTIMATIVA INFERIOR PARA O TEMPO DE CALCULO TOTAL: {round(temp * J[2], tempo.digitos)}\u0008s a {round(temp *  J[2] * erro, tempo.digitos)}\u0008s\n')


# Funcao que retorna os dois grids, modificados por funções, em um dicionario, vão definir as equações
def numgapT1s_nTs(t1, t2, t3, t4, u, delta, T, p_l):

    func1 = numpy.zeros((p_l,p_l))
    func2 = numpy.zeros((p_l,p_l))

    kx = numpy.linspace(-math.pi,math.pi,p_l)
    ky = numpy.linspace(-math.pi,math.pi,p_l)

    kx, ky = numpy.meshgrid(kx, ky)
    
    k = obterConstantes(kx, ky, t1, t2, t3, t4, u, delta, T)

    func1 = gapT1s(k['raizB'], k['g1'], k['w1'], k['tgw1'], k['w2'], k['tgw2']) / (p_l ** 2)
    func2 = nTs(k['se'], k['raizB'], k['g'], k['w1'], k['tgw1'], k['w2'], k['tgw2']) / (p_l ** 2)
    
    return ({
        'gap' : func1,
        'nts' : func2
    })



def gapT1s(raiz_b, g_1, w_1, tgw1, w_2, tgw2):
    return (1 / (4*raiz_b)) * ( ((g_1/w_1) + w_1) * tgw1 - (((g_1 / w_2) + w_2) * tgw2) )


def nTs(se, raiz_b, g, w_1, tgw1, w_2, tgw2):
    return 1 - (se / (4*raiz_b)) * ( ((g / w_1) + w_1) * tgw1 - ((g / w_2) + w_2) * tgw2 )


# Feita para obter todas constantes necessarias, como estas dependem entre si, foi posta todas em uma só função para economia de espaco e velocidade.
def obterConstantes(kx, ky, t1, t2, t3, t4, u, delta, T):
    
    cosx = numpy.cos(kx)
    cosy = numpy.cos(ky)

    cosProduct = cosx * cosy
    sinProduct = numpy.sin(kx) * numpy.sin(ky)

    e_x = -2 * (t1*cosx + t2*cosy + 2*t3*cosProduct)
    e_y = -2 * (t2*cosx + t1*cosy + 2*t3*cosProduct)

    e_xy = -4 * t4 * sinProduct
    e_xyQUAD = e_xy ** 2

    Ex = e_x - u
    Ey = e_y - u

    ExQUAD = Ex ** 2
    EyQUAD = Ey ** 2

    pe = Ex * Ey

    se = Ex + Ey
    seQUAD = ExQUAD + EyQUAD + 2 * pe

    deltaQUAD = delta ** 2

    g = e_xyQUAD - pe - deltaQUAD
    g_1 = -(deltaQUAD + EyQUAD + e_xyQUAD)

    a = ((ExQUAD+EyQUAD)/2) + e_xyQUAD + deltaQUAD
    b = (((ExQUAD-EyQUAD)/2) ** 2) + (e_xyQUAD)*(seQUAD)

    raiz_b = numpy.sqrt(numpy.abs(b))

    w_1 = numpy.sqrt(a + raiz_b)
    w_2 = numpy.sqrt(a - raiz_b)
    
    tgw1 = numpy.tanh(w_1/(2*T))
    tgw2 = numpy.tanh(w_2/(2*T))

    return {
        'raizB': raiz_b,
        'g' : g,
        'g1' : g_1,
        'se' : se,
        'w1' : w_1,
        'w2' : w_2,
        'tgw1' : tgw1,
        'tgw2' : tgw2,
    }


def perguntar(minimo = -1, fator = -1):

    while (True):
        try:
            inicial = float(input('\nValor Inicial: '))
            final = float(input('Valor Final: '))

            if (fator != -1):
                CONSTANTE = round(fator * (final - inicial))
                CONSTANTE = CONSTANTE if (CONSTANTE > minimo) else minimo
                divisoes = int(input(f'Partições entre os valores (Ex: {CONSTANTE}): '))
            else:
                divisoes = int(input(f'Partições entre os valores: '))
                
            ok = ''
            
            while (len(ok) == 0):
                ok = input('Está correto (s/n)? ')
            
            if (ok[0].lower() == 's'):
                break

        except KeyboardInterrupt:
            print('\n\nTerminando...')
            sys.exit()
        except:
            print('\nValor Inválido, tente novamente...')

    return [ inicial, final, divisoes ]


def obterConstantesIniciaisDoGrafico():

    print('\n(PRECIONE CTRL + C PARA TERMINAR O PROGRAMA A QUALQUER MOMENTO)\n\n========= ESPECIFICAÇÃO DO NÚMERO =========')
    N = perguntar(10, 27.5)

    print('\n========= ESPECIFICAÇÃO DO J (ENERGIA ATRATIVA) =========')
    J = perguntar(15, 10)

    print('\n========= ESPECIFICAÇÃO DA TEMPERATURA =========')
    T = perguntar()

    print('\n========= ESPECIFICAÇÃO DA PRECISÃO DE CADA INTEGRAL =========\n')
    while (True):
        try:
            p_l = int(input('Partições na Integral é (Ex: 100): '))
            break
        except KeyboardInterrupt:
            print('\n\nTerminando...')
            sys.exit()
        except:
            print('\nValor Inválido, tente novamente...')

    # Adicionar um \n
    print()
    return [ N, J, p_l, T ]
    

main()
