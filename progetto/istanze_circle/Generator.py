import random as r
from random   import shuffle
import sys
import json
import math
import numpy as np
import json

#parte grafica
import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.contour import ContourSet
import Tkinter as Tk

class View():
    def __init__(self,root,dim_x,dim_y):
        f = Figure()
        nticks = 10
        ax = f.add_subplot(111, aspect='1')
    	ax.set_xticks([x*(dim_x)/nticks for x in range(nticks+1)])
    	ax.set_yticks([y*(dim_y)/nticks for y in range(1,nticks+1)])
        ax.set_xlim((0, dim_x))
        ax.set_ylim((0, dim_y))
        canvas = FigureCanvasTkAgg(f, master=root)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        toolbar = NavigationToolbar2TkAgg(canvas, root)
        toolbar.update()
        self.f = f
        self.ax = ax
        self.canvas = canvas
        
    def draw(self,vec_x,vec_y):
        self.ax.plot(vec_x, vec_y, "%so" % 'w', scalex=0.0, scaley=0.0)


class Generator:

	def __init__(self, dim_x, dim_y):
		self.dim_x = dim_x
		self.dim_y = dim_y

	#generazione punti random (solo quadrante +/+)
	def random(self,points):
		vec_x = [round(r.random()*self.dim_x,2) for i in range(points)]
		vec_y = [round(r.random()*self.dim_y,2) for i in range(points)]
		return vec_x,vec_y
		
	#generazione su una matrice di punti random senza ripetizioni
	def random_matrix(self,points):
		vec_x = r.sample([i for i in range(self.dim_x)],points)
		vex_y = r.sample([i for i in range(self.dim_y)],points)

	#distribuzione in cluster (pseudo quadrati)
	#points: numero di punti totali; clusters: numero di cluster da creare;
	#cohesion: compreso tra 0 ed 1, descrive la coesione interna dei cluster
	#1-> tutti i punti del cluster sul centroide: 0-> distribuzione random
	def cluster(self,points,clusters,cohesion):
		elements = [(points/clusters)+1 if i<points%clusters else points/clusters for i in range  (clusters)]
		cohesion = min(max(cohesion,0),1)
		vec_x, vec_y = [], []
		for i in range(clusters):
			#popolo il cluster
			x,y = self.random(elements[i])
			#calcolo il centroide
			c_x = sum(x)/len(x)
			c_y = sum(y)/len(y)
			#scalo il cluster secondo la coesione interna
			# TIP: r.random() serve per concentrarne di piu sul centro e meno sull'esterno
			x = [c_x + (x[j]-c_x)*(1-cohesion)*r.random() for j in range(len(x))]
			y = [c_y + (y[j]-c_y)*(1-cohesion)*r.random() for j in range(len(y))]
			#lo sposto da qualche parte per non averli tutti piu o meno centrati
			m_x,M_x = min(x),max(x)
			m_y,M_y = min(y),max(y)
			s_x = r.random()*(m_x+self.dim_x-M_x) - m_x
			s_y = r.random()*(m_y+self.dim_y-M_y) - m_y
			x = [x[i]+s_x for i in range(len(x))]
			y = [y[i]+s_y for i in range(len(y))]
			#unisco tutti i punti, non serve averli separati
			vec_x = vec_x + x
			vec_y = vec_y + y
		return vec_x, vec_y

	#sequenze di spezzate con punti equidistanti
	#dist_min indica la distanza minima desiderata tra punti di una riga
	#puo esser modificata se vengono inseriti troppi punti per riga
	#dist_min =-1 vuol dire che fitto io le distanze
	#offset True se posso spostarli dall'inizio
	def sequence(self,points,rows,dist_min,offset):
	    elements = [(points/rows)+1 if i<points%rows else points/rows for i in range(rows)]
	    if dist_min==-1:
	        dist_min = self.dim_y*1.0 / elements[0]
	    vec_x,vec_y=[],[]
	    o = 0
	    for i in range(rows):
	        c_y = r.random()*self.dim_y
	        vec_y += [c_y for j in range(elements[i])]
	        x = [k*dist_min for k in range(elements[i])]
	        o = 0
	        if (offset):
	            o = (self.dim_x-max(x)) * r.random()
	            x = [xx+o for xx in x]
	        vec_x += x
	    return vec_x,vec_y
		
	#distribuzione dei punti su un numero ri righe orizzontali pari a rows
	def lines(self,points,rows):
		elements = [(points/rows)+1 if i<points%rows else points/rows for i in range(rows)]
		vec_x, vec_y = [], []
		for i in range(rows):
		    c_y = r.random() * self.dim_y
		    y = [c_y for j in range(elements[i])]
		    x = [r.random()*self.dim_x for j in range(elements[i])]
		    vec_x = vec_x + x
		    vec_y = vec_y + y
		return vec_x,vec_y
		
	#evoluzione del metodo precedente, stampa linee non parallele
	#coeff indica l'inclinazioni massime e minime della retta
	def lines2(self, points,rows, coeff_max, coeff_min):
	    elements = [(points/rows)+1 if i<points%rows else points/rows for i in range(rows)]
	    vec_x, vec_y = [], []
	    for i in range(rows):
		    #coefficiente angolare
		    coeff = r.random()*(coeff_max - coeff_min) + coeff_min
		    b = 0
		    x = [r.random()*self.dim_x for j in range(elements[i])]
		    if coeff>0:
		        bmax = self.dim_y-max([xx*coeff for xx in x])
		        b = r.random()*bmax
		    else:
		        bmin = min([xx*coeff for xx in x])
		        b = self.dim_y+r.random()*bmin
		    x = [r.random()*self.dim_x for j in range(elements[i])]
		    y = [xx*coeff + b for xx in x]
		    #print [yy for yy in y if yy>100 or yy<0]
		    if len([yy for yy in y if yy>100 or yy<0])>0:
		        #print 'questa soluzione fa schifo'
		        return self.lines2(points,rows,coeff_max,coeff_min)
		    vec_x += x
		    vec_y += y
	    #print vec_y
	    return vec_x,vec_y
	
	#dispone i punti su cerchi concentrici, cicles indica il numero di circonferenze
	#p = True indica che le circonferenze non sono concentriche
	# alla dimensione di essa
	def circles(self, points, circles, p):
	    elements = [(points/circles)+1 if i<points%circles else points/circles for i in range(circles)]
	    vec_x,vec_y = [],[]
	    for c in range(circles):
	        radius = r.random()*(self.dim_y/2)
	        a,b = 0,0
	        if p:
	            a = r.random()*(self.dim_x-radius*2) - (self.dim_x/2) + radius
	            b = r.random()*(self.dim_y-radius*2) - (self.dim_y/2) + radius
	        center = (self.dim_x/2 + a,self.dim_y/2 + b)
	        
	        l1,l2 = 0,elements[c]/2
	        #x = [center[0]-radius+r.random()*radius*2 for j in range(elements[c])]
	        #y = [(center[1] + abs(radius**2 - (x[i]-center[0])**2)**0.5) if i%2==1 else center[1] - abs(radius**2 - (x[i]-center[0])**2)**0.5 for i in range(elements[c])]
	        x = [center[0]-radius+r.random()*radius*2 for j in range(l2)]
	        y = [(center[1] + abs(radius**2 - (x[i]-center[0])**2)**0.5) if i%2==1 else center[1] - abs(radius**2 - (x[i]-center[0])**2)**0.5 for i in range(l2)]
	        y += [center[1]-radius+r.random()*radius*2 for j in range(l2,elements[c])]
	        x += [(center[0] + abs(radius**2 - (y[i]-center[1])**2)**0.5) if i%2==1 else center[0] - abs(radius**2 - (y[i]-center[1])**2)**0.5 for i in range(l2,elements[c])]
	        vec_x += x
	        vec_y += y
	    return vec_x,vec_y
	    
	#coeff e' un vettore di coefficienti della funzione da costruire
	# [5,1,3] -> 5 + 1*x + 3*x^2
	# error e' un bound sull'errore di ogni punto
	def function(self,points,coeff,error):
	    x = [r.random()*self.dim_x for i in range(points)]
	    #y = [  sum([xx**i * coeff[i]*1.0 for i in range(len(coeff))])      for xx in x]
	    y = [  sum([xx**i * coeff[i]*1.0 + (r.random()*error*2-error) for i in range(len(coeff))])      for xx in x]
	    return x,y
	    
	    
#true se i punti sono consistenti
def checkers(x,y):
    #if max(x)>100 or max(y)>100 or min(x)<0 or min(x)<0:
    if max(x+y)>100 or min(x+y)<0:
        print 'dati inconsistenti'
        return False
    print 'dati ok'
    return True
    
def get_matrix(x,y):
    m = [[round((abs(x[i]-x[j])**2 + abs(y[i]-y[j])**2)**0.5,4) for i in range(len(x))] for j in range(len(x))]
    return m

def main(argv):
    dim_x,dim_y = 100,100
    n_test = 4
    n_points = 100
    
    g = Generator(dim_x,dim_y)
    #root = Tk.Tk()
    #root.wm_title("MEMOC")
    #view = View(root,dim_x,dim_y)
    
    #vec_x,vec_y = g.random(200)
    #vec_x,vec_y = g.cluster(200,5,0.8)
    #vec_x,vec_y = g.cluster(200,3,0.6)
    # NO vec_x,vec_y = g.lines(200,5)
    #non funziona se metto coefficiente angolare 0
    # NO vec_x,vec_y = g.lines2(150,5,1,-1)
    # NO vec_x,vec_y = g.sequence(200,5,2,False)
    # NO vec_x,vec_y = g.sequence(50,5,-1,False)
    #vec_x,vec_y = g.sequence(200,5,2,True)
    #vec_x,vec_y = g.circles(100,2,False)
    #vec_x,vec_y = g.circles(200,3,True)
    
    # NO vec_x,vec_y = g.function(200,[60,0.8,-0.03,0.0002],2)
    vec_x,vec_y = g.function(200,[20,-0.5,0.04,-0.0003],3)
    checkers(vec_x,vec_y)
    #view.draw(vec_x,vec_y)
    
    m = get_matrix(vec_x,vec_y)
    
    data = []
    data.append([n_test,n_points])
    '''for i in range(n_test):
        vec_x,vec_y = g.random(n_points)
        #vec_x,vec_y = g.cluster(n_points,5,0.8)
        #vec_x,vec_y = g.cluster(n_points,3,0.6)
        #vec_x,vec_y = g.sequence(n_points,5,2,True)
        #vec_x,vec_y = g.circles(n_points,2,False)
        #vec_x,vec_y = g.circles(n_points,3,True)
        #data.append({'x':vec_x,'y':vec_y,'m':get_matrix(vec_x,vec_y)})
        data.append(get_matrix(vec_x,vec_y))
    '''
    with open(str(n_points)+'_circle_dataset','w') as outfile:
        json.dump([n_test,n_points],outfile)
        outfile.write('\n')
        for i in range(n_test):
            #vec_x,vec_y = g.random(n_points)
            
            #vec_x,vec_y = g.cluster(n_points,3,0.7)
			
            vec_x,vec_y = g.circles(n_points,3,False)
			
            json.dump(get_matrix(vec_x,vec_y),outfile)
            outfile.write('\n')
    
    
    #Tk.mainloop()

if __name__ == "__main__":
    main(sys.argv)


#i parametri globali di generator sono solo la dimensione che e' fissa,
#il numero di punti no perche' cosi richiamo facilmente diverse istanze.

'''
metodi che funzionano:
1) random
2) cluster coesi (5)
3) cluster sparsi (3)
4) rette orizzontali
5) rette random
6) tre tipi di sequenza
7) cerchi concentrici
8) cerchi random
9) funzioni fucking-pro

metodi che adotto:
1) random
2) 3 cluster sparsi (cohesion: 0.6)
3) 5 cluster coesi  (cohesion: 0.8)
4) sequenza da 5 spezzate (dist = 2, adcazzum = True)
5) 3 cerchi non concentrici
6) funzione fucking pro 20-0.5x+0.04x^2-0.0003x^3, err=3 

'''
