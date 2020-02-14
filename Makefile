# description file for pared_multiple, viene de pared_ramiro pero calcula varias TDP

OBJS =  pared_multiple.o integra.o tyopas.o util.o \
	opacidades_cfa.o mievo.o

paredm2: $(OBJS)
	g77 -O2 -o  paredce $(OBJS) 

pared_multiple.o: pared_multiple.f
	g77 -c pared_multiple.f	

integra.o: integra.f
	g77 -c integra.f

tyopas.o: tyopas.f
	g77 -c tyopas.f

util.o: util.f
	g77 -c util.f

opacidades_cfa.o: opacidades_cfa.f
	g77 -c opacidades_cfa.f

mievo.o: mievo.f 
	g77 -c mievo.f



