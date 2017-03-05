ARCH64=1
#turn on/off dove
WITH_DOVE = 0

#turn on/off dove
WITH_MPI = 0

#turn on/off gui
WITH_GUI = 1

#CP = g++ -pthread 
CP = g++ 

ifeq ($(WITH_DOVE),1)
EXTRA_CFLAGS = -DWITH_DOVE=0
endif

ifeq ($(WITH_MPI),1)
EXTRA_CFLAGS += -DWITH_MPI=0
CP = mpiCC
endif

ifeq ($(WITH_GUI),1)
EXTRA_CFLAGS += -DWITH_GUI=1
endif

PROJECT = mulphys

#DIRECTORIES
OBJDIR=obj
RUNDIR=run
MAINDIR=main
GUIDIR=main
IODIR=main
MODDIR=model
MESADIR=/usr/lib/x86_64-linux-gnu
#XLIBS = -L/usr/X11/lib -L/usr/X11R6/lib64 -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE
XLIBS = -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE
XMLIB = -lxml2
XMLINC = /usr/include/libxml2
BZIPINC = /usr/include
BZLIB = -lbz2
ASM_SOURCES = asm_386.S
INCDIR = -I$(MESADIR)/include -I/usr/include/GL
LIBDIR = $(MESADIR)/lib
GL_LIBS = -L$(LIBDIR) -lglut -lGLU -lGL -lm $(XLIBS)
LIB_DEP = $(LIBDIR)/$(GL_LIB) $(LIBDIR)/$(GLU_LIB) $(LIBDIR)/$(GLUT_LIB)

#COMPILERS
CPP = $(CP) -c -g $(EXTRA_CFLAGS)
#CPP = g++ -c $(EXTRA_CFLAGS) 
CC = cc -c -g $(EXTRA_CFLAGS)
#CC = cc -c $(EXTRA_CFLAGS)
FC = f77 -c $(FFLAGS) 
CPPM = $(CPP) -I./$(MAINDIR) -I./$(MODDIR) 
LD	= $(CP) -g

# HEADER FILES
MAINH = $(MAINDIR)/main.h $(MAINDIR)/geom.h $(MAINDIR)/domain.h $(MAINDIR)/var.h
GEOMH = $(MAINDIR)/main.h $(MAINDIR)/geom.h $(MAINDIR)/domain.h 
IOH = $(IODIR)/io.h $(IODIR)/output.h
FUNH = $(MAINDIR)/vecalg.h $(MAINDIR)/func.h
UTILH = $(MAINDIR)/vecalg.h $(MAINDIR)/func.h $(MAINDIR)/force.h $(MAINDIR)/tool.h
TOOLH = $(MAINDIR)/tool.h $(MAINDIR)/vecalg.h
SHAPEH = $(MAINDIR)/box.h $(MAINDIR)/cylinder.h

# OBJECT FILES
Okernel = $(OBJDIR)/main.o $(OBJDIR)/geom.o $(OBJDIR)/domain.o \
          $(OBJDIR)/var.o $(OBJDIR)/func.o $(OBJDIR)/diffalg.o
Ogui = $(OBJDIR)/gui.o
Oio = $(OBJDIR)/input.o $(OBJDIR)/output.o $(OBJDIR)/io.o
Oshape = $(OBJDIR)/box.o $(OBJDIR)/cylinder.o
#Oeispack = pythag.o tql1.o tql2.o tred1.o tred2.o #rs.o 
Oparticle = $(OBJDIR)/points.o
Otool = $(OBJDIR)/tool.o  $(OBJDIR)/force.o
Oparal = $(OBJDIR)/dove.o
Orfg = $(OBJDIR)/rfg.o $(OBJDIR)/gauss.o 
#diag.o random.o $(Oeispack)

OBJ = $(Okernel) $(Ogui) $(Oshape) $(Oio) $(Oparticle) $(Otool) $(Oparal) $(Orfg) \
	$(OBJDIR)/models.o $(OBJDIR)/mod_coral.o $(OBJDIR)/mod_crystal.o $(OBJDIR)/mod_gasflow.o \
	$(OBJDIR)/mod_particles.o $(OBJDIR)/mod_particleflow.o $(OBJDIR)/mod_incomflow.o \
	$(OBJDIR)/mod_rfg.o $(OBJDIR)/mod_prfg.o $(OBJDIR)/mod_bonding.o \
	$(OBJDIR)/mod_collapse.o $(OBJDIR)/mod_lattice.o $(OBJDIR)/mod_poisson.o \
	$(OBJDIR)/mod_membrane.o $(OBJDIR)/mod_elastic.o 
#	$(OBJDIR)/mod_pstatic.o

# EXECUTABLE:
PROG = $(RUNDIR)/$(PROJECT)

##### RULES #####

#.c: $(LIB_DEP)
#	$(CC) $(INCDIR) $(CFLAGS) $< $(GL_LIBS) -o $@

all: $(PROG)

#all: $(PROG)batch $(PROG)gui $(PROG)grid

# COMPILING EXECUTABLES:
$(PROG): $(OBJ)
	$(LD) $(OBJ) -o $(PROG) $(GL_LIBS) $(XMLIB) $(BZLIB)

#gui: $(Ogui)
#	$(LD) $(Ogui) -o gui $(GL_LIBS)

# COMPILING OBJECT FILES:

# MAIN OBJECT FILES:

$(OBJDIR)/main.o: $(MAINDIR)/main.cc $(GUIDIR)/gui.h $(MAINH) $(IOH) $(MAINDIR)/input.h
	$(CPP) -I$(XMLINC) $(MAINDIR)/main.cc -o $(OBJDIR)/main.o

$(OBJDIR)/io.o: $(MAINDIR)/io.cc $(MAINDIR)/main.h $(MAINDIR)/io.h $(MAINDIR)/input.h
	$(CPP) -I$(XMLINC) $(MAINDIR)/io.cc -o $(OBJDIR)/io.o

$(OBJDIR)/geom.o: $(MAINDIR)/geom.cc $(MAINH) $(UTIH) $(IODIR)/io.h $(MAINDIR)/tool.h  
	$(CPP) $(MAINDIR)/geom.cc -o $(OBJDIR)/geom.o

$(OBJDIR)/domain.o: $(MAINDIR)/domain.cc $(MAINH) $(TOOLH) $(SHAPEH) $(IODIR)/io.h 
	$(CPP) -I$(BZIPINC) $(MAINDIR)/domain.cc -o $(OBJDIR)/domain.o

$(OBJDIR)/dove.o: $(MAINDIR)/dove.cc $(MAINH) $(TOOLH) $(IODIR)/io.h 
	$(CPP) $(MAINDIR)/dove.cc -o $(OBJDIR)/dove.o

$(OBJDIR)/var.o: $(MAINDIR)/var.cc $(MAINH) $(TOOLH) $(IODIR)/io.h 
	$(CPP) $(MAINDIR)/var.cc -o $(OBJDIR)/var.o

$(OBJDIR)/box.o: $(MAINDIR)/box.cc $(MAINDIR)/box.h $(GEOMH) $(TOOLH) 
	$(CPP) $(MAINDIR)/box.cc -o $(OBJDIR)/box.o

$(OBJDIR)/cylinder.o: $(MAINDIR)/cylinder.cc $(MAINDIR)/cylinder.h $(GEOMH) $(TOOLH) 
	$(CPP) $(MAINDIR)/cylinder.cc -o $(OBJDIR)/cylinder.o

$(OBJDIR)/points.o: $(MAINDIR)/points.cc $(MAINH) $(TOOLH) $(MAINDIR)/force.h
	$(CPP) $(MAINDIR)/points.cc -o $(OBJDIR)/points.o

$(OBJDIR)/force.o: $(MAINDIR)/force.cc $(MAINH) $(TOOLH) $(MAINDIR)/force.h
	$(CPP) $(MAINDIR)/force.cc -o $(OBJDIR)/force.o

$(OBJDIR)/tool.o: $(MAINDIR)/tool.cc $(TOOLH) $(MAINH) 
	$(CPP) $(MAINDIR)/tool.cc -o $(OBJDIR)/tool.o

$(OBJDIR)/input.o: $(MAINDIR)/input.cc $(GEOMH) $(SHAPEH) $(MAINDIR)/tool.h $(IODIR)/eio.h $(MAINDIR)/input.h
	$(CPP) -I$(XMLINC) $(MAINDIR)/input.cc -o $(OBJDIR)/input.o

$(OBJDIR)/output.o: $(MAINDIR)/output.cc $(MAINDIR)/output.h $(GEOMH) $(MAINDIR)/tool.h  $(IODIR)/eio.h
	$(CPP) $(MAINDIR)/output.cc -o $(OBJDIR)/output.o

$(OBJDIR)/gui.o: $(GUIDIR)/gui.cc $(GUIDIR)/gui.h $(GEOMH) $(TOOLH) $(IOH) 
	$(CPP) $(INCDIR) $(GUIDIR)/gui.cc -o $(OBJDIR)/gui.o

$(OBJDIR)/diffalg.o: $(MAINDIR)/diffalg.cc $(IODIR)/io.h $(MAINH) $(TOOLH) 
	$(CPP) $(MAINDIR)/diffalg.cc -o $(OBJDIR)/diffalg.o

$(OBJDIR)/func.o: $(MAINDIR)/func.cc $(MAINDIR)/main.h $(IODIR)/io.h $(MAINDIR)/geom.h $(FUNH) 
	$(CPP) $(MAINDIR)/func.cc -o $(OBJDIR)/func.o

# COMPILING MODELS:

$(OBJDIR)/models.o: $(MODDIR)/models.cc $(MODH) $(GEOMH) $(TOOLH) $(UTILH) $(IODIR)/io.h 
	$(CPPM) $(MODDIR)/models.cc -o $(OBJDIR)/models.o

$(OBJDIR)/mod_coral.o: $(MODDIR)/Coral/main.cc $(MODH) $(GEOMH) $(UTILH) $(IODIR)/io.h $(MAINDIR)/templates.cc
	$(CPPM) $(MODDIR)/Coral/main.cc -o $(OBJDIR)/mod_coral.o

$(OBJDIR)/mod_crystal.o: $(MODDIR)/Crystal/main.cc $(MODH) $(GEOMH) $(UTILH) $(IODIR)/io.h $(MAINDIR)/templates.cc
	$(CPPM) $(MODDIR)/Crystal/main.cc -o $(OBJDIR)/mod_crystal.o

$(OBJDIR)/mod_gasflow.o: $(MODDIR)/Gasflow/main.cc $(MODDIR)/Gasflow/gasflow.h $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc 
	$(CPPM) $(MODDIR)/Gasflow/main.cc -o $(OBJDIR)/mod_gasflow.o

$(OBJDIR)/mod_particles.o: $(MODDIR)/Particles/main.cc $(MODDIR)/Particles/particles.h $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc 
	$(CPPM) -I./$(MODDIR)/Gasflow $(MODDIR)/Particles/main.cc -o $(OBJDIR)/mod_particles.o

$(OBJDIR)/mod_particleflow.o: $(MODDIR)/Particleflow/main.cc $(MODDIR)/Gasflow/gasflow.h $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc 
	$(CPPM) -I./$(MODDIR)/Gasflow $(MODDIR)/Particleflow/main.cc -o $(OBJDIR)/mod_particleflow.o

$(OBJDIR)/mod_incomflow.o: $(MODDIR)/Incomflow/main.cc $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc 
	$(CPPM) -I./$(MODDIR)/Poisson $(MODDIR)/Incomflow/main.cc -o $(OBJDIR)/mod_incomflow.o

$(OBJDIR)/mod_bonding.o: $(MODDIR)/Bonding/main.cc $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc 
	$(CPPM) $(MODDIR)/Bonding/main.cc -o $(OBJDIR)/mod_bonding.o

$(OBJDIR)/mod_collapse.o: $(MODDIR)/Collapse/main.cc $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc
	$(CPPM) $(MODDIR)/Collapse/main.cc -o $(OBJDIR)/mod_collapse.o

$(OBJDIR)/mod_lattice.o: $(MODDIR)/Lattice/main.cc $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc
	$(CPPM) $(MODDIR)/Lattice/main.cc -o $(OBJDIR)/mod_lattice.o

$(OBJDIR)/mod_elastic.o: $(MODDIR)/Elastic/main.cc $(MODDIR)/Gasflow/gasflow.h $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc 
	$(CPPM) -I./$(MODDIR)/Gasflow $(MODDIR)/Elastic/main.cc -o $(OBJDIR)/mod_elastic.o

$(OBJDIR)/mod_membrane.o: $(MODDIR)/Membrane/main.cc $(MODDIR)/Gasflow/gasflow.h $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc 
	$(CPPM) -I./$(MODDIR)/Gasflow $(MODDIR)/Membrane/main.cc -o $(OBJDIR)/mod_membrane.o

$(OBJDIR)/mod_poisson.o: $(MODDIR)/Poisson/poisson.h $(MODDIR)/Poisson/main.cc $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc 
	$(CPPM) $(MODDIR)/Poisson/main.cc -o $(OBJDIR)/mod_poisson.o

$(OBJDIR)/mod_rfg.o: $(MODDIR)/RFG/main.cc $(MODDIR)/RFG/mainrfg.h $(MODDIR)/RFG/rfg.h $(MODH) $(GEOMH) $(FUNH) $(MAINDIR)/templates.cc 
	$(CPPM) $(MODDIR)/RFG/main.cc -o $(OBJDIR)/mod_rfg.o

$(OBJDIR)/mod_prfg.o: $(MODDIR)/PRFG/main.cc $(MODDIR)/RFG/mainrfg.h $(MODDIR)/Particles/particles.h $(MODH) $(GEOMH) $(MAINDIR)/vecalg.h 
	$(CPPM) -I./$(MODDIR)/Particles -I./$(MODDIR)/RFG $(MODDIR)/PRFG/main.cc -o $(OBJDIR)/mod_prfg.o

#$(OBJDIR)/mod_pstatic.o: PStatic/main.cc RFG/main.h Particles/main.h  $(GEOMH) $(MAINDIR)/vecalg.h
#	$(CPP) PStatic/main.cc -o $(OBJDIR)/mod_pstatic.o

$(OBJDIR)/rfg.o: $(MODDIR)/RFG/rfg.c $(MODDIR)/RFG/rfg.h $(MODH) $(MAINDIR)/vecalg.h
	$(CC) -I./$(MAINDIR) $(MODDIR)/RFG/rfg.c -o $(OBJDIR)/rfg.o

$(OBJDIR)/gauss.o: $(MODDIR)/RFG/gauss.c
	$(CC) $(MODDIR)/RFG/gauss.c -o $(OBJDIR)/gauss.o

#.cc.o:
#	$(CPP) $<
#
#.c.o:
#	$(CC) $<
#
#.f.o:
#	$(FC) $<

clean:
	-rm $(PROG)
	-rm $(OBJDIR)/*.o # $(OBJDIR)/*.ti $(OBJDIR)/*.ii

tar:
	cd ..; tar cvfj mulphys-`date +%y%m%d`.tbz mulphys/ --exclude "*.o" --exclude "*.ti" --exclude "*.ii" --exclude "*.dat" --exclude $(PROG)
	cd ..; mv -i mulphys-`date +%y%m%d`.tbz arc/
