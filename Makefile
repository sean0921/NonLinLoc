# Makefile for NonLinLoc software package
#
# Invocation:
#     Solaris: make all
#     Linux:   make -R all

BINDIR=${MYBIN}
INCLUDE_DIR=./
CCFLAGS=-O3
#CCFLAGS=-O3 -pg
#CCFLAGS=-g

# Linux
TIME3D_CCFLAGS=

all : NLLoc Vel2Grid Grid2Time Grid2GMT LocSum Time2EQ PhsAssoc hypoe2hyp fpfit2hyp
distrib : NLLoc Vel2Grid Grid2Time Grid2GMT LocSum Time2EQ PhsAssoc hypoe2hyp fpfit2hyp


# --------------------------------------------------------------------------
# NLLoc
#
OBJS1=GridLib.o geo.o octtree.o util.o nrutil.o nrmatrix.o ran1.o map_project.o
PVER=1
NLLoc : ${BINDIR}/NLLoc
${BINDIR}/NLLoc : NLLoc${PVER}.o ${OBJS1}
	gcc NLLoc${PVER}.o ${OBJS1} ${CCFLAGS} -o ${BINDIR}/NLLoc -lm
NLLoc${PVER}.o : NLLoc${PVER}.c GridLib.h ${INCLUDE_DIR}nrutil.h
	gcc -c ${CCFLAGS} NLLoc${PVER}.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# Vel2Grid
#
OBJS2=GridLib.o geo.o util.o nrutil.o nrmatrix.o velmod.o map_project.o
PVER=1
Vel2Grid : ${BINDIR}/Vel2Grid
${BINDIR}/Vel2Grid : Vel2Grid${PVER}.o ${OBJS2}
	gcc Vel2Grid${PVER}.o  ${OBJS2} ${CCFLAGS} -o ${BINDIR}/Vel2Grid -lm
Vel2Grid${PVER}.o : Vel2Grid${PVER}.c GridLib.h
	gcc ${CCFLAGS} -c Vel2Grid${PVER}.c
velmod.o : ${INCLUDE_DIR}velmod.c
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}velmod.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# Grid2Time
#
OBJS3=GridLib.o geo.o util.o nrutil.o nrmatrix.o ran1.o \
		map_project.o Time_3d.o
PVER=1
Grid2Time : ${BINDIR}/Grid2Time
${BINDIR}/Grid2Time : Grid2Time${PVER}.o ${OBJS3}
	gcc Grid2Time${PVER}.o ${OBJS3} ${CCFLAGS}  \
		-o ${BINDIR}/Grid2Time -lm
Grid2Time${PVER}.o : Grid2Time${PVER}.c GridLib.h
	gcc ${CCFLAGS} -c Grid2Time${PVER}.c
Time_3d.o : Time_3d.c
	gcc ${TIME3D_CCFLAGS} -c -DNO_IEEE_PROTOCOL Time_3d.c
#	gcc ${CCFLAGS} -c -DNO_IEEE_PROTOCOL Time_3d.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# Grid2GMT
#
OBJS4=Grid2GMT.o geo.o GridLib.o GridGraphLib.o vector.o util.o nrutil.o \
	nrmatrix.o map_project.o
Grid2GMT : ${BINDIR}/Grid2GMT
${BINDIR}/Grid2GMT : ${OBJS4}
	gcc ${OBJS4} ${CCFLAGS} -o ${BINDIR}/Grid2GMT -lm
Grid2GMT.o : Grid2GMT.c GridLib.h GridGraphLib.h
	gcc ${CCFLAGS} -c Grid2GMT.c
GridGraphLib.o : GridGraphLib.c GridLib.h GridGraphLib.h
	gcc ${CCFLAGS} -c GridGraphLib.c
#
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# LocSum
#
OBJS5=LocSum.o GridLib.o geo.o util.o nrutil.o nrmatrix.o map_project.o
LocSum : ${BINDIR}/LocSum
${BINDIR}/LocSum : ${OBJS5}
	gcc ${OBJS5} ${CCFLAGS} -o ${BINDIR}/LocSum -lm
LocSum.o : LocSum.c GridLib.h
	gcc ${CCFLAGS} -c LocSum.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# Time2EQ
#
OBJS6=GridLib.o geo.o util.o nrutil.o nrmatrix.o ran1.o map_project.o
PVER=1
Time2EQ : ${BINDIR}/Time2EQ
${BINDIR}/Time2EQ : Time2EQ${PVER}.o ${OBJS6}
	gcc Time2EQ${PVER}.o ${OBJS6} ${CCFLAGS} -o ${BINDIR}/Time2EQ -lm

Time2EQ${PVER}.o : Time2EQ${PVER}.c GridLib.h ${INCLUDE_DIR}ran1.h
	gcc -c ${CCFLAGS} Time2EQ${PVER}.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# hypoe2hyp
#
OBJS7=hypoe2hyp.o GridLib.o geo.o util.o nrutil.o nrmatrix.o ran1.o map_project.o
hypoe2hyp : ${BINDIR}/hypoe2hyp
${BINDIR}/hypoe2hyp : ${OBJS7}
	gcc ${OBJS7} ${CCFLAGS} -o ${BINDIR}/hypoe2hyp -lm
hypoe2hyp.o : hypoe2hyp.c GridLib.h
	gcc ${CCFLAGS} -c hypoe2hyp.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# fpfit2hyp
#
OBJS8=fpfit2hyp.o GridLib.o geo.o util.o nrutil.o nrmatrix.o ran1.o map_project.o
fpfit2hyp : ${BINDIR}/fpfit2hyp
${BINDIR}/fpfit2hyp : ${OBJS8}
	gcc ${OBJS8} ${CCFLAGS} -o ${BINDIR}/fpfit2hyp -lm
fpfit2hyp.o : fpfit2hyp.c GridLib.h
	gcc ${CCFLAGS} -c fpfit2hyp.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# PhsAssoc
#
OBJS9=PhsAssoc.o GridLib.o geo.o util.o nrutil.o nrmatrix.o map_project.o
PhsAssoc : ${BINDIR}/PhsAssoc
${BINDIR}/PhsAssoc : ${OBJS9}
	gcc ${OBJS9} ${CCFLAGS} -o ${BINDIR}/PhsAssoc -lm
PhsAssoc.o : PhsAssoc.c GridLib.h
	gcc ${CCFLAGS} -c PhsAssoc.c
# --------------------------------------------------------------------------





# --------------------------------------------------------------------------
# Librarires
#
GridLib.o : GridLib.c GridLib.h  geometry.h nrutil.h util.h geo.h
	gcc -c ${CCFLAGS} GridLib.c

util.o : ${INCLUDE_DIR}util.c ${INCLUDE_DIR}util.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}util.c

octtree.o : ${INCLUDE_DIR}octtree.c ${INCLUDE_DIR}octtree.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}octtree.c

nrutil.o : ${INCLUDE_DIR}nrutil.c ${INCLUDE_DIR}nrutil.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}nrutil.c

nrmatrix.o : ${INCLUDE_DIR}nrmatrix.c
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}nrmatrix.c

ran1.o : ${INCLUDE_DIR}ran1.c ${INCLUDE_DIR}ran1.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}ran1.c

vector.o : ${INCLUDE_DIR}vector.c ${INCLUDE_DIR}vector.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}vector.c

geo.o : ${INCLUDE_DIR}geo.c ${INCLUDE_DIR}geo.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}geo.c

GMT_INCLUDE=./
map_project.o : ${GMT_INCLUDE}map_project.c
	gcc -c ${CCFLAGS} ${GMT_INCLUDE}map_project.c
#
# --------------------------------------------------------------------------
