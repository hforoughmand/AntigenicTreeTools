#*************************************************************************
#
#    This file is part of the software to infer antigenic trees.
#    Copyright (C) 2012  Lars Steinbrueck
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#**************************************************************************

# directories
SRC		= src
JAR		= jar
BIN		= bin

# compilers
JAVAC		= ${CONDA_PREFIX}/bin/javac
JAVAH		= ${CONDA_PREFIX}/bin/javah
CC		= ${CONDA_PREFIX}/bin/gcc

LFLAGS		= -lc -lgfortran -shared
# change the following row to match your system configuration (include directories where 'jni.h' and 'jni_md.h' are located)
LIB_LFLAGS	= -I ${CONDA_PREFIX}/lib/jvm/include/linux -I ${CONDA_PREFIX}/lib/jvm/include/

# create bin dir & compile java
$(BIN)/phyloDriver.class: $(SRC)/phyloDriver.java
	@ mkdir -p $(BIN)
	$(JAVAC) -cp $(SRC)/:$(JAR)/Jama.jar -d $(BIN) $(SRC)/phyloDriver.java
	#$(JAVAH) -bootclasspath $(BIN) -d $(BIN) -jni NNLSsolver
	$(JAVAC) -h $(BIN) -cp $(BIN):$(JAR)/Jama.jar $(SRC)/NNLSsolver.java
	$(CC) -c -g -o $(BIN)/bvls.o $(SRC)/bvls.f90 -fPIC
	$(CC) -o $(BIN)/libbvlslib.so -Wl,-soname,$(BIN)/libbvlslib.so $(LIB_LFLAGS) $(SRC)/bvls.c $(BIN)/bvls.o $(LFLAGS) -fPIC

$(BIN)/AntigenicTreeTools.jar: $(BIN)/phyloDriver.class
	jar vcfe $(BIN)/AntigenicTreeTools.jar phyloDriver $(BIN)/phyloDriver.class $(BIN)/*.class $(BIN)/*.o $(BIN)/*.so jar/Jama.jar