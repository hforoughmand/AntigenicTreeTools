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

ifdef CONDA_PREFIX
USR := ${CONDA_PREFIX}
else
USR := /usr
endif

# compilers
JAVAC		= $(USR)/bin/javac
JAVAH		= $(USR)/bin/javah
CC			= $(USR)/bin/gcc

LFLAGS		= -lc -lgfortran -shared
# change the following row to match your system configuration (include directories where 'jni.h' and 'jni_md.h' are located)
ifdef CONDA_PREFIX
LIB_LFLAGS	= -I ${CONDA_PREFIX}/lib/jvm/include/linux -I ${CONDA_PREFIX}/lib/jvm/include/
else 
LIB_LFLAGS	= -I /usr/lib/jvm/java-*/include/ -I /usr/lib/jvm/java-*/include/linux/
endif

# create bin dir & compile java
$(BIN)/libbvlslib.so: $(SRC)/phyloDriver.java
	@ mkdir -p $(BIN)
	$(JAVAC) -cp $(SRC)/:$(JAR)/Jama.jar -d $(BIN) $(SRC)/phyloDriver.java
	#$(JAVAH) -bootclasspath $(BIN) -d $(BIN) -jni NNLSsolver
	$(JAVAC) -h $(BIN) -cp $(BIN):$(JAR)/Jama.jar $(SRC)/NNLSsolver.java
	$(CC) -c -g -o $(BIN)/bvls.o $(SRC)/bvls.f90 -fPIC
	$(CC) -o $(BIN)/libbvlslib.so -Wl,-soname,$(BIN)/libbvlslib.so $(LIB_LFLAGS) $(SRC)/bvls.c $(BIN)/bvls.o $(LFLAGS) -fPIC

$(BIN)/AntigenicTreeTools.jar: $(BIN)/phyloDriver.class
	jar vcfe $(BIN)/AntigenicTreeTools.jar phyloDriver $(BIN)/phyloDriver.class $(BIN)/*.class $(BIN)/*.o $(BIN)/*.so jar/Jama.jar