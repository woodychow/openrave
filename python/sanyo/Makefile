#!/bin/bash
CXX = g++
CXXFLAGS = -O3 -fPIC
LDFLAGS = -llapack
SOFLAGS = -shared

SO = ikfast0x1000004a.TranslationDirection5D.1_2_3_4_5_f0.so
EXEC = $(basename $(SO))

all: $(SO) $(EXEC)

$(SO): $(SO:.so=.cpp)
	$(CXX) -o $@ $? $(CXXFLAGS) $(SOFLAGS) $(LDFLAGS)

$(EXEC): % : %.cpp
	$(CXX) -o $@ $? $(CXXFLAGS) $(LDFLAGS)

cp:
	-@DEST=~/mujin/checkoutroot/openrave/python/sanyo; \
	rm -rf $$DEST; \
	mkdir $$DEST; \
	cp -r `pwd`/* $$DEST; \
	rm $$DEST/$(SO) $$DEST/$(EXEC)

clean:
	rm -rf $(EXEC) a.out *.so *~ *# ./.~

re:
	make clean && make && make cp
