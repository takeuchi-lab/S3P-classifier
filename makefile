AR       = ar
ARFLAGS  = rus

CXX      = g++

CXXFLAGS = -g -MMD -MP -Wall -Wextra -Winit-self -Wno-unused-parameter -fopenmp -std=c++11 -O3

RM       = rm -f
LDFLAGS  = -fopenmp
LIBS     =
INCLUDE  = -I../include

TARGET   = ./train

OBJDIR   = ./obj

SOURCES  = $(wildcard *.cc)

OBJECTS  = $(addprefix $(OBJDIR)/, $(SOURCES:.cc=.o))

DEPENDS  = $(OBJECTS:.o=.d)

$(TARGET): $(OBJECTS) $(LIBS)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: %.cc
	@if [ ! -d $(OBJDIR) ];\
	then echo "mkdir -p $(OBJDIR)";mkdir -p $(OBJDIR);\
	fi
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

#clean and build
.PHONY: all
all: clean $(TARGET)

#clean
.PHONY:clean
clean:
	$(RM) $(OBJECTS) $(DEPENDS) $(TARGET)

-include $(DEPENDS)
