
CXX = g++

CXXFLAGS = -DNDEBUG -Wno-deprecated -finline-functions -I../eigen/include/eigen3/Eigen/
LINKFLAGS = -lm

SRCS = \
	Main.cpp \
	Lcm.cpp \
	Database.cpp \
	OccurenceDeriver.cpp \
	plars.cpp \

OBJS = $(SRCS:%.cpp=%.o)

SRCS2 = \
	Predict.cpp \
	
OBJS2 = $(SRCS2:%.cpp=%.o)

all: plars predict

plars: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) $(OBJS) $(LINKFLAGS) -o plars

predict: $(OBJS2)
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) $(OBJS2) $(LINKFLAGS) -o predict

debug:
	make all CXXFLAGS="-ggdb -W -Wall -pedantic"

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) -c $<

clean:
	rm -f plars predict *.o *~
