INCLUDE = -Iinclude -I`root-config --incdir`
LIBS = `root-config --libs`  -l RooFit -lMathMore
FLAGS = -g

HFILES = include/RooBsTimeAngle.h include/RooAngle.h include/RooErrPdf.h
OFILES =     obj/RooBsTimeAngle.o    obj/RooAngle.o      obj/RooErrPdf.o

all: WsBs GenBs DataBs FitBs FitAcceptance

obj/%.o:src/%.cxx
	g++ -fPIC $(FLAGS) $(INCLUDE) $(LIBS) -o $@ -c $<

lib/libBFitter.so: $(HFILES) $(OFILES)
	rootcint -f dict.cxx -c $(HFILES)
	g++ -fPIC $(FLAGS) $(INCLUDE) $(LIBS) -o obj/dict.o -c dict.cxx
	rm dict.cxx dict.h
	g++ $(FLAGS) -shared -o $@ $(INCLUDE) $(LIBS) $(OFILES) obj/dict.o 
#	ar rcs lib/libBFitter.a $(OFILES) obj/dict.o

WsBs: bs_work.cpp lib/libBFitter.so
	g++ $(FLAGS) -o $@ $(INCLUDE) $(LIBS) $(OFILES) obj/dict.o $<
GenBs: bs_gen.cpp lib/libBFitter.so
	g++ $(FLAGS) -o $@ $(INCLUDE) $(LIBS) $(OFILES) obj/dict.o $<
DataBs: bs_data.cpp lib/libBFitter.so
	g++ $(FLAGS) -o $@ $(INCLUDE) $(LIBS) $(OFILES) obj/dict.o $<
FitBs: bs_fit.cpp lib/libBFitter.so
	g++ $(FLAGS) -o $@ $(INCLUDE) $(LIBS) $(OFILES) obj/dict.o $<
FitAcceptance: fitAcceptance.cpp lib/libBFitter.so
	g++ $(FLAGS) -o $@ $(INCLUDE) $(LIBS) $(OFILES) obj/dict.o $<

clean:
	rm -v obj/*.o lib/lib* WsBs GenBs DataBs FitBs

#%: %.cpp lib/libBFitter.so
#	g++ $(FLAGS) -o $@ $(INCLUDE) $(LIBS) $(OFILES) obj/dict.o $<

