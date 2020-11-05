all:
	g++ -o rcdaq2hits RCDAQ_TO_HITS.cpp -I${ONLINE_MAIN}/include -I{OFFLINE_MAIN}/include -L${ONLINE_MAIN}/lib -lSpectrum  -lMinuit -L${OFFLINE_MAIN}/lib -lpmonitor -lEvent -lNoRootEvent -lmessage `root-config --cflags --glibs` -fpermissive
	g++ -o hits2centroid HITS_TO_CENTROID.cpp `root-config --cflags --glibs`
	#g++ -o example connectionExample.C `root-config --cflags --glibs`
