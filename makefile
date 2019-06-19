RenderSpherical	:	IDLinterface.o Render.o
			g++ -shared -o RenderSpherical.so IDLInterface.o Render.o
IDLInterface.o	:	IDLInterface.cpp Render.h
			g++ -c -O3 -fPIC -D LINUX IDLInterface.cpp
Render.o	:	Render.cpp
			g++ -c -O3 -fPIC -D LINUX Render.cpp			
