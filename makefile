RenderSpherical	:	IDLinterface.o Render.o
			g++ -shared -o RenderSpherical.so IDLinterface.o Render.o
IDLinterface.o	:	IDLinterface.cpp Render.h
			g++ -c -O3 -fPIC -D LINUX IDLinterface.cpp
Render.o	:	Render.cpp
			g++ -c -O3 -fPIC -D LINUX Render.cpp			
