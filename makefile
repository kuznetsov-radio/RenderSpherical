RenderSpherical	:	IDLinterface.o Render.o Trace.o
			g++ -shared -o RenderSpherical.so IDLinterface.o Render.o Trace.o
IDLinterface.o	:	IDLinterface.cpp Render.h Trace.h
			g++ -c -O3 -fPIC -D LINUX IDLinterface.cpp
Render.o	:	Render.cpp
			g++ -c -O3 -fPIC -D LINUX Render.cpp			
Trace.o		:	Trace.cpp
			g++ -c -O3 -fPIC -D LINUX Trace.cpp