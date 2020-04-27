#include <stdio.h>
#include <string.h>

#include "interface.h"


int main()
{
	printf("running example....\nminimap2 version = ");

	enum { kMaxArgs = 64 };
	int argc = 0;
	char *argv[kMaxArgs];

	 char commandLine[250] = "minimap2 --version";
	char *p2 = strtok(commandLine, " ");
	
	while (p2 && argc < kMaxArgs-1)
	  {
	    argv[argc++] = p2;
	    p2 = strtok(0, " ");
	  }
	argv[argc] = 0;
	
	init_minimap2(argc,argv);

	return 0;
}

