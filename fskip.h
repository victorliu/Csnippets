#ifndef FSKIP_H_INCLUDED
#define FSKIP_H_INCLUDED

#include <stdio.h>

// Skip whitespace
// Skip shell script style comments starting with #
// Skip C style single line comments, and multi-line comments
void fskip(FILE *fp);

#endif
