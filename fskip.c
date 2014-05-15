#include <stdio.h>

// Skip whitespace
// Skip shell script style comments starting with #
// Skip C style single line comments, and multi-line comments
void fskip(FILE *fp){
	//while(isspace(is.peek())){ is.ignore(1); }
	//is >> std::ws;
	// parser inspired by JSON_parser.c
	enum states{
		GO,//S_DEFAULT,           // default state, not in comment
		SL,//S_GOT_SLASH,         // got the first slash, possibly starting a comment
		LN,//S_IN_LINE_COMMENT,   // in a single line comment
		BK,//S_IN_BLOCK_COMMENT,  // in a block (multi-line) comment
		ST,//S_GOT_STAR           // got (possibly) ending star of block comment,
		// these two are not real states; reaching them exits the function
		//S_ERROR_SLASH,       // non space, need to put back a slash
		//S_ERROR,             // non space
		NR_STATES
	};
	enum actions{ // these must be negative (to be different from states)
		RET = -1, // return
		RPS = -2  // return and put back a slash character
	};
	enum classes{
		C_SLASH,
		C_HASH,
		C_STAR,
		C_EOL,
		C_SPACE,
		C_ETC,
		NR_CLASSES
	};
		
	static int ascii_class[64] = {
	// This array maps the lower 64 ASCII characters into character classes.
	// The remaining characters should be mapped to C_ETC.
		C_ETC,   C_ETC,   C_ETC,  C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
		C_ETC,   C_SPACE, C_EOL,  C_ETC,   C_ETC,   C_EOL,   C_ETC,   C_ETC,
		C_ETC,   C_ETC,   C_ETC,  C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
		C_ETC,   C_ETC,   C_ETC,  C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,

		C_SPACE, C_ETC,   C_ETC,  C_HASH,  C_ETC,   C_ETC,   C_ETC,   C_ETC,
		C_ETC,   C_ETC,   C_STAR, C_ETC,   C_ETC,   C_ETC,   C_ETC, C_SLASH,
		C_ETC,   C_ETC,   C_ETC,  C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
		C_ETC,   C_ETC,   C_ETC,  C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC
	};

	
	/*
	State transition table:
	State\Char |  /   |  #   |  *   |  \n  | space (non \n) | anything else
	-----------------------------------------------------------------------
	       DEF | GOT/ | LINE | ERR  | DEF  | DEF            | ERR
	      GOT/ | LINE | ERR/ | BLOK | ERR/ | ERR/           | ERR/
	      LINE | LINE | LINE | LINE | DEF  | LINE           | LINE
	      BLOK | BLOK | BLOK | GOT* | BLOK | BLOK           | BLOK
	      GOT* | DEF  | BLOK | GOT* | BLOK | BLOK           | BLOK
	      ERR/ | ERR  | ERR  | ERR  | ERR  | ERR            | ERR
	       ERR | ERR  | ERR  | ERR  | ERR  | ERR            | ERR
	When hitting the error state, we may have to put back some chars into the stream
	*/
	
		
	static int state_transition_table[NR_STATES][NR_CLASSES] = {
	// The state transition table takes the current state and the current symbol,
	// and returns either a new state or an action.

	//       /    #    *   EOL space etc
	/*GO*/ {GO , LN , RET, GO , GO , RET},
	/*SL*/ {LN , RPS, BK , RPS, RPS, RPS},
	/*LN*/ {LN , LN , LN , GO , LN , LN },
	/*BK*/ {BK , BK , ST , BK , BK , BK },
	/*ST*/ {GO , BK , ST , BK , BK , BK }
	};

	int state = GO;
	while(!feof(fp)){
		int cclass;
		int c = fgetc(fp); ungetc(c, fp);
		cclass = (c < 64) ? ascii_class[c] : C_ETC;
		state = state_transition_table[state][cclass];
		if(state < 0){
			if(RET == state){
				return;
			}else if(RPS == state){
				ungetc('/', fp);
				return;
			}
		}else{
			fgetc(fp);
			if(C_EOL == cclass){
				if('\r' == c){ // possibly need to eat a \n
					if((c = fgetc(fp)) == '\n'){
						// ignore
					}else{
						ungetc(c, fp);
					}
				}
			}
		}
	}
}
