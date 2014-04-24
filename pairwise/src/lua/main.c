
/**
  * This is just a demonstration and test platform for development of
  * Lua extensions. 
  * It translates command line arguments into a call to a named Lua function
  * defined in the provided Lua script.
  */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <err.h>

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

static bool opt_call_as_coroutine = false;
static int  opt_verbosity         = 0;

static int _load_script( const char *source, lua_State *state ) {

	struct stat info;
	if( (access( source, R_OK ) == 0) && (stat( source, &info ) == 0) ) {
		return luaL_dofile( state, source );
	} else
		return luaL_dostring( state, source );
}

static const char *USAGE = 
"%s -c <file|script> <fxn name> [args ...]\n"
"Output options:\n"
"  --coroutine | -c   call <fxn name> as a coroutine [%s]\n"
"  --verbosity | -v verbosity [%d]\n"
"Notes:\n"
#ifdef _DEBUG
"This is a _DEBUG build.\n"
#endif
"Written by rkramer@systemsbiology.org.\n";


static void _print_usage( const char *exename, FILE *fp ) {
	static const char *_T = "yes";
	static const char *_F = "no";
	fprintf( fp, USAGE, 
		exename,
		opt_call_as_coroutine ? _T : _F,
		opt_verbosity );
}


int main( int argc, char *argv[] ) {

	struct stat info;
	double result;
	char *program = NULL;
	const char *func = NULL;
	const char *fsrc = NULL;

	int nargs   = 0;
	float *args = NULL;

	if( argc < 2 ) {
		_print_usage( argv[0], stdout );
		exit(0);
	}

	lua_State *L = luaL_newstate();
	luaL_openlibs(L);

	do {
		static const char *CHAR_OPTIONS 
			= "cv:?";
		static struct option LONG_OPTIONS[] = {

			{"coroutine",  0,0,'c'},
			{"verbosity",  1,0,'v'},
			{ NULL,        0,0, 0 }
		};

		int opt_offset = 0;
		const int c = getopt_long( argc, argv, CHAR_OPTIONS, LONG_OPTIONS, &opt_offset );
		switch (c) {

		case 'c': opt_call_as_coroutine = true; break;

		case 'v':
			opt_verbosity = atoi( optarg );
			break;

		case '?':
			_print_usage( argv[0], stdout );
			exit(0);
			break;
		case -1: // ...signals no more options.
			break;
		default:
			printf ("error: unknown option: %c\n", c );
			exit(-1);
		}
		if( -1 == c ) break;

	} while( true );

	fsrc = argv[optind++];
	func = argv[optind++];

	if( optind < argc ) {
		nargs = argc-optind;
		int i = 0;
		args = calloc( nargs, sizeof(float) );
		while( optind < argc )
			args[i++] = atoi( argv[optind++] );
	}


	if( _load_script( fsrc, L ) == LUA_OK ) {

		int isnum;
		lua_getglobal( L, func ); // "pushes onto stack value of the global name"
		if( lua_isnil( L, -1 ) )
			errx( -1, "%s not defined (in Lua's global namespace)", func );

		if( opt_call_as_coroutine ) {
			int lua_status = LUA_YIELD;
			while( lua_status == LUA_YIELD ) {
				lua_status = lua_resume( L, NULL, 0 );
				if( lua_status <= LUA_YIELD /* OK == 0, YIELD == 1*/ ) {
					const int b = lua_tonumberx( L, -1, &isnum );
					const int a = lua_tonumberx( L, -2, &isnum );
					lua_pop( L, 2 );
					fprintf( stdout, "%d %d\n", a, b );
				} else
					fputs( lua_tostring(L,-1), stderr );
			}

		} else {

			int i = nargs;
#if 0
			while( i-- > 0 ) {
				lua_pushnumber( L, args[i] );
			}
#else
			for(i = 0; i < nargs; i++ ) {
				lua_pushnumber( L, args[i] );
			}
#endif
			if( lua_pcall( L, 2, 1, 0 ) == LUA_OK ) {
				result = lua_tonumberx( L, -1, &isnum );
				if( isnum )
					fprintf( stdout, "result = %f\n", result );
				else
					fprintf( stdout, "result was not a number!\n" );
			} else
				fprintf( stderr, "lua_pcall failed: %s\n", lua_tostring(L,-1) );
			lua_pop( L, 1 );

		}
	} else
		errx( -1, "failed executing %s: %s", fsrc, lua_tostring(L,-1) );

	if( args )
		free( args );

	lua_close(L);

	return 0;
}


