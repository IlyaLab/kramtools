
#ifndef _varfmt_h_
#define _varfmt_h_

#define MAX_OUTPUT_COLUMNS (24)

struct feature_pair;

/**
  * This (admittedly odd) macro is used to insure that all
  * emitter methods have the exactly same functional signature.
  */
#define EMITTER_SIG const struct feature_pair *pair, const struct CovariateAnalysis *covan, FILE *fp

typedef void (*EMITTER_FXN)( EMITTER_SIG );

#define FORMAT_TABULAR (0)
#define FORMAT_JSON    (1)

const char *emit_config( const char *specifier, int format );
void emit_exec( EMITTER_SIG );

#endif

