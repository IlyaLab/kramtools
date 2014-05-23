
#include <Python.h>

static PyObject * _run( PyObject *self, PyObject *args ) {
	FILE *fp = NULL;
	PyObject *fobj = NULL;
	if( ! PyArg_ParseTuple( args, "O", &fobj ) )
		return NULL;
	if( ! PyFile_Check( fobj ) )
		return NULL;
	fp = PyFile_AsFile( fobj );
	if( fp ) {
		char *line = NULL;
		size_t n   = 0;
		PyFile_IncUseCount(fobj);
		Py_BEGIN_ALLOW_THREADS
		while( getline( &line, &n, fp ) > 0 ) {
			fputs( line, stdout );
		}
		if( line )
			free( line );
		Py_END_ALLOW_THREADS
		PyFile_DecUseCount(fobj);
	}
	//fprintf( stdout, "Hello from _run (%p)\n", fobj );
	Py_RETURN_NONE;
}

static PyMethodDef methods[] = {
	{"run",_run,METH_VARARGS},
	{NULL,NULL},
};

PyMODINIT_FUNC initpairwise(void) {
	PyObject *m
		= Py_InitModule( "pairwise", methods );
}

