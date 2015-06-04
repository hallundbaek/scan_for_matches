#include <Python.h>
#include "../../scan_for_matches.h"
#include "../../parser.h"
#include "../../scanner.h"
#include <stdio.h>
#include <string.h>

//static char module_docstring[] = "The scan_for_matches module provides an interface for finding patterns in a file!";

static char sfm_docstring[] = "Finds patterns in a file!";

static PyObject* scan_for_matches(PyObject *self, PyObject *args);

static PyObject *
scan_for_matches(self, args)
  PyObject *self, *args;
{
  static PyObject *ScanError;

  ScanError = PyErr_NewException("scan.error", NULL, NULL);
  Py_INCREF(ScanError);

  match_t *matches, *tmp;
  char *fasta_file, *line, *ignore_file;
  int protein, complements, show_overlaps, max_hits, stop_after, verbose, count;
  
  if (!PyArg_ParseTuple(args, "sssiiiiii", &fasta_file, &line, &ignore_file, &protein, &complements, &show_overlaps, &max_hits, &stop_after, &verbose)){
    PyErr_SetString(ScanError, "Parsing arguments failed");
    return NULL;
  }
  
  matches = find_matches(fasta_file, line, ignore_file, protein, complements, show_overlaps, max_hits, stop_after, verbose);
  tmp = matches;
  
  count = 1;
  while(tmp->next){
    count++;
    tmp = tmp->next;
  }
  
  PyObject *ret_list, *ret_list_elem, *start_end_list;
  
  ret_list = PyList_New(count); 
  ret_list_elem = PyList_New(2);
  start_end_list = PyList_New(2);

  int i; 
  for(i = 0; i < count; i++){
    //start_end_list: [start,end]
    PyList_Insert(start_end_list, 0, Py_BuildValue("i",matches->start));
    PyList_Insert(start_end_list, 1, Py_BuildValue("i",matches->end));

    //ret_elem_list: [[start,end],match]
    PyList_Insert(ret_list_elem, 0, start_end_list);
    PyList_Insert(ret_list_elem, 1, Py_BuildValue("s",matches->match));

    //ret_list: [[[start,end],match],...]
    PyList_Insert(ret_list, i, ret_list_elem);
    matches = matches->next;
  }
  
  Py_DECREF(ret_list_elem);
  Py_DECREF(start_end_list);
  
  return ret_list;

}

static PyMethodDef
module_methods[] =
{
  {"scan_for_matches",scan_for_matches, METH_VARARGS, sfm_docstring},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC 
initscan_for_matches(void)
{
  Py_InitModule("scan_for_matches", module_methods);
}
