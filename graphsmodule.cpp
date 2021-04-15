#include <Python.h>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <set>

typedef struct {
    PyObject_HEAD;
    uint64_t Vertices;
    std::vector<int> EdgesList[64];
} AdjacencyList;

static PyObject* AdjacencyList__new__(
    PyTypeObject* type, PyObject* args) {
    return type->tp_alloc(type, 0);
}
static void AdjacencyList__del__(AdjacencyList* self) {
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int AdjacencyList__init__(AdjacencyList* self, PyObject* args) {
    char* text;
    // self->EdgesList = new std::vector<int>[64];
    if (PyArg_ParseTuple(args, "s", &text)) {
        int c = int(text[0]);
        int k = 0;
        int i;
        int j = 0;
        self->Vertices = 0;
        if (c < 126) {
            i = 1;
            for (int x = 0; x < (c - 63); x++) {
                self->Vertices += pow(2, x);
                j++;
            }
        }
        else {
            i = 4;
            for (int x = 0; x < (int(text[2])); x++) {
                self->Vertices += pow(2, x);
                j++;
            }
        }
        for (int x = 1; x < j; x++) {
            for (int y = 0; y < x; y++) {
                if (k == 0) {
                    c = int(text[i]) - 63;
                    i++;
                    k = 6;
                }
                k--;
                int64_t bit = 1;
                if ((c & (bit << k)) != 0) {
                    if (x < y) {
                        self->EdgesList[x].push_back(y);
                    }
                    else {
                        self->EdgesList[y].push_back(x);
                    }
                }
            }
        }
        for (int x = 0; x < j; x++) {
            sort(self->EdgesList[x].begin(), self->EdgesList[x].end());
        }
        return 0;
    }
    return -1;
}



static PyObject* number_of_vertices(AdjacencyList* self) {
    uint64_t count = 0;
    uint64_t copy = self->Vertices;
	while (copy) {
		count += copy & 1;
		copy >>= 1;
	}

    return Py_BuildValue("i", count);
}

bool _is_bit_set(AdjacencyList* self, uint64_t n) {
	uint64_t bit = 1;
	return self->Vertices & (bit << n);
}

void _set_bit(AdjacencyList* self, uint64_t n) {
	uint64_t bit = 1;
	self->Vertices |= bit << n;
}

void _clear_bit(AdjacencyList* self, uint64_t n) {
	uint64_t bit = 1;
	self->Vertices &= ~(bit << n);
}

void _print_graph(AdjacencyList* self) {
    std::cout << "self->Vertices: " << self->Vertices << std::endl;
    for (int v = 0; v < 64; v++) {
        if (_is_bit_set(self, v)) {
            std::cout << "[" << v << "]";
            for (std::vector<int>::iterator it = self->EdgesList[v].begin(); it != self->EdgesList[v].end(); it++) {
                std::cout << "->" <<  *it; 
            }
            std::cout << std::endl;
        }
    }
}

static PyObject* vertices(AdjacencyList* self) {
    PyObject *set = PySet_New(0);

    for (uint64_t v = 0; v < 64; v++) {
		if (_is_bit_set(self, v)) {
            PyObject *t = Py_BuildValue("i", v);
            PySet_Add(set, t);
            Py_DECREF(t);
		}
	}
    return set;
}

int _vertex_degree(AdjacencyList* self, uint64_t vertex) {
	return self->EdgesList[vertex].size();
}

static PyObject* number_of_edges(AdjacencyList* self) {
    uint64_t sum = 0;

    for (uint64_t v = 0; v < 64; v++) {
        if (_is_bit_set(self, v)) {
            sum += _vertex_degree(self, v);
        }
    }

    return Py_BuildValue("i", sum);
}

static PyObject* edges(AdjacencyList* self) {
   PyObject* set = PySet_New(0);
    for (int v = 0; v < 64; v++) {
        if (_is_bit_set(self, v)) {
            for (std::vector<int>::iterator it = self->EdgesList[v].begin(); it != self->EdgesList[v].end(); it++) {
                if (v < *it) {
                    PyObject* tup = Py_BuildValue("(ii)", v, *it);
                    PySet_Add(set, tup);
                    Py_DECREF(tup);
                }
                else {
                    PyObject* tup = Py_BuildValue("(ii)", *it, v);
                    PySet_Add(set, tup);
                    Py_DECREF(tup);
                }
            }
        }
    }
    return set;
}

static PyObject* is_edge(AdjacencyList* self, PyObject* args) {
    uint64_t x = 0;
    uint64_t y = 0;
    uint64_t temp = 0;
    PyArg_ParseTuple(args, "ii", &x, &y);
    if (y < x) {
        temp = x;
        x = y;
        y = temp;
    }
    for (std::vector<int>::iterator it = self->EdgesList[x].begin(); it != self->EdgesList[x].end(); it++) {
        if (*it == y) {
            Py_RETURN_TRUE;
        }
    }
    Py_RETURN_FALSE;
}

static PyObject* vertex_degree(AdjacencyList* self, PyObject *args) {
    uint64_t vertex = 0;
    PyArg_ParseTuple(args, "i", &vertex);
    uint64_t degree = 0;
    degree = self->EdgesList[vertex].size();

    for (int v = 0; v < vertex; v++) {
        for (std::vector<int>::iterator it = self->EdgesList[v].begin(); it != self->EdgesList[v].end(); it++) {
            if (*it == vertex) {
                degree += 1;
            }
        }
    }

    return Py_BuildValue("i", degree);
}

static PyObject* vertex_neighbors(AdjacencyList* self, PyObject* args) {
    uint64_t x = 0;
    PyArg_ParseTuple(args, "i", &x);
    PyObject* set = PySet_New(0);
    for (std::vector<int>::iterator it = self->EdgesList[x].begin(); it != self->EdgesList[x].end(); it++) {
        PyObject* t = Py_BuildValue("i", *it);
        PySet_Add(set, t);
        Py_DECREF(t);
    }
 
    for (int v = 0; v < x; v++) {
        for (std::vector<int>::iterator it = self->EdgesList[v].begin(); it != self->EdgesList[v].end(); it++) {
            if (*it == x) {
                PyObject* t = Py_BuildValue("i", v);
                PySet_Add(set, t);
                Py_DECREF(t);
            }
        }
    }
    return set;
}

static PyObject* add_vertex(AdjacencyList* self, PyObject* args) {
    uint64_t x = 0;
    PyArg_ParseTuple(args, "i", &x);
    _set_bit(self, x);

    Py_RETURN_NONE;
}



static PyObject* delete_vertex(AdjacencyList* self, PyObject* args) {
    uint64_t x = 0;
    PyArg_ParseTuple(args, "i", &x);
    if (_is_bit_set(self, x)) {
        self->EdgesList[x].clear();
        for (int v = 0; v < 64; v++) {
            if (_is_bit_set(self, v)) {
                for (std::vector<int>::iterator it = self->EdgesList[v].begin(); it != self->EdgesList[v].end(); it++) {
                    if (*it == x) {
                        self->EdgesList[v].erase(it);
                    }
                }
            }
        }
    }
    _clear_bit(self, x);

    Py_RETURN_NONE;
}

static PyObject* add_edge(AdjacencyList* self, PyObject* args) {

    uint64_t u = 0;
    uint64_t v = 0;
    uint64_t temp = 0;
    PyArg_ParseTuple(args, "ii", &u, &v);

    if (u > v) {
        temp = u;
        u = v;
        v = temp;
    }

    self->EdgesList[u].push_back(v);
    sort(self->EdgesList[u].begin(), self->EdgesList[u].end());

    Py_RETURN_NONE;
}

static PyObject* delete_edge(AdjacencyList* self, PyObject* args) {
    uint64_t u = 0;
    uint64_t v = 0;
    uint64_t temp = 0;
    PyArg_ParseTuple(args, "|ii", &u, &v);

    if (u > v) {
        temp = u;
        u = v;
        v = temp;
    }
    for (std::vector<int>::iterator it = self->EdgesList[u].begin(); it != self->EdgesList[u].end(); it++) {
        if (*it == v){
            self->EdgesList[u].erase(it);
            break;
        }
    }

    Py_RETURN_NONE;
}

static PyObject* square(AdjacencyList*);

static PyObject* AdjacencyList__str__(AdjacencyList* self) {
    PyObject* x = PyFloat_FromDouble(self->Vertices);
    PyObject* r =
        PyUnicode_FromFormat("(%S)", x);
    Py_XDECREF(x);
    return r;
}

static PyMethodDef AdjacencyListMethods[] = {
    {"number_of_vertices",(PyCFunction)number_of_vertices, METH_NOARGS, "Return number of vertices"},
    {"vertices",(PyCFunction)vertices, METH_NOARGS, "Return vertices"},
    {"number_of_edges",(PyCFunction)number_of_edges, METH_NOARGS, "Return number of edges"},
    {"edges",(PyCFunction)edges, METH_NOARGS, "Return edges"},
    {"is_edge",(PyCFunction)is_edge, METH_VARARGS, "Check if edge"},
    {"vertex_degree",(PyCFunction)vertex_degree, METH_VARARGS, "Return vertex degree"},
    {"vertex_neighbors",(PyCFunction)vertex_neighbors, METH_VARARGS, "Return vertex neighbours"},
    {"add_vertex",(PyCFunction)add_vertex, METH_VARARGS, "Add vertex"},
    {"delete_vertex",(PyCFunction)delete_vertex, METH_VARARGS, "Delete vertex"},
    {"add_edge",(PyCFunction)add_edge, METH_VARARGS, "Add edge"},
    {"delete_edge",(PyCFunction)delete_edge, METH_VARARGS, "Delete edge"},
    {"square",(PyCFunction)square, METH_NOARGS, "Return square of input graph"},
    {NULL}
};
static PyTypeObject AdjacencyListType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "simple_graphs.AdjacencyList",
    sizeof(AdjacencyList),
    0,
    (destructor)AdjacencyList__del__,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    (reprfunc)AdjacencyList__str__,
    0, 0, 0,
    Py_TPFLAGS_DEFAULT,
    "AdjacencyList", 
    0, 0, 0, 0, 0, 0,
    AdjacencyListMethods,
    0, 0, 0, 0, 0, 0, 0, 
    (initproc)AdjacencyList__init__,
    0,
    (newfunc)AdjacencyList__new__ 
};

static PyObject* square(AdjacencyList* self) {
    PyObject *argList = Py_BuildValue("(s)", "?");
    PyObject *squareG = PyObject_CallObject((PyObject*)&AdjacencyListType, argList);
	for (int v = 0; v < 64; v++) {
        if (_is_bit_set(self, v)) {
            PyObject* tV =  Py_BuildValue("(i)", v);
	        add_vertex((AdjacencyList*)squareG, tV);
		    for (int u = 0; u < 64; u++) {
                if (_is_bit_set(self, u)) {
                    PyObject *edgesUV = Py_BuildValue("(ii)", u, v);
                    if (is_edge(self, edgesUV) == Py_True) {
                        add_edge((AdjacencyList*)squareG, edgesUV);
                    }
                    else {
	        	        for (int u0 = 0; u0 < 64; u0++) {
                            if (_is_bit_set(self, u0)) {
                                PyObject *edgesU0V = Py_BuildValue("(ii)", u0, v);
                                PyObject *edgesU0U = Py_BuildValue("(ii)", u0, u);
                                if (is_edge(self, edgesU0V) == Py_True && is_edge(self, edgesU0U) == Py_True) {
                                    add_edge((AdjacencyList*)squareG, edgesUV);
                                }
                                Py_DECREF(edgesU0V);
                                Py_DECREF(edgesU0U);
                            }
                        }
                    }
                    Py_DECREF(edgesUV);
                }
            }
            Py_DECREF(tV);
        }
	}
    Py_DECREF(argList);

    return squareG;
}

static PyModuleDef graph_module = {
PyModuleDef_HEAD_INIT,
"simple_graphs",
"Basic graph operations",
-1,
NULL, NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC PyInit_simple_graphs(void) {
    if (PyType_Ready(&AdjacencyListType) < 0) return NULL;
    PyObject* m = PyModule_Create(&graph_module);
    if (m == NULL) return NULL;
    Py_INCREF(&AdjacencyListType);
    PyModule_AddObject(m, "AdjacencyList",
        (PyObject*)&AdjacencyListType);
    return m;
}