
#include "DawesOzone.hpp"
#include "RynTypes.hpp"


        Real_t DawesOzone_calc_potential_(
            FlatCoordinates coords,
            Names atoms,
            ExtraBools extra_bools,
            ExtraInts extra_ints,
            ExtraFloats extra_floats
            ) {
        
            // Load extra args (if necessary)
            
            int nwalkers = extra_ints[0];
            
            
            // Get data as raw array
            RawWalkerBuffer raw_coords = coords.data();
        
            
            double masses[3]; masses[0] = 15.99994; masses[1] = masses[0]; masses[2] = masses[1];  // feed in masses
            
            
            Real_t energy = -100000;
            calc_potential_(raw_coords,&nwalkers,masses,&energy);
            return energy;
            ;
        }
        
        static PyObject* DawesOzone_calc_potential_Wrapper = PyCapsule_New((void *)DawesOzone_calc_potential_, "_potential", NULL);
        

bool _AttachCapsuleToModule(PyObject* module, PyObject* capsule, const char* name) {
    bool i_did_good = (PyModule_AddObject(module, name, capsule) == 0);
    if (!i_did_good) {
        Py_XDECREF(capsule);
        Py_DECREF(module);
    }

    return i_did_good;
}

static PyMethodDef DawesOzoneMethods[] = {
        {NULL, NULL, 0, NULL}
};

// TODO: ADD IN SOMETHING THAT LETS US GET THE ARGUMENT NAMES DIRECTLY FROM THE LIB

#if PY_MAJOR_VERSION > 2

const char DawesOzone_doc[] = "DawesOzone provides a potential";
static struct PyModuleDef DawesOzoneModule = {
    PyModuleDef_HEAD_INIT,
    "DawesOzone",   /* name of module */
    DawesOzone_doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    DawesOzoneMethods
};

PyMODINIT_FUNC PyInit_DawesOzone(void)
{
    PyObject *m;
    m = PyModule_Create(&DawesOzoneModule);
    if (m == NULL) {
        return NULL;
    }

    if (!_AttachCapsuleToModule(m, DawesOzone_calc_potential_Wrapper, "_potential")) { return NULL; }

    return m;
}

#else

PyMODINIT_FUNC initDawesOzone(void)
{
    PyObject *m;
    m = Py_InitModule("DawesOzone", DawesOzoneMethods);
    if (m == NULL) {
    return NULL
    }

    _AttachCapsuleToModule(m, DawesOzone_PotentialWrapper);
}

#endif