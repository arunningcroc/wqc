cimport scqo
cdef class ClassicalField:
    cdef double omega
    cdef double fstrength
    def __cinit__(self, omega = 1, fstrength = 1):
        self.omega = omega
        self.fstrength = fstrength
    def set_omega(self, omega):
        self.omega = omega
    def get_fstrength(self):
        return self.fstrength
    def set_fstrength(self, fst):
        self.fstrength = fst
    def get_omega(self):
        return self.omega
cdef class NLevelAtom:
    cdef double omega
    cdef int nlevels
    cdef bool RWA
    cdef ClassicalField field
    cdef char *datafile
    cdef int counter
    cdef double fstrength
    def __cinit__(self, omega = 1, nlevels = 2, RWA = True):
        self.omega = omega
        self.nlevels = nlevels
        self.RWA = RWA
        self.counter = 0
    def set_field(self, field):
        self.field = field
    def set_name(self, name):
        kek = bytes(name, 'UTF8')
        self.datafile = kek
    def set_rwa(self, RWA):
        self.RWA = RWA
    def set_omega(self, omega):
        self.omega = omega
    def get_detuning(self):
        return self.field.get_omega()-self.omega
    def calculate(self):
        if self.nlevels == 2 and self.RWA == True:
            scqo.solve_ode(scqo.RABI, self.field.get_omega(),self.omega, self.field.get_fstrength(), self.nlevels, self.counter)
            self.counter = self.counter + 1
        elif self.nlevels == 2 and self.RWA == False:
            scqo.solve_ode(scqo.NORWA, self.field.get_omega(),self.omega, self.field.get_fstrength(), self.nlevels, self.counter)
            self.counter = self.counter+1
        else:
            print("Only n=2 implemented so far")
        