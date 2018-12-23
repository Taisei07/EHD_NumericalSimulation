# -*- coding: utf-8 -*-

import numpy

a = numpy.random.rand(100)
b = numpy.random.rand(100)

print(numpy.linalg.norm(a - b))

a = numpy.random.rand(100)
B = numpy.random.rand(1000000, 100)

print([numpy.linalg.norm(a - b) for b in B])
