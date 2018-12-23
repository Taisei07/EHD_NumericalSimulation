# -*- coding: utf-8 -*-

import numpy
import cupy
import time
start_time = time.time()

a = numpy.random.rand(100)
b = numpy.random.rand(100)

numpy.linalg.norm(a - b)

end1 = time.time()
print(end1-start_time)

a = numpy.random.rand(100)
B = numpy.random.rand(1000000, 100)

[numpy.linalg.norm(a - b) for b in B]

end2 = time.time()
print(end2-end1)

numpy.linalg.norm(a - B, axis=1)

end3 = time.time()
print(end3-end2)

a = cupy.random.rand(100)
B = cupy.random.rand(1000000, 100)

cupy.linalg.norm(a - B, axis=1)

end4 = time.time()
print(end4-end3)
