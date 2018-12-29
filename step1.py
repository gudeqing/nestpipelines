import os
import time

n = 5
while n > 0:
    n -= 1
    for i in os.listdir('.'):
        print(i*5)
    time.sleep(1)
print('this is step1')
