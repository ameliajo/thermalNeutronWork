import unittest
from convolution import interpolate


class MyFirstTests(unittest.TestCase):
    def test_hello(self):
        xVals = [1,2,3,5]
        yVals = [1,3,2,6]
        for i in range(4):
            self.assertAlmostEqual(interpolate(xVals,yVals,xVals[i]),yVals[i])
        xTests = [1.5,3.2,4.0,4.1,4.7]
        yTests = [2.0,2.4,4.0,4.2,5.4]
        for i in range(len(xTests)):
            self.assertAlmostEqual(interpolate(xVals,yVals,xTests[i]),yTests[i])
        xTests = [0.9,0.999,5.001,5.1]
        for i in range(len(xTests)):
            self.assertAlmostEqual(interpolate(xVals,yVals,xTests[i]),0.0)


if __name__ == "__main__":
     unittest.main()
