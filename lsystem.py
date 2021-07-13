import math
from typing import Tuple
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from numpy import double

class LSystem:
  forwardSymbols = ('A', 'B', 'F', 'G')
  rotationSymbols = ('+', '-')
  def __init__(self, variables : Tuple[str, ...], constants : Tuple[str, ...], rules : dict[str, str], angles : Tuple[double, ...], initialString : str):
    self.variables = variables
    self.constants = constants
    self.rules = rules
    self.angles = angles
    self.initialString = initialString
    self.rotationMatrices = []
    for angle in self.angles:
      self.rotationMatrices.append([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])

  def rotate(self, v : Tuple[double, ...], nth : int):
    R = self.rotationMatrices[nth]
    return (R[0][0] * v[0] + R[0][1] * v[1], R[1][0] * v[0] + R[1][1] * v[1])

  def getFinalString(self, iter : int):
    currentString = self.initialString
    for _ in range(iter):
      nextString = ''
      for symbol in currentString:
        if symbol in self.variables:
          nextString += self.rules[symbol]
        else:
          nextString += symbol
      currentString = nextString
    return currentString

  def getSegments(self, finalString : str):
    currentPoint = (0, 0)
    currentDirection = (1, 0) # TODO: compute initial direction
    segments = []
    stack = []
    for symbol in finalString:
      if symbol in self.forwardSymbols:
        nextPoint = (currentPoint[0] + currentDirection[0], currentPoint[1] + currentDirection[1])
        segments.append([currentPoint, nextPoint])
        currentPoint = nextPoint
      elif symbol in self.rotationSymbols:
        currentDirection = self.rotate(currentDirection, self.rotationSymbols.index(symbol))
      elif symbol == '[':
        stack.append(currentPoint)
      elif symbol == ']':
        currentPoint = stack.pop()
    return segments

class KochCurve(LSystem):
  def __init__(self):
    LSystem.__init__(
      self, 
      ('A'),
      ('+', '-'),
      {'A' : 'A+A-A-A+A'},
      (math.pi / 2, -math.pi / 2),
      'A'
      )

class SerpenskiArrowheadCurve(LSystem):
  def __init__(self):
    LSystem.__init__(
      self,
      ('A', 'B'),
      ('+', '-'),
      {'A' : 'B+A+B', 'B' : 'A-B-A'},
      (math.pi / 3, -math.pi / 3),
      'A'
      )

class SerpenskiTriangle(LSystem):
  def __init__(self):
    LSystem.__init__(
      self,
      ('A', 'B'),
      ('+', '-'),
      {'A' : 'A-B+A+B-A', 'B' : 'BB'},
      (2 * math.pi / 3, -2 * math.pi / 3),
      'A-B-B'
      )

class SerpenskiSquare(LSystem):
  def __init__(self):
    LSystem.__init__(
      self,
      ('X'),
      ('F', '+', '-'),
      {'X' : 'XF−F+F−XF+F+XF−F+F−X'}, # incorrect?
      (math.pi / 2, -math.pi / 2),
      'F+XF+F+XF'
      )

class SerpenskiCurve(LSystem):
  def __init__(self):
    LSystem.__init__(
      self,
      ('X'),
      ('F', 'G', '+', '-'),
      {'X' : 'XF+G+XF−−F−−XF+G+X'},
      (math.pi / 4, -math.pi / 4),
      'F−−XF−−F−−XF'
      )

class DragonCurve(LSystem):
  def __init__(self):
    LSystem.__init__(
      self,
      ('A', 'B'),
      ('+', '-'),
      {'A' : 'A+B', 'B' : 'A-B'},
      (math.pi / 2, -math.pi / 2),
      'A'
      )

class FractalPlant(LSystem):
  def __init__(self):
    LSystem.__init__(
      self,
      ('X', 'F'),
      ('+', '-', '[', ']'),
      {'X' : 'F+[[X]-X]-F[-FX]+X', 'F' : 'FF'},
      (math.pi * 25 / 180, -math.pi * 25 / 180),
      'X'
      )

def plot(segments):
  lineSegments = LineCollection(segments)
  _, ax = plt.subplots()
  ax.add_collection(lineSegments)
  ax.autoscale()
  ax.set_aspect(1)
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  plt.tight_layout()
  plt.show()

def run():
  system = FractalPlant()
  iter = 6
  finalString = system.getFinalString(iter)
  print(finalString)
  segments = system.getSegments(finalString)
  print(segments)
  plot(segments)

if __name__ == '__main__':
  run()