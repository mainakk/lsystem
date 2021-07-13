import math
from typing import Tuple
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from numpy import double

class LSystem:
  forwardDrawSymbols = ('F', 'G')
  forwardJumpSymbols = ('f', 'g')
  def __init__(self, variables : Tuple[str, ...], constants : Tuple[str, ...], rules : dict[str, str], angle : double, initialString : str):
    self.variables = variables
    self.constants = constants
    self.rules = rules
    self.angle = angle
    self.initialString = initialString
    self.rotatedVectors = {0 : (1, 0)} # cached

  def getKthRotatedUnitVector(self, k : int):
    if k not in self.rotatedVectors: # TODO: modulo 2 * pi / theta
      self.rotatedVectors[k] = math.cos(self.angle * k), math.sin(self.angle * k)
    return self.rotatedVectors[k]

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
    rotationIndex = 0
    for symbol in finalString:
      if symbol in self.forwardDrawSymbols:
        nextPoint = (currentPoint[0] + currentDirection[0], currentPoint[1] + currentDirection[1])
        segments.append((currentPoint, nextPoint))
        currentPoint = nextPoint
      if symbol in self.forwardJumpSymbols:
        nextPoint = (currentPoint[0] + currentDirection[0], currentPoint[1] + currentDirection[1])
        currentPoint = nextPoint
      elif symbol == '+':
        rotationIndex += 1
        currentDirection = self.getKthRotatedUnitVector(rotationIndex)
      elif symbol == '-':
        rotationIndex -= 1
        currentDirection = self.getKthRotatedUnitVector(rotationIndex)
      elif symbol == '[':
        stack.append(currentPoint)
      elif symbol == ']':
        currentPoint = stack.pop()
    return segments

class KochCurve(LSystem):
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 9
  def __init__(self):
    LSystem.__init__(
      self, 
      ('F'),
      ('+', '-'),
      {'F' : 'F+F-F-F+F'},
      math.pi / 2,
      'F'
      )

class QuadraticKochIsland(LSystem):
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 9
  def __init__(self):
    LSystem.__init__(
      self, 
      ('F'),
      ('+', '-'),
      {'F' : 'F+FF-FF-F-F+F+FF-F-F+F+FF+FF-F'},
      math.pi / 2,
      'F-F-F-F'
      )

class IslandSnakeCombination(LSystem):
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 9
  def __init__(self):
    LSystem.__init__(
      self, 
      ('F', 'f'),
      ('+', '-'),
      {'F' : 'F+f-FF+F+FF+Ff+FF-f+FF-F-FF-Ff-FFF', 'f' : 'ffffff'},
      math.pi / 2,
      'F+F+F+F'
      )

class SerpenskiTriangle(LSystem):
  # https://en.wikipedia.org/wiki/L-system#Example_5:_Sierpinski_triangle
  def __init__(self):
    LSystem.__init__(
      self,
      ('F', 'G'),
      ('+', '-'),
      {'F' : 'F-G+F+G-F', 'G' : 'GG'},
      2 * math.pi / 3,
      'F-G-G'
      )

class SerpenskiSquare(LSystem):
  # https://en.wikipedia.org/wiki/Sierpi%C5%84ski_curve#Representation_as_Lindenmayer_system
  def __init__(self):
    LSystem.__init__(
      self,
      ('X'),
      ('F', '+', '-'),
      {'X' : 'XF−F+F−XF+F+XF−F+F−X'}, # incorrect?
      math.pi / 2,
      'F+XF+F+XF'
      )

class SerpenskiCurve(LSystem):
  # https://en.wikipedia.org/wiki/Sierpi%C5%84ski_curve#Representation_as_Lindenmayer_system
  def __init__(self):
    LSystem.__init__(
      self,
      ('X'),
      ('F', 'G', '+', '-'),
      {'X' : 'XF+G+XF−−F−−XF+G+X'},
      math.pi / 4,
      'F−−XF−−F−−XF'
      )

class DragonCurve(LSystem):
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 11
  def __init__(self):
    LSystem.__init__(
      self,
      ('F', 'G'),
      ('+', '-'),
      {'F' : 'F+G+', 'G' : '-F-G'},
      math.pi / 2,
      'F'
      )

class SerpenskiGasket(LSystem):
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 11
  def __init__(self):
    LSystem.__init__(
      self,
      ('F', 'G'),
      ('+', '-'),
      {'F' : 'G+F+G', 'G' : 'F-G-F'},
      math.pi / 3,
      'G'
      )

class FractalPlant(LSystem):
  # https://en.wikipedia.org/wiki/L-system#Example_7:_Fractal_plant
  def __init__(self):
    LSystem.__init__(
      self,
      ('X', 'F'),
      ('+', '-', '[', ']'),
      {'X' : 'F+[[X]-X]-F[-FX]+X', 'F' : 'FF'},
      math.pi * 25 / 180,
      'X'
      )

class FractalPlantA(LSystem):
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 25
  def __init__(self):
    LSystem.__init__(
      self,
      ('F'),
      ('+', '-', '[', ']'),
      {'F' : 'F[+F]F[-F]F'},
      math.pi * 25.7 / 180,
      'F'
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
  iter = 2
  finalString = system.getFinalString(iter)
  print(finalString)
  segments = system.getSegments(finalString)
  plot(segments)

if __name__ == '__main__':
  run()