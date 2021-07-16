import math
from typing import Tuple
from operator import add, neg
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation as R

class LSystem:
  forwardDrawSymbols = ('F', 'G')
  forwardJumpSymbols = ('f', 'g')
  def __init__(self, variables : Tuple[str, ...], constants : Tuple[str, ...], rules : dict[str, str], angle : float, initialString : str, initialDirection = (1, 0)):
    self.variables = variables
    self.constants = constants
    self.rules = rules
    self.angle = angle
    self.initialString = initialString
    self.initialDirection = initialDirection
    self.initRotatedVectors()

  def initRotatedVectors(self):
    self.rotatedVectors = {0 : self.initialDirection} # cached

  def getKthRotatedUnitVector(self, k : int):
    if k not in self.rotatedVectors: # TODO: modulo 2 * pi / theta
      idx0, idx1 = self.initialDirection
      coskth = math.cos(self.angle * k)
      sinkth = math.sin(self.angle * k)
      self.rotatedVectors[k] = coskth * idx0 - sinkth * idx1, sinkth * idx0 + coskth * idx1
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
    currentDirection = self.initialDirection
    segments = []
    stack = []
    rotationIndex = 0
    for symbol in finalString:
      if symbol in self.forwardDrawSymbols:
        nextPoint = tuple(map(add, currentPoint, currentDirection))
        segments.append((currentPoint, nextPoint))
        currentPoint = nextPoint
      if symbol in self.forwardJumpSymbols:
        nextPoint = tuple(map(add, currentPoint, currentDirection))
        currentPoint = nextPoint
      elif symbol == '+':
        rotationIndex += 1
        currentDirection = self.getKthRotatedUnitVector(rotationIndex)
      elif symbol == '-':
        rotationIndex -= 1
        currentDirection = self.getKthRotatedUnitVector(rotationIndex)
      elif symbol == '[':
        stack.append((currentPoint, currentDirection, rotationIndex))
      elif symbol == ']':
        currentPoint, currentDirection, rotationIndex = stack.pop()
    return segments

class LSystem3D(LSystem):
  def __init__(self, variables : Tuple[str, ...], constants : Tuple[str, ...], rules : dict[str, str], angle : float, initialString : str,
               initialHeadDirection = np.array([0, 0, 1]), initialLeftArmDirection = np.array([-1, 0, 0])):
    LSystem.__init__(self, variables, constants, rules, angle, initialString)
    self.initialHeadDirection = initialHeadDirection
    self.initialLeftArmDirection = initialLeftArmDirection

  def getSegments(self, finalString : str):
    currentPoint = np.array([0, 0, 0])
    currentHeadDirection = self.initialHeadDirection
    currentLeftArmDirection = self.initialLeftArmDirection
    segments = []
    stack = []
    for symbol in finalString:
      if symbol in self.forwardDrawSymbols:
        nextPoint = currentPoint + currentHeadDirection
        segments.append((currentPoint, nextPoint))
        currentPoint = nextPoint
      if symbol in self.forwardJumpSymbols:
        currentPoint = currentPoint + currentHeadDirection
      elif symbol in '+-&^\\/':
        if symbol == '+':
          rotationAxis = np.cross(currentHeadDirection, currentLeftArmDirection)
        elif symbol == '-':
          rotationAxis = np.cross(currentLeftArmDirection, currentHeadDirection)
        elif symbol == '&':
          rotationAxis = currentLeftArmDirection
        elif symbol == '^':
          rotationAxis = -currentLeftArmDirection
        elif symbol == '\\':
          rotationAxis = currentHeadDirection
        elif symbol == '/':
          rotationAxis = -currentHeadDirection
        r = R.from_rotvec(self.angle * rotationAxis)
        currentHeadDirection = r.apply(currentHeadDirection)
        currentLeftArmDirection = r.apply(currentLeftArmDirection)
      elif symbol == '|':
        currentHeadDirection = -currentHeadDirection
        currentLeftArmDirection = -currentLeftArmDirection
      elif symbol == '[':
        stack.append((currentPoint, currentHeadDirection, currentLeftArmDirection))
      elif symbol == ']':
        currentPoint, currentHeadDirection, currentLeftArmDirection = stack.pop()
    return segments

class KochCurve(LSystem):
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 9
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
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 9
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
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 9
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
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 11
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
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 11
  def __init__(self):
    LSystem.__init__(
      self,
      ('F', 'G'),
      ('+', '-'),
      {'F' : 'G+F+G', 'G' : 'F-G-F'},
      math.pi / 3,
      'G'
      )

class FractalPlantA(LSystem):
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 25
  def __init__(self):
    LSystem.__init__(
      self,
      ('F'),
      ('+', '-', '[', ']'),
      {'F' : 'F[+F]F[-F]F'},
      math.pi * 25.7 / 180,
      'F',
      (0, 1)
      )

class FractalPlantB(LSystem):
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 25
  def __init__(self):
    LSystem.__init__(
      self,
      ('F'),
      ('+', '-', '[', ']'),
      {'F' : 'F[+F]F[-F][F]'},
      math.pi * 20 / 180,
      'F',
      (0, 1)
      )

class FractalPlantC(LSystem):
  # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 25
  def __init__(self):
    LSystem.__init__(
      self,
      ('F'),
      ('+', '-', '[', ']'),
      {'F' : 'FF-[-F+F+F]+[+F-F-F]'},
      math.pi * 22.5 / 180,
      'F',
      (0, 1)
      )

class FractalPlantD(LSystem):
    # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 25
  def __init__(self):
    LSystem.__init__(
      self,
      ('X', 'F'),
      ('+', '-', '[', ']'),
      {'X' : 'F[+X]F[-X]+X', 'F' : 'FF'},
      math.pi * 20 / 180,
      'X',
      (0, 1)
      )

class FractalPlantE(LSystem):
    # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 25
  def __init__(self):
    LSystem.__init__(
      self,
      ('X', 'F'),
      ('+', '-', '[', ']'),
      {'X' : 'F[+X][-X]FX', 'F' : 'FF'},
      math.pi * 25.7 / 180,
      'X',
      (0, 1)
      )

class FractalPlantF(LSystem):
    # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 25
  def __init__(self):
    LSystem.__init__(
      self,
      ('X', 'F'),
      ('+', '-', '[', ']'),
      {'X' : 'F-[[X]+X]+F[+FX]-X', 'F' : 'FF'},
      math.pi * 22.5 / 180,
      'X',
      (0, 1)
      )

class HilberCurve3D(LSystem3D):
    # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 20
  def __init__(self):
    LSystem3D.__init__(
      self,
      ('A', 'B', 'C', 'D'),
      ('F', '+', '-', '&', '^', '\\', '/', '|'),
      {
        'A' : 'B-F+CFC+F-D&F^D-F+&&CFC+F+B//',
        'B' : 'A&F^CFB^F^D^^-F-D^|F^B|FC^F^A//',
        'C' : '|D^|F^B-F+C^F^A&&FA&F^C+F+B^F^D//',
        'D' : '|CFB-F+B|FA&F^A&&FB-F+B|FC//'
      },
      math.pi / 2,
      'A'
      )

class Bush3D(LSystem3D):
    # Lindenmayer, A., Prusinkiewicz, P. (2012). The Algorithmic Beauty of Plants. United States: Springer New York. p. 26
  def __init__(self):
    LSystem3D.__init__(
      self,
      ('A', 'F', 'S', 'L'),
      ('f', '+', '-', '&', '^', '\\', '/', '|', '[', ']', '!', '\''),
      {
        'A' : '[&FL!A]/////\'[&FL!A]///////\'[&FL!A]',
        'F' : 'S/////F',
        'S' : 'FL',
        'L' : '[\'\'\'^^{-f+f+f-|-f+f+f}]'
      },
      math.pi / 8,
      'A'
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

def getAxesLimits(segments):
  coordinates = [c for s in segments for p in s for c in p]
  return (min(coordinates), max(coordinates))

def plot3D(segments):
  lineSegments = Line3DCollection(segments)
  ax = plt.figure().add_subplot(projection='3d')
  ax.add_collection3d(lineSegments)
  axesLimits = getAxesLimits(segments)
  ax.set_xlim3d(axesLimits[0], axesLimits[1])
  ax.set_ylim3d(axesLimits[0], axesLimits[1])
  ax.set_zlim3d(axesLimits[0], axesLimits[1])
  ax.set_aspect('auto')
  plt.show()

def run():
  system = Bush3D()
  iter = 7
  finalString = system.getFinalString(iter)
  print(finalString)
  segments = system.getSegments(finalString)
  print(segments)
  plot3D(segments)

if __name__ == '__main__':
  run()