# Import Required Libraries
import pygame as pg
import pymunk as pm
import pymunk.pygame_util as pmg
import math
import numpy as np
from scipy.interpolate import interp1d

# Initialize PyGame
pg.init()

# Set Up Window
WIDTH, HEIGHT = 1280, 720
window = pg.display.set_mode((WIDTH, HEIGHT))

# Create pymunk space
space = pm.Space()

# Draw The Simulation


def draw(space, window, draw_options, wave):
    window.fill("black")

    # ====================================
    # Create Smooth Curve
    # ====================================

    # ====================================
    # FFT
    # ====================================
    # Extract x and y coordinates from the points
    x_coords = [point.body.position[0] for point in wave]
    y_coords = [point.body.position[1] for point in wave]

    # Compute the FFT of the y-coordinates
    y_fft = np.fft.fft(y_coords)

    # Apply a low-pass filter by zeroing out higher frequency components
    cutoff = 6  # Modify this value to adjust the smoothness of the filtered line
    y_fft[cutoff:-cutoff] = 0

    # Compute the inverse FFT to obtain the filtered y-coordinates
    y_filtered = np.fft.ifft(y_fft)

    # Draw lines between the filtered points
    for i in range(len(x_coords) - 1):
        pg.draw.line(window, (255, 255, 255, 100),
                     (x_coords[i], y_filtered[i].real), (x_coords[i + 1], y_filtered[i + 1].real), 3)

    # ====================================
    # FFT
    # ====================================

    # ====================================
    # SPLINE (Slower than FFT? + Freaking out at certain points)
    # ====================================

    # Extract x and y coordinates from the points
    # xCoords = [point.body.position[0] for point in wave]
    # yCoords = [point.body.position[1] for point in wave]

    # # Interpolation function
    # interp = interp1d(xCoords, yCoords, kind='cubic', fill_value="extrapolate")

    # # Generate new interpolated points
    # xInterp = np.linspace(xCoords[0], yCoords[-1], num=len(wave))
    # yInterp = interp(xInterp)

    # # Draw lines between the interpolated points
    # for i in range(len(xInterp) - 1):
    #     pg.draw.line(window, (255, 255, 255, 100),
    #                  (xInterp[i], yInterp[i]), (xInterp[i + 1], yInterp[i + 1]), 3)

    # for object in range(len(wave)-1):
    #     # Draw Lines Connecting Adjacent Points
    #     # if object != 0 and object != len(wave)-1:
    #     pg.draw.line(window, (255, 255, 255, 100),
    #                  (wave[object].body.position), (wave[object + 1].body.position), 3)

    # ====================================
    # SPLINE
    # ====================================

    # ====================================
    # Create Smooth Curve
    # ====================================

    space.debug_draw(draw_options)

    pg.display.update()


# Get Distance Between Two Points
def distance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

# ====================================
# CREATION FUNCTIONS
# ====================================

# Create Boundaries for our Simulation


def create_boundaries(space, width, height):
    rectRadius = 3
    borderColor = (0, 0, 0, 100)
    rects = [
        # Ground
        [(width/2, height-rectRadius), (width, rectRadius)],
        # Ceiling
        [(width/2, rectRadius/2), (width, rectRadius)],
        # Walls
        [(rectRadius/2, height/2), (rectRadius, height)],
        [(width-rectRadius, height/2), (rectRadius, height)],
    ]

    for pos, size in rects:
        body = pm.Body(body_type=pm.Body.STATIC)
        body.position = pos
        shape = pm.Poly.create_box(body, size)
        shape.elasticity = 0.4
        shape.friction = 0.5
        shape.color = borderColor
        space.add(body, shape)


# ====================================
# END CREATION FUNCTIONS
# ====================================

# ====================================
# Object Classes
# ====================================

# Create Spring Points Class
class SpringPoints:
    def __init__(self, x=0, y=0, height=None):
        self.dampening = 0.1
        self.tension = 0.1
        self.height = height
        self.velocity = 0
        self.x = x
        self.y = y

        self.dragging = False

        body = pm.Body(body_type=pm.Body.DYNAMIC)
        body.position = (self.x, self.y)
        # shape = pm.Circle(body, 4.5)
        shape = pm.Circle(body, 0.01)
        shape.mass = 1
        shape.elasticity = 0.5
        shape.friction = 0.2
        shape.color = (255, 255, 255, 100)

        # Spring
        b0 = pm.Body(body_type=pm.Body.KINEMATIC)
        b0.position = (int(self.x), int(HEIGHT))
        p0 = self.x, HEIGHT
        joint = pm.constraints.DampedSpring(
            # b0, body, (self.x, HEIGHT), (self.x, HEIGHT/2), HEIGHT, 50, 5)
            # b0, body, (0, 0), (0, 0), (HEIGHT/3), 10, 30)
            b0, body, (0, 0), (0, 0), (HEIGHT/3), 0.5, 0.8)
        space.add(b0, body, shape, joint)
        # space.add(body, shape)

        # joint.b.gravity_scale = (0, 0)
        # joint.a.gravity_scale = (0, 0)

        # Slide joint constraint to constrain the motion of the spring along the x-axis
        # slide_joint = pm.constraints.SlideJoint(
        # body, b0, (0, 0), (0, 0), -5, 5)
        # space.add(slide_joint)

        self.body = body
        self.shape = shape
        self.joint = joint
        self.anchor = b0
        # self.slide_joint = slide_joint

    # def draw(self, surface):
        # pg.draw.circle(surface, "white", (self.x, self.y), 25)
        # pg.draw.circle(surface, "white", (self.x, self.height), 25)

    # ====================================
    # END Object Classes
    # ====================================


def main(window, width, height):
    # Main Loop Vars

    # Drawing Options
    draw_options = pmg.DrawOptions(window)
    draw_options.flags = pm.SpaceDebugDrawOptions.DRAW_SHAPES

    run = True
    paused = False
    clock = pg.time.Clock()
    fps = 60
    dt = 1 / fps

    # Check for Mouse Drag
    dragging = False

    # Objects to Draw
    objects = []

    # Setup Gravity for the Space
    # space.gravity = (0, 9.81 * 100)
    space.gravity = (0, 9.81 * 100 * 0)

    # Call Creation Functions
    create_boundaries(space, width, height)

    # springPoint = SpringPoints(width/2, height/2)

    # Create Wave
    # wave = [SpringPoints(loc, height/2) for loc in range(width//5)]
    intervals = np.linspace(20, WIDTH-20, 100)
    # intervals = np.linspace(20, WIDTH-20, 5)
    print(intervals)
    wave = [SpringPoints(loc, height/2) for loc in intervals]
    wavePointsPos = [(point.x, point.y) for point in wave]
    # wave = [SpringPoints(width/2, height/2)]

    # Current Wave Height (change if want to account for displacement)
    WAVEHEIGHT = 480
    stiffness = 60
    damping = 5
    # Rest Length Scale
    SCALE = 1/4
    SCALE = 1/3
    SCALE = 1/2
    # SCALE = 1

    # Link Adjacent Wave Points
    for ind in range(len(wave)):
        # Join Farmost Left Point to Left Boundary
        if ind == 0:
            print(f'{wave[ind].body.position=}')
            leftBoundary = pm.Body(body_type=pm.Body.KINEMATIC)
            leftBoundary.position = (20, WAVEHEIGHT)
            print(f'{leftBoundary.position=}')
            joint = pm.constraints.DampedSpring(
                leftBoundary, wave[ind].body, (0, 0), (0, 0),
                wave[ind].x, 10, 10
            )
            joint2 = pm.constraints.DampedSpring(
                wave[ind].body, wave[ind + 1].body, (0, 0), (0, 0), distance(
                    wave[ind].body.position, wave[ind + 1].body.position)*SCALE, stiffness, damping
            )
            space.add(joint2)

            space.add(joint)

        # Join Farmost Right Point to Right Boundary
        elif ind == len(wave)-1:
            rightBoundary = pm.Body(body_type=pm.Body.KINEMATIC)
            rightBoundary.position = (WIDTH-20, WAVEHEIGHT)
            joint = pm.constraints.DampedSpring(
                wave[ind].body, rightBoundary, (0, 0), (0, 0), distance(
                    wave[ind].body.position, rightBoundary.position)*SCALE, 10, 10
            )

        # Join To Adjacent Points
        else:
            joint = pm.constraints.DampedSpring(
                wave[ind].body, wave[ind + 1].body, (0, 0), (0, 0), distance(
                    wave[ind].body.position, wave[ind + 1].body.position)*SCALE, stiffness, damping
            )
            space.add(joint)

    decreasingS = False
    increasingS = False
    decreasingD = False
    increasingD = False

    # Main Simulation Loop
    while run:
        # Get Events

        # For each point
        for springPoint in wave:
            # print(f'{springPoint.body.position=}')

            if springPoint.dragging:
                # Stop that point from moving if dragging
                springPoint.body.velocity = (0, 0)
                springPoint.joint.b.position = (
                    springPoint.joint.b.position[0], pg.mouse.get_pos()[1])

            for event in pg.event.get():
                # Handle Quite Event
                if event.type == pg.QUIT:
                    run = False
                    break

                # Handle Pausing/Playing The Simulation
                if event.type == pg.KEYDOWN and event.key == pg.K_SPACE:
                    if paused:
                        paused = False
                    else:
                        paused = True

                # Adjust stiffness and damping
                if event.type == pg.KEYDOWN and event.key == pg.K_p or increasingS:
                    springPoint.joint.stiffness += 1
                    increasingS = True
                    print(f'{springPoint.joint.stiffness=}')

                if event.type == pg.KEYUP and event.key == pg.K_p:
                    increasingS = False

                if event.type == pg.KEYDOWN and event.key == pg.K_o or decreasingS:
                    springPoint.joint.stiffness -= 1
                    decreasingS = True
                    print(f'{springPoint.joint.stiffness=}')

                if event.type == pg.KEYUP and event.key == pg.K_o:
                    decreasingS = False

                if event.type == pg.KEYDOWN and event.key == pg.K_SEMICOLON or increasingD:
                    springPoint.joint.damping += 1
                    increasingD = True
                    print(f'{springPoint.joint.damping=}')

                if event.type == pg.KEYUP and event.key == pg.K_SEMICOLON:
                    increasingD = False

                if event.type == pg.KEYDOWN and event.key == pg.K_l or decreasingD:
                    springPoint.joint.springPoint.joint.damping -= 1
                    decreasingD = True
                    print(f'{springPoint.joint.damping=}')

                if event.type == pg.KEYUP and event.key == pg.K_l:
                    decreasingD = False

                # Moving Ball On Mouseclick

                # Check for drag
                if (event.type == pg.MOUSEBUTTONDOWN and event.button == 1) or dragging:
                    if springPoint.shape.point_query(pg.mouse.get_pos()):

                        # Get closest Spring Point
                        min_distance = float("inf")
                        closest_spring = None

                        for springPoint in wave:
                            dist = distance(
                                springPoint.body.position, event.pos)
                            if dist < min_distance:
                                min_distance = dist
                                closest_spring = springPoint

                        if closest_spring is not None:
                            closest_spring.dragging = True
                    print(f'{pg.mouse.get_pos()=}')

                if event.type == pg.MOUSEBUTTONUP:
                    for springPoint in wave:
                        springPoint.dragging = False

                # for springPoint in wave:
                    # if springPoint.dragging:
                        # springPoint.body.velocity = (0, 0)
                        # springPoint.joint.b.position = pg.mouse.get_pos()

        # Check for Pause/Play State
        if not paused:
            draw(space, window, draw_options, wave)
            space.step(dt)
            clock.tick(fps)

    pg.quit()


if __name__ == "__main__":
    main(window, WIDTH, HEIGHT)
