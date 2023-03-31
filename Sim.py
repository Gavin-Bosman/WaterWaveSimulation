# Import Required Libraries
import pygame as pg
import pymunk as pm
import pymunk.pygame_util as pmg
import math
import numpy as np

# Initialize PyGame
pg.init()

# Set Up Window
WIDTH, HEIGHT = 1280, 720
window = pg.display.set_mode((WIDTH, HEIGHT))

# Create pymunk space
space = pm.Space()

# Draw The Simulation


def draw(space, window, draw_options, objects):
    window.fill("black")

    # for object in objects:
    # object.draw(window)

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
        shape = pm.Circle(body, 1)
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
            b0, body, (0, 0), (0, 0), (HEIGHT/3), 50, 5)
        space.add(b0, body, shape, joint)
        # space.add(body, shape)

        # joint.b.gravity_scale = (0, 0)
        # joint.a.gravity_scale = (0, 0)

        self.body = body
        self.shape = shape
        self.joint = joint
        self.anchor = b0

    # def draw(self, surface):
        # pg.draw.circle(surface, "white", (self.x, self.y), 25)
        # pg.draw.circle(surface, "white", (self.x, self.height), 25)

    # ====================================
    # END Object Classes
    # ====================================


def main(window, width, height):
    # Main Loop Vars
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
    intervals = np.linspace(10, WIDTH-10, 500)
    print(intervals)
    wave = [SpringPoints(loc, height/2) for loc in intervals]
    # wave = [SpringPoints(WIDTH/2, height/2)]

    # Drawing Options
    draw_options = pmg.DrawOptions(window)
    draw_options.flags = pm.SpaceDebugDrawOptions.DRAW_SHAPES

    # Main Simulation Loop
    while run:
        # Get Events

        # For each point
        for springPoint in wave:

            if springPoint.dragging:
                # Stop that point from moving if dragging
                springPoint.body.velocity = (0, 0)
                springPoint.joint.b.position = pg.mouse.get_pos()

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

                if event.type == pg.MOUSEBUTTONUP:
                    for springPoint in wave:
                        springPoint.dragging = False

                if dragging:
                    springPoint.body.velocity = (0, 0)
                    springPoint.joint.b.position = pg.mouse.get_pos()

            objects.append(springPoint)

        # Check for Pause/Play State
        if not paused:
            draw(space, window, draw_options, objects)
            space.step(dt)
            clock.tick(fps)

    pg.quit()


if __name__ == "__main__":
    main(window, WIDTH, HEIGHT)
