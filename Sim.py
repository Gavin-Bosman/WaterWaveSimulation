# library for visualizing our simulation
import sys

import pygame
# library for 2D physics simulation
import pymunk

# pymunk defines a space, which is an area within we are calculating physics
    # we define the gravity of this space

# spaces contain bodies, which by default are particles
# we can apply a shape around this body that can collide with other matter

def create_obj(space):
    shape_radius = 25
    density = 3

    # static bodies do not move, but can collide with other bodies
    # dynamic bodies can be moved by physics
    # kinematic bodies can be moved by the user
    pymunk.Body()
    # setting mass and intertia to 0 allows chipmunk to calculate those values based on the shapes density
    body = pymunk.Body(0, 0, body_type=pymunk.Body.DYNAMIC) 
    body.position=(400,shape_radius)

    shape = pymunk.Circle(body, shape_radius)
    shape.density = density
    space.add(body,shape)
    return shape


def draw_shapes(shapes):
    for shape in shapes:
        # pygame expects int, pymunk uses float
        pos_x = int(shape.body.position.x)
        pos_y = int(shape.body.position.y)

        pygame.draw.circle(screen, (255,255,255), (pos_x,pos_y), 25)

pygame.init()
screen = pygame.display.set_mode((800, 800)) # our display surface
clock = pygame.time.Clock()
space = pymunk.Space()
space.gravity = (0, 500)
shapes = []
shapes.append(create_obj(space))

# main display loop
while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()
    
    screen.fill((0,0,0))
    draw_shapes(shapes)
    space.step(1/30)
    
    for body in space.bodies:

        if (body.position.y >= 800):
            pygame.quit()
            sys.exit()

    pygame.display.update()
    clock.tick(30)