# Import Required Libraries
import pygame as pg
import pygame.gfxdraw as gfx
import pymunk as pm
import pymunk.pygame_util as pmg
import math
import numpy as np
from scipy.interpolate import interp1d
import pygame.mask


# Initialize PyGame
pg.init()

# Load the icon image
icon = pygame.image.load("icon.png")


# Set Up Window
WIDTH, HEIGHT = 1280, 720
window = pg.display.set_mode(
    (WIDTH, HEIGHT), pg.HWSURFACE | pg.DOUBLEBUF | pg.HWACCEL)
pygame.display.set_caption('Ocean Wave Simulation')

# Set the window icon
pygame.display.set_icon(icon)

# Background Image
background_image = pg.image.load('background.png')
background_image = pg.transform.scale(background_image, (WIDTH, HEIGHT))

# Surface for Background
background_surface = pygame.Surface((WIDTH, HEIGHT))
background_surface.blit(background_image, (0, 0))

# Surface for Wave and Objects
wave_surface = pg.Surface((WIDTH, HEIGHT), pg.SRCALPHA)
rock_surface = pg.Surface((WIDTH, HEIGHT), pg.SRCALPHA)
rock_surface.blit(pg.image.load('Rock.png'), (0, 0))

imageRock = pg.image.load("Rock.png").convert_alpha()
imageRockSubmerged = pg.image.load("Rock.png")


# If WaveMode False, Object Mode is on
WaveMode = True
# WindMode
WindMode = True
WINDMAX = 1.1

# Text
HELP = False
pygame.font.init()
font = pygame.font.Font('Avenir LT Std 65 Medium.otf', 15)
textLine1 = font.render("Toggle Help: Press H", True, (37, 37, 37))
windLine = font.render("Wind Intensity:", True, (97,97,97))
windIntensity = font.render("Wind Intensity:", True, (97,97,97))
windIntensityOFF = font.render("OFF", True, (97,97,97))

textLine2 = font.render(
    "Toggle Wind Mode Press E", True, (37, 37, 37))

textLine3 = font.render(
    "Toggle Mode (Object/Wave): Press W", True, (37, 37, 37))

# font = pygame.font.Font('Avenir LT 65 Medium Bold.ttf', 15)
# ModeTitleColor = (42,93,224)
ModeTitleColor = (225,92,17)
WaveModeTitle = font.render("Wave", True, ModeTitleColor)
ObjectModeTitle = font.render("Object", True, ModeTitleColor)
WindModeTitle = font.render("Wind", True, ModeTitleColor)

textLine4 = font.render(
    "Adjust Wind (Wave Mode): (Ctrl + Scroll)", True, (37, 37, 37))
textLine5 = font.render(
    "Adjust Object Size (Object Mode): (Ctrl + Scroll)", True, (37, 37, 37))

# ====================================
# SPACES
# ====================================
# Create pymunk space
# Wave Space
space = pm.Space()
# Object Space
spaceObj = pm.Space()

handler = space.add_collision_handler(1, 2)


def begin(space1, arbiter, space2):
    for c in arbiter.contact_point_set:
        handler.pre_solve(arbiter, space1)
        c.normal_impulse, c.tangent_impulse = arbiter.total_impulses
        handler.post_solve(arbiter, space1)
    return True


# Texture
image = pg.image.load("texture.png")
# Blit the image onto the surface
image_surface = pg.Surface(
    (image.get_width()+95, image.get_height()), pg.SRCALPHA)
image_surface.blit(image, (0, 0))
polygon_surface = pg.Surface((WIDTH+40, HEIGHT), pg.SRCALPHA)
polygon_surface.blit(image_surface, (0, 0), special_flags=pg.BLEND_RGBA_MULT)

# Draw The Simulation


def draw(space, spaceObj, window, draw_options, wave, objects, WAVEHEIGHT, RADIUS, INDICATOR_RADIUS, objectResize, windIntensity):

    wave_surface.fill((0, 0, 0, 0))

    # Draw resizing shape indicator
    if objectResize:
        pg.draw.circle(wave_surface, (255, 255, 255, 100),pg.mouse.get_pos(), INDICATOR_RADIUS)

    # ====================================
    # Create Smooth Curve
    # ====================================

    # FFT

    # Numpy for speed
    x_coords = np.array([point.body.position[0] for point in wave])
    y_coords = np.array([point.body.position[1] for point in wave])
    # Numpy for speed

    # Compute the FFT of the y-coordinates
    y_fft = np.fft.fft(y_coords)

    # Apply a low-pass filter by zeroing out higher frequency components
    cutoff = 6  # Modify this value to adjust the smoothness of the filtered line
    y_fft[cutoff:-cutoff] = 0

    # Compute the inverse FFT to obtain the filtered y-coordinates
    y_filtered = np.fft.ifft(y_fft)

    # Numpy for speed
    vertices = np.array([(x_coords[i], y_filtered[i].real)
                        for i in range(len(x_coords))])
    vertices = np.append(
        vertices, [(x_coords[-1], HEIGHT), (x_coords[0], HEIGHT)], axis=0)

    # Create a new surface for the polygon
    polygon_surface = pg.Surface((WIDTH+40, HEIGHT), pg.SRCALPHA)

    # Draw the polygon onto the surface
    pg.gfxdraw.filled_polygon(polygon_surface, vertices, (255, 255, 255, 255))
    # pg.gfxdraw.aapolygon(polygon_surface, vertices, (14,49,81, 255))

    # Blit the image surface onto the polygon surface
    polygon_surface.blit(image_surface, (0, 0),
                         special_flags=pg.BLEND_RGBA_MULT)

    # Blit the polygon surface onto the main surface
    pg.gfxdraw.aapolygon(wave_surface, vertices, (124, 166, 203, 255))
    wave_surface.blit(polygon_surface, (0, 0))

    # Blit Renders
    window.blit(background_surface, (0, 0))
    window.blit(polygon_surface, (0, 0))
    window.blit(wave_surface, (0, 0))
    
    # Blit Text
    window.blit(textLine1, (10, 10))
    
    # Wind Intensity
    windIntensity = font.render(str(int(abs(windIntensity/0.6*5) - 11) * -1), True, (97,97,97))
    window.blit(windLine, (WIDTH - 145, 10))

    if HELP:
        window.blit(textLine2, (10, 35))
        window.blit(textLine3, (10, 60))
        window.blit(textLine4, (10, 85))
        window.blit(textLine5, (10, 110))
        if WaveMode:
            window.blit(WaveModeTitle, (160.5, 60))
        else:
            window.blit(ObjectModeTitle, (108, 60))
        if WindMode:
            window.blit(WindModeTitle, (60.5, 35))
    if WindMode:
        window.blit(windIntensity, (WIDTH - 30, 10.2))
    else:
        window.blit(windIntensityOFF, (WIDTH - 39, 10.2))


    # ====================================
    # Create Smooth Curve
    # ====================================

    # space.debug_draw(draw_options)

    global imageRock, imageRockSubmerged
    for object in objects:
        # object.body.radius = RADIUS
        # print(object.body.radius)

        if (object.body.position[1] > HEIGHT-WAVEHEIGHT):
            pass
        
        if not object.body.image:
            object.body.image = pygame.transform.scale(
                imageRock, (object.body.radius*2, object.body.radius*2))
        
        rockSurface = pg.Surface(
            (object.body.radius*2, object.body.radius*2), pg.SRCALPHA)
        pg.draw.circle(rockSurface, (255, 255, 255),
                    pg.mouse.get_pos(), RADIUS)
        image_rect = object.body.image.get_rect(center=(object.body.radius, object.body.radius))
        rockSurface.blit(object.body.image, image_rect)
        window.blit(
            rockSurface, (object.body.position[0]-object.body.radius, object.body.position[1] - object.body.radius), special_flags=pg.BLENDMODE_ADD)

    pygame.display.flip()


# Get Distance Between Two Points
def distance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

# ====================================
# CREATION FUNCTIONS
# ====================================

# Create Boundaries for our Simulation

def create_boundaries(space, width, height):
    rectRadius = 3
    borderColor = (255, 255, 255, 20)
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
        shape.filled = False
        space.add(body, shape)


def createObject(space, pos, radius=20, mass=100000, elasticity=0, friction=100):
    body = pm.Body(body_type=pm.Body.DYNAMIC)
    body.position = pos
    body.velocity = (0, 981/2)
    body.radius = radius
    body.image = None
    body.submerged = False
    shape = pm.Circle(body, radius)
    shape.mass = mass
    # shape.elasticity = elasticity
    shape.elasticity = 0
    shape.friction = friction
    shape.color = (255, 255, 255, 100)

    space.add(body, shape)
    return shape

# ====================================
# END CREATION FUNCTIONS
# ====================================


def apply_wind_force(wave, wind_strength, target_point, wind_force_duration, elapsed_wind_time):
    if target_point is not None:
        # Calculate the current wind force based on elapsed time and duration
        current_wind_force = wind_strength * (elapsed_wind_time / wind_force_duration)

        # If the elapsed time has passed half the duration, start reducing the force
        if elapsed_wind_time >= wind_force_duration / 2:
            current_wind_force = wind_strength - current_wind_force

        # Apply the force to the target point
        target_point.body.apply_force_at_local_point((0, current_wind_force), (0, 0))


def is_submerged(obj, wave):
    x_pos = obj.body.position[0]
    y_pos_obj = obj.body.position[1]
    
    # Scale x_pos to the length of the wave array
    x_index = int(x_pos / WIDTH * len(wave))
    
    # Ensure the index is within the bounds of the wave array
    x_index = max(0, min(x_index, len(wave) - 1))
    
    y_pos_wave = wave[x_index].body.position[1]

    # Get closest Spring Point
    min_distance = float("inf")
    closest_spring = None

    for springPoint in wave:
        dist = distance(
            springPoint.body.position, obj.body.position)
        if dist < min_distance:
            min_distance = dist
            closest_spring = springPoint

    return y_pos_obj > closest_spring.body.position[1]


# ====================================
# Object Classes
# ====================================

# Create Spring Points Class


class SpringPoints:
    def __init__(self, space, x=0, y=0, height=HEIGHT/3):
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
            b0, body, (0, 0), (0, 0), (height), 0.5, 0.8)
        space.add(b0, body, shape, joint)

        self.body = body
        self.shape = shape
        self.joint = joint
        self.anchor = b0
        # space.gravity = (0,0)

    # ====================================
    # END Object Classes
    # ====================================


def main(window, width, height):
    # Main Loop Vars

    # Drawing Options
    draw_options = pmg.DrawOptions(window)
    draw_options.flags = pm.SpaceDebugDrawOptions.DRAW_SHAPES

    global HELP

    # If WaveMode False, Object Mode is on
    global WaveMode, WindMode

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
    spaceObj.gravity = (0, 9.81 * 100)

    # Call Creation Functions
    create_boundaries(space, width, height)


    # Create Wave
    WaveHeight = HEIGHT/2
    intervals = np.linspace(-40, WIDTH+40, 100)
    wave = [SpringPoints(space, loc, WaveHeight, HEIGHT//2.4)
            for loc in intervals]
    wave = np.array(wave)

    # Save X Coordinates
    x_coords = np.array([point.body.position[0] for point in wave])

    # Rest Length Scale
    SCALE = 1/4
    SCALE = 1/3
    SCALE = 1/2
    # SCALE = 1
    LEFTjoint = None
    RIGHTjoint = None


    WAVEHEIGHT = WaveHeight + WaveHeight * (SCALE /2.7)
    # WAVEHEIGHT = HEIGHT - HEIGHT * 0.6667
    stiffness = 144
    damping = 5


    # Link Adjacent Wave Points
    for ind in range(len(wave)):
        # Join Farmost Left Point to Left Boundary
        if ind == 0:
            # print(f'{wave[ind].body.position=}')
            leftBoundary = pm.Body(body_type=pm.Body.KINEMATIC)
            leftBoundary.position = (0, WAVEHEIGHT)
            # print(f'{leftBoundary.position=}')
            LEFTjoint = pm.constraints.DampedSpring(
                leftBoundary, wave[ind].body, (0, 0), (0, 0),
                wave[ind].x, 1000000000, 1000000000
            )
            joint2 = pm.constraints.DampedSpring(
                wave[ind].body, wave[ind + 1].body, (0, 0), (0, 0), distance(
                    wave[ind].body.position, wave[ind + 1].body.position)*SCALE, stiffness, damping
            )
            space.add(joint2)

            space.add(LEFTjoint)

        # Join Farmost Right Point to Right Boundary
        elif ind == len(wave)-1:
            rightBoundary = pm.Body(body_type=pm.Body.KINEMATIC)
            rightBoundary.position = (WIDTH, WAVEHEIGHT)
            joint = pm.constraints.DampedSpring(
                wave[ind].body, rightBoundary, (0, 0), (0, 0), distance(
                    wave[ind].body.position, rightBoundary.position)*SCALE, 1000000000, 1000000000
            )

        # Join To Adjacent Points
        else:
            RIGHTjoint = pm.constraints.DampedSpring(
                wave[ind].body, wave[ind + 1].body, (0, 0), (0, 0), distance(
                    wave[ind].body.position, wave[ind + 1].body.position)*SCALE, stiffness, damping
            )
            space.add(RIGHTjoint)

    # Key Hold Check
    scrolling = False #Check for Ctrl hold for use with scrolling
    objectResize = False

    RADIUS = 15
    INDICATOR_RADIUS = 15
    
    # Wind Variables
    target_point = None
    wind_force_duration = 0
    elapsed_wind_time = 0
    timer = 0
    wind_incrementer = 0.899

    # Main Simulation Loop
    while run:
        # Get Events

        # For each point

        LEFTjoint.position = (-50,WAVEHEIGHT)
        RIGHTjoint.position = (WIDTH+50,WAVEHEIGHT)

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

                # Toggle Help
                if event.type == pg.KEYDOWN and event.key == pg.K_h:
                    if HELP:
                        HELP = False
                    else:
                        HELP = True

                # Toggle Modes
                if event.type == pg.KEYDOWN and event.key == pg.K_w:
                    if WaveMode:
                        WaveMode = False
                    else:
                        WaveMode = True
                    print(f'{WaveMode=}')
                if event.type == pg.KEYDOWN and event.key == pg.K_e:
                    if WindMode:
                        WindMode = False
                    else:
                        WindMode = True
                    print(f'{WindMode=}')

                # Drop Object on Mouse Click When In Object Mode
                if not WaveMode:
                    if event.type == pg.MOUSEBUTTONDOWN and event.button == 1:
                        # createObject(spaceObj, event.pos, currentRadius, currentMass, currentElasticity, currentFriction)
                        RADIUS = INDICATOR_RADIUS
                        ball = createObject(space, event.pos, RADIUS, 6670 * RADIUS)
                        objects.append(ball)

                    # Adjust Density on Ctrl+Scroll When In Object Mode
                    if event.type == pg.KEYDOWN and event.key == pg.K_LCTRL or scrolling:
                        scrolling = True
                        objectResize = True
                        if event.type == pg.MOUSEWHEEL:
                            if scrolling:
                                if INDICATOR_RADIUS < 21:
                                    if event.y > 0:
                                        INDICATOR_RADIUS += 2
                                        print("Density Increase")
                                        print(INDICATOR_RADIUS)
                                if INDICATOR_RADIUS > 7:
                                    if event.y < 0:
                                        print("Density Decrease")
                                        INDICATOR_RADIUS -= 2
                                        print(INDICATOR_RADIUS)

                elif WaveMode:
                    if event.type == pg.KEYDOWN and event.key == pg.K_LCTRL or scrolling:
                        scrolling = True
                        if event.type == pg.MOUSEWHEEL:
                            if scrolling:
                                    if event.y < 0:
                                        if wind_incrementer < WINDMAX:
                                            wind_incrementer += 0.1
                                            print("Wave Intensity Decrease")
                                            print(wind_incrementer)
                                    if event.y > 0:
                                        if wind_incrementer > 0.6:
                                            print("Wave Intensity Increase")
                                            wind_incrementer -= 0.1
                                            print(wind_incrementer)
                                        

                if event.type == pg.KEYUP and event.key == pg.K_LCTRL:
                    scrolling = False
                    objectResize = False
                    print(f'{scrolling=}')


                # Adjust Wind Speed on Ctrl+Scroll When In Wave Mode

                # Handle Pausing/Playing The Simulation
                if event.type == pg.KEYDOWN and event.key == pg.K_SPACE:
                    if paused:
                        paused = False
                    else:
                        paused = True


                # Check for Wave Drag in Wave Mode
                if WaveMode:
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

        # Check for Pause/Play State
        if not paused:

            # Apply damping force to submerged objects
            for obj in objects:
                if is_submerged(obj, wave):
                    if obj.body.velocity[1] > 50 or obj.body.velocity[1] < 0:
                        obj.body.velocity *= 0.96

                vx,vy = obj.body.velocity
                obj.body.velocity = pm.Vec2d(vx*0.96, vy)
                    


            # Add Wind Force
            if WindMode:
                hist, bins  = np.histogram(x_coords, bins=20)

                if timer >= len(bins) - 2:
                    timer = -1

                if target_point is None or elapsed_wind_time >= wind_force_duration:
                    
                    timer += 1

                    # Choose target point sequentially from random bins in according to where the clock is
                    # to simulate wind in one direction
                    target_points = np.random.uniform(bins[timer], bins[timer+1])
                    index = np.abs(x_coords - target_points).argmin()
                    target_point = wave[index]

                    wind_force_duration = 15  # Adjust the duration range (2 to 5 seconds) as needed
                    elapsed_wind_time = 0

                # Apply wind force to the target point
                wind_strength = np.random.choice(np.linspace(-100 ,-10000))
                apply_wind_force(wave, wind_strength, target_point=target_point,
                                wind_force_duration=wind_force_duration, elapsed_wind_time=elapsed_wind_time)
            

            # Updates
            elapsed_wind_time += wind_incrementer
            draw(space, spaceObj, window, draw_options, wave, objects, WaveHeight, RADIUS, INDICATOR_RADIUS, objectResize, wind_incrementer)
            space.step(dt)
            spaceObj.step(dt)
            clock.tick(fps)

    pg.quit()


if __name__ == "__main__":
    main(window, WIDTH, HEIGHT)
