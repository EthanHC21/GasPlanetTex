import java.util.LinkedList;

PGraphics pg;

// simulation properties
float noiseRad = 5; // radius used to calculate noise (essentially noise scale)
float dynRad = 50; // radius used by numerical integrator
float inc = .01; // w incremement
float fieldStrength = 2; // strength of the velocity field
int numParticles = 1000000; // figure it out
float colorAttenuation = .5; // factor by which the difference between a particle's color and the color of its current cell is multiplied
// 0 means they keep their previous color completely, and 1 makes it not do anything
float outerAttenuation = 0; // factor by which the difference between a particle's color and the color of an outer cell is multiplied

// band parameters
int bandVel = 25; // velocity of the bands
int numBands = 7; // number of bands
int bandExp = 1; // exponent of sin function used in band calculation
float minBands = 0.0; // minimum band blend factor
float minNoise = .5; // minimum noise blend factor
float minBandFac = 0.0; // minimum absolute output of the band sinusoidal function

float blend;

// grid discretization settings
int scl = 1;
int rows, cols;

// store frames so we don't print every iteration
int frames = 0;

// particle position and previous position array
PVector[] pos, prevPos;

// array holding grid velocities (note these are in what I'm coining "Surface Coordinates™"
PVector[][] gridVels;

// simplex wrapper for noise generation
SimplexWrapper noise = new SimplexWrapper();

// time coordinate for noise
float w = 0;

// timestep
float dt = .01;

// boolean to store whether we draw particles or not
boolean drawParticles = false;
// boolean to store whether we draw lines or not
boolean drawLines = false;
// boolean to store whether we draw from pixel containers
boolean drawContainers = true;

// reference image for bands
PImage img;
// array to hold the particle colors
color[] pColors = new color[numParticles];

// array of pixelContainers for drawing
PixelContainer[][] containers;

void setup() {
  
  //size(400, 400, P3D);
  pg = createGraphics(4096, 2048);
  size(4096, 2048);
  
  background(0);
  
  containers = new PixelContainer[height][width];
  
  // create particle position array
  pos = new PVector[numParticles];
  // create particle previous position array
  prevPos = new PVector[numParticles];
  
  // check grid size is even multiple of dimensions
  if (!(width % scl == 0 && rows % scl == 0)) {
    print("Error: Discretization not square.");
    exit();
  }
  
  // check that there is actually room for variation between bands and noise
  if (minBands + minNoise <= 1) {
    blend = 1 - (minBands + minNoise);
    //println("Blend is " + blend);
  } else {
    print("Error: Minimums add to greater than 1.");
    exit();
  }
  
  // get the reference image
  img = loadImage("ripplesbrightblur.png");
  // resize it
  img.resize(width/2, height);
  // made a new one to hold the mirrored image
  PImage fullImage = createImage(width, height, RGB);
  // now create a new array to hold all the pixels of the mirrored image
  int[] fullPixels = new int[width*height];
  for (int i = 0; i < width*height; i++) {
    int row = floor(i / width);
    if (i - row*width < width/2) {
      fullImage.pixels[i] = img.pixels[row  * width/2 + i % width];
    } else {
      fullImage.pixels[i] = img.pixels[row * width/2 + width/2 - (i % width - width/2) - 1];
    }
  }
  // save the pixels
  fullImage.updatePixels();
  
  // set this as background image
  image(fullImage, 0, 0);
  
  // update the display window's pixels
  loadPixels();
  
  // loop through pixels all pixels and draw the appropriate rectangles
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      containers[i][j] = new PixelContainer(pixels[i * width + j]);
      // test it by making it black
      //containers[i][j] = new PixelContainer(color(0, 0, 0));
      color temp = containers[i][j].drawPixel(j, i);
      noStroke();
      fill(temp);
      rect(j, i, 1, 1);
    }
  }
  
  //println("Length of thing is " + img.pixels.length);
  
  // initialize the particles and their velocities
  for (int i = 0; i < pos.length; i++) {
    pos[i] = new PVector(dynRad, random(-PI, PI), random(0, PI));
    PVector posEqui = sphere2Equi(pos[i]);
    int index = int(floor(posEqui.y * height) * width + floor(posEqui.x * width));
    //println("i is " + i);
    pColors[i] = pixels[index];
  }
  
  // initialize rows and columns
  rows = int(height/scl);
  cols = int(width/scl);
  // initialize array of grid velocities
  gridVels = new PVector[rows][cols];
  updateVelocityField();
  addBands();
  
  println("Velocity field generated");
  
  println("Setup complete");
  
  for (int j = 0; j < 50; j++) {
    removeDivergence();
  }
  
  
  colorMode(RGB);
  
}

void draw() {
  
  if (drawParticles) {
    background(255);
  }
  
  // calculate pair forces
  for (int i = 0; i < pos.length; i++) {
    
    // get position of particle i in rect
    PVector posRect = sphere2Rect(pos[i]);
    
    // save previous position for plotting purposes
    prevPos[i] = pos[i].copy();
    
    // perform verlet algorithm on particle i because we have completed its calculation
    
    // get particle's tangential velocity
    PVector tanVel = lookupTanVelocity(pos[i]).mult(fieldStrength);
    
    // calculate the radial velocity necessary to keep the particle on the surface of the sphere
    PVector radVel = getRadUnitVec(pos[i]).mult(-1 * tanVel.magSq()/this.dynRad);
    // add these two for total velocity
    PVector velTot = tanVel.copy().add(radVel);
    
    // multiply velocity by timestep and update the position
    posRect.add(velTot.mult(dt));
    
    // store this new position
    pos[i] = rect2SphereFix(posRect, dynRad);
    
    // draw the particle if that's what we decided
    if (drawParticles) {
      
      //println("Color is " + pColors[i]);
      
      stroke(pColors[i]);
      strokeWeight(2);
      
      // get equi coords
      PVector equiPos = sphere2Equi(pos[i]);
      point(equiPos.x * width, equiPos.y * height);
    }
    
    // draw the lines if that's what we decide
    if (drawLines) {
      float thetaDisp = pos[i].y - prevPos[i].y;
      float thetaVel = getThetaVel(pos[i], tanVel);
      
      stroke(pColors[i], 5);
      strokeWeight(1);
      
      // get positions
      PVector pos1 = sphere2Equi(prevPos[i]);
      PVector pos2 = sphere2Equi(pos[i]);
      
      // check theta velocity and displacement are of the same sign
      if ((thetaDisp > 0 && thetaVel > 0) || (thetaDisp < 0 && thetaVel < 0)) {
        
        line(pos1.x * width, pos1.y * height, pos2.x * width, pos2.y * height);
        
      } else {
        
        //stroke(255, 0, 0);
        float yMid;
        
        // make sure none are at the borders exactly
        if (!(pos1.x == 0 || pos1.x == 1 || pos2.x == 0 || pos2.x == 1)) {
          if (thetaVel > 0) { // right side
            
            // get y coord of crossing
            yMid = pos1.y + (pos2.y - pos1.y)/(pos2.x + 1 - pos1.x) * (1 - pos1.x);
            
            line(pos1.x * width, pos1.y * height, width, yMid * height);
            line(0, yMid * height, pos2.x * width, pos2.y * height);
            
          } else { // left side
          
            // get y coord of crossing
            yMid = pos1.y + (pos2.y - pos1.y)/(pos2.x - 1 - pos1.x) * (-1 * pos1.x);
            
            line(pos1.x * width, pos1.y * height, 0, yMid * height);
            line(width, yMid * height, pos2.x * width, pos2.y * height);
            
          }
        }
      }
    }
    
    
    // need equi
    if (drawContainers) {
      // get equi coords
    PVector equiPos = sphere2Equi(pos[i]);
      // edit the color of the particles to be similar to the pixel they're in
      color cellCol = containers[int(equiPos.y * (height - 1))][int(equiPos.x * (width - 1))].getCurrentColor(); //<>//
      float r_c = red(cellCol);
      float g_c = green(cellCol);
      float b_c = blue(cellCol);
      color partCol = pColors[i];
      float r_p = red(partCol);
      float g_p = green(partCol);
      float b_p = blue(partCol);
      // calculate new color
      int r_new = int(r_p + (r_c - r_p) * colorAttenuation);
      int g_new = int(g_p + (g_c - g_p) * colorAttenuation);
      int b_new = int(b_p + (b_c - b_p) * colorAttenuation);
      // set new color
      pColors[i] = color(r_new, g_new, b_new);
      
      // get the center cell that we're drawing in
      int y = int(equiPos.y * (height - 1));
      int x = int(equiPos.x * (width - 1));
      // indices
      int xr, xl, xd, xu, xur, xul, xdr, xdl, yr, yl, yd, yu, yur, yul, ydr, ydl;
      // check if we're on the corners
      if (x == 0 && y == 0) {
        
        // x coords
        xr = x + 1;
        xl = width - 1;
        xd = x;
        xu = width/2;
        xur = width/2 - 1;
        xul = width/2 + 1;
        xdr = x + 1;
        xdl = width - 1;
        
        // y coords
        yr = y;
        yl = y;
        yd = y + 1;
        yu = y;
        yur = y;
        yul = y;
        ydr = y + 1;
        ydl = y + 1;
        
      } else if (x == width - 1 && y == 0) {
        
        // x coords
        xr = 0;
        xl = x - 1;
        xd = x;
        xu = width/2 - 1;
        xur = width/2 - 2;
        xul = width/2;
        xdr = 0;
        xdl = x - 1;
        
        // y coords
        yr = y;
        yl = y;
        yd = y + 1;
        yu = y;
        yur = y;
        yul = y;
        ydr = y + 1;
        ydl = y + 1;
        
      } else if (x == 0 && y == height - 1) {
        
        // x coords
        xr = x + 1;
        xl = width - 1;
        xd = width/2;
        xu = x;
        xur = x + 1;
        xul = width - 1;
        xdr = width/2 - 1;
        xdl = width/2 + 1;
        
        // y coords
        yr = y;
        yl = y;
        yd = y;
        yu = y - 1;
        yur = y - 1;
        yul = y - 1;
        ydr = y;
        ydl = y;
        
      } else if (x == width - 1 && y == height - 1) {
        
        // x coords
        xr = 0;
        xl = x - 1;
        xd = width/2 - 1;
        xu = x;
        xur = 0;
        xul = x - 1;
        xdr = width/2 - 2;
        xdl = width/2 - 1;
        
        // y coords
        yr = y;
        yl = y;
        yd = y;
        yu = y - 1;
        yur = y - 1;
        yul = y - 1;
        ydr = y;
        ydl = y;
        
      } else if (y == 0) { // top edge
        
        // x coords
        xr = x + 1;
        xl = x - 1;
        xd = x;
        xu = x + width/2;
        xur = xu - 1;
        xul = xu + 1;
        xdr = x + 1;
        xdl = x - 1;
        if (xu >= width) xu -= width;
        if (xur >= width) xur -= width;
        if (xul >= width) xul -= width;
        
        // y coords
        yr = y;
        yl = y;
        yd = y + 1;
        yu = y;
        yur = y;
        yul = y;
        ydr = y + 1;
        ydl = y + 1;
        
      } else if (y == height - 1) { // bottom edge
        
        // x coords
        xr = x + 1;
        xl = x - 1;
        xd = x + width/2;
        xu = x;
        xur = x + 1;
        xul = x - 1;
        xdr = xd - 1;
        xdl = xd + 1;
        if (xd >= width) xd -= width; 
        if (xdr >= width) xdr -= width;
        if (xdl >= width) xdl -= width;
        
        // y coords
        yr = y;
        yl = y;
        yd = y;
        yu = y - 1;
        yur = y - 1;
        yul = y - 1;
        ydr = y;
        ydl = y;
        
      } else {
        
        if (x == 0) {
          
          // x coords
          xr = x + 1;
          xl = width - 1;
          xd = x;
          xu = x;
          xur = x + 1;
          xul = width - 1;
          xdr = x + 1;
          xdl = width - 1;
          
        } else if (x == width - 1) {
          
          // x coords
          xr = 0;
          xl = x - 1;
          xd = x;
          xu = x;
          xur = 0;
          xul = x - 1;
          xdr = 0;
          xdl = x - 1;
          
        } else {
          
          // x coords
          xr = x + 1;
          xl = x - 1;
          xd = x;
          xu = x;
          xur = x + 1;
          xul = x - 1;
          xdr = x + 1;
          xdl = x - 1;
          
        }
        
        // y coords
        yr = y;
        yl = y;
        yd = y + 1;
        yu = y - 1;
        yur = y - 1;
        yul = y - 1;
        ydr = y + 1;
        ydl = y + 1;
        
      }
      
      // add the particle's new color to the new cell
      containers[y][x].addColor(pColors[i]);
      /*
      // set the outer color
      int r_o = int(r_p + (r_c - r_p) * outerAttenuation);
      int g_o = int(g_p + (g_c - g_p) * outerAttenuation);
      int b_o = int(b_p + (b_c - b_p) * outerAttenuation);
      color outerColor = color(r_o, g_o, b_o);
      
      containers[yl][xl].addColor(outerColor);
      containers[yr][xr].addColor(outerColor);
      containers[yu][xu].addColor(outerColor);
      containers[yd][xd].addColor(outerColor);
      containers[yul][xul].addColor(outerColor);
      containers[yur][xur].addColor(outerColor);
      containers[ydl][xdl].addColor(outerColor);
      containers[ydr][xdr].addColor(outerColor);
      */
    }
    
    
  }
  
  // advance noise field time
  //w += inc;
  
  // update the velocity field
  //updateVelocityField();
  
  // loop through pixels all pixels and draw the appropriate rectangles
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      color temp = containers[i][j].drawPixel(j, i);
      noStroke();
      fill(temp);
      rect(j, i, 1, 1);
    }
  }
  
}

// convert from sphere to rectangular
PVector sphere2Rect(PVector sphere) {
  
  float x = sphere.x * sin(sphere.z) * cos(sphere.y);
  float y = sphere.x * sin(sphere.z) * sin(sphere.y);
  float z = sphere.x * cos(sphere.z);
  
  return new PVector(x, y, z);
}

// convert from rect to sphere but fix the radius
PVector rect2SphereFix(PVector rect, float newRad) {
  return new PVector(newRad, atan2(rect.y,rect.x), acos(rect.z/sqrt(sq(rect.x) + sq(rect.y) + sq(rect.z))));
}

// convert spherical to equirectangular
PVector sphere2Equi(PVector sphere) {

  // map the coords so they can be multiplied later
  return new PVector(sphere.y/TWO_PI + .5, sphere.z/PI);
}

// get the unit vector in the radial direction
PVector getRadUnitVec(PVector spherePos) {
  
  return new PVector(sin(spherePos.z) * cos(spherePos.y), sin(spherePos.z) * sin(spherePos.y), cos(spherePos.z)).normalize();
}

// get velocity in terms of spherical
float getThetaVel(PVector spherePos, PVector rectVel) {
  
  // get theta unit vector
  float theta = spherePos.y;
  // normalize the vector to decrease floating point error
  PVector e_theta = new PVector(-1 * sin(theta), cos(theta), 0).normalize();
  
  return e_theta.dot(rectVel);
  
}

PVector equi2NoiseSphere(float x, float y) {
  
  float theta = TWO_PI * (x - .5);
  float phi = PI * y;
  
  return new PVector(noiseRad, theta, phi);
  
}

PVector equi2Sphere(float x, float y) {
  
  float theta = TWO_PI * (x - .5);
  float phi = PI * y;
  
  return new PVector(dynRad, theta, phi);
  
}

PVector rectVec2SurfaceVec(PVector rectVec, PVector spherePos) {
  
  // get spherical unit vectors
  float phi = spherePos.z;
  float theta = spherePos.y;
  PVector e_theta = new PVector(-1 * sin(theta), cos(theta), 0).normalize();
  PVector e_phi = new PVector(cos(phi) * cos(theta), cos(phi) * sin(theta), -1 * sin(phi)).normalize();
  
  // get the dot product of each. that's the definition of Surface Coordinates™
  return new PVector(e_theta.dot(rectVec), e_phi.dot(rectVec));
  
}

PVector surfaceVec2RectVec(PVector surfVec, float theta, float phi) {
  
  // get spherical unit vectors
  PVector e_theta = new PVector(-1 * sin(theta), cos(theta), 0).normalize();
  PVector e_phi = new PVector(cos(phi) * cos(theta), cos(phi) * sin(theta), -1 * sin(phi)).normalize();
  
  // just multiply each unit vector by its respective component in surface
  return e_theta.mult(surfVec.x).add(e_phi.mult(surfVec.y));
  
}

// update the velocity field with new values from the noise
void updateVelocityField() {
  
  // loop through grid velocities array
  for (int y = 0;  y < rows; y++){
    for (int x = 0; x < cols; x++){
      
      // get normalized equirectangular coordinate of the center of each cell
      float xEqui = (float(x * scl) + float(scl)/2)/width;
      float yEqui = (float(y * scl) + float(scl)/2)/height;
      
      // convert these to coordinates on the noise sphere
      PVector noisePos = equi2NoiseSphere(xEqui, yEqui);
      // convert this to rectangular coordinates for the noise function
      PVector noiseRect = sphere2Rect(noisePos);
      
      // get the noise in rectangular
      PVector vel = noise.getUCCurl4D(noiseRect.x, noiseRect.y, noiseRect.z, w).mult(fieldStrength);
      
      // get phi to attenuate velocity
      float phi = PI * yEqui;
      
      // convert this to Surface Coordinates™ and add to the array
      gridVels[y][x] = rectVec2SurfaceVec(vel, noisePos).mult(sin(phi));
      //gridVels[y][x] = new PVector(0, 0, 0);
    }
  }
  
}

// return the velocity from the lookup grid velocities based on the particle's position
PVector lookupTanVelocity(PVector spherePos) {
  
  // get position of particle i in equi
  PVector posEqui = sphere2Equi(spherePos);
  // get the cell that the particle is in
  int x = floor(posEqui.x * width/scl);
  int y = floor(posEqui.y * height/scl);
  
  if (x == cols) {
    x = 0;
  } else if (y == rows) {
    y = 0;
  }
  
  // get the velocity in rectangular coordinates
  return surfaceVec2RectVec(gridVels[y][x], spherePos.y, spherePos.z);
  
}

// remove divergence from velocity field
void removeDivergence() {
  
  // store the new tangential velocities
  PVector[][] reducedGridVels = new PVector[rows][cols];
  
  PVector[] tempConvolution = new PVector[9];
  
  // get pixel scale
  //float p = height/PI;
  // get epsilon in the y (phi)
  //float epsilon_y = dynRad * p;
  
  //println("removeDivergence begin");
  
  // loop through each cell
  for (int y = 0; y < rows; y++){
    
    // get the epsilon in the x (theta)
    //float epsilon_x = p * dynRad * sin((float(y * scl) + float(scl)/2.0) * PI);
    
    for(int x = 0; x < cols; x++){
      
      //fill(255, 0, 0);
      //rect(x * scl, y * scl, scl, scl);
      
      // indices
      int xr, xl, xd, xu, xur, xul, xdr, xdl, yr, yl, yd, yu, yur, yul, ydr, ydl;
      
      // define and calculate offset variables for x and y
      // u = up, d = down, r = right, l = left, combinations are acceptable
      
      // check if we're on the corners
      if (x == 0 && y == 0) {
        
        // x coords
        xr = x + 1;
        xl = cols - 1;
        xd = x;
        xu = cols/2;
        xur = cols/2 - 1;
        xul = cols/2 + 1;
        xdr = x + 1;
        xdl = cols - 1;
        
        // y coords
        yr = y;
        yl = y;
        yd = y + 1;
        yu = y;
        yur = y;
        yul = y;
        ydr = y + 1;
        ydl = y + 1;
        
        // get vectors
        tempConvolution[0] = gridVels[yul][xul].copy();
        tempConvolution[1] = gridVels[yu][xu].copy();
        tempConvolution[2] = gridVels[yur][xur].copy();
        tempConvolution[3] = gridVels[yl][xl].copy();
        tempConvolution[4] = gridVels[y][x].copy();
        tempConvolution[5] = gridVels[yr][xr].copy();
        tempConvolution[6] = gridVels[ydl][xdl].copy();
        tempConvolution[7] = gridVels[yd][xd].copy();
        tempConvolution[8] = gridVels[ydr][xdr].copy();
        
        // mirror the ones that flipped over the pole
        tempConvolution[0].mult(-1);
        tempConvolution[1].mult(-1);
        tempConvolution[2].mult(-1);
        
      } else if (x == cols - 1 && y == 0) {
        
        // x coords
        xr = 0;
        xl = x - 1;
        xd = x;
        xu = cols/2 - 1;
        xur = cols/2 - 2;
        xul = cols/2;
        xdr = 0;
        xdl = x - 1;
        
        // y coords
        yr = y;
        yl = y;
        yd = y + 1;
        yu = y;
        yur = y;
        yul = y;
        ydr = y + 1;
        ydl = y + 1;
        
        // get vectors
        tempConvolution[0] = gridVels[yul][xul].copy();
        tempConvolution[1] = gridVels[yu][xu].copy();
        tempConvolution[2] = gridVels[yur][xur].copy();
        tempConvolution[3] = gridVels[yl][xl].copy();
        tempConvolution[4] = gridVels[y][x].copy();
        tempConvolution[5] = gridVels[yr][xr].copy();
        tempConvolution[6] = gridVels[ydl][xdl].copy();
        tempConvolution[7] = gridVels[yd][xd].copy();
        tempConvolution[8] = gridVels[ydr][xdr].copy();
        
        // mirror the ones that flipped over the pole
        tempConvolution[0].mult(-1);
        tempConvolution[1].mult(-1);
        tempConvolution[2].mult(-1);
        
      } else if (x == 0 && y == rows - 1) {
        
        // x coords
        xr = x + 1;
        xl = cols - 1;
        xd = cols/2;
        xu = x;
        xur = x + 1;
        xul = cols - 1;
        xdr = cols/2 - 1;
        xdl = cols/2 + 1;
        
        // y coords
        yr = y;
        yl = y;
        yd = y;
        yu = y - 1;
        yur = y - 1;
        yul = y - 1;
        ydr = y;
        ydl = y;
        
        // get vectors
        tempConvolution[0] = gridVels[yul][xul].copy();
        tempConvolution[1] = gridVels[yu][xu].copy();
        tempConvolution[2] = gridVels[yur][xur].copy();
        tempConvolution[3] = gridVels[yl][xl].copy();
        tempConvolution[4] = gridVels[y][x].copy();
        tempConvolution[5] = gridVels[yr][xr].copy();
        tempConvolution[6] = gridVels[ydl][xdl].copy();
        tempConvolution[7] = gridVels[yd][xd].copy();
        tempConvolution[8] = gridVels[ydr][xdr].copy();
        
        // mirror the ones that flipped over the pole
        tempConvolution[6].mult(-1);
        tempConvolution[7].mult(-1);
        tempConvolution[8].mult(-1);
        
      } else if (x == cols - 1 && y == rows - 1) {
        
        // x coords
        xr = 0;
        xl = x - 1;
        xd = cols/2 - 1;
        xu = x;
        xur = 0;
        xul = x - 1;
        xdr = cols/2 - 2;
        xdl = cols/2 - 1;
        
        // y coords
        yr = y;
        yl = y;
        yd = y;
        yu = y - 1;
        yur = y - 1;
        yul = y - 1;
        ydr = y;
        ydl = y;
        
        // get vectors
        tempConvolution[0] = gridVels[yul][xul].copy();
        tempConvolution[1] = gridVels[yu][xu].copy();
        tempConvolution[2] = gridVels[yur][xur].copy();
        tempConvolution[3] = gridVels[yl][xl].copy();
        tempConvolution[4] = gridVels[y][x].copy();
        tempConvolution[5] = gridVels[yr][xr].copy();
        tempConvolution[6] = gridVels[ydl][xdl].copy();
        tempConvolution[7] = gridVels[yd][xd].copy();
        tempConvolution[8] = gridVels[ydr][xdr].copy();
        
        // mirror the ones that flipped over the pole
        tempConvolution[6].mult(-1);
        tempConvolution[7].mult(-1);
        tempConvolution[8].mult(-1);
        
      } else if (y == 0) { // top edge
        
        // x coords
        xr = x + 1;
        xl = x - 1;
        xd = x;
        xu = x + cols/2;
        xur = xu - 1;
        xul = xu + 1;
        xdr = x + 1;
        xdl = x - 1;
        if (xu >= cols) xu -= cols;
        if (xur >= cols) xur -= cols;
        if (xul >= cols) xul -= cols;
        
        // y coords
        yr = y;
        yl = y;
        yd = y + 1;
        yu = y;
        yur = y;
        yul = y;
        ydr = y + 1;
        ydl = y + 1;
        
        // get vectors
        tempConvolution[0] = gridVels[yul][xul].copy();
        tempConvolution[1] = gridVels[yu][xu].copy();
        tempConvolution[2] = gridVels[yur][xur].copy();
        tempConvolution[3] = gridVels[yl][xl].copy();
        tempConvolution[4] = gridVels[y][x].copy();
        tempConvolution[5] = gridVels[yr][xr].copy();
        tempConvolution[6] = gridVels[ydl][xdl].copy();
        tempConvolution[7] = gridVels[yd][xd].copy();
        tempConvolution[8] = gridVels[ydr][xdr].copy();
        
        // mirror the ones that flipped over the pole
        tempConvolution[0].mult(-1);
        tempConvolution[1].mult(-1);
        tempConvolution[2].mult(-1);
        
      } else if (y == rows - 1) { // bottom edge
        
        // x coords
        xr = x + 1;
        xl = x - 1;
        xd = x + cols/2;
        xu = x;
        xur = x + 1;
        xul = x - 1;
        xdr = xd - 1;
        xdl = xd + 1;
        if (xd >= cols) xd -= cols; 
        if (xdr >= cols) xdr -= cols;
        if (xdl >= cols) xdl -= cols;
        
        // y coords
        yr = y;
        yl = y;
        yd = y;
        yu = y - 1;
        yur = y - 1;
        yul = y - 1;
        ydr = y;
        ydl = y;
        
        // get vectors
        tempConvolution[0] = gridVels[yul][xul].copy();
        tempConvolution[1] = gridVels[yu][xu].copy();
        tempConvolution[2] = gridVels[yur][xur].copy();
        tempConvolution[3] = gridVels[yl][xl].copy();
        tempConvolution[4] = gridVels[y][x].copy();
        tempConvolution[5] = gridVels[yr][xr].copy();
        tempConvolution[6] = gridVels[ydl][xdl].copy();
        tempConvolution[7] = gridVels[yd][xd].copy();
        tempConvolution[8] = gridVels[ydr][xdr].copy();
        
        // mirror the ones that flipped over the pole
        tempConvolution[6].mult(-1);
        tempConvolution[7].mult(-1);
        tempConvolution[8].mult(-1);
        
      } else {
        
        if (x == 0) {
          
          // x coords
          xr = x + 1;
          xl = cols - 1;
          xd = x;
          xu = x;
          xur = x + 1;
          xul = cols - 1;
          xdr = x + 1;
          xdl = cols - 1;
          
        } else if (x == cols - 1) {
          
          // x coords
          xr = 0;
          xl = x - 1;
          xd = x;
          xu = x;
          xur = 0;
          xul = x - 1;
          xdr = 0;
          xdl = x - 1;
          
        } else {
          
          // x coords
          xr = x + 1;
          xl = x - 1;
          xd = x;
          xu = x;
          xur = x + 1;
          xul = x - 1;
          xdr = x + 1;
          xdl = x - 1;
          
        }
        
        // y coords
        yr = y;
        yl = y;
        yd = y + 1;
        yu = y - 1;
        yur = y - 1;
        yul = y - 1;
        ydr = y + 1;
        ydl = y + 1;
        
        tempConvolution[0] = gridVels[yul][xul].copy();
        tempConvolution[1] = gridVels[yu][xu].copy();
        tempConvolution[2] = gridVels[yur][xur].copy();
        tempConvolution[3] = gridVels[yl][xl].copy();
        tempConvolution[4] = gridVels[y][x].copy();
        tempConvolution[5] = gridVels[yr][xr].copy();
        tempConvolution[6] = gridVels[ydl][xdl].copy();
        tempConvolution[7] = gridVels[yd][xd].copy();
        tempConvolution[8] = gridVels[ydr][xdr].copy();
        
      }
      
      // flip the top and bottom row to make the y direction consistent
      PVector[] convolution = new PVector[9];
      convolution[0] = tempConvolution[6];
      convolution[1] = tempConvolution[7];
      convolution[2] = tempConvolution[8];
      convolution[3] = tempConvolution[3];
      convolution[4] = tempConvolution[4];
      convolution[5] = tempConvolution[5];
      convolution[6] = tempConvolution[0];
      convolution[7] = tempConvolution[1];
      convolution[8] = tempConvolution[2];
      
      // apply the reduction algorithm
      
      PVector term1 = new PVector(1, 1).mult(new PVector(1, 1).dot(convolution[2].copy().add(convolution[6])));
      PVector term2 = new PVector(1, -1).mult(new PVector(1, -1).dot(convolution[0].copy().add(convolution[8])));
      PVector term3 = convolution[3].copy().add(convolution[5]).sub(convolution[1]).sub(convolution[7]).mult(2);
      term3.y *= -1;
      PVector term4 = convolution[4].mult(-4);
      
      PVector incDiv = term1.add(term2).add(term3).add(term4);
      incDiv.mult(1.0/16);
      //incDiv.x *= 1.0/epsilon_x;
      //incDiv.y *= 1.0/epsilon_y;
      
      // save this velocity with reduced divergence
      reducedGridVels[y][x] = incDiv.copy().add(gridVels[y][x]);
      
    }
  }
  
  //println("removeDivergence complete");
  
  // update gridvels
  gridVels = reducedGridVels;
  
}

void addBands() {
  
  float bandFac;
  float theta;
  float yEqui, xEqui;
  PVector e_theta;
  
  // loop through rows
  for (int y = 0; y < rows; y++) {
    
    // get normalized equirectangular coordinate of the center of each cell
    yEqui = (float(y * scl) + float(scl)/2)/height;
    
    if (bandExp % 2 == 0) {
      
      // get percentage of band velocity at given height
      bandFac = pow(-1, floor(numBands * yEqui)) * pow(sin(numBands * PI * yEqui), bandExp);
      
      // check it's within the allowed range
      if (abs(bandFac) < minBandFac) bandFac = bandFac/abs(bandFac) * minBandFac;
      
    } else {
      
      // get percentage of band velocity at given height
      bandFac = pow(sin(numBands * PI * yEqui), bandExp);
      
      // check it's within the allowed range
      if (abs(bandFac) < minBandFac) bandFac = bandFac/abs(bandFac) * minBandFac;
      
    }
    
    // now we get loop through all of the x values and add the velocity in accordance with the percentage
    for (int x = 0; x < cols; x++) {
      
      // get theta from normalized equi coordinate of the center of each cell
      xEqui = (float(x * scl) + float(scl)/2)/width;
      theta = TWO_PI * (xEqui - .5);
      // get theta unit vector
      e_theta = new PVector(-1 * sin(theta), cos(theta), 0).normalize();
      
      // get the band vector - multiply unit vector by magnitude
      PVector bandVec = e_theta.mult(bandVel);
      
      PVector surfaceVec = rectVec2SurfaceVec(bandVec, equi2Sphere(xEqui, yEqui)).mult(abs(bandFac)/bandFac);
      
      // perform the weighted average in the cell
      gridVels[y][x].mult(minNoise + blend * (1 - abs(bandFac))).add(surfaceVec.mult(minBands + (blend * abs(bandFac))));
    }
    
  }
  
}

// function to stop and save the image
void keyPressed() {
  if (key == ' ') {
    save("out.tiff");
    exit();
  }
}
