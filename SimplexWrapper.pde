class SimplexWrapper {
  
  OpenSimplexNoise simplex;
  OpenSimplexNoise simplexX;
  OpenSimplexNoise simplexY;
  OpenSimplexNoise simplexZ;
  
  float delta = .01; // TODO: test different values of this for gradient and curl
  float octaves = 1;
  
  public SimplexWrapper() {
    simplex = new OpenSimplexNoise();
    
    long xSeed = floor(2147483647L * random(1));
    simplexX = new OpenSimplexNoise(xSeed);
    long ySeed = floor(2147483647L * random(1));
    simplexY = new OpenSimplexNoise(ySeed);
    long zSeed = floor(2147483647L * random(1));
    simplexZ = new OpenSimplexNoise(zSeed);
    
  }
  
  float getNoise3D(float x, float y, float z){
    float total = 0;
    float frequency = 1;
    float amplitude = 1;
    float lacunarity = 1.5;
    float persistence = .5;
    float totalAmplitude = 0;
    
    for (int i = 0; i < octaves; i++) {
      total += (float) this.simplex.eval(x * frequency, y * frequency, z * frequency) * amplitude;
      totalAmplitude += amplitude;
      amplitude *= persistence;
      frequency *= lacunarity;
    }
    
    return total/totalAmplitude;
  }
  
  // get 4D noise
  float getNoise4D(float x, float y, float z, float w){
    float total = 0;
    float frequency = 1;
    float amplitude = 1;
    float lacunarity = 1.5;
    float persistence = .5;
    float totalAmplitude = 0;
    
    for (int i = 0; i < octaves; i++) {
      total += (float) this.simplex.eval(x * frequency, y * frequency, z * frequency, w * frequency) * amplitude;
      totalAmplitude += amplitude;
      amplitude *= persistence;
      frequency *= lacunarity;
    }
    
    return total/totalAmplitude;
  }
  
  // get uncorrelated 4D noise
  float getNoise4DX(float x, float y, float z, float w){
    float total = 0;
    float frequency = 1;
    float amplitude = 1;
    float lacunarity = 1.5;
    float persistence = .5;
    float totalAmplitude = 0;
    
    for (int i = 0; i < octaves; i++) {
      total += (float) this.simplexX.eval(x * frequency, y * frequency, z * frequency, w * frequency) * amplitude;
      totalAmplitude += amplitude;
      amplitude *= persistence;
      frequency *= lacunarity;
    }
    
    return total/totalAmplitude;
  }
  
  float getNoise4DY(float x, float y, float z, float w){
    float total = 0;
    float frequency = 1;
    float amplitude = 1;
    float lacunarity = 1.5;
    float persistence = .5;
    float totalAmplitude = 0;
    
    for (int i = 0; i < octaves; i++) {
      total += (float) this.simplexY.eval(x * frequency, y * frequency, z * frequency, w * frequency) * amplitude;
      totalAmplitude += amplitude;
      amplitude *= persistence;
      frequency *= lacunarity;
    }
    
    return total/totalAmplitude;
  }
  
  float getNoise4DZ(float x, float y, float z, float w){
    float total = 0;
    float frequency = 1;
    float amplitude = 1;
    float lacunarity = 1.5;
    float persistence = .5;
    float totalAmplitude = 0;
    
    for (int i = 0; i < octaves; i++) {
      total += (float) this.simplexZ.eval(x * frequency, y * frequency, z * frequency, w * frequency) * amplitude;
      totalAmplitude += amplitude;
      amplitude *= persistence;
      frequency *= lacunarity;
    }
    
    return total/totalAmplitude;
  }
  
  // gets a correlated gradient
  PVector getGrad4D(float x, float y, float z, float w){
    
    float gradX = (getNoise4D(x + delta, y, z, w) - getNoise4D(x - delta, y, z, w))/(2 * delta);
    float gradY = (getNoise4D(x, y + delta, z, w) - getNoise4D(x, y - delta, z, w))/(2 * delta);
    float gradZ = (getNoise4D(x, y, z + delta, w) - getNoise4D(x, y, z - delta, w))/(2 * delta);
    
    return new PVector(gradX, gradY, gradZ);
    
  }
  
  PVector getUCCurl4D(float x, float y, float z, float w) {
    
    float curlX = (getUCGrad4DZ(x, y + delta, z, w).z - getUCGrad4DZ(x, y - delta, z, w).z)/(2 * delta) - (getUCGrad4DY(x, y, z + delta, w).y - getUCGrad4DY(x, y, z - delta, w).y)/(2 * delta);
    float curlY = (getUCGrad4DX(x, y, z + delta, w).x - getUCGrad4DX(x, y, z - delta, w).x)/(2 * delta) - (getUCGrad4DZ(x + delta, y, z, w).z - getUCGrad4DZ(x - delta, y, z, w).z)/(2 * delta);
    float curlZ = (getUCGrad4DY(x + delta, y, z, w).y - getUCGrad4DY(x - delta, y, z, w).y)/(2 * delta) - (getUCGrad4DX(x, y + delta, z, w).x - getUCGrad4DX(x, y - delta, z, w).x)/(2 * delta);
    
    return new PVector(curlX, curlY, curlZ);
    
  }
  
  // get an uncorrelated gradients
  PVector getUCGrad4DX(float x, float y, float z, float w) {
    
    float gradX = (getNoise4DX(x + delta, y, z, w) - getNoise4DX(x - delta, y, z, w))/(2 * delta);
    float gradY = (getNoise4DX(x, y + delta, z, w) - getNoise4DX(x, y - delta, z, w))/(2 * delta);
    float gradZ = (getNoise4DX(x, y, z + delta, w) - getNoise4DX(x, y, z - delta, w))/(2 * delta);
    
    return new PVector(gradX, gradY, gradZ);
    
  }
  
  PVector getUCGrad4DY(float x, float y, float z, float w) {
    
    float gradX = (getNoise4DY(x + delta, y, z, w) - getNoise4DY(x - delta, y, z, w))/(2 * delta);
    float gradY = (getNoise4DY(x, y + delta, z, w) - getNoise4DY(x, y - delta, z, w))/(2 * delta);
    float gradZ = (getNoise4DY(x, y, z + delta, w) - getNoise4DY(x, y, z - delta, w))/(2 * delta);
    
    return new PVector(gradX, gradY, gradZ);
    
  }
  
  PVector getUCGrad4DZ(float x, float y, float z, float w) {
    
    float gradX = (getNoise4DZ(x + delta, y, z, w) - getNoise4DZ(x - delta, y, z, w))/(2 * delta);
    float gradY = (getNoise4DZ(x, y + delta, z, w) - getNoise4DZ(x, y - delta, z, w))/(2 * delta);
    float gradZ = (getNoise4DZ(x, y, z + delta, w) - getNoise4DZ(x, y, z - delta, w))/(2 * delta);
    
    return new PVector(gradX, gradY, gradZ);
    
  }
  
}
