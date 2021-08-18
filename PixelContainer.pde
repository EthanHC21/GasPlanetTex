import java.util.LinkedList;

class PixelContainer {
  
  LinkedList<Integer> particleColors = new LinkedList<Integer>();
  
  color currentColor;
  
  public PixelContainer(color c) {
    currentColor = c;
  }
  
  void addColor(int pColor) {
    particleColors.add(pColor);
  }
  
  color drawPixel(int x, int y) {
    
    if (particleColors.size() > 0) {
      
      int r = 0;
      int g = 0;
      int b = 0;
      for (int i : particleColors) {
        r += red(i);
        g += green(i);
        b += blue(i);
      }
      
      currentColor = color(r/particleColors.size(), g/particleColors.size(), b/particleColors.size());
      
      // delete the LinkedList
      particleColors = new LinkedList<Integer>();
      
      return currentColor;
      
    } else {
      return currentColor;
    }
  }
  
  color getCurrentColor() {
     return currentColor;
  }
}
