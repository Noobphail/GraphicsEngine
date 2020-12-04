//
// Created by ralph on 26/11/2020.
//
#include <CanvasPoint.h>
#include <Colour.h>
#include <DrawingWindow.h>

//returns vector of points along a line between from and to
std::vector<CanvasPoint> findLine(CanvasPoint from, CanvasPoint to)
{
    std::vector<CanvasPoint> ret;
    float diffX = to.x - from.x;
    float diffY = to.y - from.y;
    float stepNum = fmax(1, fmax(abs(diffX), abs(diffY)));
    float stepSizeX = diffX / stepNum;
    float stepSizeY = diffY / stepNum;
    for (float i = 0.0; i <= stepNum; i++)
    {
        float x = from.x + (i * stepSizeX);
        float y = from.y + (i * stepSizeY);
        CanvasPoint pt = CanvasPoint(round(x), round(y));
        //could be more efficient with an early exit, but would need to check the line is going off screen not coming from offscreen onto the screen
        if (round(x) >= 0 && round(y) >= 0 && round(x) < WIDTH && round(y) < HEIGHT)
        {
            ret.push_back(pt);
        }
    }

    return ret;
}

//draws line in 2D
void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour col)
{
    std::vector<CanvasPoint> line = findLine(from, to);
    for (int i = 0; i < line.size(); i++)
    {
        uint32_t packed = colToInt32(col);
        window.setPixelColour(line[i].x, line[i].y, packed);
    }
}

//idk if this is necessary i was tired ok (it isn't ;-; )
std::vector<TexturePoint> findTextureLine(TexturePoint from, TexturePoint to)
{
    std::vector<TexturePoint> ret;
    float diffX = to.x - from.x;
    float diffY = to.y - from.y;
    float stepNum = fmax(abs(diffX), abs(diffY));
    float stepSizeX = diffX / stepNum;
    float stepSizeY = diffY / stepNum;
    for (float i = 0.0; i < stepNum; i++)
    {
        float x = from.x + (i * stepSizeX);
        float y = from.y + (i * stepSizeY);
        TexturePoint pt = TexturePoint(x, y);
        ret.push_back(pt);
    }
    return ret;
}

//draws textured line between from and to
void drawTexture(DrawingWindow &window, CanvasPoint from, CanvasPoint to, TextureMap m)
{
    std::vector<CanvasPoint> line = findLine(from, to);
    std::vector<float> textLineX = interpolateSingleFloats(from.texturePoint.x, to.texturePoint.x, line.size());
    std::vector<float> textLineY = interpolateSingleFloats(from.texturePoint.y, to.texturePoint.y, line.size());
    for (int i = 0; i < line.size(); i++)
    {
        TexturePoint temp = TexturePoint(textLineX[i], textLineY[i]);
        uint32_t packed = m.pixels[texturePointToIndex(temp, m)];
        window.setPixelColour(line[i].x, line[i].y, packed);
    }
}
