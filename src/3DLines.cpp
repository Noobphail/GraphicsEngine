//
// Created by ralph on 26/11/2020. ok
//

//draws line accounting for depth
void drawLine3D(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour col, float zBuffer[WIDTH][HEIGHT], float inverseZL, float inverseZR)
{
    std::vector<CanvasPoint> line = findLine(from, to);
    std::vector<float> inverseZVals = interpolateSingleFloats(inverseZL, inverseZR, line.size());
    assert(from.y == to.y);
    for (int i = 0; i < line.size(); i++)
    {
        if ((line[i].x >= 0 && line[i].x < window.width) && (line[i].y >= 0 && line[i].y < window.height))
        {
            if (inverseZVals[i] > zBuffer[(int)round(line[i].x)][(int)round(line[i].y)])
            {
                zBuffer[(int)round(line[i].x)][(int)round(line[i].y)] = inverseZVals[i];
                uint32_t packed = colToInt32(col);
                window.setPixelColour(line[i].x, line[i].y, packed);
            }
        }
    }
}

//draws textured line in 3D space (uses z-buffer)
void drawTexture3D(DrawingWindow &window, CanvasPoint from, CanvasPoint to, TextureMap m, float zBuffer[WIDTH][HEIGHT], float zL, float zR)
{
    std::vector<CanvasPoint> line = findLine(from, to);
    std::vector<float> textLineX = interpolateSingleFloats(from.texturePoint.x, to.texturePoint.x, line.size());
    std::vector<float> textLineY = interpolateSingleFloats(from.texturePoint.y, to.texturePoint.y, line.size());
    std::vector<float> zVals = interpolateSingleFloats(zL, zR, line.size());
    assert(from.y == to.y);
    for (int i = 0; i < line.size(); i++)
    {
        if (1 / zVals[i] > zBuffer[(int)round(line[i].x)][(int)round(line[i].y)])
        {
            zBuffer[(int)round(line[i].x)][(int)round(line[i].y)] = 1 / zVals[i];
            TexturePoint temp = TexturePoint(textLineX[i], textLineY[i]);
            uint32_t packed = m.pixels[texturePointToIndex(temp, m)];
            window.setPixelColour(line[i].x, line[i].y, packed);
        }
    }
}

//draws perspective corrected textured line in 3D space (uses z-buffer)
void drawPerspectiveTexture3D(DrawingWindow &window, CanvasPoint from, CanvasPoint to, TextureMap m, float zBuffer[WIDTH][HEIGHT], float zL, float zR)
{
    std::vector<CanvasPoint> line = findLine(from, to);
    std::vector<float> cValuesX = interpolateSingleFloats(from.texturePoint.x / zL, to.texturePoint.x / zR, line.size());
    std::vector<float> cValuesY = interpolateSingleFloats(from.texturePoint.y / zL, to.texturePoint.y / zR, line.size());

    std::vector<float> inverseZVals = interpolateSingleFloats(1 / zL, 1 / zR, line.size());
    for (int i = 0; i < line.size(); i++)
    {
        if (inverseZVals[i] > zBuffer[(int)round(line[i].x)][(int)round(line[i].y)])
        {
            zBuffer[(int)round(line[i].x)][(int)round(line[i].y)] = inverseZVals[i];
            TexturePoint temp = TexturePoint((1 / inverseZVals[i]) * cValuesX[i], (1 / inverseZVals[i]) * cValuesY[i]);
            uint32_t packed = m.pixels[texturePointToIndex(temp, m)];
            window.setPixelColour(line[i].x, line[i].y, packed);
        }
    }
}