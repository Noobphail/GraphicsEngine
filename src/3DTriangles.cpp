//
// Created by ralph on 26/11/2020.
//

//helper function for fillTriangle3D. Fills bottom-flat triangle in 3D space
void fillBottomFlatTriangle3D(DrawingWindow &window, CanvasTriangle tri, Colour col, float zBuffer[WIDTH][HEIGHT])
{
    float leftSlope = (tri.v1().x - tri.v0().x) / (tri.v1().y - tri.v0().y);
    float rightSlope = (tri.v2().x - tri.v0().x) / (tri.v2().y - tri.v0().y);
    float xL = tri.v0().x;
    float xR = tri.v0().x;
    std::vector<float> inverseZEdgesL = interpolateSingleFloats(tri.v0().inverseDepth, tri.v1().inverseDepth, tri.v1().y - tri.v0().y + 1);
    std::vector<float> inverseZEdgesR = interpolateSingleFloats(tri.v0().inverseDepth, tri.v2().inverseDepth, tri.v1().y - tri.v0().y + 1);
    for (int i = tri.v0().y; i <= tri.v1().y; i++)
    {
        CanvasPoint from = CanvasPoint(round(xL), i);
        CanvasPoint to = CanvasPoint(round(xR), i);
        drawLine3D(window, from, to, col, zBuffer, inverseZEdgesL[i - tri.v0().y], inverseZEdgesR[i - tri.v0().y]);
        xL += leftSlope;
        xR += rightSlope;
    }
}

//helper function for fillTriangle3D. Fills top-flat triangle in 3D space
void fillTopFlatTriangle3D(DrawingWindow &window, CanvasTriangle tri, Colour col, float zBuffer[WIDTH][HEIGHT])
{
    float leftSlope = (tri.v2().x - tri.v0().x) / (tri.v2().y - tri.v0().y);
    float rightSlope = (tri.v2().x - tri.v1().x) / (tri.v2().y - tri.v1().y);
    float xL = tri.v2().x;
    float xR = tri.v2().x;
    std::vector<float> inverseZEdgesL = interpolateSingleFloats(tri.v0().inverseDepth, tri.v2().inverseDepth, tri.v2().y - tri.v0().y);
    std::vector<float> inverseZEdgesR = interpolateSingleFloats(tri.v1().inverseDepth, tri.v2().inverseDepth, tri.v2().y - tri.v0().y);
    for (int i = tri.v2().y; i > tri.v0().y; i--)
    {
        CanvasPoint from = CanvasPoint(round(xL), i);
        CanvasPoint to = CanvasPoint(round(xR), i);
        drawLine3D(window, from, to, col, zBuffer, inverseZEdgesL[i - (int)round(tri.v0().y) - 1], inverseZEdgesR[i - (int)round(tri.v0().y) - 1]);
        xL -= leftSlope;
        xR -= rightSlope;
    }
}

//fills triangle with solid colour, taking depth into account
void fillTriangle3D(DrawingWindow &window, CanvasTriangle tri, Colour col, float zBuffer[WIDTH][HEIGHT])
{
    CanvasPoint points[] = {tri.v0(), tri.v1(), tri.v2()};
    sortPointsByY(points, 3);
    if (points[1].y == points[2].y)
    {
        CanvasTriangle t = CanvasTriangle(points[0], points[1], points[2]);
        fillBottomFlatTriangle3D(window, t, col, zBuffer);
    }
    else if (points[0].y == points[1].y)
    {
        CanvasTriangle t = CanvasTriangle(points[0], points[1], points[2]);
        fillTopFlatTriangle3D(window, t, col, zBuffer);
    }
    else
    {
        int extraX = round(points[0].x + (points[2].x - points[0].x) * (points[1].y - points[0].y) / (points[2].y - points[0].y));
        CanvasPoint extraPoint = CanvasPoint(extraX, points[1].y);
        extraPoint.inverseDepth = points[0].inverseDepth + (points[2].inverseDepth - points[0].inverseDepth) * (points[1].y - points[0].y) / (points[2].y - points[0].y);
        CanvasTriangle upper = CanvasTriangle(points[0], extraPoint, points[1]);
        CanvasTriangle lower = CanvasTriangle(extraPoint, points[1], points[2]);
        fillBottomFlatTriangle3D(window, upper, col, zBuffer);
        fillTopFlatTriangle3D(window, lower, col, zBuffer);
        Colour white = Colour(255, 255, 255);
    }
}

//helper function for textureTriangle3D. Textures bottom-flat triangles, taking depth into account (scuffed)
void textureBottomFlatTriangle3D(DrawingWindow &window, CanvasTriangle tri, TextureMap m, float zBuffer[WIDTH][HEIGHT])
{

    float leftSlope = (tri.v1().x - tri.v0().x) / (tri.v1().y - tri.v0().y);
    float rightSlope = (tri.v2().x - tri.v0().x) / (tri.v2().y - tri.v0().y);
    float xL = tri.v0().x;
    float xR = tri.v0().x;
    std::vector<float> zEdgesL = interpolateSingleFloats(tri.v0().depth, tri.v1().depth, tri.v1().y - tri.v0().y + 1);
    std::vector<float> zEdgesR = interpolateSingleFloats(tri.v0().depth, tri.v2().depth, tri.v1().y - tri.v0().y + 1);
    for (int i = tri.v0().y; i <= tri.v1().y; i++)
    {
        CanvasPoint from = CanvasPoint(round(xL), i);
        float fromFraction = (i - tri.v0().y) / (tri.v1().y - tri.v0().y);
        float fTPX = tri.v0().texturePoint.x + fromFraction * (tri.v1().texturePoint.x - tri.v0().texturePoint.x);
        float fTPY = tri.v0().texturePoint.y + fromFraction * (tri.v1().texturePoint.y - tri.v0().texturePoint.y);
        TexturePoint fromTP = TexturePoint(fTPX, fTPY);
        from.texturePoint = fromTP;
        CanvasPoint to = CanvasPoint(round(xR), i);
        float toFraction = (i - tri.v0().y) / (tri.v2().y - tri.v0().y);
        float tTPX = tri.v0().texturePoint.x + toFraction * (tri.v2().texturePoint.x - tri.v0().texturePoint.x);
        float tTPY = tri.v0().texturePoint.y + toFraction * (tri.v2().texturePoint.y - tri.v0().texturePoint.y);
        TexturePoint toTP = TexturePoint(tTPX, tTPY);
        to.texturePoint = toTP;
        drawTexture3D(window, from, to, m, zBuffer, zEdgesL[i - (int)round(tri.v0().y)], zEdgesR[i - (int)round(tri.v0().y)]);
        xL += leftSlope;
        xR += rightSlope;
    }
}

//helper function for textureTriangle3D. Textures top-flat triangles, taking depth into account (scuffed)
void textureTopFlatTriangle3D(DrawingWindow &window, CanvasTriangle tri, TextureMap m, float zBuffer[WIDTH][HEIGHT])
{
    float leftSlope = (tri.v2().x - tri.v0().x) / (tri.v2().y - tri.v0().y);
    float rightSlope = (tri.v2().x - tri.v1().x) / (tri.v2().y - tri.v1().y);
    float xL = tri.v2().x;
    float xR = tri.v2().x;
    std::vector<float> zEdgesL = interpolateSingleFloats(tri.v0().depth, tri.v2().depth, tri.v2().y - tri.v0().y);
    std::vector<float> zEdgesR = interpolateSingleFloats(tri.v1().depth, tri.v2().depth, tri.v2().y - tri.v0().y);
    for (int i = tri.v2().y; i > tri.v0().y; i--)
    {
        CanvasPoint from = CanvasPoint(round(xL), i);
        float fromFraction = (tri.v2().y - i) / (tri.v2().y - tri.v0().y); //(tri.v2().y - tri.v0().y)/i;
        float fTPX = tri.v2().texturePoint.x - fromFraction * (tri.v2().texturePoint.x - tri.v0().texturePoint.x);
        float fTPY = tri.v2().texturePoint.y - fromFraction * (tri.v2().texturePoint.y - tri.v0().texturePoint.y);
        TexturePoint fromTP = TexturePoint(fTPX, fTPY);
        from.texturePoint = fromTP;
        CanvasPoint to = CanvasPoint(round(xR), i);
        float toFraction = (tri.v2().y - i) / (tri.v2().y - tri.v1().y); //(tri.v2().y - tri.v1().y)/i;
        assert(toFraction <= 1 && toFraction >= 0);
        float tTPX = tri.v2().texturePoint.x - toFraction * (tri.v2().texturePoint.x - tri.v1().texturePoint.x);
        float tTPY = tri.v2().texturePoint.y - toFraction * (tri.v2().texturePoint.y - tri.v1().texturePoint.y);
        TexturePoint toTP = TexturePoint(tTPX, tTPY);
        to.texturePoint = toTP;
        drawTexture3D(window, from, to, m, zBuffer, zEdgesL[i - (int)round(tri.v0().y) - 1], zEdgesR[i - (int)round(tri.v0().y) - 1]);
        xL -= leftSlope;
        xR -= rightSlope;
    }
}

//draws a (scuffed) 3d textured triangle
void textureTriangle3D(DrawingWindow &window, CanvasTriangle tri, TextureMap m, float zBuffer[WIDTH][HEIGHT])
{
    CanvasPoint points[] = {tri.v0(), tri.v1(), tri.v2()};
    sortPointsByY(points, 3);
    if (points[1].y == points[2].y)
    {
        CanvasTriangle t = CanvasTriangle(points[0], points[1], points[2]);
        textureBottomFlatTriangle3D(window, t, m, zBuffer);
    }
    else if (points[0].y == points[1].y)
    {
        CanvasTriangle t = CanvasTriangle(points[0], points[1], points[2]);
        textureTopFlatTriangle3D(window, t, m, zBuffer);
    }
    else
    {
        int extraX = round(points[0].x + (points[2].x - points[0].x) * (points[1].y - points[0].y) / (points[2].y - points[0].y));
        CanvasPoint extraPoint = CanvasPoint(extraX, points[1].y);
        extraPoint.depth = points[0].depth + (points[2].depth - points[0].depth) * (points[1].y - points[0].y) / (points[2].y - points[0].y);
        float EPFrac = (points[1].y-points[0].y) / (points[2].y - points[0].y);
        float ePTPX = points[0].texturePoint.x + EPFrac * (points[2].texturePoint.x - points[0].texturePoint.x);
        float ePTPY = points[0].texturePoint.y + EPFrac * (points[2].texturePoint.y - points[0].texturePoint.y);
        TexturePoint EPTP = TexturePoint(ePTPX, ePTPY);
        extraPoint.texturePoint = EPTP;
        CanvasTriangle upper = CanvasTriangle(points[0], extraPoint, points[1]);
        CanvasTriangle lower = CanvasTriangle(extraPoint, points[1], points[2]);
        textureBottomFlatTriangle3D(window, upper, m, zBuffer);
        textureTopFlatTriangle3D(window, lower, m, zBuffer);
        Colour white = Colour(255, 255, 255);
    }
}

void texturePerspectiveBottomFlatTriangle3D(DrawingWindow &window, CanvasTriangle tri, TextureMap m, float zBuffer[WIDTH][HEIGHT])
{
    float leftSlope = (tri.v1().x - tri.v0().x) / (tri.v1().y - tri.v0().y);
    float rightSlope = (tri.v2().x - tri.v0().x) / (tri.v2().y - tri.v0().y);
    float xL = tri.v0().x;
    float xR = tri.v0().x;
    float sharedPointCOverZX = tri.v0().texturePoint.x/tri.v0().depth;
    float sharedPointCOverZY = tri.v0().texturePoint.y/tri.v0().depth;
    float cLOverzX = tri.v1().texturePoint.x/tri.v1().depth;
    float cLOverzY = tri.v1().texturePoint.y/tri.v1().depth;
    float cROverzX = tri.v2().texturePoint.x/tri.v2().depth;
    float cROverzY = tri.v2().texturePoint.y/tri.v2().depth;
    float inverseSharedDepth = 1/tri.v0().depth;
    float inverseLDepth = 1/tri.v1().depth;
    float inverseRDepth = 1/tri.v2().depth;
    int rowCount = 0;
    for (int i = tri.v0().y; i <= tri.v1().y; i++)
    {
        CanvasPoint from = CanvasPoint(round(xL), i);
        float q = (tri.v2().y - i) / (tri.v2().y - tri.v0().y);
        float currentInverseZL = (inverseSharedDepth * q + (1-q) * inverseLDepth);
        float currentInverseZR = (inverseSharedDepth * q + (1-q) * inverseRDepth);
        TexturePoint fromTP = TexturePoint((q * sharedPointCOverZX + (1-q) *cLOverzX) / currentInverseZL, (q * sharedPointCOverZY + (1-q) *cLOverzY) / currentInverseZL);
        from.texturePoint = fromTP;
        CanvasPoint to = CanvasPoint(round(xR), i);
        TexturePoint toTP = TexturePoint((q * sharedPointCOverZX + (1-q) *cROverzX) / currentInverseZR, (q * sharedPointCOverZY + (1-q) *cROverzY) / currentInverseZR);
        to.texturePoint = toTP;
        drawPerspectiveTexture3D(window, from, to, m, zBuffer, 1 / currentInverseZL, 1 / currentInverseZR);
        xL += leftSlope;
        xR += rightSlope;
        rowCount++;
    }
}

//helper function for textureTriangle3D. Textures top-flat triangles, taking depth into account (scuffed)
void texturePerspectiveTopFlatTriangle3D(DrawingWindow &window, CanvasTriangle tri, TextureMap m, float zBuffer[WIDTH][HEIGHT])
{
    float leftSlope = (tri.v2().x - tri.v0().x) / (tri.v2().y - tri.v0().y);
    float rightSlope = (tri.v2().x - tri.v1().x) / (tri.v2().y - tri.v1().y);
    float xL = tri.v2().x;
    float xR = tri.v2().x;

    float sharedPointCOverZX = tri.v2().texturePoint.x/tri.v2().depth;
    float sharedPointCOverZY = tri.v2().texturePoint.y/tri.v2().depth;
    float cLOverzX = tri.v0().texturePoint.x/tri.v0().depth;
    float cLOverzY = tri.v0().texturePoint.y/tri.v0().depth;
    float cROverzX = tri.v1().texturePoint.x/tri.v1().depth;
    float cROverzY = tri.v1().texturePoint.y/tri.v1().depth;
    float inverseSharedDepth = 1/tri.v2().depth;
    float inverseLDepth = 1/tri.v0().depth;
    float inverseRDepth = 1/tri.v1().depth;
    int rowCount = 0;
    for (int i = tri.v2().y; i > tri.v0().y; i--)
    {
        float q = (tri.v2().y - i) / (tri.v2().y - tri.v0().y);
        float currentInverseZL = (inverseSharedDepth * (1-q) + q * inverseLDepth);
        float currentInverseZR = (inverseSharedDepth * (1-q) + q * inverseRDepth);
        CanvasPoint from = CanvasPoint(round(xL), i);
        TexturePoint fromTP = TexturePoint(((1-q) * sharedPointCOverZX + q *cLOverzX) / currentInverseZL, ((1-q) * sharedPointCOverZY + q *cLOverzY) / currentInverseZL);
        from.texturePoint = fromTP;
        CanvasPoint to = CanvasPoint(round(xR), i);
        TexturePoint toTP = TexturePoint(((1-q) * sharedPointCOverZX + q *cROverzX) / currentInverseZR, ((1-q) * sharedPointCOverZY + q *cROverzY) / currentInverseZR);
        to.texturePoint = toTP;
        drawPerspectiveTexture3D(window, from, to, m, zBuffer, 1/currentInverseZL, 1/currentInverseZR);
        xL -= leftSlope;
        xR -= rightSlope;
        rowCount++;
    }
}

//draws a 3d textured triangle
void texturePerspetiveTriangle3D(DrawingWindow &window, CanvasTriangle tri, TextureMap m, float zBuffer[WIDTH][HEIGHT])
{
    CanvasPoint points[] = {tri.v0(), tri.v1(), tri.v2()};
    sortPointsByY(points, 3);
    if (points[1].y == points[2].y)
    {
        CanvasTriangle t = CanvasTriangle(points[0], points[1], points[2]);
        texturePerspectiveBottomFlatTriangle3D(window, t, m, zBuffer);
    }
    else if (points[0].y == points[1].y)
    {
        CanvasTriangle t = CanvasTriangle(points[0], points[1], points[2]);
        texturePerspectiveTopFlatTriangle3D(window, t, m, zBuffer);
    }
    else
    {
        int extraX = round(points[0].x + (points[2].x - points[0].x) * (points[1].y - points[0].y) / (points[2].y - points[0].y));
        CanvasPoint extraPoint = CanvasPoint(extraX, points[1].y);
        float fracAlong = (points[1].y -points[0].y) / (points[2].y - points[0].y);
        extraPoint.depth = 1/(1/points[0].depth + (1/points[2].depth - 1/points[0].depth) * fracAlong);
        float extraTextureX = extraPoint.depth * (points[0].texturePoint.x/points[0].depth + fracAlong * (points[2].texturePoint.x/points[2].depth - points[0].texturePoint.x/points[0].depth));
        float extraTextureY = extraPoint.depth * (points[0].texturePoint.y/points[0].depth + fracAlong * (points[2].texturePoint.y/points[2].depth - points[0].texturePoint.y/points[0].depth));
        TexturePoint newTexturePoint = TexturePoint(extraTextureX, extraTextureY);
        extraPoint.texturePoint = newTexturePoint;
        CanvasTriangle upper = CanvasTriangle(points[0], extraPoint, points[1]);
        CanvasTriangle lower = CanvasTriangle(extraPoint, points[1], points[2]);
        texturePerspectiveBottomFlatTriangle3D(window, upper, m, zBuffer);
        texturePerspectiveTopFlatTriangle3D(window, lower, m, zBuffer);
        Colour white = Colour(255, 255, 255);
    }
}