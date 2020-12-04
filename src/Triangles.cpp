//
// Created by ralph on 26/11/2020.
//

//draws a stroked 2d triangle
void drawStrokedTriangle(DrawingWindow &window, CanvasPoint p1, CanvasPoint p2, CanvasPoint p3, Colour col)
{
    CanvasTriangle tri = CanvasTriangle(p1, p2, p3);
    drawLine(window, tri.v0(), tri.v1(), col);
    drawLine(window, tri.v1(), tri.v2(), col);
    drawLine(window, tri.v2(), tri.v0(), col);
}

//draws a stroked 2d triangle
void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle t, Colour col)
{
    drawStrokedTriangle(window, t.v0(), t.v1(), t.v2(), col);
}

//fill bottom-flat 2d triangle (solid col)
void fillBottomFlatTriangle(DrawingWindow &window, CanvasTriangle tri, Colour col)
{
    float leftSlope = (tri.v1().x - tri.v0().x) / (tri.v1().y - tri.v0().y);
    float rightSlope = (tri.v2().x - tri.v0().x) / (tri.v2().y - tri.v0().y);
    float xL = tri.v0().x;
    float xR = tri.v0().x;
    for (int i = tri.v0().y; i <= tri.v1().y; i++)
    {
        CanvasPoint from = CanvasPoint(round(xL), i);
        CanvasPoint to = CanvasPoint(round(xR), i);
        drawLine(window, from, to, col);
        xL += leftSlope;
        xR += rightSlope;
    }
}

//fill top-flat 2d triangle (solid col)
void fillTopFlatTriangle(DrawingWindow &window, CanvasTriangle tri, Colour col)
{
    float leftSlope = (tri.v2().x - tri.v0().x) / (tri.v2().y - tri.v0().y);
    float rightSlope = (tri.v2().x - tri.v1().x) / (tri.v2().y - tri.v1().y);
    float xL = tri.v2().x;
    float xR = tri.v2().x;
    for (int i = tri.v2().y; i > tri.v0().y; i--)
    {
        CanvasPoint from = CanvasPoint(round(xL), i);
        CanvasPoint to = CanvasPoint(round(xR), i);
        drawLine(window, from, to, col);
        xL -= leftSlope;
        xR -= rightSlope;
    }
}

// fill 2d triangle with solid colour
void fillTriangle(DrawingWindow &window, CanvasTriangle tri, Colour col)
{
    CanvasPoint points[] = {tri.v0(), tri.v1(), tri.v2()};
    sortPointsByY(points, 3);
    if (points[1].y == points[2].y)
    {
        CanvasTriangle t = CanvasTriangle(points[0], points[1], points[2]);
        fillBottomFlatTriangle(window, t, col);
    }
    else if (points[0].y == points[1].y)
    {
        CanvasTriangle t = CanvasTriangle(points[0], points[1], points[2]);
        fillBottomFlatTriangle(window, t, col);
    }
    else
    {
        int extraX = round(points[0].x + (points[2].x - points[0].x) * (points[1].y - points[0].y) / (points[2].y - points[0].y));
        CanvasPoint extraPoint = CanvasPoint(extraX, points[1].y);
        CanvasTriangle upper = CanvasTriangle(points[0], extraPoint, points[1]);
        CanvasTriangle lower = CanvasTriangle(extraPoint, points[1], points[2]);
        fillBottomFlatTriangle(window, upper, col);
        fillTopFlatTriangle(window, lower, col);
    }
}

//helper function for textureTriangle. Textures bottom-flat triangle in 2D space.
void textureBottomFlatTriangle(DrawingWindow &window, CanvasTriangle tri, TextureMap m)
{
    float leftSlope = (tri.v1().x - tri.v0().x) / (tri.v1().y - tri.v0().y);
    float rightSlope = (tri.v2().x - tri.v0().x) / (tri.v2().y - tri.v0().y);
    float xL = tri.v0().x;
    float xR = tri.v0().x;
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
        drawTexture(window, from, to, m);
        xL += leftSlope;
        xR += rightSlope;
    }
}

//helper function for textureTriangle. Textures top-flat triangle in 2D space.
void textureTopFlatTriangle(DrawingWindow &window, CanvasTriangle tri, TextureMap m)
{
    float leftSlope = (tri.v2().x - tri.v0().x) / (tri.v2().y - tri.v0().y);
    float rightSlope = (tri.v2().x - tri.v1().x) / (tri.v2().y - tri.v1().y);
    float xL = tri.v2().x;
    float xR = tri.v2().x;
    for (int i = tri.v2().y; i > tri.v0().y; i--)
    {
        CanvasPoint from = CanvasPoint(round(xL), i);
        float fromFraction = (tri.v2().y - i) / (tri.v2().y - tri.v0().y); //(tri.v2().y - tri.v0().y)/i;
        assert(fromFraction <= 1 && fromFraction >= 0);
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
        drawTexture(window, from, to, m);
        xL -= leftSlope;
        xR -= rightSlope;
    }
}

//Draws 2D textured triangle
void textureTriangle(DrawingWindow &window, CanvasTriangle tri)
{
    const std::string fileName = "src/texture.ppm";
    TextureMap m = TextureMap(fileName);
    CanvasPoint points[] = {tri.v0(), tri.v1(), tri.v2()};
    sortPointsByY(points, 3);
    int extraX = round(points[0].x + (points[2].x - points[0].x) * (points[1].y - points[0].y) / (points[2].y - points[0].y));
    CanvasPoint extraPoint = CanvasPoint(extraX, points[1].y);
    float EPFrac = points[1].y / (points[2].y - points[0].y);
    float ePTPX = points[0].texturePoint.x + EPFrac * (points[2].texturePoint.x - points[0].texturePoint.x);
    float ePTPY = points[0].texturePoint.y + EPFrac * (points[2].texturePoint.y - points[0].texturePoint.y);
    TexturePoint EPTP = TexturePoint(ePTPX, ePTPY);
    extraPoint.texturePoint = EPTP;
    CanvasTriangle upper = CanvasTriangle(points[0], extraPoint, points[1]);
    CanvasTriangle lower = CanvasTriangle(extraPoint, points[1], points[2]);
    Colour white = Colour(255, 255, 255);
    textureTopFlatTriangle(window, lower, m);
    textureBottomFlatTriangle(window, upper, m);
}

