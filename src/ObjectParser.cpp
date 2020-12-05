//
// Created by ralph on 26/11/2020.
//bruh i wrote this not you :angry:
//
//Ns - specular exponent
//Ni - optical density
//Ka - reflectivity

//
//loads object file, returning triangles contained within the file coloured corresponding to a given palette (map)
std::vector<ModelTriangle> loadOBJ(const std::string fileName, float scalingFactor, std::map<std::string, Material> m) //, std::map<std::string, TextureMap> t)
{
    std::ifstream inputStream(fileName, std::ifstream::in);
    std::string nextLine;
    std::vector<ModelTriangle> ret;
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> vertexNorms;
    std::vector<glm::vec2> textureVertices;
    Colour currentCol = Colour(255, 255, 255);
    TextureMap currentText;
    //bool isTexture = false;
    bool usingNorms = false;
    while (!inputStream.eof())
    {
        std::getline(inputStream, nextLine);
        if (nextLine.length() > 0)
        {
            if (nextLine.substr(0, 6) == "usemtl")
            {
                std::string colStr = nextLine.substr(7);
                //std::map<std::string, Material>::iterator it = m.find(colStr);
                std::cout << "string: " << colStr << std::endl;
                currentCol = m[colStr].diffuseColour;
                currentCol.name = colStr;
                if(m[colStr].usingTexture){
                    currentText = m[colStr].textureMap;
                    std::cout << "TEXTURES" << currentText << std::endl;
                }

                // if (m[colStr].usingTexture)
                // if (it != m.end()) //if colour
                // {
                //     currentCol = m[colStr].diffuseColour;
                //     currentCol.name = colStr;
                //     isTexture = false;
                // }
                // else //if texture
                // {
                //     currentCol = Colour(-1, -1, -1);
                //     currentCol.name = colStr;
                //     currentText = m[colStr].textureMap;
                //     isTexture = true;
                // }
            }
            if (nextLine.at(0) == 'v')
            {

                std::string workingString = nextLine.substr(nextLine.find(' ') + 1);
                //ngl this do be kinda long
                float x = stof(workingString.substr(0, workingString.find(' ')));
                x *= scalingFactor;
                workingString = workingString.substr(workingString.find(' ') + 1);
                float y = stof(workingString.substr(0, workingString.find(' ')));
                y *= -1 * scalingFactor;
                if (nextLine.at(1) == 't')
                {
                    std::cout << "wwwwww: " << currentText.height << std::endl;
                    float xReal = x * currentText.width;
                    std::cout<<currentText.height<<"|||| ok"<<std::endl;
                    float yReal = -y * currentText.height;
                    glm::vec2 point = glm::vec2(xReal, yReal);
                    textureVertices.push_back(point);
                }
                else
                {
                    workingString = workingString.substr(workingString.find(' ') + 1);
                    float z = stof(workingString.substr(0, workingString.find(' ')));
                    z *= scalingFactor;
                    glm::vec3 point = glm::vec3(x, y, z);
                    if (nextLine.at(1) == 'n')
                    {
                        usingNorms = true;
                        //std::cout << "yo" << std::endl;
                        vertexNorms.push_back(point);
                    }
                    else
                    {
                        vertices.push_back(point);
                    }
                }
            }
            else if (nextLine.at(0) == 'f')
            {
                std::string workingString = nextLine.substr(2);
                std::vector<std::string> faceChunks = split(workingString, ' ');
                int index1 = 0, index2 = 0, index3 = 0;
                std::array<int, 3> vArray = {index1, index2, index3};
                TexturePoint tp1, tp2, tp3;
                std::array<TexturePoint, 3> tArray = {tp1, tp2, tp3};
                int nIndex1 = 0, nIndex2 = 0, nIndex3 = 0;
                std::array<int, 3> nArray = {nIndex1, nIndex2, nIndex3};
                for (int i = 0; i < 3; i++)
                {
                    //std::cout << faceChunks[i] << std::endl;
                    std::vector<std::string> splitIndices = split(faceChunks[i], '/');
                    vArray[i] = stoi(splitIndices[0]) - 1;
                    //std::cout << splitIndices.size() << std::endl;
                    if (splitIndices.size() > 1 && splitIndices[1] != "")
                    {
                        //std::cout<<"yeehaw"<<std::endl;
                        //why must you bring such pain
                        //std::cout<< textureVertices.size()<< std::endl;
                        print(splitIndices[1]);
                        //std::cout<<textureVertices[stoi(splitIndices[1])-1].x<<std::endl;
                        tArray[i] = TexturePoint(textureVertices[stoi(splitIndices[1])-1].x, textureVertices[stoi(splitIndices[1])-1].y);
                        //std::cout<<"womp"<<std::endl;
                    }
                    if (splitIndices.size() == 3 && splitIndices[2] != "")
                    {
                        nArray[i] = stoi(splitIndices[2]) - 1;
                    }
                }
                // if (tArray[0].x != 0){
                //     std::cout << "tarraywswwwwwww"<< std::endl;
                //     std::cout << tArray[0] << std::endl;
                // }
                std::array<TexturePoint, 3> tPoints = {tArray[0], tArray[1], tArray[2]};
                std::cout << "I1: " << vArray[0] << ", I2: " << vArray[1] << ", I3: " << vArray[2] << std::endl;
                ModelTriangle tri = ModelTriangle(vertices[vArray[0]], vertices[vArray[1]], vertices[vArray[2]], currentCol);
                tri.mat = m[tri.colour.name];
                tri.normal = glm::normalize(glm::cross(vertices[vArray[1]] - vertices[vArray[0]], vertices[vArray[2]] - vertices[vArray[0]]));
                if (usingNorms)
                {
                    std::array<glm::vec3, 3> vNorms = {vertexNorms[nArray[0]], vertexNorms[nArray[1]], vertexNorms[nArray[2]]};
                    tri.vertexNorms = vNorms;
                    if (vertexNorms[nArray[0]] == vertexNorms[nArray[1]] && vertexNorms[nArray[1]] == vertexNorms[nArray[2]] && vertexNorms[nArray[2]] == tri.normal)
                    {
                        tri.isPhongShaded = false;
                    }
                    else
                    {
                        tri.isPhongShaded = true;
                    }
                }
                tri.texturePoints = tPoints;
                if (tri.texturePoints[0].x != 0){
                //std::cout << "parse"<< std::endl;
                std::cout << tri.texturePoints[0] << std::endl;}
                ret.push_back(tri);
                std::cout << "ooooo" << std::endl; //printVector(tri.normal) << "-----------" << printVector(tri.vertexNorms[0]) << "----------" << printVector(tri.vertexNorms[1]) << "----------" << printVector(tri.vertexNorms[2]) << std::endl;
            }
        }
    }
    return ret;
}

//load material file into passed map (palette)
void loadMTL(const std::string fileName, std::map<std::string, Material> *m) //, std::map<std::string, TextureMap> *t)
{
    std::ifstream inputStream(fileName, std::ifstream::in);
    std::string nextLine;

    bool firstMtl = true;
    Material newMat;
    while (!inputStream.eof())
    {
        std::getline(inputStream, nextLine);
        //std::cout << "xd" << std::endl;
        if (nextLine.length() > 0)
        {
            if (nextLine.substr(0, 6) == "newmtl")
            {
                std::cout << "new material" << std::endl;
                if (!firstMtl)
                {
                    std::cout << "adding material" << std::endl;
                    m->insert(std::make_pair(newMat.name, newMat));
                }
                firstMtl = false;
                std::string name = nextLine.substr(7);
                newMat = Material(name);
                newMat.specularColour = Colour(255,255,255);
            }
            else if (nextLine.substr(0, 6) == "map_Kd")
            {
                std::string workingString = "src/" + nextLine.substr(7);
                TextureMap tm = TextureMap(workingString);
                //t->insert(std::make_pair(name, tm));
                newMat.textureMap = tm;
                newMat.usingTexture = true;
            }
            else if (nextLine.substr(0, 2) == "Kd")
            {
                std::string workingString = nextLine.substr(3);
                std::vector<std::string> splitVals = split(workingString, ' ');
                float r = round(stof(splitVals[0]) * 255);
                float g = round(stof(splitVals[1]) * 255);
                float b = round(stof(splitVals[2]) * 255);
                Colour col = Colour(r, g, b);
                newMat.diffuseColour = col;
            }
            else if (nextLine.substr(0, 2) == "Ns")
            {
                std::string workingString = nextLine.substr(3);
                float exponent = stof(workingString);
                newMat.specularExponent = exponent;
            }
            else if (nextLine.substr(0, 2) == "Ni")
            {
                std::string workingString = nextLine.substr(3);
                float optDens = stof(workingString);
                newMat.opticalDensity = optDens;
            }
            else if (nextLine.substr(0, 2) == "Tr")
            {
                std::string workingString = nextLine.substr(3);
                float transparency = stof(workingString);
                newMat.transparency = transparency;
            }
            else if (nextLine.substr(0, 2) == "Rf")
            {
                std::cout<<"peeeeeen"<<std::endl;
                std::string workingString = nextLine.substr(3);
                float reflectivity = stof(workingString);
                std::cout<<reflectivity<<std::endl;
                newMat.reflectivity = reflectivity;
                std::cout<<newMat.name<<std::endl;
            }
        }
    }
    m->insert(std::make_pair(newMat.name, newMat));
    std::cout << "peepee" << std::endl;
}
