//
// Created by gian on 12/19/25.
#pragma once

class Mesh;

class Face {
public:
    const int a,b,c;
    Face(int a, int b, int c);
friend class Mesh;
};
