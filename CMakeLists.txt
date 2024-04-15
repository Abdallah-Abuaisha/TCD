cmake_minimum_required(VERSION 3.23)
project(TCDGitHub)

set(CMAKE_CXX_STANDARD 17)

add_executable(TCDGitHub main.cpp main.h Station.cpp Station.h Connection.cpp Connection.h Trip.cpp Trip.h Label.cpp
        Label.h Leg.cpp Leg.h Transfer.cpp Transfer.h Stop.cpp Stop.h StopEvent.cpp StopEvent.h Service.cpp
        Service.h Temporary_Label.cpp Temporary_Label.h)
