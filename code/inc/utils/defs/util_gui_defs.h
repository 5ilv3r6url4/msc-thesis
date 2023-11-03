#ifndef UTIL_GUI_DEFS_H
#define UTIL_GUI_DEFS_H

#define GRID_WIDTH 1000
#define GRID_HEIGHT 800
#define PIXEL_WIDTH  1.0 / GRID_WIDTH
#define PIXEL_HEIGHT 1.0 / GRID_HEIGHT

enum editor_tool_t  { TOOL_DRAW_FREE,
                      TOOL_DRAW_POLYGON,
                      TOOL_SELECT,
                      TOOL_INTERSECTION_NORMAL_FLIP,
                      TOOL_INTERSECTION_NORMAL_FLOP,
                      TOOL_REMOVE_INTERSECTION,
                      TOOL_INSPECTOR };

enum mouse_action_t { MOUSE_PRESS, MOUSE_MOVE, MOUSE_RELEASE };
enum mouse_button_t { NO_BUTTON, LEFT_BUTTON, RIGHT_BUTTON };

enum display_mode_t { DISPLAY_2D, DISPLAY_3D };

#endif // UTIL_GUI_DEFS_H
