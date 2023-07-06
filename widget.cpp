void draw_ball(mat4 P, mat4 V, vec3 s, vec3 color = monokai.white, real scale = 1.0) { library.meshes.sphere.draw(P, V, M4_Translation(s) * M4_Scaling(scale * 0.01), color); };
void draw_pipe(mat4 P, mat4 V, vec3 s, vec3 t, vec3 color = monokai.white, real scale = 1.0) {
    real r = scale * 0.007;
    vec3 e = t - s;
    real mag_e = norm(e);
    vec3 E = { 0.0, 1.0, 0.0 };
    vec3 e_hat = e / mag_e;
    vec3 u = cross(e_hat, E);
    real norm_u = norm(u);
    mat4 R = (norm_u < 1.0e-7) ? globals.Identity : M4_RotationAxisAngle(u / norm_u, ((dot(e_hat, E) < 0) ? 1.0 : -1.0) * asin(norm_u));
    library.meshes.cylinder.draw(P, V, M4_Translation(0.5 * (s + t)) * R * M4_Scaling(r, mag_e, r) * M4_Translation(0.0, -0.5, 0.0), color);
};

// TODO: MOUSE_OWNER stuff
// TODO: returns whether the central targetPosition was just pressed
// FOREVER: because we're calling this widget one at a time (instead of on an array) multiple widget's handles can show up as hot


struct WidgetResult {
    bool mouseClickConsumed;
    bool pleaseDisableHandle;
    bool mouseHotConsumed;
    bool recastFeaturePoint;
};

WidgetResult widget(mat4 P, mat4 V, int i, vec3 *targetPosition, vec3 featurePointPosition, vec3 widgetColor = { 1.0, 0.0, 1.0 }) {
    mat4 PV = P * V;


    // statics
    enum State { Normal, DraggingHandle, DraggingFeaturePoint };
    static State state;
    static int selected_i = -1; // dragging_i
    static int selected_handle = -1;

    // params
    double tol = .02;

    // building the widget
    double L_handle; {
        L_handle = .07;
        // double _L_handle_NDC = .07;
        //norm(*targetPosition - transformPoint(inverse(PV), transformPoint(PV, *targetPosition) + V3(_L_handle_NDC, 0, 0)));
    }
    vec3 handle_positions[3] = { *targetPosition + V3(L_handle, 0, 0), *targetPosition + V3(0, L_handle, 0), *targetPosition + V3(0, 0, L_handle) };

    // hots
    bool targetPositionIsHot = (norm(transformPoint(PV, *targetPosition).xy - globals.mouse_position_NDC) < tol);
    bool featurePointPositionIsHot = (norm(transformPoint(PV, featurePointPosition).xy - globals.mouse_position_NDC) < tol);
    int hot_handle; {
        hot_handle = -1;
        if (!targetPositionIsHot && (state == Normal)) {
            for_(d, 3) {
                vec3 tmp = transformPoint(PV, handle_positions[d]);
                tmp.z = 0;
                if (norm(V2(tmp.x, tmp.y) - globals.mouse_position_NDC) < tol) {
                    hot_handle = d;
                }
            }
        }
    }

    { // draw (incorporating infor about hots)

        draw_ball(P, V, *targetPosition, widgetColor);
        draw_ball(P, V, featurePointPosition, widgetColor, ((state == DraggingFeaturePoint) && (selected_i == i)) ? 0.8 : (featurePointPositionIsHot ? 1.2 : 1.0));
        if (!targetPositionIsHot) { // translate widget handles
            draw_pipe(P, V, *targetPosition, featurePointPosition, widgetColor, 0.7);
            for_(d, 3) {
                bool selected = (selected_i == i) && (selected_handle == d);
                bool hot = (hot_handle == d);
                real scale = selected ? 0.6 : hot ? 1.0 : 0.8;
                vec3 color = selected ? AVG(monokai.white, AVG(monokai.white, widgetColor)) : hot ? AVG(monokai.white, widgetColor) : widgetColor;
                draw_pipe(P, V, handle_positions[d], *targetPosition, color, scale);
                draw_ball(P, V, handle_positions[d], color, scale);
            }
        } else {
            draw_pipe(P, V, *targetPosition, featurePointPosition, AVG(monokai.black, widgetColor), 0.4);
        }
    }


    // bail conditions
    if (globals._mouse_owner != COW_MOUSE_OWNER_NONE && globals._mouse_owner != COW_MOUSE_OWNER_WIDGET) return {};
    if ((selected_i != -1) && (selected_i != i)) return {};

    if (state == Normal) {
        if (globals.mouse_left_pressed) { // *
            if (targetPositionIsHot) return { true, true };
            if (featurePointPositionIsHot) {
                globals._mouse_owner = COW_MOUSE_OWNER_WIDGET;
                state = DraggingFeaturePoint;
                selected_i = i;
                return { true };
            }
            if ((hot_handle != -1) && (selected_i == -1)) {
                globals._mouse_owner = COW_MOUSE_OWNER_WIDGET;
                state = DraggingHandle;
                selected_i = i;
                selected_handle = hot_handle;
                return { true, false };
            }
        }
    } else if (state == DraggingHandle) {
        if (globals.mouse_left_held && (selected_handle != -1)) {
            mat4 World_from_NDC = inverse(PV);
            vec3 s = transformPoint(World_from_NDC, V3(globals.mouse_position_NDC, -1.0));
            vec3 t = transformPoint(World_from_NDC, V3(globals.mouse_position_NDC,  1.0));
            vec3 new_handle_position;
            _line_line_closest_points(*targetPosition, handle_positions[selected_handle], s, t, &new_handle_position, 0);
            *targetPosition += new_handle_position - handle_positions[selected_handle];

        } else if (globals.mouse_left_released) {
            globals._mouse_owner = COW_MOUSE_OWNER_NONE;
            state = Normal;
            selected_i = -1;
            selected_handle = -1;
        }
    } else if (state == DraggingFeaturePoint) {
        if (globals.mouse_left_released) {
            globals._mouse_owner = COW_MOUSE_OWNER_NONE;
            selected_i = -1;
            state = Normal;
        }
    }

    return { false, false, (targetPositionIsHot || featurePointPositionIsHot || (hot_handle != -1)), (state == DraggingFeaturePoint)};
}
