"""Built-in shader sources for the rendering pipeline."""

# --- Blinn-Phong shading with per-vertex colour ---

VERTEX_SHADER = """
#version 330 core

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

in vec3 in_position;
in vec3 in_normal;
in vec3 in_color;

out vec3 frag_pos;
out vec3 frag_normal;
out vec3 frag_color_v;

void main() {
    vec4 world_pos = model * vec4(in_position, 1.0);
    frag_pos = world_pos.xyz;
    frag_normal = mat3(transpose(inverse(model))) * in_normal;
    frag_color_v = in_color;
    gl_Position = projection * view * world_pos;
}
"""

FRAGMENT_SHADER = """
#version 330 core

uniform vec3 light_pos;
uniform vec3 light_color;
uniform vec3 view_pos;
uniform float alpha;

in vec3 frag_pos;
in vec3 frag_normal;
in vec3 frag_color_v;

out vec4 frag_color;

void main() {
    // Ambient
    float ambient_strength = 0.2;
    vec3 ambient = ambient_strength * light_color;

    // Diffuse
    vec3 norm = normalize(frag_normal);
    vec3 light_dir = normalize(light_pos - frag_pos);
    float diff = max(dot(norm, light_dir), 0.0);
    vec3 diffuse = diff * light_color;

    // Specular (Blinn-Phong)
    float specular_strength = 0.4;
    vec3 view_dir = normalize(view_pos - frag_pos);
    vec3 halfway = normalize(light_dir + view_dir);
    float spec = pow(max(dot(norm, halfway), 0.0), 64.0);
    vec3 specular = specular_strength * spec * light_color;

    vec3 result = (ambient + diffuse + specular) * frag_color_v;
    frag_color = vec4(result, alpha);
}
"""

# --- Simple line shader for paths and wireframe ---

LINE_VERTEX_SHADER = """
#version 330 core

uniform mat4 view;
uniform mat4 projection;

in vec3 in_position;
in vec3 in_color;

out vec3 frag_color_v;

void main() {
    frag_color_v = in_color;
    gl_Position = projection * view * vec4(in_position, 1.0);
}
"""

LINE_FRAGMENT_SHADER = """
#version 330 core

in vec3 frag_color_v;
out vec4 frag_color;

void main() {
    frag_color = vec4(frag_color_v, 1.0);
}
"""
