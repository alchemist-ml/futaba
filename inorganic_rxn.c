#define MAX 4

typedef enum {
    METAL,
    NONMETAL,
    METTALOID
} ElementType;

typedef struct {
    const char * symbol;
    ElementType type;
    int oxi_states[MAX];
    int oxi_count;
    int group;
    int period;
    int valance_electrons;
    int electronegativity;
} Element;

typedef struct {
    const Element elements[MAX];
    const int counts[MAX];
    int element_count;
    int oxi_states[MAX];
    int charge;
} Species;

