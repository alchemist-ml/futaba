typedef enum {
    ALKALI_METALS,
    ALKALINE_EARTH_METALS,
    TRANSITION_METALS,
    BORON_GROUP,
    CARBON_GROUP,
    NITROGEN_GROUP,
    OXYGEN_GROUP,
    HALOGENS,
    NOBLE_GASES
} ElementGroup;

typedef struct {
    const char* symbol;
    const char* name;
    ElementGroup group;
    int atomic_number;
    int group_number;
    int oxi_states[4];
    int electronegativity; // Pauling Electronegativity (0.7 to 4.0)
} Element;


typedef struct {
    const char* formula;
    const char* name;
    int net_charge;
    char charge_type;
    Element* element[4];
    int element_count[4]; // count of corresponding elements
} Ion;

// Compound Type Bitmasks
// Acid, Base
#define ACID (1u << 0)
#define BASE (1u << 1)
#define AMPHOTERIC (1u << 2)

// SALTS
#define SALT (1u << 3)
#define HALIDE (1u << 4)
#define OXYSALT (1u << 5)

// Oxides
#define OXIDE (1u << 6)
#define HYDROXIDE (1u << 7)
#define PEROXIDE (1u << 8)
#define SUPEROXIDE (1u << 9)

// Others
#define HYDRIDE (1u << 10)
#define HYDRATE (1u << 11)
#define AMMONIATE (1u << 12)
#define SULPHIDE (1u << 13)
#define NITRIDE (1u << 14)
#define CARBIDE (1u << 15)
#define PHOSPHIDE (1u << 16)
#define BORIDE (1u << 17)
#define SILICIDE (1u << 18)
#define COORDINATION (1u << 19)

#define HAS_FLAG(p, c) (((p) & (c)) != 0)
#define HAS_FLAG_ALL(p, c) (((p) & (c)) == (c))

typedef struct {
    const char* formula;
    const char * name;
    unsigned int classification;
    Ion ion[4];
    int ion_count[4]; // count of corresponding ions
    float pH;
} Compound;