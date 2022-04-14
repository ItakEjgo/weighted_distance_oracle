# weighted_distance_oracle
EAR-Oralce implementation

### prameters

generate arbitrary point-to-arbitrary point query or not:
bool generate_flag = getarg(0, "--generate")

string input = getarg("", "--input"),

output = getarg("", "--output");

unsigned grid_num = getarg(4, "--grid-num");

unsigned q_num = getarg(100, "--query-num");

bool weighted_flag = getarg(0, "--weighted");

string method_type = getarg("", "--method");
float err = getarg(0.2, "--eps");

unsigned sp_num = getarg(5, "--sp-num");

unsigned parallel_num = getarg(1, "--parallel-num");

unsigned parallel_id = getarg(0, "--parallel-id");