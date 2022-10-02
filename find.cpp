// Include files
#include "find.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
namespace coder {
    void eml_find(const ::coder::array<boolean_T, 1U> &x, ::coder::array<int, 1U>
    &i) {
        int counter = 0;
        int index = 0;
        int nx = x.size(0);
        i.set_size(nx);
        while (index <= nx - 1) {
            if (x[index]) {
                i[counter] = index + 1;
                counter++;
                index++;
            } else {
                index++;
            }
        }

        i.set_size(counter);
    }
}
