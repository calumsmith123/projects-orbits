/* checks character string is a positive integer */

static int number(const char *num) {
    int i = 0;
    for (; num[i] != 0; ++i) {
        if (!isdigit(num[i])) {
            i = 0;
            break;
        }
    }
    return i;
}
