#include <stdio.h>
#include "orthrus/orthrus.h"

int main(void) {
    printf("[*] Testing Orthrus C-ABI from pure C...\n");
    printf("[+] Linked with liborthrus version: %s\n", orthrus_version());

    orthrus_quat_t result;
    orthrus_status_t status = orthrus_klpt_navigate("431", "97", 32, &result);

    if (status == ORTHRUS_OK) {
        printf("[SUCCESS] C-API returned solution:\n");
        printf("    gamma = (%s + %s*i + %s*j + %s*k) / 2\n",
               result.a0, result.a1, result.a2, result.a3);
        return 0;
    } else {
        printf("[ERROR] C-API failed with status code: %d\n", status);
        return 1;
    }
}
