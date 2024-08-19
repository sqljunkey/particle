#include <stdio.h>
#include <stdlib.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

int main(void) {
    // Initialize the library and check potential ABI mismatches
    LIBXML_TEST_VERSION

    // Read the entire stdin into a buffer
    char *buffer = NULL;
    size_t size = 0;
    size_t capacity = 1024;
    buffer = (char *)malloc(capacity);
    if (!buffer) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    size_t read_size;
    while ((read_size = fread(buffer + size, 1, capacity - size, stdin)) > 0) {
        size += read_size;
        if (size >= capacity) {
            capacity *= 2;
            buffer = (char *)realloc(buffer, capacity);
            if (!buffer) {
                fprintf(stderr, "Memory reallocation failed\n");
                return 1;
            }
        }
    }
    buffer[size] = '\0'; // Null-terminate the buffer

    // Parse the XML from the buffer
    xmlDoc *doc = xmlReadMemory(buffer, size, "noname.xml", NULL, 0);
    if (doc == NULL) {
        fprintf(stderr, "Failed to parse XML\n");
        free(buffer);
        return 1;
    }
    free(buffer);

    // Print the XML document to stdout
    xmlSaveFormatFileEnc("-", doc, "UTF-8", 1);

    // Free the document
    xmlFreeDoc(doc);

    // Cleanup function for the XML library
    xmlCleanupParser();

    return 0;
}
