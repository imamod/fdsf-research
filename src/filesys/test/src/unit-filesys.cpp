#include "Common.h"
#include "FileSys.h"

TEST_CASE("FolderCreation") {
    for (const auto& it : { 3, 5, 7, 9 }) {
        std::string absFilename = filesys::createDirectory(12, it);
        filesys::writeFile(absFilename + "y.txt", { 1,2,3,4,5 });
        filesys::writeFile(absFilename + "I.txt", { 1,2,3,4,5 });
        filesys::writeFile(absFilename + "I_.txt", { 1,2,3,4,5 });
        filesys::writeFile(absFilename + "Y.txt", { 1,2,3,4,5 });
    }

    for (const auto& it : { 3, 5, 7, 9 }) {
        std::string absFilename = filesys::createDirectory(12, it, "test/");
        filesys::writeFile(absFilename + "y.txt", { 1,2,3,4,5 });
        filesys::writeFile(absFilename + "I.txt", { 1,2,3,4,5 });
        filesys::writeFile(absFilename + "I_.txt", { 1,2,3,4,5 });
        filesys::writeFile(absFilename + "Y.txt", { 1,2,3,4,5 });
    }

}
