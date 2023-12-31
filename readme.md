### Описание библиотеки
Библиотека помехозащищенного кодирования кодом Рида-Соломона (Reed–Solomon codes).

### Включение библиотки
Библиотека реализована в виде одного заголовочного файла rscode_headeronly.h.
```cpp
#include "rscode_headeronly.h"
```
### Пример использования
Кодируем последовательность слов размером 8 бит { 240, 229, 234, 224, 213 } с возможностью исправления до 2 кодовых слов, используя статические таблицы кодирования.
```cpp
std::vector<uint8_t> data{ 240, 229, 234, 224, 213 };

auto coder{ ZLib::ErrCorr::RSCoder<uint8_t, 2, true>() };
auto encoded_data{ coder.encode(data) };

//искажение данных
encoded_data[0] = 203; encoded_data[1] = 236; //{ 203, 229, 234, 224, 213 }

auto decoded_data = coder.decode(encoded_data); //{ 240, 229, 234, 224, 213 }
```
