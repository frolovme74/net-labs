"""
ДЗ по дисциплине «Сети и телекоммуникации», вариант 22.

Циклический код Ц[15,11], образующий полином g(x) = x^4 + x + 1.
Информационный вектор: 11111010011.
Определяется корректирующая способность кода Ck для всех кратностей.

Модель: кодер -> канал (e = 0...01 .. 1...11) -> декодер (синдром)
"""

from math import comb

N = 15                  # длина кодового слова
K = 11                  # число информационных разрядов
R = N - K               # число проверочных разрядов
G = 0b10011             # g(x) = x^4 + x + 1
INFO = "11111010011"    # информационный вектор


# ---------- Арифметика полиномов над GF(2) ----------
# Полином хранится как целое число: бит j - коэффициент при x^j.

def degree(p: int) -> int:
    return p.bit_length() - 1


def poly_mod(a: int, g: int) -> int:
    """Остаток от деления a(x) на g(x) по модулю 2."""
    dg = degree(g)
    while a and degree(a) >= dg:
        a ^= g << (degree(a) - dg)
    return a


def poly_str(p: int) -> str:
    """Запись полинома в виде x^6 + x^5 + 1."""
    if p == 0:
        return "0"
    terms = []
    for j in range(degree(p), -1, -1):
        if p >> j & 1:
            terms.append("1" if j == 0 else "x" if j == 1 else f"x^{j}")
    return " + ".join(terms)


def bits(p: int, width: int = N) -> str:
    return format(p, f"0{width}b")


# ---------- Кодер ----------

def encode(m: int) -> int:
    """v(x) = x^(n-k)*m(x) + [x^(n-k)*m(x) mod g(x)]."""
    shifted = m << R
    return shifted | poly_mod(shifted, G)


# ---------- Декодер ----------

def build_syndrome_table() -> dict[int, int]:
    """Синдром -> однократная ошибка e(x) = x^j, j = 0..n-1."""
    return {poly_mod(1 << j, G): 1 << j for j in range(N)}


SYNDROMES = build_syndrome_table()


def syndrome(r: int) -> int:
    return poly_mod(r, G)


def decode(r: int) -> tuple[int, int]:
    """
    Возвращает (исправленное слово, синдром).
    При ненулевом синдроме инвертируется разряд,
    указанный таблицей синдромов однократных ошибок.
    """
    s = syndrome(r)
    if s == 0:
        return r, s
    return r ^ SYNDROMES[s], s


# ---------- Модель канала и подсчёт Ck ----------

def channel(v: int, e: int) -> int:
    """Канал: к переданному слову по mod 2 прибавляется ошибка."""
    return v ^ e


def analyse(v: int):
    """
    Перебор всех векторов ошибки от 0...01 до 1...11.
    Для каждой кратности i считается:
      Nk - исправленные ошибки (восстановлено исходное слово),
      No - обнаруженные ошибки (синдром не равен нулю).
    Возвращает строки (i, C(n,i), Nk, Ck, No).
    """
    nk = [0] * (N + 1)
    no = [0] * (N + 1)
    for e in range(1, 1 << N):
        i = e.bit_count()
        corrected, s = decode(channel(v, e))
        if s != 0:
            no[i] += 1
        if corrected == v:
            nk[i] += 1
    rows = []
    for i in range(1, N + 1):
        c = comb(N, i)
        rows.append((i, c, nk[i], nk[i] / c, no[i]))
    return rows


def note(i: int, c: int, nk: int, no: int) -> str:
    """Текст столбца «Примечание» таблицы результатов."""
    if nk == c:
        return "все ошибки исправлены"
    if no == c:
        return "все обнаружены, исправление ложное"
    if no == 0:
        return "e(x) - кодовое слово, не обнаружена"
    return f"не обнаружено {c - no}, остальные исправлены ложно"


# ---------- Вывод ----------

def print_encoding(m: int, v: int) -> None:
    shifted = m << R
    p = poly_mod(shifted, G)
    print(f"g(x) = {poly_str(G)}  ({bits(G, R + 1)})")
    print(f"m    = {bits(m, K)},  m(x) = {poly_str(m)}")
    print(f"x^{R}*m(x) = {bits(shifted)}")
    print(f"p(x) = {bits(p, R)},  p(x) = {poly_str(p)}")
    print(f"v    = {bits(v)}  ({bits(m, K)}.{bits(p, R)})")
    print(f"Проверка: v(x) mod g(x) = {syndrome(v)}")
    print()


def print_syndrome_table() -> None:
    print("Синдромы однократных ошибок")
    print(f"{'e(x)':>6} | {'s(x)':<18} | s3s2s1s0")
    print("-" * 40)
    for j in range(N):
        s = poly_mod(1 << j, G)
        print(f"{poly_str(1 << j):>6} | {poly_str(s):<18} | {bits(s, R)}")
    print()


def print_example(v: int, e: int) -> None:
    r = channel(v, e)
    corrected, s = decode(r)
    found = poly_str(SYNDROMES[s]) if s else "нет"
    result = "верно" if corrected == v else "НЕВЕРНО"
    print(f"e = {bits(e)}  (e(x) = {poly_str(e)})")
    print(f"  r = {bits(r)}, s = {bits(s, R)}, ошибка: {found}")
    print(f"  исправлено: {bits(corrected)} - {result}")
    print()


def print_results(rows) -> None:
    print("Таблица 1 - Корректирующая способность кода Ц[15,11]")
    print(f"{'i':>3} | {'C(n,i)':>6} | {'Nk':>3} | {'Ck':>6} | Примечание")
    print("-" * 72)
    for i, c, nk, ck, no in rows:
        print(f"{i:>3} | {c:>6} | {nk:>3} | {ck:>6.4f} | {note(i, c, nk, no)}")
    print("-" * 72)
    total = sum(row[1] for row in rows)
    fixed = sum(row[2] for row in rows)
    print(f"Всего векторов ошибки: {total}, исправлено: {fixed}")


def main() -> None:
    m = int(INFO, 2)
    v = encode(m)
    print_encoding(m, v)
    print_syndrome_table()
    print_example(v, 1 << 6)
    print_example(v, (1 << 6) | (1 << 2))
    print_results(analyse(v))


if __name__ == "__main__":
    main()
