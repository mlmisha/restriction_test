import argparse

import Bio
import Bio.Restriction
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dataclasses import dataclass

#переводит последовательность из FASTA-файла в объект Seq() Bio.Seq
def read_fasta(path):
    record = SeqIO.read(path, "fasta")
    return record.seq

#описание аргументов
def get_arguments():
    parser = argparse.ArgumentParser(
        description = "Симуляция рестрикционного клонирования"
    )

    parser.add_argument(
        "-fv", "--fasta-vector",
        metavar = ("vector.fasta"),
        help = "Путь к FASTA-файлу с последовательностью векторной плазмиды"
    )

    parser.add_argument(
        "-fd", "--fasta-donor",
        metavar = ("donor.fasta"),
        help = "Путь к FASTA-файлу с последовательностью донора"
    )

    parser.add_argument(
        "-vr", "--vector-restriction",
        nargs = 2,
        metavar = ("R1", "R2"),
        help = "Две рестриктазы для вектора через пробел (например EcoRI XhoI)"
    )

    parser.add_argument(
        "-dr", "--donor-restriction",
        nargs = 2,
        metavar = ("R1", "R2"),
        help = "Две рестриктазы для донора через пробел (например EcoRI XhoI)"
    )

    parser.add_argument(
        "-o", "--output",
        metavar = ("output.fasta"),
        help = "Путь к выходному FASTA-файлу. Иначе вывод в консоль"
    )

    return parser.parse_args()

#парсинг аргументов, для неиспользованных аргументов вызывается ввод с клавиатуры (кроме --output - он опционален)
def get_data(args):

    if args.fasta_vector:
        vector = read_fasta(args.fasta_vector)
    else:
        vector = Seq(input("Последовательность вектора: "))


    if args.fasta_donor:
        donor = read_fasta(args.fasta_donor)
    else:
        donor = Seq(input("Последовательность донора: "))


    if args.vector_restriction:
        v_restrictase = args.vector_restriction
    else:
        v_restrictase = input(
            "Название двух рестриктаз для вектора через пробел (например EcoRI XhoI): "
        ).split()


    if args.donor_restriction:
        d_restrictase = args.donor_restriction
    else:
        d_restrictase = input(
            "Название двух рестриктаз для донора через пробел (например EcoRI XhoI): "
        ).split()


    return vector, donor, v_restrictase, d_restrictase, args.output


#оборачиваю рестриктазы в отдельные классы, чтобы хранить вместе с ними индексы начала их сайтов и понятно обращаться
#к белку или позиции
@dataclass
class Restrictase:
    enzyme: object
    site_ind: int

#проверка наличия сайтов рестрикции для каждой рестриктазы по-отдельности
#на вход принимает последовательность (донор/вектор) и рестриктазу, vec - вектор ли (нужно для вывода сообщения)
def check_restrictase (seq, restrictase, vec = True):

    if vec: seq_type = "векторе" #это просто для вывода сообщения
    else: seq_type = "доноре"
    failure = True #есть ли проблема с сайтами рестриктазы

    res = restrictase.search(seq, linear = False) #поиск сайта рестриктазы в последовательности

    if len(res) == 0: #нет сайта
        print(f"Сайт рестриктазы {restrictase.__name__} в {seq_type} отсутсвует")
    elif len(res) > 1:#несколько сайтов
        print(f"В {seq_type} содержится более одного сайта рестриктазы {restrictase.__name__}: {res}")
    else:#1 сайт
        print(f"Сайт {restrictase.__name__} в {seq_type}: {res}")
        failure = False

    return failure, res

#для каждой рестриктазы проверяет, есть ли у нее сайт в соотв. последовательности и один ли,
#а потом возвращает объекты класса Restrictase вместе с позициями их сайтов
#vec/don_seq - последовательность, v_1/2 и d_1/2 - рестриктазы для вектора и донора, соотв. 
def check_all_restrictases (vec_seq, don_seq, v_1, v_2, d_1, d_2):
    v_1_failure, v_1_ind = check_restrictase(vec_seq, v_1)
    v_2_failure,v_2_ind = check_restrictase(vec_seq, v_2)

    d_1_failure, d_1_ind = check_restrictase(don_seq, d_1,vec = False)
    d_2_failure, d_2_ind = check_restrictase(don_seq, d_2,vec = False)

    if v_1_failure or v_2_failure or d_1_failure or d_2_failure:
        raise ValueError("Количество сайтов одной или нескольких рестриктаз не равно 1")

    vec_1 = Restrictase(v_1, v_1_ind[0])
    vec_2 = Restrictase(v_2, v_2_ind[0])
    don_1 = Restrictase(d_1, d_1_ind[0])
    don_2 = Restrictase(d_2, d_2_ind[0])

    return vec_1, vec_2, don_1, don_2

#проверка совместимости. Чтобы вставилось в той же ориентации левая рестриктаза должна давать совместимые сайты с левой,
#а правая с правой. Для начала определяю левые и правые рестриктазы
def define_left_and_right_restrictase (rest1, rest2):
    if rest1.site_ind > rest2.site_ind:
        return (rest2, rest1)
    else:
        return (rest1, rest2)

#для проверки совместимости использую оператор % из BioPython, проверяет можно ли лигировать концы получаемые парой рестриктаз
def check_compatibility (v_1, v_2, d_1, d_2):
    v_l, v_r = define_left_and_right_restrictase(v_1, v_2)
    d_l, d_r = define_left_and_right_restrictase(d_1, d_2)

    #контроль расстояния между сайтами
    if (v_r.site_ind - v_l.site_ind-len(v_l.enzyme.site)) <5 or (d_r.site_ind - d_l.site_ind-len(d_l.enzyme.site)) <5:
        raise ValueError("Сайты некоторых рестриктаз расположены слишком близко друг к другу")

    if v_r.enzyme % d_r.enzyme and v_l.enzyme % d_l.enzyme:
        print()
        if v_r.enzyme % v_l.enzyme: #дополнительно проверяю, не образуют ли вообще все рестриктазы совместимые концы
        print("ВНИМАНИЕ: Выбранные рестриктазы являются изошизомерами или изокаудомерами друг друга. Высока вероятность залипания вектора/вставки"+
                        " самих на себя (самолигирования) или встраивания вставки в вектор в обратной ориентации. Рекомендуется выбрать другие рестриктазы")
        print("Рестриктазы подходят")
        return v_l, v_r, d_l, d_r
    else:
        raise ValueError("Рестриктазы не образуют лигируемые концы")

#непосредственно клонирование
#vec/don_seq - последовательности вектора/донора, v_l/r и d_l/r - левые и правые рестриктазы вектора и донора, соотв. 
def clone (vec_seq, don_seq, v_l, v_r, d_l, d_r):
    d_st1 = d_l.enzyme.catalyze(don_seq, linear=False) #сначала режем левой рестриктазой
    d_st2 = d_r.enzyme.catalyze(d_st1[0]) #потом правой
    to_clone = d_st2[0] #последовательность вставки

    v_st1 = v_l.enzyme.catalyze(vec_seq, linear = False)
    v_st2 = v_r.enzyme.catalyze(v_st1[0]) 

    result = to_clone+v_st2[1] # вставка + вектор с вырезанным кусочком

    return result # полученная последовательность будет начинаться со вставки, а заканчиваться будет 5' концом сайта рестрикции левой рестриктазы


#последовательный запуск всех описанных выше функций, анализирующих введенные данные 
def process (vec_seq, don_seq, v_1, v_2, d_1, d_2):

    vec_r1, vec_r2, don_r1, don_r2 = check_all_restrictases(vec_seq, don_seq, v_1, v_2, d_1, d_2)

    v_l, v_r, d_l, d_r = check_compatibility(vec_r1, vec_r2, don_r1, don_r2)

    final_seq = clone(vec_seq, don_seq, v_l, v_r, d_l, d_r)

    return final_seq

#запуск скрипта
def main():
    args = get_arguments()

    vector, donor, vr, dr, output = get_data(args) #данные

    restr1_vector = getattr(Bio.Restriction, vr[0]) #получение объектов Bio.Restriction
    restr2_vector = getattr(Bio.Restriction, vr[1])

    restr1_donor = getattr(Bio.Restriction, dr[0])
    restr2_donor = getattr(Bio.Restriction, dr[1])  

    res = process(vector, donor, restr1_vector, restr2_vector, restr1_donor, restr2_donor)

    if output: #запись в FASTA файл полученной плазмиды, если указан флажок -o или --output
        SeqIO.write(
            SeqRecord(res, id="cloned_plasmid"),
            output,
            "fasta"
        )
    else:
        print()
        print(res)

if __name__ == "__main__":
    main()