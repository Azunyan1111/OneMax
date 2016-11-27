# -*- coding: utf-8 -*-
# 著作 Azunyan
# コーディング規約は基本的にPEP8に準じて作成します。
# 本コードはOneMax問題を遺伝的アルゴリズムを用いて解くプログラムコードである。

import GeneticAlgorithm as ga
import random
from decimal import Decimal

# 遺伝子情報の長さ
GENOM_LENGTH = 100
# 遺伝子集団の大きさ
MAX_GENOM_LIST = 100
# 遺伝子選択数
SELECT_GENOM = 20
# 個体突然変異確率
INDIVIDUAL_MUTATION = 0.1
# 遺伝子突然変異確率
GENOM_MUTATION = 0.1
# 繰り返す世代数
MAX_GENERATION = 40


def create_genom(length):
    """
    引数で指定された桁のランダムな遺伝子情報を生成、格納したgenomClassで返します。
    :param length: 遺伝子情報の長さ
    :return: 生成した個体集団genomClass
    """
    genome_list = []
    for i in range(length):
        genome_list.append(random.randint(0, 1))
    return ga.genom(genome_list, 0)


def evaluation(ga):
    """評価関数です。今回は全ての遺伝子が1となれば最適解となるので、
    合計して遺伝子と同じ長さの数となった場合を1として0.00~1.00で評価します
    :param ga: 評価を行うgenomClass
    :return: 評価処理をしたgenomClassを返す
    """
    genom_total = sum(ga.getGenom())
    return Decimal(genom_total) / Decimal(GENOM_LENGTH)


def select(ga, elite):
    """選択関数です。エリート選択を行います
    評価が高い順番にソートを行った後、一定以上
    :param ga: 選択を行うgenomClassの配列
    :return: 選択処理をした一定のエリート、genomClassを返す
    """
    # 現行世代個体集団の評価を高い順番にソートする
    sort_result = sorted(ga, reverse=True, key=lambda u: u.evaluation)
    # 一定の上位を抽出する
    result = [sort_result.pop(0) for i in range(elite)]
    return result


def crossover(ga_one, ga_second):
    """交叉関数です。二点交叉を行います。
    :param ga: 交叉させるgenomClassの配列
    :param ga_one:
    :param ga_second:
    :return: 二つの子孫genomClassを格納したリスト返す
    """
    # 子孫を格納するリストを生成します
    genom_list = []
    # 入れ替える二点の点を設定します→[1:25]
    cross_one = random.randint(0, GENOM_LENGTH)
    cross_second = random.randint(cross_one, GENOM_LENGTH)
    # 遺伝子を取り出します
    one = ga_one.getGenom()
    second = ga_second.getGenom()
    # 交叉させます
    progeny_one = one[:cross_one] + second[cross_one:cross_second] + one[cross_second:]
    progeny_second = second[:cross_one] + one[cross_one:cross_second] + second[cross_second:]
    # genomClassインスタンスを生成して子孫をリストに格納する
    genom_list.append(ga.genom(progeny_one, 0))
    genom_list.append(ga.genom(progeny_second, 0))
    return genom_list


def next_generation_gene_create(ga, ga_elite, ga_progeny):
    """
    世代交代処理を行います
    :param ga: 現行世代個体集団
    :param ga_elite: 現行世代エリート集団
    :param ga_progeny: 現行世代子孫集団
    :return: 次世代個体集団
    """
    # 現行世代個体集団の評価を低い順番にソートする
    next_generation_geno = sorted(ga, reverse=False, key=lambda u: u.evaluation)
    # 追加するエリート集団と子孫集団の合計ぶんを取り除く
    for i in range(0, len(ga_elite) + len(ga_progeny)):
        next_generation_geno.pop(0)
    # エリート集団と子孫集団を次世代集団を次世代へ追加します
    next_generation_geno.extend(ga_elite)
    next_generation_geno.extend(ga_progeny)
    return next_generation_geno


def mutation(ga, induvidual_mutation, genom_mutation):
    """突然変異関数です。
    :param ga: genomClass
    :return: 突然変異処理をしたgenomClassを返す"""
    ga_list = []
    for i in ga:
        # 個体に対して一定の確率で突然変異が起きる
        if induvidual_mutation > (random.randint(0, 100) / Decimal(100)):
            genom_list = []
            for i_ in i.getGenom():
                # 個体の遺伝子情報一つ一つに対して突然変異がおこる
                if genom_mutation > (random.randint(0, 100) / Decimal(100)):
                    genom_list.append(random.randint(0, 1))
                else:
                    genom_list.append(i_)
            i.setGenom(genom_list)
            ga_list.append(i)
        else:
            ga_list.append(i)
    return ga_list


if __name__ == '__main__':

    # 一番最初の現行世代個体集団を生成します。
    current_generation_individual_group = []
    for i in range(MAX_GENOM_LIST):
        current_generation_individual_group.append(create_genom(GENOM_LENGTH))

    for count_ in range(1, MAX_GENERATION + 1):
        # 現行世代個体集団の遺伝子を評価し、genomClassに代入します
        for i in range(MAX_GENOM_LIST):
            evaluation_result = evaluation(current_generation_individual_group[i])
            current_generation_individual_group[i].setEvaluation(evaluation_result)
        # エリート個体を選択します
        elite_genes = select(current_generation_individual_group,SELECT_GENOM)
        # エリート遺伝子を交叉させ、リストに格納します
        progeny_gene = []
        for i in range(0, SELECT_GENOM):
            progeny_gene.extend(crossover(elite_genes[i - 1], elite_genes[i]))
        # 次世代個体集団を現行世代、エリート集団、子孫集団から作成します
        next_generation_individual_group = next_generation_gene_create(current_generation_individual_group,
                                                                       elite_genes, progeny_gene)
        # 次世代個体集団全ての個体に突然変異を施します。
        next_generation_individual_group = mutation(next_generation_individual_group,INDIVIDUAL_MUTATION,GENOM_MUTATION)

        # 1世代の進化的計算終了。評価に移ります

        # 各個体適用度を配列化します。
        fits = [i.getEvaluation() for i in current_generation_individual_group]

        # 進化結果を評価します
        min_ = min(fits)
        max_ = max(fits)
        avg_ = sum(fits) / Decimal(len(fits))

        # 現行世代の進化結果を出力します
        print "-----第{}世代の結果-----".format(count_)
        print "  Min:{}".format(min_)
        print "  Max:{}".format(max_)
        print "  Avg:{}".format(avg_)

        # 現行世代と次世代を入れ替えます
        current_generation_individual_group = next_generation_individual_group

    # 最終結果出力
    print "最も優れた個体は{}".format(elite_genes[0].getGenom())