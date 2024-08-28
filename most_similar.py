import click as ck
from gensim import models
from gensim import utils

@ck.command()
def main():
    model = models.Word2Vec.load('data/hpo.wv')
    print(model.wv.most_similar(['HP:0001596', 'HP:0000951'], topn=10))
    

if __name__ == '__main__':
    main()
