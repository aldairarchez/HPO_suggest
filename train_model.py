import click as ck
from gensim import models
from gensim import utils

class HPOCorpus:
    """An iterator that yields sentences (lists of str)."""

    def __init__(self, data_file):
        self.data_file = data_file
        
    def __iter__(self):
        with open(self.data_file) as f:
            for line in f:
                # assume there's one document per line, tokens separated by whitespace
                yield line.strip().split(' ')


@ck.command()
def main():
    sentences = HPOCorpus('data/corpus.txt')
    model = models.Word2Vec(
        sentences=sentences,
        vector_size=128,
        min_count=1,
        window=20,
        compute_loss=True,
        epochs=50,
        workers=4)
    print(model.wv.most_similar(['HP:0001596', 'HP:0000951'], topn=10))
    model.save('data/hpo.wv')


if __name__ == '__main__':
    main()
