import pandas as pd
from wordcloud import WordCloud
from wordcloud import ImageColorGenerator
from wordcloud import STOPWORDS
import matplotlib.pyplot as plt


df = pd.read_pickle('oxphos_merged_aids_cids_clnsmi_murcko_assaydesc_ETCassay_pmid.pkl')

# get abstract text for different sets of AIDs
text1 = " ".join( df['Abstract'].drop_duplicates( keep='first' ) )
text2 = " ".join( df.loc[ df['ETC_linked_AID'] == 1, 'Abstract' ].drop_duplicates( keep='first' ) )
text3 = " ".join( df.loc[ df['ETC_linked_PMID'] == 1, 'Abstract' ].drop_duplicates( keep='first' ) )
textx = " ".join( df.loc[ df['ETC_linked_AID'] == 0, 'Abstract' ].drop_duplicates( keep='first' ) )

stopwords = set(STOPWORDS)

plt.figure( figsize=(15,10) )
wordcloud1 = WordCloud( width=1500, height=1000, stopwords=stopwords, background_color='white').generate(text1)
plt.imshow( wordcloud1, interpolation='bilinear' )
plt.axis("off")
plt.savefig( 'wordcloud_all_abstracts.png', dpi=1200 )
plt.close()

plt.figure( figsize=(15,10) )
wordcloud2 = WordCloud( width=1500, height=1000, stopwords=stopwords, background_color='white').generate(text2)
plt.imshow( wordcloud2, interpolation='bilinear' )
plt.axis("off")
plt.savefig( 'wordcloud_ETC-linked_abstracts.png', dpi=1200 )
plt.close()

plt.figure( figsize=(15,10) )
wordcloud3 = WordCloud( width=1500, height=1000, stopwords=stopwords, background_color='white').generate(text3)
plt.imshow( wordcloud3, interpolation='bilinear' )
plt.axis('off')
plt.savefig( 'wordcloud_ETC-linked_abstracts_pmids.png', dpi=1200 )
plt.close()

plt.figure( figsize=(15,10) )
wordcloudx = WordCloud( width=1500, height=1000, stopwords=stopwords, background_color='white').generate(textx)
plt.imshow( wordcloudx, interpolation='bilinear' )
plt.axis('off')
plt.savefig( 'wordcloud_non-ETC-linked_abstracts.png', dpi=1200 )
plt.close()

