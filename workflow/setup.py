
df=pd.read_table("config/samples.tsv")
sm=pd.read_table("config/sm.txt", header=None)

sam=df[df['sample'].str.startswith("SAMN")].reset_index(drop=True)
bgi=df[~df['sample'].str.startswith("SAMN")]

sm.columns=['SAM_TAG']
sam['SAM_TAG']=sm['SAM_TAG']

bgi['SAM_TAG']=bgi['sample']
df=pd.concat([bgi, sam])

df=df[['sample','SAM_TAG','bam_path']]

df.to_csv("config/samples.tsv", index=False, sep='\t')

#df=df.drop(columns='bai_path')
#df['bam_path']="/lustre/home/bovo/RABBIT_BGI/BAM/BAM_RMDUP/"
#df['bam_path']=df['bam_path'].where(~df['sample'].str.startswith("SAM"), "/lustre/home/bovo/BAM_RABBIT_RMDUP/")
#df['bam_path']=df['bam_path']+df['sample']+".rmdup.bam"
#df.to_csv("config/samples.tsv", index=False, sep='\t')
