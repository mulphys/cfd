#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "main.h"
#include "io.h"
#include "input.h"

//int	getrank
//(
//	char	*word
//)
//{	int	rank;
//	char	*p;
//	//count variable's rank
//	for (p=word;!isspace(*p)&&*p!='['&&*p!='='&&*p!='\0';p++);
//	if (*p!='[')return 0;
//	rank=1;
//	for(;!isspace(*p)&&*p!=']'&&*p!='(';p++)
//		if (*p==',')rank++;
//	if (*p!=']') ERROR1("No closing ']' in",word);	
//	return rank;
//}
//void	setrank
//(
//	char	*word,
//	int	*Rank,
//	int	**ind
//)
//{	int	rank;
//	char	*p,*s;
//	//count variable's rank
//	for (p=word;!isspace(*p)&&*p!='['&&*p!='='&&*p!='\0';p++);
//	if (*p!='[')
//	{	*Rank=0;
//		*p='\0';
//		return;
//	}
//	rank=1;
//	*p++='\0';
//	s=p;
//	while(!isspace(*p)&&*p!=']'&&*p!='(')p++;
//	if (*p!=']') ERROR1("No closing ']' in",word);	
//	*p='\0';
//	p=s;
//	while ((p=strchr(++p,','))!=NULL) rank++;
//	//store indexes
//	ALLOC(*ind,int,rank);
//	for (int i=0; i<rank; i++)
//	{
//		*(*ind+i)=(int)*s;
//		while(*s!=','&&*s!='\0')s++;s++;
//	}
//	*Rank=rank;
//}
//void	setind
//(
//	char	*word,
//	int	rank,
//	int	**ind
//)
//{	int	i;
//	char	*p,*s;
//	//count variable's rank
//	if (rank==0) 
//	{	*ind=NULL;
//		return;
//	}
//	for (p=word;!isspace(*p)&&*p!='['&&*p!='='&&*p!='\0';p++);
//	if (*p!='[')
//	{
//		fprintf
//		(	stderr,
//			"Specified rank=%d to big for the word='%s'\n",
//			rank,word
//		);exit(1);
//		//*ind==NULL;
//		//*p='\0';
//		return;
//	}
//	p++;
//	s=p;//right after '['
//	i=1;
//	for (;!isspace(*p)&&*p!=']'&&*p!='(';p++)
//		if (*p==',') i++;
//	if (i!=rank) ERROR1("Rank argument in setind does not match one in",word);
//	if (*p!=']') ERROR1("No closing ']' in",word);	
//	//store indexes
//	ALLOC(*ind,int,rank);
//	for (int i=0; i<rank; i++)
//	{
//		*(*ind+i)=(int)*s;
//		while(*s!=','&&*s!=']')s++;s++;
//	}
//}
namespace Input
{
int getCharAttr(xmlNodePtr cur, char *attr, char *result) 
{
	xmlChar *key;
	key = xmlGetProp(cur, (const xmlChar *)attr);
	if(key!=NULL)
	{	SCOPY(result,(char*)key);
		xmlFree(key);
		return 1;
	}
	if(option.verbose)
	{	fprintf
		(	stderr,
			"\tAttribute %s not found\n",
			attr
		);
	}
	return 0;
}
int getIntAttr(xmlNodePtr cur, char *attr) 
{
	xmlChar *key;
	key = xmlGetProp(cur, (const xmlChar *)attr);
	if(key!=NULL)
	{	int result;
		result=atoi((char*)key);
		xmlFree(key);
		return result;
	}
	if(option.verbose)
	{	fprintf
		(	stderr,
			"\tAttribute %s not found\n",
			attr
		);
	}
	return -1;
}
int getIntAttr
(	xmlDocPtr doc, 
	char *tag, 
	char *keyname, 
	char *key, 
	char *attr
) 
{
	xmlNodePtr cur;
///	xmlChar *key;
	cur = xmlDocGetRootElement(doc);
	cur = cur->xmlChildrenNode;
	if(option.verbose) 
	{	printf
		(	"getIntAttr:tag=%s, keyname=%s, key=%s, attr=%s\n",
			tag,keyname,key,attr
		);	fflush(stdout);
	}
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)tag)))
		{	char keyval[MAXLINLEN];
			if(getCharAttr(cur,keyname,keyval)!=0)
			{	if(strcmp(keyval,key)==0)
				{	int	value=0;
					///val=getIntAttr(cur,attr);
					value=parseInt(doc,cur,attr);
					return value;
				}
			}
		}
		cur=cur->next;
	}
	return -1;
}
int parseWord(xmlDocPtr doc, xmlNodePtr cur, char *keyword, char *result) 
{
	xmlChar *key;
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)keyword))) 
		{	key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			SCOPY(result,(char*)key);
			xmlFree(key);
			return 1;
		}
		cur = cur->next;
	}
	if(option.verbose)
	{	fprintf
		(	stderr,
			"\tTag '%s' not found\n",
			keyword
		);
	}
	return 0;
}
int parseInt(xmlDocPtr doc, xmlNodePtr cur, char *keyword) 
{
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)keyword))) 
		{	xmlChar *key;
			int keyval;
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			keyval=(int)atoi((char*)key);
			xmlFree(key);
			return keyval;
		}
		cur = cur->next;
	}
	if(option.verbose)
	{
		fprintf(stderr,"Integer data for tag '%s' not found\n",keyword);
		fflush(stderr);
	}
	return -1;
}
int parseFloat(xmlDocPtr doc, xmlNodePtr cur, char *keyword, double &result) 
{
	xmlChar *key;
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)keyword))) 
		{
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			//printf("%s: %s\n", keyword,key);///DDD
			if(key!=NULL) result=(double)atof((char*)key);
			xmlFree(key);
			return 1;
		}
		cur = cur->next;
	}
	return 0;
}
double parseFloat(xmlDocPtr doc, xmlNodePtr cur, char *keyword) 
{
	xmlChar *key;
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)keyword))) 
		{	double result;
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			//printf("%s: %s\n", keyword,key);///DDD
			result=(double)atof((char*)key);
			xmlFree(key);
			return result;
		}
		cur = cur->next;
	}
	return 0.0;
}
void	getTimeCFG
(
	char inpfilename[]
)
{	int	level=0,timefound=0;
	char	*p,word[MAXLINLEN];
	FILE	*inp;
	//Default times
	runtime.start=0.0;
	runtime.end=1.0;
	runtime.step0=runtime.step0=0.01;
	//get time
	OPENREAD(inpfilename,inp);
	while (!feof(inp)&&timefound==0&&level>=0)
	{	GETWORD(word,inp);
		if (*word=='!') {SKIPLINE(inp);continue;}
		if (*word=='{') level++;
		if (*word=='}') level--;
		if (level==0&&strcmp(word,TIME_KEYWORD)==0)
		{	timefound=1;
			GETWORD(word,inp);
			while (!feof(inp)&&*word!='{')
				GETWORD(word,inp);
			if (*word!='{') 
			{	fprintf
				(	stderr,
					"Keyword %d should be followed by {",TIME_KEYWORD
				);exit(1);
			}
			level++;
			while (!feof(inp)&&level>=1)
			{	GETWORD(word,inp);
				if (*word=='!') {SKIPLINE(inp);continue;}
				if (*word=='{') level++;
				if (*word=='}') level--; 
				if (level==0) break;
				if (level==1)
				{	if (strcmp(word,TIME_START_KEYWORD)==0)
						GETREAL(runtime.start,inp);
					if (strcmp(word,TIME_END_KEYWORD)==0)
						GETREAL(runtime.end,inp);
					if (strcmp(word,TIME_STEP_KEYWORD)==0)
					{	GETREAL(runtime.step,inp);
						runtime.step0=runtime.step;
					}
				}
			}
		}
	}
	runtime.current=runtime.prev=runtime.start;
	if (!timefound)
		printf("Keyword %s not found in %d - using defaults\n");
	if (option.verbose)
		printf
		(	"Start time = %g\nEnd time = %g\nTime step = %g\n",
			runtime.start,runtime.end,runtime.step
		);
	fclose(inp);
}
void	getTimeXML
(
	char inpfilename[]
)
{//?	extern double parseFloat(xmlDocPtr doc, xmlNodePtr cur, char *keyword);
	int	level=0,timefound=0;
	char	*p,word[MAXLINLEN];
	FILE	*inp;
	const char *doctype=DOCTYPE;
	xmlDocPtr doc;
	xmlNodePtr cur;
	//Default times
	runtime.start=0.0;
	runtime.end=1.0;
	runtime.step0=runtime.step0=0.01;
	doc = xmlParseFile(inpfilename);
	if (doc == NULL ) 
	{	fprintf(stderr,"XML parser failed in %s\n",inpfilename);
		return;
	}
	cur = xmlDocGetRootElement(doc);
	if (cur == NULL) 
	{	fprintf(stderr,"Empty document: %s\n",inpfilename);
		xmlFreeDoc(doc);
		return;
	}
	if (xmlStrcmp(cur->name, (const xmlChar *) doctype)) 
	{	fprintf(stderr,"document of the wrong type, root node != %s in %s\n",doctype,inpfilename);
		xmlFreeDoc(doc);
		return;
	}
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)TIME_KEYWORD)))
		{	runtime.start=parseFloat(doc,cur,(char*)TIME_START_KEYWORD);
			runtime.end=parseFloat(doc,cur,(char*)TIME_END_KEYWORD);
			runtime.step=parseFloat(doc,cur,(char*)TIME_STEP_KEYWORD);
		}
		cur = cur->next;
	}
	xmlFreeDoc(doc);
}
void	getTime
(
	char inpfilename[]
)
{
	if(strstr(inpfilename,".xml\0")!=NULL)
		getTimeXML(inpfilename);
	else
		getTimeCFG(inpfilename);
	if (option.verbose)
		printf("Time: start=%g, end=%g, step=%g\n",runtime.start,runtime.end,runtime.step);
}
int	countObjectsXML
(	char	*keyword,
	char	*filename
)
{	int	count;
	xmlDocPtr doc;
	xmlNodePtr cur;
	doc = xmlParseFile(filename);
	if (doc == NULL ) 
	{	fprintf(stderr,"Document not parsed successfully in %s\n",filename);
		return 0;
	}
	cur = xmlDocGetRootElement(doc);
	if (cur == NULL) 
	{	fprintf(stderr,"empty document in %s\n",filename);
		xmlFreeDoc(doc);
		return 0;
	}
	if (xmlStrcmp(cur->name, (const xmlChar *) DOCTYPE)) 
	{	fprintf
		(	stderr,"document of the wrong type, root node != %s in %s\n",
			DOCTYPE,filename
		);
		xmlFreeDoc(doc);
		return 0;
	}
	cur = cur->xmlChildrenNode;
	count=0;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)keyword)))
			count++;
		cur = cur->next;
	}
	xmlFreeDoc(doc);
	return count;
}
int	countObjectsCFG
(	char	*keyword,
	char	*filename
)
{	int	count,wordlen=strlen(keyword);
	char	*s,line[MAXLINLEN];
	FILE	*inp;
	OPENREAD(filename,inp);
	count=0;
	while ((s=fgets(line,MAXLINLEN,inp))!=NULL)
	{
		while(isspace(*s))s++;
		if (memcmp(s,keyword,wordlen)==0&&isspace(*(s+wordlen))) count++;
	}
	fclose(inp);
	return count;
}
int	countObjects
(	char	*keyword,
	char	*filename
)
{	int count;
	if(strstr(filename,".xml\0")!=NULL)
		count=countObjectsXML(keyword,filename);
	else
		count=countObjectsCFG(keyword,filename);
	if (option.verbose)
		printf("Number of domains: count=%d\n",count);
	return count;
}
} //END namespace Input
