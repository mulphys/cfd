
namespace Input
{
int getCharAttr(xmlNodePtr cur, char *attr, char *result);
int getIntAttr(xmlNodePtr cur, char *attr);
int getIntAttr
(	xmlDocPtr doc, 
	char *tag, 
	char *keyname, 
	char *key, 
	char *attr
); 
int parseWord(xmlDocPtr doc, xmlNodePtr cur, char *keyword, char *result); 
int parseInt(xmlDocPtr doc, xmlNodePtr cur, char *keyword);
double parseFloat(xmlDocPtr doc, xmlNodePtr cur, char *keyword);
int parseFloat(xmlDocPtr doc, xmlNodePtr cur, char *keyword, double &result);
void	getTimeCFG(char inpfilename[]);
void	getTimeXML(char inpfilename[]);
void	getTime(char inpfilename[]);
int	countObjectsXML(char	*keyword, char	*filename);
int	countObjectsCFG(char	*keyword, char	*filename);
int	countObjects(char	*keyword, char	*filename);
}

