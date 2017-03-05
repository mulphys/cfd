template	<class Element>
void	Domain::alloc
(	int	type,
	Element	*&element
)
{	int	ne;
	if (type<0||type>=maxelements)
	{	fprintf(stderr,"Unknown element type no. %d\n",type);
		exit(1);
	}
	ne=Ne[type];
	if(option.debug)
	{	printf
		(	"Allocating %d*%d=%d bytes for %s\n",
			ne,sizeof(Element),
			ne*sizeof(Element),
			elementype[type]
		);FLUSH;
	}
	if (ne>=0)
	{	ALLOC(element,Element,ne);}
	else
		element=NULL;
}
template	<class Element>
void	Domain::alloc
(	int	ne,
	int	type,
	Element	*&root
)
{	if (type<0||type>=maxelements)
	{	fprintf(stderr,"Unknown element type no. %d\n",type);
		exit(1);
	}
//	ne=Ne[type];
	if(option.debug)
	{	printf
		(	"Allocating %d*%d=%d bytes for %s\n",
			ne,sizeof(Element),
			ne*sizeof(Element),
			elementype[type]
		);FLUSH;
	}
	if(ne>0)
	{	Element	*element;
		if(root==NULL)
		{	root=new Element;
			root->next=root->prev=root;
			ne--;
		}
		element=root;
		for(int i=0;i<ne;i++)
		{	Element *next=new Element;
			insert(element,next);
			element=next;
		}
	}
}
template <class Element>
void setZero	
(	int	loc,
 	int	size,
 	Element	*root
)
{
	if (root!=NULL) 
	{
		Element
			*element=root;
		if (option.verbose)
		{	printf
			(	"Initializing %d*%d bytes for variables at each element ... ",
				size,sizeof(double)
			);FLUSH;
		}
		do
		{	double	*var=element->var+loc;
			for (int i=0; i<size; i++)
				var[i]=0.0;
			element=element->next;
		}	while(element!=root);
		if (option.verbose)printf("done\n");
	}
}
template	<class Element>
void	Domain::insert
(	Element *&current,
	Element	*&inserted
//	Element	*&first
	//Element	*&last
)
{	Element *p,*pp;
	p=current->next;
	current->next=inserted;
	inserted->prev=current;
	inserted->next=p;
	p->prev=inserted;
	//last=first->prev;
}
template	<class Element>
void	Domain::insert_before
(	Element *&current,
	Element	*&inserted
//	Element	*&first
	//Element	*&last
)
{	Element *p,*pp;
	p=current->prev;
	current->prev=inserted;
	inserted->prev=p;
	inserted->next=current;
	p->next=inserted;
	//last=first->prev;
}
template	<class Element>
void	Domain::deleteList
(	Element *&head
)
{	if(head==NULL) return;
	if(head->next!=NULL) deleteList(head->next);
	delete head;
}
