#ifndef SAVE_H
#define SAVE_H

extern void DeleteFile(char* fname);
/*
	Function to delete a file, if it exists and is allowed to be deleted.
	
	Input: name of a file (fname)
	Output: Nothing -- File with name fname is deleted using the system command "rm -f".
	Algorithm: Uses access(fname, R_OK|F_OK) of unistd.h to test various perperties of the file.
*/

extern void Save(struct tiling, char* fname);
/*
	Input: pointer to a Tiling (pG)
		   the name of the file to which the tiling needs to be saved
	Output: Nothing
	Algorithm: Write the edges, faces and the indices of open edges of a tiling, into a file.
			   The file names are: Edges ---  (pG->type)_edges.txt
			   					   Faces ---  (pG->type)_faces.txt
			   					   Open Edges ---  (pG->type)_doEdges.txt
*/

#endif	/* SAVE_H */
