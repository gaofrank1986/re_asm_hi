
# for each node on free surface
		# assign this node as glbl src point
		# calculate solid angle
		#  A_MATRIX(INODE,INODE,IS)  = solid_angle
				# IS is symetric parts
		# for each element on free surface
			# if this src node is not on element
				#  do normal integration get amatrix,bmatrix
			# else 
				# do singular integration

			#  get normal vector and pos on each node
				# get partial (potentail)/partial(noraml vector


# for each node on body surface
		# assign this node as glbl src point
		# calculate solid angle
		#  A_MATRIX(INODE,INODE,IS)  = solid_angle
				# IS is symetric parts

		# for each element on free surface
			# call norm0 integration
			
