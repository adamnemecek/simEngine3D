#allow the addition of points to bodies through the construction of joints,
# so the joint constructor would take a point, determine if that body has the point
#already, and if it does, then it register the ID appropriately, but if it doesn't
#then it adds the point to the list. This relieves the user of th burden of keeping
#track of the ID's of the points that have been added to that body. points should belong
#to the body though, conceptually that is true.
