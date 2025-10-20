module global_kdtree

    use m_KdTree, only: KdTree

    implicit none

    type(KdTree), save :: tree
    logical :: is_built = .false.

end module global_kdtree