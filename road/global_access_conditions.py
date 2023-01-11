def is_owner(request, view, action: str) -> bool:
    return request.user == view.get_object().owner
