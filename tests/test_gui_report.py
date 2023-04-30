from rnalysis.gui.gui_report import *


def test_Node():
    # Test creating a node object
    predecessors = [1, 2, 3]
    node = Node(4, 'Node 4', predecessors, 'popup', 'Count matrix', 'filename.txt')
    assert node.node_id == 4
    assert node.node_name == 'Node 4 (#4)'
    assert node.predecessors == set(predecessors)
    assert node.popup_element == 'popup'
    assert node.node_type == 'Count matrix'
    assert node.is_active is True
    assert node.filename == 'filename.txt'

    # Test setting and getting is_active property
    node.set_active(False)
    assert node.is_active is False
    node.set_active(True)
    assert node.is_active is True

    # Test adding a predecessor
    node.add_predecessor(5)
    assert node.predecessors == {1, 2, 3, 5}
