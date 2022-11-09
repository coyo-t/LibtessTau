from node import *

def __main_node ():
	class TN(Node):
		def __init__ (self, key = 'None'):
			super().__init__()
			self.key: str = key

	thing = DxDList[TN]()
	thing.start = TN('A haha')
	thing.start >> TN('B wow!!!') >> TN('C') >> TN('D')

	head = thing.head

	for node in reversed(+head):
		print(node.key)

	# for node in thing:
	# 	print(node.key)


if __name__ == '__main__':
	__main_node()

