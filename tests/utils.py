class MockDdbClient:
    def __init__(self, table_name, items, key_name):
        self.items = items
        self.table_name = table_name
        self.key_name = key_name

    def batch_get_item(self, RequestItems):
        if self.table_name not in RequestItems:
            return {}
        keys = RequestItems[self.table_name]["Keys"]
        result_list = []
        for item in self.items:
            item_key = {self.key_name: item[self.key_name]}
            if item_key in keys:
                result_list.append(item)
        result = {"Responses": {self.table_name: result_list}}
        return result
