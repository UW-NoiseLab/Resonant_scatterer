import lumapi
import numpy as np

def is_small_value(val, max_len=10):
    """只保留小规模数据"""
    if isinstance(val, (int, float, bool, str)):
        return True
    if isinstance(val, np.ndarray):
        return val.size <= max_len
    return False


def dump_project(fdtd):
    fdtd.eval("names = getobjects;")
    names = fdtd.getv("names")

    data = {}

    for name in names:
        fdtd.eval(f'props = getproperties("{name}");')
        props = fdtd.getv("props")

        obj_data = {}

        for p in props:
            # ❌ 过滤明显没用的字段
            if p.lower() in ["vertices", "points", "data", "dataset"]:
                continue

            try:
                val = fdtd.getObjectById(name)[p]

                # ❌ 跳过大数据
                if not is_small_value(val):
                    continue

                # ❌ 跳过超长字符串
                if isinstance(val, str) and len(val) > 50:
                    continue

                obj_data[p] = val

            except:
                pass

        # ❗ 如果这个 object 没有有效参数就跳过
        if len(obj_data) > 0:
            data[name] = obj_data

    return data


if __name__ == "__main__":
    fdtd = lumapi.FDTD()
    fdtd.newproject()
    print(fdtd.ge)
    # data = dump_project(fdtd)

    # import json
    # with open("Scatterer_lsf_gen.json", "w") as f:
    #     json.dump(data, f, indent=2)