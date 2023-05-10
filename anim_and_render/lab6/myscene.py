import maya.cmds as cmds

# Создание кубика
cube = cmds.polyCube()[0]

# Назначение имени объекту
cmds.rename(cube, 'MyCube')

cmds.move(-0.5, 0.5, -0.5, 'MyCube')

# Установка ключевых кадров
cmds.currentTime(0)
cmds.setKeyframe('MyCube', attribute='rotateX', value=0)
cmds.setKeyframe('MyCube', attribute='rotateY', value=0)
cmds.setKeyframe('MyCube', attribute='rotateZ', value=0)
cmds.xform('MyCube', piv=(0.5, -0.5, 0))

cmds.currentTime(30)
cmds.xform('MyCube', piv=(0.5, -0.5, 0))
cmds.setKeyframe('MyCube', attribute='rotateX', value=0)
cmds.setKeyframe('MyCube', attribute='rotateY', value=0)
cmds.setKeyframe('MyCube', attribute='rotateZ', value=-90)

cmds.currentTime(60)
cmds.xform('MyCube', piv=(0, 0.5, 0))
cmds.setKeyframe('MyCube', attribute='rotateX', value=-0)
cmds.setKeyframe('MyCube', attribute='rotateY', value=90)
cmds.setKeyframe('MyCube', attribute='rotateZ', value=-90)

cmds.currentTime(90)
cmds.xform('MyCube', piv=(0, -0.5, 0.5))
cmds.setKeyframe('MyCube', attribute='rotateX', value=-90)
cmds.setKeyframe('MyCube', attribute='rotateY', value=90)
cmds.setKeyframe('MyCube', attribute='rotateZ', value=-90)

cmds.currentTime(120)
cmds.xform('MyCube', piv=(-0.5, -0.5, 0))
cmds.setKeyframe('MyCube', attribute='rotateX', value=-180)
cmds.setKeyframe('MyCube', attribute='rotateY', value=90)
cmds.setKeyframe('MyCube', attribute='rotateZ', value=-90)

cmds.currentTime(150)
cmds.xform('MyCube', piv=(-0.5, 0, -0.5))
cmds.setKeyframe('MyCube', attribute='rotateX', value=-270)
cmds.setKeyframe('MyCube', attribute='rotateY', value=90)
cmds.setKeyframe('MyCube', attribute='rotateZ', value=-90)

cmds.currentTime(180)
cmds.xform('MyCube', piv=(0, 0, -0.5))
cmds.setKeyframe('MyCube', attribute='rotateX', value=-270)
cmds.setKeyframe('MyCube', attribute='rotateY', value=180)
cmds.setKeyframe('MyCube', attribute='rotateZ', value=-90)

cmds.currentTime(210)
cmds.xform('MyCube', piv=(0, -0.5, -0.5))
cmds.setKeyframe('MyCube', attribute='rotateX', value=-270)
cmds.setKeyframe('MyCube', attribute='rotateY', value=270)
cmds.setKeyframe('MyCube', attribute='rotateZ', value=-90)

cmds.currentTime(240)
cmds.xform('MyCube', piv=(0.5, -0.5, -0.5))
cmds.setKeyframe('MyCube', attribute='rotateX', value=-270)
cmds.setKeyframe('MyCube', attribute='rotateY', value=360)
cmds.setKeyframe('MyCube', attribute='rotateZ', value=-90)

cmds.currentTime(270)
cmds.xform('MyCube', piv=(0.5, 0.5, 0.5))
cmds.setKeyframe('MyCube', attribute='rotateX', value=-360)
cmds.setKeyframe('MyCube', attribute='rotateY', value=360)
cmds.setKeyframe('MyCube', attribute='rotateZ', value=-90)
