import maya.cmds as cmds

selectionList = cmds.ls(selection=True, type='transform')

if len(selectionList) >= 1:
    startTime = cmds.playbackOptions(query=True, minTime=True)
    endTime = cmds.playbackOptions(query=True, maxTime=True)
    for objectName in selectionList:
        cmds.cutKey(objectName, time=(startTime, endTime), attribute='rotate')
        cmds.setKeyframe(objectName, time=startTime, attribute='rotateY', value=0)
        cmds.setKeyframe(objectName, time=endTime, attribute='rotateY', value=360)
else:
    print('Please select at least one object')
