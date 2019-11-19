function triplesave(handle,filename)
if isgraphics(handle,'Axis')
    saveas(handle.Parent,strcat(filename,'.fig'));
    saveas(handle.Parent,strcat(filename,'.png'));
    saveas(handle.Parent,strcat(filename,'.jpg'));
elseif isgraphics(handle,'Figure')
    saveas(handle,strcat(filename,'.fig'));
    saveas(handle,strcat(filename,'.png'));
    saveas(handle,strcat(filename,'.jpg'));
else
    error('handle argument must be either a ''Axis'' or a ''Figure'' object');
end 
end