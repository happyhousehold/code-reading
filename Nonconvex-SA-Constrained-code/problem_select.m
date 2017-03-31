mtitle = 'Solving semi-supervised supprot vector machine problem';
problem_selection = menu (mtitle, '                  Start                ', '                 Quit                 ');
if problem_selection == 1
    data.problem ='S3VM';
else
    data.problem ='';
end