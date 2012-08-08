function tail = htmltail(addline)
tail{1} = ' </body>';
tail{2} = ' </html>';
if nargin == 1
    tail = horzcat(addline, tail);
end
    