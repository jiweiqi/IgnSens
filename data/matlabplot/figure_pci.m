function H = figure_pci(PW, PL)

if nargin == 0
    PW=6.7;
    PL=6;
end

H=figure();
set(H,'PaperPositionMode','manual');
set(H,'PaperUnits','centimeters');
set(H,'PaperPosition',[0 0 PW PL]);
set(H,'Units','centimeters');
set(H,'Position',[10 10 PW PL]);

end