function IncludedCells = grabFromFilter(str, varargin)

h = batchAnalyses_OKR_Physiology;
IncludedCells = h.outsideCall(str, varargin);
close(h.UIFigure)
end