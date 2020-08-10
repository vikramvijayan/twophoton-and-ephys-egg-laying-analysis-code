function [out] = patch_errorbar(data_mean, std_or_sem, x, color)

hold on;

lo = data_mean - std_or_sem;
hi = data_mean + std_or_sem;

keepIndex = isnan(lo);
lo(keepIndex) = 0;

keepIndex = isnan(hi);
hi(keepIndex) = 0;


hp = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)],color,'FaceAlpha',.2,'EdgeColor','none');

%hp = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)]);
end
