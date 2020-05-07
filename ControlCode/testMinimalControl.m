
xt = allControlTrajectories_emotionrec.stateTrajectories_persistence{10};
ut = allControlTrajectories_emotionrec.controlInputs_persistence{10};
xf = allControlTrajectories_emotionrec.xf{10};

figure; hold on;
for i = 1:1001
    firstTerm = sumsqr((xt(i, :)-xf'));
    secondTerm = sumsqr(ut(i, :));
    plot(i, firstTerm, 'r.', 'MarkerSize', 20);
    plot(i, secondTerm, 'b.', 'MarkerSize', 20);
end
