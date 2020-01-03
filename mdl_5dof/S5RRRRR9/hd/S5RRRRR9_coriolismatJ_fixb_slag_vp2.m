% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:50
% EndTime: 2019-12-31 22:28:26
% DurationCPUTime: 18.29s
% Computational Cost: add. (28120->699), mult. (59384->955), div. (0->0), fcn. (62214->8), ass. (0->381)
t777 = qJD(3) + qJD(4);
t410 = sin(qJ(4));
t411 = sin(qJ(3));
t414 = cos(qJ(4));
t415 = cos(qJ(3));
t534 = t414 * t415;
t366 = -t410 * t411 + t534;
t367 = -t410 * t415 - t411 * t414;
t409 = sin(qJ(5));
t413 = cos(qJ(5));
t492 = t366 * t413 + t367 * t409;
t710 = Ifges(6,5) * t492;
t279 = t366 * t409 - t367 * t413;
t730 = Ifges(6,6) * t279;
t526 = t710 - t730;
t681 = -pkin(8) - pkin(7);
t385 = t681 * t411;
t386 = t681 * t415;
t306 = t385 * t410 - t386 * t414;
t257 = pkin(9) * t366 + t306;
t705 = t385 * t414 + t386 * t410;
t725 = pkin(9) * t367 + t705;
t758 = -t257 * t409 + t413 * t725;
t763 = t758 * mrSges(6,2);
t140 = t257 * t413 + t409 * t725;
t766 = t140 * mrSges(6,1);
t769 = -t766 / 0.2e1 - t763 / 0.2e1;
t775 = 0.2e1 * t769 + t526;
t776 = t775 * qJD(5);
t585 = t279 * Ifges(6,4);
t159 = Ifges(6,2) * t492 + t585;
t638 = pkin(3) * t415;
t400 = -pkin(2) - t638;
t614 = Ifges(5,4) * t367;
t615 = Ifges(5,4) * t366;
t664 = t279 / 0.2e1;
t668 = -t279 / 0.2e1;
t715 = -Ifges(5,2) + Ifges(5,1);
t717 = t492 / 0.2e1;
t274 = Ifges(6,4) * t492;
t161 = t279 * Ifges(6,1) + t274;
t720 = -Ifges(6,2) * t279 + t274;
t738 = t720 + t161;
t740 = Ifges(6,1) * t492 - t585;
t322 = -pkin(4) * t366 + t400;
t722 = mrSges(6,1) * t279 + mrSges(6,2) * t492;
t745 = t322 * t722;
t774 = (-t400 * mrSges(5,1) - t614) * t367 + (mrSges(5,2) * t400 - t367 * t715 + t615) * t366 + t159 * t668 + t745 + t740 * t664 + t738 * t717;
t679 = -t140 / 0.2e1;
t157 = -mrSges(6,1) * t492 + mrSges(6,2) * t279;
t773 = m(6) * t322 + t157;
t360 = Ifges(5,6) * t367;
t361 = Ifges(5,5) * t366;
t489 = t361 + t360 + t526;
t753 = -t306 * mrSges(5,1) - mrSges(5,2) * t705 + t489;
t759 = -t763 - t766;
t772 = t753 + t759;
t416 = cos(qJ(2));
t412 = sin(qJ(2));
t340 = t367 * t412;
t579 = t340 * mrSges(5,3);
t307 = mrSges(5,2) * t416 + t579;
t541 = t411 * t412;
t339 = t410 * t541 - t412 * t534;
t580 = t339 * mrSges(5,3);
t309 = -mrSges(5,1) * t416 + t580;
t663 = t306 / 0.2e1;
t493 = t339 * t409 + t340 * t413;
t589 = t493 * mrSges(6,3);
t207 = mrSges(6,2) * t416 + t589;
t678 = -t207 / 0.2e1;
t245 = t339 * t413 - t340 * t409;
t590 = t245 * mrSges(6,3);
t209 = -mrSges(6,1) * t416 + t590;
t754 = t209 * t679;
t430 = -t754 + t758 * t678 - t705 * t307 / 0.2e1 + t309 * t663;
t680 = -t758 / 0.2e1;
t749 = t758 / 0.2e1;
t764 = (t680 + t749) * mrSges(6,2) + (t679 + t140 / 0.2e1) * mrSges(6,1);
t771 = qJD(2) * t764;
t770 = qJD(5) * t764;
t708 = t492 * mrSges(6,3);
t765 = -t140 * t413 + t409 * t758;
t613 = Ifges(6,4) * t245;
t117 = Ifges(6,2) * t493 - Ifges(6,6) * t416 - t613;
t407 = t412 * pkin(6);
t377 = pkin(3) * t541 + t407;
t281 = -pkin(4) * t340 + t377;
t723 = -mrSges(6,1) * t245 + mrSges(6,2) * t493;
t739 = Ifges(6,1) * t493 + t613;
t757 = -(t739 / 0.4e1 - t117 / 0.4e1) * t279 - t281 * t722 / 0.2e1 - t322 * t723 / 0.2e1;
t752 = t740 / 0.4e1 - t159 / 0.4e1;
t646 = t412 / 0.2e1;
t746 = t281 * t723;
t405 = Ifges(4,4) * t415;
t479 = -Ifges(4,1) * t411 - t405;
t643 = t415 / 0.2e1;
t742 = t479 * t643;
t709 = Ifges(6,5) * t493;
t731 = Ifges(6,6) * t245;
t528 = t709 + t731;
t236 = Ifges(6,4) * t493;
t119 = -Ifges(6,1) * t245 - Ifges(6,5) * t416 + t236;
t721 = Ifges(6,2) * t245 + t236;
t737 = t721 + t119;
t497 = t731 / 0.2e1 + t709 / 0.2e1;
t639 = pkin(3) * t414;
t399 = pkin(4) + t639;
t548 = t409 * t410;
t352 = -pkin(3) * t548 + t399 * t413;
t545 = t410 * t413;
t353 = pkin(3) * t545 + t399 * t409;
t358 = (-t409 * t414 - t545) * pkin(3);
t359 = (t413 * t414 - t548) * pkin(3);
t632 = t412 * pkin(7);
t380 = -pkin(2) * t416 - pkin(1) - t632;
t364 = t415 * t380;
t539 = t412 * t415;
t488 = -pkin(8) * t539 + t364;
t280 = (-pkin(6) * t411 - pkin(3)) * t416 + t488;
t531 = t415 * t416;
t321 = pkin(6) * t531 + t380 * t411;
t298 = -pkin(8) * t541 + t321;
t284 = t410 * t298;
t169 = t280 * t414 - t284;
t331 = t339 * pkin(9);
t142 = t169 + t331;
t121 = -pkin(4) * t416 + t142;
t286 = t414 * t298;
t170 = t280 * t410 + t286;
t635 = pkin(9) * t340;
t143 = t170 + t635;
t568 = t143 * t409;
t52 = t121 * t413 - t568;
t567 = t143 * t413;
t53 = t121 * t409 + t567;
t62 = -t142 * t409 - t567;
t63 = t142 * t413 - t568;
t626 = t63 * mrSges(6,2);
t627 = t62 * mrSges(6,1);
t622 = t627 / 0.2e1 - t626 / 0.2e1;
t653 = -t358 / 0.2e1;
t691 = -m(6) / 0.2e1;
t727 = (t352 * t62 + t353 * t63 + t358 * t52 + t359 * t53) * t691 + t209 * t653 + t359 * t678 - t622;
t341 = t367 * t416;
t342 = t366 * t416;
t247 = t341 * t413 - t342 * t409;
t250 = t341 * t409 + t342 * t413;
t461 = -Ifges(6,5) * t250 / 0.2e1 - Ifges(6,6) * t247 / 0.2e1;
t633 = t412 * pkin(2);
t387 = -pkin(7) * t416 + t633;
t323 = pkin(6) * t541 + t387 * t415;
t283 = pkin(3) * t412 - pkin(8) * t531 + t323;
t324 = -pkin(6) * t539 + t387 * t411;
t540 = t411 * t416;
t300 = -pkin(8) * t540 + t324;
t175 = t283 * t414 - t300 * t410;
t134 = pkin(4) * t412 - pkin(9) * t342 + t175;
t176 = t283 * t410 + t300 * t414;
t144 = pkin(9) * t341 + t176;
t60 = t134 * t413 - t144 * t409;
t61 = t134 * t409 + t144 * t413;
t486 = Ifges(6,3) * t646 - t61 * mrSges(6,2) / 0.2e1 + t60 * mrSges(6,1) / 0.2e1 - t461;
t686 = m(6) * pkin(4);
t724 = -t765 * t686 / 0.2e1 - t769;
t716 = t493 / 0.2e1;
t607 = Ifges(4,6) * t411;
t612 = Ifges(4,5) * t415;
t462 = t612 / 0.2e1 - t607 / 0.2e1;
t707 = Ifges(3,4) - t462;
t704 = -Ifges(4,2) * t411 + t405;
t490 = Ifges(5,5) * t340 + Ifges(5,6) * t339 + t528;
t703 = -t323 * t411 + t324 * t415;
t702 = -mrSges(4,1) * t415 + mrSges(4,2) * t411;
t520 = pkin(6) * t540;
t297 = t488 - t520;
t185 = t297 * t414 - t284;
t700 = t185 / 0.2e1 - t169 / 0.2e1;
t354 = t358 * mrSges(6,1);
t578 = t359 * mrSges(6,2);
t698 = (mrSges(5,1) * t410 + mrSges(5,2) * t414) * pkin(3) - t354 + t578;
t630 = t53 * mrSges(6,1);
t631 = t52 * mrSges(6,2);
t697 = -t630 / 0.2e1 - t631 / 0.2e1 + t497;
t210 = mrSges(6,1) * t412 - mrSges(6,3) * t250;
t655 = t352 / 0.2e1;
t208 = -mrSges(6,2) * t412 + mrSges(6,3) * t247;
t677 = t208 / 0.2e1;
t684 = -mrSges(4,2) / 0.2e1;
t685 = mrSges(4,1) / 0.2e1;
t690 = m(6) / 0.2e1;
t695 = (t352 * t60 + t353 * t61) * t690 + t353 * t677 + t210 * t655 + t324 * t684 + t323 * t685;
t694 = -t400 * (-mrSges(5,1) * t339 + mrSges(5,2) * t340) / 0.2e1 - t377 * (-mrSges(5,1) * t367 + mrSges(5,2) * t366) / 0.2e1 - t737 * t492 / 0.4e1 - t738 * t493 / 0.4e1 + t752 * t245 + t757;
t693 = -m(5) / 0.2e1;
t692 = m(5) / 0.2e1;
t689 = pkin(3) / 0.2e1;
t683 = -t52 / 0.2e1;
t682 = -t53 / 0.2e1;
t672 = t247 / 0.2e1;
t670 = -t245 / 0.2e1;
t669 = t250 / 0.2e1;
t667 = -t492 / 0.2e1;
t661 = -t339 / 0.2e1;
t660 = t340 / 0.2e1;
t658 = t341 / 0.2e1;
t657 = t342 / 0.2e1;
t656 = -t352 / 0.2e1;
t654 = -t353 / 0.2e1;
t652 = t359 / 0.2e1;
t651 = t366 / 0.2e1;
t649 = -t409 / 0.2e1;
t648 = -t411 / 0.2e1;
t647 = t411 / 0.2e1;
t645 = -t413 / 0.2e1;
t644 = t413 / 0.2e1;
t642 = -t416 / 0.2e1;
t641 = -t416 / 0.4e1;
t640 = pkin(3) * t411;
t637 = pkin(4) * t367;
t636 = pkin(4) * t409;
t408 = t416 * pkin(6);
t184 = -t297 * t410 - t286;
t148 = t184 - t635;
t149 = t331 + t185;
t66 = t148 * t413 - t149 * t409;
t625 = t66 * mrSges(6,1);
t67 = t148 * t409 + t149 * t413;
t624 = t67 * mrSges(6,2);
t621 = t625 / 0.2e1 - t624 / 0.2e1;
t618 = Ifges(4,4) * t411;
t617 = Ifges(5,4) * t339;
t616 = Ifges(5,4) * t340;
t611 = Ifges(5,5) * t342;
t610 = Ifges(5,5) * t416;
t606 = Ifges(5,6) * t341;
t605 = Ifges(5,6) * t416;
t602 = pkin(3) * qJD(3);
t601 = pkin(4) * qJD(4);
t594 = t169 * mrSges(5,2);
t593 = t170 * mrSges(5,1);
t592 = t184 * mrSges(5,1);
t591 = t185 * mrSges(5,2);
t586 = t279 * mrSges(6,3);
t577 = t366 * mrSges(5,3);
t576 = t367 * mrSges(5,3);
t118 = Ifges(6,4) * t250 + Ifges(6,2) * t247 + Ifges(6,6) * t412;
t120 = Ifges(6,1) * t250 + Ifges(6,4) * t247 + Ifges(6,5) * t412;
t123 = -mrSges(6,1) * t493 - mrSges(6,2) * t245;
t124 = -mrSges(6,1) * t247 + mrSges(6,2) * t250;
t238 = Ifges(5,4) * t342 + Ifges(5,2) * t341 + Ifges(5,6) * t412;
t239 = Ifges(5,1) * t342 + Ifges(5,4) * t341 + Ifges(5,5) * t412;
t251 = -mrSges(5,1) * t340 - mrSges(5,2) * t339;
t252 = -mrSges(5,1) * t341 + mrSges(5,2) * t342;
t378 = pkin(3) * t540 + t408;
t282 = -pkin(4) * t341 + t378;
t308 = -mrSges(5,2) * t412 + mrSges(5,3) * t341;
t310 = mrSges(5,1) * t412 - mrSges(5,3) * t342;
t320 = t364 - t520;
t335 = Ifges(4,6) * t412 + t416 * t704;
t480 = Ifges(4,1) * t415 - t618;
t337 = Ifges(4,5) * t412 + t416 * t480;
t484 = mrSges(4,1) * t411 + mrSges(4,2) * t415;
t357 = t484 * t416;
t373 = mrSges(4,2) * t416 - mrSges(4,3) * t541;
t374 = -mrSges(4,2) * t412 - mrSges(4,3) * t540;
t375 = -mrSges(4,1) * t416 - mrSges(4,3) * t539;
t376 = mrSges(4,1) * t412 - mrSges(4,3) * t531;
t473 = Ifges(5,2) * t340 - t617;
t478 = -Ifges(5,1) * t339 + t616;
t446 = t480 * t412;
t336 = -Ifges(4,5) * t416 + t446;
t533 = t415 * t336;
t334 = -Ifges(4,6) * t416 + t412 * t704;
t543 = t411 * t334;
t5 = (t533 / 0.2e1 - t543 / 0.2e1 - t606 - t611 - pkin(1) * mrSges(3,2) + (-Ifges(6,3) + Ifges(3,1) - Ifges(3,2) - Ifges(5,3) - Ifges(4,3) + (m(4) * pkin(6) + t484) * pkin(6)) * t412 + t461 + t707 * t416) * t416 + t176 * t307 + t170 * t308 + t175 * t309 + t169 * t310 + t281 * t124 + t282 * t123 + t52 * t210 + t61 * t207 + t53 * t208 + t60 * t209 + (-pkin(1) * mrSges(3,1) + Ifges(5,5) * t661 + Ifges(6,5) * t670 + Ifges(5,6) * t660 + Ifges(6,6) * t716 + pkin(6) * t357 + t335 * t648 + t337 * t643 - t707 * t412) * t412 + t118 * t716 + t117 * t672 + t478 * t657 + t473 * t658 + t238 * t660 + t239 * t661 + t119 * t669 + t120 * t670 + m(4) * (t320 * t323 + t321 * t324) + t324 * t373 + t321 * t374 + t323 * t375 + t320 * t376 + t377 * t252 + t378 * t251 + m(5) * (t169 * t175 + t170 * t176 + t377 * t378) + m(6) * (t281 * t282 + t52 * t60 + t53 * t61);
t573 = t5 * qJD(1);
t521 = pkin(3) * t539;
t299 = -pkin(4) * t339 + t521;
t381 = Ifges(4,2) * t415 + t618;
t420 = t746 + t53 * t590 - t169 * t579 - t52 * t589 + t245 * t117 / 0.2e1 + t739 * t670 + t490 * t642 + t737 * t716;
t436 = t170 * mrSges(5,3) - t617 - t605 / 0.2e1 - t377 * mrSges(5,1);
t441 = t377 * mrSges(5,2) - t610 / 0.2e1 + t616;
t468 = Ifges(4,5) * t411 + Ifges(4,6) * t415;
t6 = t420 + t185 * t307 + t184 * t309 + t299 * t123 + t67 * t207 + t66 * t209 + (-t340 * t715 + t436) * t339 + t441 * t340 + m(5) * (t169 * t184 + t170 * t185) + m(6) * (t281 * t299 + t52 * t66 + t53 * t67) + (-t415 * t334 / 0.2e1 + t336 * t648 + t416 * t468 / 0.2e1 + (-pkin(6) * t702 + t381 * t647 + t742) * t412 + (m(5) * t377 + t251) * t638 + (t320 * t411 - t321 * t415) * mrSges(4,3)) * t412 + t320 * t373 - t321 * t375;
t572 = t6 * qJD(1);
t9 = t420 + t169 * t307 - t170 * t309 + t63 * t207 + t62 * t209 + m(6) * (t52 * t62 + t53 * t63) + ((-m(6) * t281 - t123) * pkin(4) + t436) * t339 + (-t339 * t715 + t441) * t340;
t571 = t9 * qJD(1);
t12 = t52 * t207 - t53 * t209 + t746 + t528 * t642 - (-t53 * mrSges(6,3) + t739 / 0.2e1 - t117 / 0.2e1) * t245 + (-t52 * mrSges(6,3) + t119 / 0.2e1 + t721 / 0.2e1) * t493;
t570 = t12 * qJD(1);
t558 = t339 * t157;
t557 = t352 * t493;
t556 = t352 * t492;
t555 = t353 * t245;
t554 = t353 * t279;
t553 = t367 * t123;
t552 = t377 * t411;
t550 = t409 * t245;
t549 = t409 * t279;
t547 = t410 * t309;
t546 = t410 * t339;
t544 = t411 * t251;
t542 = t411 * t373;
t538 = t413 * t493;
t537 = t413 * t492;
t536 = t414 * t307;
t535 = t414 * t340;
t532 = t415 * t375;
t107 = -mrSges(6,1) * t353 - mrSges(6,2) * t352;
t522 = qJD(5) * t107;
t519 = pkin(7) * mrSges(4,3) / 0.2e1;
t516 = -Ifges(4,1) / 0.4e1 + Ifges(4,2) / 0.4e1;
t515 = Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t514 = Ifges(5,3) * t646;
t513 = t590 / 0.2e1;
t512 = -t589 / 0.2e1;
t511 = mrSges(6,3) * t668;
t510 = t708 / 0.2e1;
t509 = -t708 / 0.2e1;
t508 = -t586 / 0.2e1;
t505 = t580 / 0.2e1;
t504 = -t579 / 0.2e1;
t503 = t576 / 0.2e1;
t499 = t184 / 0.2e1 + t170 / 0.2e1;
t491 = -mrSges(4,3) * t632 / 0.2e1;
t487 = pkin(4) * t645 + t656;
t477 = -Ifges(5,1) * t367 + t615;
t472 = Ifges(5,2) * t366 - t614;
t469 = -t607 + t612;
t292 = -mrSges(5,1) * t366 - mrSges(5,2) * t367;
t426 = t611 / 0.2e1 + t606 / 0.2e1 + t175 * mrSges(5,1) / 0.2e1 - t176 * mrSges(5,2) / 0.2e1 + t486;
t417 = t694 + t426 + (t339 * t367 - t340 * t366) * Ifges(5,4) + t715 * (t339 * t651 + t367 * t660);
t338 = -t637 + t640;
t421 = (t281 * t338 + t299 * t322 + (t53 + t66) * t758 + (-t52 + t67) * t140) * t691 - t299 * t157 / 0.2e1 - t338 * t123 / 0.2e1 + t543 / 0.4e1 - t533 / 0.4e1 + t430;
t431 = t361 / 0.2e1 + t360 / 0.2e1 + t710 / 0.4e1 - t730 / 0.4e1;
t432 = t245 * t679 - t279 * t682 + t52 * t717 + t716 * t758;
t435 = (t170 + t184) * t705 + (-t169 + t185) * t306;
t442 = (t175 * t414 + t176 * t410) * t692;
t457 = t410 * t308 / 0.2e1 + t414 * t310 / 0.2e1;
t460 = t306 * t661 + t660 * t705;
t1 = (0.3e1 / 0.4e1 * t612 - 0.3e1 / 0.4e1 * t607 + t431) * t416 + t421 + (t532 / 0.2e1 + t542 / 0.2e1) * pkin(7) + (-t544 / 0.2e1 + t457) * pkin(3) + (t66 * t664 + t667 * t67 + t432) * mrSges(6,3) + pkin(3) * t442 + ((-pkin(6) * mrSges(4,1) / 0.2e1 - t479 / 0.4e1 + t405 / 0.4e1 + pkin(2) * t684 + (t519 - t516) * t411) * t411 + (pkin(6) * t684 + 0.3e1 / 0.4e1 * t618 + t381 / 0.4e1 + pkin(2) * t685 + (t519 + t516) * t415 + (-t292 / 0.2e1 + t400 * t693) * pkin(3)) * t415 + t515) * t412 + (pkin(3) * t552 + t435) * t693 + t417 + (-t366 * t700 - t367 * t499 + t460) * mrSges(5,3) + t695;
t456 = t381 * t648 - t742;
t13 = t480 * t647 - pkin(2) * t484 + t704 * t643 + t456 + (m(5) * t400 + t292) * t640 + t773 * t338 + t774;
t467 = -t1 * qJD(1) + t13 * qJD(2);
t16 = -t637 * t773 + t774;
t428 = t210 * t644 + (t409 * t61 + t413 * t60) * t690 + t409 * t677;
t437 = (t53 + t62) * t758 + (-t52 + t63) * t140;
t464 = -t281 * t367 - t322 * t339;
t4 = t431 * t416 + t437 * t691 + (t553 / 0.2e1 + t558 / 0.2e1 + t464 * t691 + t428) * pkin(4) + (t62 * t664 + t63 * t667 + t432) * mrSges(6,3) + t460 * mrSges(5,3) + t417 + t514 + t430;
t466 = -qJD(1) * t4 + qJD(2) * t16;
t23 = t745 + (t740 / 0.2e1 - t159 / 0.2e1) * t279 + (t161 / 0.2e1 + t720 / 0.2e1) * t492;
t418 = (t119 / 0.4e1 + t721 / 0.4e1) * t492 + (mrSges(6,3) * t680 + t161 / 0.4e1 + t720 / 0.4e1) * t493 - (mrSges(6,3) * t679 + t752) * t245 + t207 * t749 + t754 + t526 * t641 - t757;
t8 = t418 - t486;
t465 = qJD(1) * t8 + qJD(2) * t23;
t463 = (-t636 / 0.2e1 + t654) * mrSges(6,1);
t458 = t550 / 0.2e1 - t538 / 0.2e1;
t443 = (t409 * t67 + t413 * t66) * t690;
t105 = -m(6) * (t352 * t358 + t353 * t359) + t698;
t11 = -t700 * mrSges(5,2) + t499 * mrSges(5,1) + pkin(4) * t443 + (t547 / 0.2e1 - t536 / 0.2e1 + (t535 / 0.2e1 - t546 / 0.2e1) * mrSges(5,3)) * pkin(3) + (-t555 / 0.2e1 + t557 / 0.2e1 + t458 * pkin(4)) * mrSges(6,3) + t621 + t727;
t427 = ((t353 + t358) * t758 + (-t352 + t359) * t140) * t690 + t769;
t429 = -t554 / 0.2e1 - t556 / 0.2e1 + t492 * t652 + t279 * t653;
t20 = (t663 - t306 / 0.2e1) * mrSges(5,1) + ((t549 / 0.2e1 + t537 / 0.2e1) * pkin(4) + t429) * mrSges(6,3) + t427 + t724;
t440 = -qJD(1) * t11 + qJD(2) * t20 - qJD(3) * t105;
t425 = (-t245 * t654 + t493 * t656) * mrSges(6,3) + t207 * t655 + t209 * t654 + t497;
t15 = (t683 + t67 / 0.2e1) * mrSges(6,2) + (t682 - t66 / 0.2e1) * mrSges(6,1) + t425 - t497;
t439 = qJD(1) * t15 + qJD(3) * t107 + t771;
t423 = (t207 * t644 + t209 * t649 + (-t245 * t649 + t493 * t645) * mrSges(6,3)) * pkin(4) + t497;
t19 = (t683 + t63 / 0.2e1) * mrSges(6,2) + (t682 - t62 / 0.2e1) * mrSges(6,1) + t423 - t497;
t379 = (mrSges(6,1) * t409 + mrSges(6,2) * t413) * pkin(4);
t85 = -t354 / 0.2e1 + t463 + (t652 + t487) * mrSges(6,2);
t434 = -qJD(1) * t19 - qJD(3) * t85 + qJD(4) * t379 - t771;
t419 = -t694 + t426 + t52 * t509 + t53 * t511 + t758 * t512 + t140 * t513 - t339 * (Ifges(5,1) * t366 + t614) / 0.4e1 - t367 * (Ifges(5,1) * t340 + t617) / 0.4e1 + t367 * (t473 - t605) / 0.4e1 + t170 * t503 + t705 * t504 + t306 * t505 + t339 * t472 / 0.4e1 + (Ifges(5,2) * t367 + t477 + t615) * t340 / 0.4e1 + (Ifges(5,2) * t339 + t478 - t610 + t616) * t366 / 0.4e1 + t489 * t641;
t365 = t379 * qJD(5);
t84 = t354 / 0.2e1 - t578 / 0.2e1 + t487 * mrSges(6,2) + t463;
t18 = t423 + t622 + t697;
t17 = pkin(4) * t413 * t509 + mrSges(6,3) * t429 + t511 * t636 + t427 - t724 + t753;
t14 = t425 + t621 + t697;
t10 = t490 + t592 / 0.2e1 - t591 / 0.2e1 - t593 / 0.2e1 - t594 / 0.2e1 + t504 * t639 + t536 * t689 + t352 * t512 + t353 * t513 + (mrSges(6,3) * t458 + t443) * pkin(4) + t621 + (t410 * t505 - t547 / 0.2e1) * pkin(3) - t727;
t7 = t418 + t486;
t3 = t419 + (pkin(4) * t464 + t437) * t690 + t428 * pkin(4) - t170 * t576 / 0.2e1 + t63 * t510 + t62 * t508 + t514 - t430 - (t553 + t558) * pkin(4) / 0.2e1;
t2 = t695 - t421 + t419 + t66 * t508 + t67 * t510 + t484 * t407 / 0.2e1 + t702 * t633 / 0.2e1 - t381 * t539 / 0.2e1 + t292 * t521 / 0.2e1 + t515 * t412 + t544 * t689 + ((t400 * t539 + t552) * pkin(3) + t435) * t692 + t700 * t577 + t184 * t503 + (t446 / 0.4e1 + t491 * t415) * t415 + t462 * t416 + t411 ^ 2 * t491 - (t542 + t532) * pkin(7) / 0.2e1 + t469 * t641 + (t442 + t457) * pkin(3) + (t479 / 0.2e1 - t704 / 0.4e1) * t541;
t21 = [qJD(2) * t5 + qJD(3) * t6 + qJD(4) * t9 + qJD(5) * t12, t2 * qJD(3) + t3 * qJD(4) + t7 * qJD(5) + t573 + ((m(4) * t703 + t415 * t374 - t411 * t376) * pkin(7) + t703 * mrSges(4,3) + t322 * t124 + t306 * t308 + t282 * t157 + t140 * t208 + (Ifges(3,5) + (-mrSges(3,1) + t702) * pkin(6) + t456) * t416 + (-m(4) * t408 - t357) * pkin(2) + t61 * t708 + t118 * t717 + t705 * t310 + 0.2e1 * (t175 * t705 + t176 * t306 + t378 * t400) * t692 - t60 * t586 + t159 * t672 + t477 * t657 + t472 * t658 + t120 * t664 + t161 * t669 + t335 * t643 + t337 * t647 + t238 * t651 - t367 * t239 / 0.2e1 + t378 * t292 + t400 * t252 - Ifges(3,6) * t412 + (-Ifges(5,5) * t367 + Ifges(6,5) * t279 + Ifges(5,6) * t366 + Ifges(6,6) * t492 + t468) * t646 + t175 * t576 + t176 * t577 + mrSges(3,2) * t407 + t758 * t210 + 0.2e1 * (t140 * t61 + t282 * t322 + t60 * t758) * t690) * qJD(2), t572 + t2 * qJD(2) + (m(6) * (t352 * t66 + t353 * t67) - t624 + t625 + t592 - t591 - t320 * mrSges(4,2) - t321 * mrSges(4,1) - Ifges(4,5) * t541 - Ifges(4,6) * t539 + t490 + (t555 - t557) * mrSges(6,3)) * qJD(3) + t10 * qJD(4) + t14 * qJD(5) + (m(5) * (t184 * t414 + t185 * t410) + (-t535 + t546) * mrSges(5,3)) * t602, t571 + t3 * qJD(2) + t10 * qJD(3) + (t490 - t593 - t594 - t626 + t627) * qJD(4) + t18 * qJD(5) + (m(6) * (t409 * t63 + t413 * t62) + (-t538 + t550) * mrSges(6,3)) * t601, t570 + t7 * qJD(2) + t14 * qJD(3) + t18 * qJD(4) + (t528 - t630 - t631) * qJD(5); -qJD(3) * t1 - qJD(4) * t4 + qJD(5) * t8 - t573, qJD(3) * t13 + qJD(4) * t16 + qJD(5) * t23, (m(6) * (-t140 * t352 + t353 * t758) + t469 + t702 * pkin(7) + (-t554 - t556) * mrSges(6,3) + t772) * qJD(3) + t17 * qJD(4) + t776 + (m(5) * (-t306 * t414 + t410 * t705) + (-t366 * t414 + t367 * t410) * mrSges(5,3)) * t602 + t467, t17 * qJD(3) + t772 * qJD(4) + t776 + (m(6) * t765 + (-t537 - t549) * mrSges(6,3)) * t601 + t466, (t526 + t759) * qJD(5) + t465 + t777 * t775; qJD(2) * t1 - qJD(4) * t11 + qJD(5) * t15 - t572, qJD(4) * t20 - t467 + t770, -qJD(4) * t105 + t522, ((t358 * t413 + t359 * t409) * t686 - t698) * qJD(4) + t84 * qJD(5) + t440, qJD(4) * t84 + t439 + t522; qJD(2) * t4 + qJD(3) * t11 + qJD(5) * t19 - t571, -qJD(3) * t20 - t466 + t770, qJD(5) * t85 - t440, -t365, -t365 - t434; -qJD(2) * t8 - qJD(3) * t15 - qJD(4) * t19 - t570, -t764 * t777 - t465, -qJD(4) * t85 - t439, t434, 0;];
Cq = t21;
