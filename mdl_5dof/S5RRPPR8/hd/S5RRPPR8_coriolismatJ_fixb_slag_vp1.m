% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:09
% EndTime: 2019-12-31 19:38:36
% DurationCPUTime: 18.48s
% Computational Cost: add. (28456->658), mult. (44379->848), div. (0->0), fcn. (48494->8), ass. (0->391)
t446 = sin(qJ(2));
t448 = cos(qJ(2));
t571 = pkin(8) + qJ(5);
t537 = sin(t571);
t538 = cos(t571);
t361 = t446 * t538 - t448 * t537;
t455 = t446 * t537 + t448 * t538;
t449 = cos(qJ(1));
t320 = t455 * t449;
t307 = Icges(6,4) * t320;
t319 = t361 * t449;
t447 = sin(qJ(1));
t218 = Icges(6,2) * t319 - Icges(6,6) * t447 + t307;
t605 = -Icges(6,1) * t319 + t218 + t307;
t306 = Icges(6,4) * t319;
t220 = Icges(6,1) * t320 - Icges(6,5) * t447 + t306;
t725 = -Icges(6,2) * t320 + t220 + t306;
t100 = t361 * t605 + t455 * t725;
t689 = -t447 / 0.4e1;
t263 = -Icges(6,5) * t455 - Icges(6,6) * t361;
t345 = Icges(6,4) * t361;
t267 = -Icges(6,2) * t455 + t345;
t271 = Icges(6,1) * t455 + t345;
t597 = t267 + t271;
t344 = Icges(6,4) * t455;
t268 = Icges(6,2) * t361 + t344;
t270 = Icges(6,1) * t361 - t344;
t599 = t268 - t270;
t89 = t263 * t447 + t319 * t599 + t320 * t597;
t800 = (t89 + t100) * t689;
t317 = t361 * t447;
t318 = t455 * t447;
t87 = -t263 * t449 + t317 * t599 + t318 * t597;
t799 = -t87 / 0.2e1;
t798 = -t87 / 0.4e1;
t796 = t100 / 0.2e1 + t89 / 0.2e1;
t305 = Icges(6,4) * t318;
t640 = Icges(6,2) * t317;
t497 = t640 + t305;
t637 = Icges(6,6) * t449;
t465 = t637 + t497;
t304 = Icges(6,4) * t317;
t649 = Icges(6,1) * t318;
t501 = t304 + t649;
t642 = Icges(6,5) * t449;
t466 = t642 + t501;
t486 = -t317 * t218 - t318 * t220;
t493 = Icges(6,5) * t318 + Icges(6,6) * t317;
t463 = Icges(6,3) * t449 + t493;
t731 = t463 - t493;
t443 = t449 ^ 2;
t569 = 0.2e1 * t305;
t496 = t569 + 0.2e1 * t637;
t568 = 0.2e1 * t642;
t722 = -t317 * (t496 + t640) - t318 * (t568 + t649) - Icges(6,3) * t443;
t780 = t722 * t449;
t16 = -t780 + (t731 * t447 + (-t466 + t501) * t320 - (t465 - t497) * t319 + t486) * t447;
t462 = Icges(6,5) * t320 + Icges(6,6) * t319 - Icges(6,3) * t447;
t459 = t449 * t462;
t113 = t459 - t486;
t114 = t319 * t218 + t320 * t220 - t447 * t462;
t17 = (-t319 * t497 - t320 * t501 + t113 - 0.2e1 * t459 + t486) * t449 + (-t317 * t465 - t318 * t466 - t731 * t449 + t114 - t722) * t447;
t57 = t113 * t447 + t780;
t58 = t114 * t447 - (t319 * t465 + t320 * t466 - t447 * t463) * t449;
t684 = t449 / 0.4e1;
t685 = -t449 / 0.4e1;
t795 = 0.2e1 * t17 * t684 + 0.2e1 * t58 * t685 + 0.2e1 * (t16 + t57) * t689 + (t89 / 0.4e1 + t100 / 0.4e1) * t447;
t433 = Icges(4,5) * t446;
t650 = Icges(4,1) * t448;
t503 = t433 + t650;
t339 = Icges(4,4) * t447 + t449 * t503;
t646 = Icges(3,4) * t446;
t398 = Icges(3,1) * t448 - t646;
t341 = Icges(3,5) * t447 + t398 * t449;
t792 = -t339 - t341;
t273 = rSges(6,1) * t361 - rSges(6,2) * t455;
t633 = qJ(3) * t448;
t401 = pkin(2) * t446 - t633;
t664 = pkin(3) * t446;
t543 = t401 + t664;
t444 = sin(pkin(8));
t622 = t448 * t444;
t652 = cos(pkin(8));
t427 = pkin(4) * t652 + pkin(3);
t662 = -pkin(3) + t427;
t482 = -pkin(4) * t622 + t446 * t662 + t273 + t543;
t174 = t482 * t447;
t176 = t482 * t449;
t540 = t446 * t652;
t383 = t540 - t622;
t627 = t446 * t444;
t472 = t448 * t652 + t627;
t513 = rSges(5,1) * t383 - rSges(5,2) * t472 + t543;
t212 = t513 * t447;
t214 = t513 * t449;
t660 = rSges(4,1) * t446;
t402 = -rSges(4,3) * t448 + t660;
t580 = t401 + t402;
t291 = t580 * t447;
t293 = t580 * t449;
t440 = t449 * rSges(4,2);
t441 = t449 * pkin(6);
t653 = rSges(4,3) + qJ(3);
t683 = rSges(4,1) + pkin(2);
t723 = t446 * t653 + t448 * t683 + pkin(1);
t260 = -t723 * t447 + t440 + t441;
t261 = (rSges(4,2) + pkin(6)) * t447 + t723 * t449;
t621 = t448 * t449;
t623 = t447 * t448;
t601 = t260 * t621 + t261 * t623;
t634 = qJ(3) * t446;
t713 = pkin(2) + pkin(3);
t724 = t448 * t713 + pkin(1) + t634;
t557 = t447 * t622;
t346 = -t447 * t540 + t557;
t347 = t472 * t447;
t729 = -t347 * rSges(5,1) + t346 * rSges(5,2);
t188 = t441 + (-rSges(5,3) - qJ(4)) * t449 - t724 * t447 + t729;
t428 = t447 * qJ(4);
t556 = t444 * t621;
t348 = -t449 * t540 + t556;
t349 = t472 * t449;
t479 = t349 * rSges(5,1) - t348 * rSges(5,2) - t447 * rSges(5,3);
t189 = t447 * pkin(6) + t724 * t449 - t428 + t479;
t611 = t188 * t621 + t189 * t623;
t405 = pkin(2) * t448 + t634;
t445 = -pkin(7) - qJ(4);
t564 = pkin(4) * t627;
t581 = t427 * t623 + t447 * t564;
t730 = -t318 * rSges(6,1) - t317 * rSges(6,2);
t167 = t441 + (-rSges(6,3) + t445) * t449 + (-pkin(1) - t405) * t447 - t581 + t730;
t223 = t320 * rSges(6,1) + t319 * rSges(6,2) - t447 * rSges(6,3);
t541 = pkin(4) * t444 + qJ(3);
t663 = pkin(2) + t427;
t168 = (pkin(6) + t445) * t447 + (t446 * t541 + t448 * t663 + pkin(1)) * t449 + t223;
t614 = t167 * t621 + t168 * t623;
t624 = t446 * t449;
t626 = t446 * t447;
t714 = m(6) / 0.2e1;
t715 = m(5) / 0.2e1;
t716 = m(4) / 0.2e1;
t561 = (-t291 * t624 + t293 * t626 + t601) * t716 + (-t174 * t624 + t176 * t626 + t614) * t714 + (-t212 * t624 + t214 * t626 + t611) * t715;
t424 = pkin(2) * t626;
t508 = rSges(5,1) * t346 + rSges(5,2) * t347;
t207 = t424 + (-t633 + t664) * t447 - t508;
t416 = qJ(3) * t621;
t590 = t348 * rSges(5,1) + t349 * rSges(5,2);
t208 = -t624 * t713 + t416 + t590;
t287 = t424 + (-t448 * t653 + t660) * t447;
t423 = rSges(4,3) * t621;
t288 = -t624 * t683 + t416 + t423;
t231 = -rSges(6,1) * t317 + rSges(6,2) * t318;
t182 = t424 + (t427 * t446 - t448 * t541) * t447 - t231;
t234 = t319 * rSges(6,1) - t320 * rSges(6,2);
t727 = pkin(4) * t556 - t234;
t183 = -t624 * t663 + t416 + t727;
t487 = t182 * t449 + t183 * t447;
t562 = ((t287 * t449 + t288 * t447) * t446 + t601) * t716 + (t446 * t487 + t614) * t714 + ((t207 * t449 + t208 * t447) * t446 + t611) * t715;
t756 = t561 - t562;
t791 = t756 * qJD(1);
t618 = Icges(6,1) - Icges(6,2);
t536 = t618 * t318;
t470 = t536 + t642;
t749 = -0.2e1 * t304;
t778 = t470 - t749;
t790 = t778 * t319 + t447 * (-Icges(6,5) * t317 + Icges(6,6) * t318);
t789 = -t271 / 0.2e1;
t221 = rSges(6,3) * t449 - t730;
t158 = -t447 * t221 - t449 * t223;
t775 = -t447 * t231 + t234 * t449;
t788 = t158 * t775;
t615 = (-t447 * t182 + t183 * t449) * t714 + (-t447 * t207 + t208 * t449) * t715;
t123 = t447 * t174 + t176 * t449;
t616 = t123 * t714 + (t447 * t212 + t214 * t449) * t715;
t25 = t616 - t615;
t787 = t25 * qJD(1);
t225 = -Icges(6,5) * t319 + Icges(6,6) * t320;
t786 = (-t225 * t447 - t319 * t725 + t320 * t605) * t447;
t785 = (t225 * t449 - t317 * t725 + t318 * t605) * t447;
t784 = t775 * t448;
t783 = (t270 / 0.2e1 - t268 / 0.2e1) * t455;
t782 = (-Icges(3,6) + Icges(4,6)) * t448 + (-Icges(4,4) - Icges(3,5)) * t446;
t393 = Icges(3,2) * t448 + t646;
t636 = Icges(4,3) * t448;
t494 = t636 - t433;
t781 = -t792 + (-t393 - t494) * t449;
t272 = -rSges(6,1) * t455 - rSges(6,2) * t361;
t442 = t447 ^ 2;
t576 = t442 + t443;
t735 = t576 * t446;
t667 = m(6) * (t272 * t735 + t784);
t419 = Icges(4,5) * t621;
t331 = Icges(4,6) * t447 + Icges(4,3) * t624 + t419;
t391 = Icges(3,5) * t448 - Icges(3,6) * t446;
t333 = Icges(3,3) * t447 + t391 * t449;
t392 = Icges(4,4) * t448 + Icges(4,6) * t446;
t335 = Icges(4,2) * t447 + t392 * t449;
t779 = -t331 * t624 + t792 * t621 + (-t333 - t335) * t447;
t332 = Icges(3,5) * t623 - Icges(3,6) * t626 - Icges(3,3) * t449;
t420 = Icges(3,4) * t626;
t340 = Icges(3,1) * t623 - Icges(3,5) * t449 - t420;
t324 = Icges(5,4) * t347;
t242 = Icges(5,2) * t346 - Icges(5,6) * t449 - t324;
t323 = Icges(5,4) * t346;
t244 = Icges(5,1) * t347 + Icges(5,5) * t449 - t323;
t609 = -t242 * t348 - t349 * t244;
t777 = -t447 * t332 - t340 * t621 + t609;
t772 = t346 * t242 + t347 * t244;
t694 = t267 / 0.2e1;
t771 = t599 * t455 / 0.2e1 + (-t694 + t789) * t361;
t747 = -t383 / 0.2e1;
t688 = t447 / 0.2e1;
t745 = -t449 / 0.2e1;
t489 = t167 * t449 + t168 * t447;
t765 = t272 * t489;
t526 = t576 * t273;
t436 = Icges(3,4) * t448;
t641 = Icges(3,2) * t446;
t337 = Icges(3,6) * t447 + (t436 - t641) * t449;
t762 = -t337 * t624 - t779;
t761 = t782 * t447;
t629 = (-Icges(4,2) * t449 + t447 * t392) * t449;
t759 = -t629 + t779;
t758 = t781 * t447;
t535 = t618 * t317;
t757 = t535 - t637;
t238 = Icges(5,5) * t347 - Icges(5,6) * t346 + Icges(5,3) * t449;
t336 = Icges(3,4) * t623 - Icges(3,2) * t626 - Icges(3,6) * t449;
t297 = t341 * t623;
t521 = t333 * t449 - t297;
t755 = -t238 * t447 - t336 * t624 - t337 * t626 - t521 - t777;
t643 = Icges(4,5) * t448;
t390 = Icges(4,3) * t446 + t643;
t395 = Icges(4,1) * t446 - t643;
t549 = Icges(5,2) * t383 / 0.2e1 + Icges(5,4) * t472 + Icges(5,1) * t747;
t651 = Icges(3,1) * t446;
t754 = t549 * t472 - (t436 + t651 / 0.2e1 - t641 / 0.2e1 + t395 / 0.2e1 - t390 / 0.2e1) * t448 - (t398 / 0.2e1 - t393 / 0.2e1 + t433 + t650 / 0.2e1 - t636 / 0.2e1) * t446;
t504 = -t436 - t651;
t531 = (t504 * t449 - t337) * t447;
t532 = (-t504 * t447 + t336) * t449;
t533 = (-Icges(4,1) * t624 + t331 + t419) * t447;
t330 = -Icges(4,6) * t449 + t390 * t447;
t534 = (-t395 * t447 + t330) * t449;
t750 = (-t534 + t533 + t532 + t531) * t448;
t746 = -t447 / 0.2e1;
t338 = -Icges(4,4) * t449 + t447 * t503;
t741 = (t330 * t446 + t338 * t448) * t447;
t566 = m(6) / 0.4e1 + m(5) / 0.4e1;
t625 = t446 * t448;
t578 = t576 * t625;
t737 = (m(4) / 0.4e1 + t566) * (t578 - t625);
t124 = t238 * t449 + t772;
t240 = Icges(5,5) * t349 - Icges(5,6) * t348 - Icges(5,3) * t447;
t325 = Icges(5,4) * t349;
t243 = -Icges(5,2) * t348 - Icges(5,6) * t447 + t325;
t645 = Icges(5,4) * t348;
t246 = Icges(5,1) * t349 - Icges(5,5) * t447 - t645;
t608 = -t348 * t243 + t349 * t246;
t127 = -t447 * t240 + t608;
t733 = t124 + t127;
t728 = t782 * t449;
t726 = qJD(2) * t449;
t522 = -0.2e1 * t305;
t456 = t522 + t757;
t477 = (-t786 - (t456 * t320 + t790) * t449) * t688 + (-t785 - ((t522 - 0.2e1 * t637) * t318 - (-0.2e1 * t536 + t749 - 0.2e1 * t642) * t317) * t449) * t745;
t458 = t569 - t757;
t478 = (t786 - (t458 * t320 - t790) * t449) * t688 + (t785 - (-(-t749 + t568) * t317 + (t496 - 0.2e1 * t535) * t318) * t449) * t745;
t385 = pkin(3) * t623 + qJ(4) * t449;
t721 = 0.2e1 * t735;
t719 = 0.4e1 * qJD(1);
t718 = 0.2e1 * qJD(2);
t582 = t576 * t405;
t518 = t447 * t385 + t449 * (pkin(3) * t621 - t428) + t582;
t134 = t447 * (rSges(5,3) * t449 - t729) + t449 * t479 + t518;
t607 = -t212 * t623 - t214 * t621;
t710 = m(5) * (t134 * t735 + t607);
t707 = m(5) * (t188 * t207 + t189 * t208);
t704 = m(6) * (-t174 * t234 - t176 * t231 - t765);
t703 = m(6) * (t273 * t487 - t765);
t322 = t448 * t662 + t564;
t102 = (t221 - t385 + t581) * t447 + (t322 * t449 + t223 + t428) * t449 + t518;
t612 = -t174 * t623 - t176 * t621;
t702 = m(6) * (t102 * t735 + t612);
t524 = t576 * t272;
t600 = t448 * t526;
t699 = m(6) * (-t784 + (t158 - t524) * t446 + t600);
t697 = m(6) * (t167 * t182 + t168 * t183);
t696 = m(6) * (t167 * t231 + t168 * t234);
t695 = t244 / 0.2e1;
t661 = rSges(3,1) * t448;
t547 = pkin(1) + t661;
t577 = rSges(3,2) * t626 + t449 * rSges(3,3);
t289 = -t447 * t547 + t441 + t577;
t422 = rSges(3,2) * t624;
t290 = -t422 + t547 * t449 + (rSges(3,3) + pkin(6)) * t447;
t403 = rSges(3,1) * t446 + rSges(3,2) * t448;
t379 = t403 * t447;
t381 = t403 * t449;
t682 = m(3) * (t289 * t379 - t290 * t381);
t406 = rSges(4,1) * t448 + rSges(4,3) * t446;
t185 = t447 * (t406 * t447 - t440) + (t447 * rSges(4,2) + t406 * t449) * t449 + t582;
t595 = -t291 * t623 - t293 * t621;
t678 = m(4) * (t185 * t735 + t595);
t676 = m(4) * (t260 * t287 + t261 * t288);
t675 = m(4) * (-t260 * t626 + t261 * t624);
t674 = m(5) * (-t188 * t626 + t189 * t624);
t673 = m(5) * (-t188 * t449 - t189 * t447);
t670 = m(6) * (-t167 * t626 + t168 * t624);
t669 = m(6) * (t158 * t735 + t600);
t668 = m(6) * t489;
t628 = t336 * t446;
t587 = -t494 * t447 + t338;
t585 = -Icges(3,2) * t623 + t340 - t420;
t583 = t447 * (qJ(3) * t623 - t424) + t449 * (-pkin(2) * t624 + t416);
t579 = -t405 - t406;
t471 = (t231 * t449 + t234 * t447) * t446;
t104 = -t471 * m(6) / 0.2e1;
t574 = t104 * qJD(3);
t468 = t775 * t714;
t516 = m(6) * t526;
t117 = t468 + t516 / 0.2e1;
t573 = t117 * qJD(1);
t567 = t714 + t715;
t199 = t567 * t721;
t572 = t199 * qJD(1);
t559 = t16 / 0.2e1 + t57 / 0.2e1;
t558 = t58 / 0.2e1 - t17 / 0.2e1;
t366 = Icges(5,4) * t383;
t281 = -Icges(5,2) * t472 + t366;
t284 = Icges(5,1) * t472 + t366;
t550 = t281 / 0.2e1 + t284 / 0.2e1;
t542 = -pkin(3) * t448 - t405;
t539 = -t478 - t477;
t530 = t587 * t449;
t528 = t585 * t449;
t520 = t337 * t446 - t332;
t97 = t456 * t361 - (t470 + 0.2e1 * t304) * t455;
t515 = t704 / 0.2e1 + t800 + (-t87 + t97) * t685;
t99 = t458 * t361 + t455 * t778;
t514 = t703 / 0.2e1 + t800 + (t87 + t99) * t684;
t512 = -rSges(5,1) * t472 - rSges(5,2) * t383 + t542;
t511 = -t331 * t626 + t335 * t449 - t339 * t623;
t510 = t392 / 0.2e1 + t391 / 0.2e1 - Icges(5,5) * t472 / 0.2e1 + Icges(5,6) * t747;
t481 = t272 - t322 + t542;
t175 = t481 * t447;
t177 = t481 * t449;
t488 = t175 * t447 + t177 * t449;
t485 = t449 * (Icges(5,5) * t346 + Icges(5,6) * t347) - t447 * (Icges(5,5) * t348 + Icges(5,6) * t349);
t258 = Icges(5,1) * t346 + t324;
t259 = Icges(5,1) * t348 + t325;
t461 = (t243 + t259) * t447 - (-t242 + t258) * t449;
t256 = Icges(5,2) * t347 + t323;
t257 = Icges(5,2) * t349 + t645;
t460 = (t246 - t257) * t447 - (t244 - t256) * t449;
t457 = t447 * t559 + t449 * t558;
t190 = -t629 + t741;
t454 = t558 + (t190 + t733 - t741 - t759 - t772) * t746 + (t127 + t762) * t688 + (-t297 + (t333 + t628) * t449 + t755 + t777) * t745;
t453 = t511 * t746 + t559 + (-t447 * (-t340 * t448 + t628) + t124 + t190 - t332 * t449) * t745 + (t449 * t520 - t608 + t733 + t759 + t762) * t449 / 0.2e1 + (t449 * t240 + t330 * t624 + t338 * t621 + t511 + t521 + (t238 + t520) * t447 + t755 + t609) * t688;
t407 = -rSges(3,2) * t446 + t661;
t294 = t579 * t449;
t292 = t579 * t447;
t215 = t512 * t449;
t213 = t512 * t447;
t200 = t566 * t721 - (m(5) + m(6)) * t735 / 0.2e1;
t198 = t449 * (-rSges(4,1) * t624 + t423) - t402 * t442 + t583;
t184 = 0.4e1 * t737;
t149 = -pkin(3) * t735 + t447 * t508 + t449 * t590 + t583;
t116 = t468 - t516 / 0.2e1;
t115 = -t667 / 0.2e1;
t107 = t727 * t449 + (pkin(4) * t557 + t231) * t447 - t427 * t735 + t583;
t106 = t669 / 0.2e1;
t105 = t471 * t714;
t68 = -t273 * t524 + t788;
t67 = -t668 + t673;
t62 = t699 / 0.2e1;
t45 = t670 + t674 + t675;
t44 = t696 + (t789 - t267 / 0.2e1) * t361 - t783;
t29 = t106 + t62 + t667 / 0.2e1;
t28 = t115 + t106 - t699 / 0.2e1;
t27 = t115 + t62 - t669 / 0.2e1;
t24 = t615 + t616;
t22 = t678 + t702 + t710;
t11 = t550 * t383 + (t271 / 0.2e1 + t694) * t361 + t783 + t682 + t676 + t707 + t697 - t754;
t10 = m(6) * t68 + t478;
t9 = m(6) * (t102 * t775 + t123 * t272) + t477;
t7 = t561 + t562;
t4 = t447 * t453 + t449 * t454;
t3 = -t704 / 0.2e1 + (t798 + t97 / 0.4e1) * t449 + t514 + t795;
t2 = t457 + t514 + t515;
t1 = -t703 / 0.2e1 + (t798 - t99 / 0.4e1) * t449 + t515 + t795;
t5 = [t11 * qJD(2) + t45 * qJD(3) + t67 * qJD(4) + t44 * qJD(5), t11 * qJD(1) + t7 * qJD(3) + t24 * qJD(4) + t2 * qJD(5) + ((t188 * t215 + t189 * t213 - t207 * t214 - t208 * t212) * t715 + (t167 * t177 + t168 * t175 - t174 * t183 - t176 * t182) * t714 + (t260 * t294 + t261 * t292 - t287 * t293 - t288 * t291) * t716) * t718 + (t799 + m(3) * (-t289 * t407 - t379 * t403) + (t242 / 0.2e1 - t258 / 0.2e1) * t383 - (t695 - t256 / 0.2e1) * t472 - t550 * t347 + t549 * t346 + t510 * t449 - t454 - t99 / 0.2e1) * t726 + ((t510 * t447 + m(3) * (-t290 * t407 + t381 * t403) - (-t246 / 0.2e1 + t257 / 0.2e1) * t472 - t549 * t348 + (t243 / 0.2e1 + t259 / 0.2e1) * t383 - t453 + t550 * t349 + t796) * t447 + (-t534 / 0.2e1 + t532 / 0.2e1 + t533 / 0.2e1 + t531 / 0.2e1) * t446 + (-t530 / 0.2e1 - t528 / 0.2e1 + t781 * t688) * t448) * qJD(2), qJD(1) * t45 + qJD(2) * t7 + qJD(4) * t200 - qJD(5) * t104, qJD(1) * t67 + qJD(2) * t24 + qJD(3) * t200 + qJD(5) * t116, t44 * qJD(1) + t2 * qJD(2) + t116 * qJD(4) - t574 + ((m(6) * (t167 * t272 + t231 * t273) + t799 + t97 / 0.2e1 - t558) * t449 + (m(6) * (t168 * t272 + t234 * t273) - t559 + t796) * t447) * qJD(5); t4 * qJD(2) + t756 * qJD(3) + t25 * qJD(4) + t1 * qJD(5) + (-t682 / 0.4e1 - t676 / 0.4e1 - t697 / 0.4e1 - t707 / 0.4e1) * t719 + ((t695 - t244 / 0.2e1) * t383 * t447 + (t281 + t284) * t747 + t754 + t771) * qJD(1), t4 * qJD(1) + t22 * qJD(3) + t9 * qJD(5) - (t346 * t460 + t347 * t461 - t449 * t485 + (-t728 * t449 + t750 + ((t585 + t587) * t449 - t758) * t446) * t447 + t761 * t443) * t726 / 0.2e1 + (m(6) * (t102 * t107 - t174 * t175 - t176 * t177) + m(5) * (t134 * t149 - t212 * t213 - t214 * t215) + m(4) * (t185 * t198 - t291 * t292 - t293 * t294) + m(3) * ((t447 * (rSges(3,1) * t623 - t577) + t449 * (rSges(3,1) * t621 + t447 * rSges(3,3) - t422)) * (-t447 * t379 - t381 * t449) + t576 * t407 * t403) + t478 + (t348 * t460 + t349 * t461 + t728 * t442 + t447 * t485 + (t750 - t761 * t447 + (t530 + t528 - t758) * t446) * t449) * t688) * qJD(2), t791 + t22 * qJD(2) + t28 * qJD(5) + (-0.4e1 * t737 + 0.2e1 * (t716 + t567) * (-t448 * t735 + t578)) * qJD(3), t787, t1 * qJD(1) + t9 * qJD(2) + t28 * qJD(3) + ((-t68 + (-t102 + t158) * t775 + (-t526 - t123) * t272) * m(6) - t477 + t539) * qJD(5); -t756 * qJD(2) - t199 * qJD(4) + t105 * qJD(5) + (-t670 / 0.4e1 - t675 / 0.4e1 - t674 / 0.4e1) * t719, -t791 + t184 * qJD(3) + t27 * qJD(5) + 0.4e1 * (-t702 / 0.4e1 - t710 / 0.4e1 - t678 / 0.4e1) * qJD(2) + ((-t448 * t107 + t612) * t714 + (-t448 * t149 + t607) * t715 + (-t448 * t198 + t595) * t716 + ((t102 + t488) * t714 + (t213 * t447 + t215 * t449 + t134) * t715 + (t292 * t447 + t294 * t449 + t185) * t716) * t446) * t718, t184 * qJD(2), -t572, t105 * qJD(1) + t27 * qJD(2) + qJD(5) * t667; -t25 * qJD(2) + t199 * qJD(3) + t117 * qJD(5) + (t668 / 0.4e1 - t673 / 0.4e1) * t719, -t787 + ((t175 * t449 - t447 * t177) * t714 + (t213 * t449 - t447 * t215) * t715) * t718, t572, 0, t573; (-t696 - t771) * qJD(1) + t3 * qJD(2) + t574 - t117 * qJD(4) + t457 * qJD(5), t3 * qJD(1) + ((t158 * t107 + t273 * t488) * m(6) - t478 + t539) * qJD(2) + t29 * qJD(3) + t10 * qJD(5), qJD(1) * t104 + qJD(2) * t29, -t573, t457 * qJD(1) + t10 * qJD(2) + (m(6) * (t272 * t526 - t788) + t477) * qJD(5);];
Cq = t5;
