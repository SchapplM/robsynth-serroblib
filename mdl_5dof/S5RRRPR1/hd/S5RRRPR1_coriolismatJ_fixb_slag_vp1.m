% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:49
% EndTime: 2019-12-05 18:38:19
% DurationCPUTime: 19.29s
% Computational Cost: add. (64262->657), mult. (49461->848), div. (0->0), fcn. (46270->10), ass. (0->403)
t550 = qJ(2) + qJ(3);
t525 = pkin(9) + t550;
t519 = sin(t525);
t520 = cos(t525);
t526 = sin(t550);
t527 = cos(t550);
t893 = Icges(4,5) * t527 + Icges(5,5) * t520 - Icges(4,6) * t526 - Icges(5,6) * t519;
t892 = Icges(4,3) + Icges(5,3);
t554 = cos(qJ(1));
t891 = t893 * t554;
t552 = sin(qJ(1));
t711 = t527 * t552;
t715 = t526 * t552;
t719 = t520 * t552;
t723 = t519 * t552;
t885 = Icges(4,5) * t711 + Icges(5,5) * t719 - Icges(4,6) * t715 - Icges(5,6) * t723 - t892 * t554;
t890 = t892 * t552 + t891;
t756 = Icges(5,4) * t519;
t456 = Icges(5,1) * t520 - t756;
t364 = Icges(5,5) * t552 + t456 * t554;
t757 = Icges(4,4) * t526;
t473 = Icges(4,1) * t527 - t757;
t405 = Icges(4,5) * t552 + t473 * t554;
t889 = -t364 * t719 - t405 * t711;
t714 = t526 * t554;
t624 = -rSges(4,2) * t714 + rSges(4,3) * t552;
t642 = rSges(4,2) * t715 + t554 * rSges(4,3);
t710 = t527 * t554;
t262 = t552 * (rSges(4,1) * t711 - t642) + t554 * (rSges(4,1) * t710 + t624);
t555 = -pkin(7) - pkin(6);
t524 = t552 * t555;
t546 = t554 * pkin(6);
t689 = t554 * t555;
t553 = cos(qJ(2));
t545 = t553 * pkin(2);
t523 = t545 + pkin(1);
t774 = pkin(1) - t523;
t663 = -t552 * (t552 * t774 - t546 - t689) + t554 * (-pkin(6) * t552 - t554 * t774 - t524);
t181 = t262 + t663;
t474 = rSges(4,1) * t526 + rSges(4,2) * t527;
t449 = t474 * t552;
t450 = t474 * t554;
t279 = -t552 * t449 - t554 * t450;
t771 = rSges(4,1) * t527;
t475 = -rSges(4,2) * t526 + t771;
t551 = sin(qJ(2));
t777 = pkin(2) * t551;
t569 = t474 + t777;
t854 = t569 * t554;
t855 = t569 * t552;
t114 = t181 * t279 + (t552 * t855 + t554 * t854) * t475;
t776 = pkin(3) * t526;
t598 = pkin(4) * t519 + t776;
t451 = t554 * t598;
t548 = t552 ^ 2;
t522 = qJ(5) + t525;
t511 = cos(t522);
t727 = t511 * t552;
t510 = sin(t522);
t731 = t510 * t552;
t393 = rSges(6,1) * t731 + rSges(6,2) * t727;
t446 = rSges(6,1) * t510 + rSges(6,2) * t511;
t396 = t446 * t554;
t850 = t552 * t393 + t554 * t396;
t159 = -t554 * t451 - t548 * t598 - t850;
t448 = -t598 - t777;
t567 = t446 - t448;
t267 = t567 * t552;
t269 = t567 * t554;
t769 = rSges(6,1) * t511;
t447 = -rSges(6,2) * t510 + t769;
t521 = pkin(3) * t527;
t775 = pkin(4) * t520;
t576 = -t447 - t775 - t521;
t281 = t576 * t552;
t283 = t576 * t554;
t681 = -t267 * t281 - t269 * t283;
t348 = rSges(6,1) * t727 - rSges(6,2) * t731 - t554 * rSges(6,3);
t730 = t510 * t554;
t622 = -rSges(6,2) * t730 + rSges(6,3) * t552;
t726 = t511 * t554;
t244 = t552 * t348 + t554 * (rSges(6,1) * t726 + t622);
t476 = t521 + t523;
t427 = t476 + t775;
t547 = -qJ(4) + t555;
t631 = -pkin(8) + t547;
t508 = t552 * t631;
t517 = t552 * t547;
t650 = -t552 * t427 - t554 * t631;
t645 = t476 - t523;
t675 = -t552 * ((-t547 + t555) * t554 - t645 * t552) + t554 * (t554 * t645 - t517 + t524);
t690 = t554 * t547;
t107 = -t552 * (t552 * t476 + t650 + t690) + t244 + t554 * (-t508 + t517 + (t427 - t476) * t554) + t675;
t89 = t107 + t663;
t46 = t89 * t159 + t681;
t722 = t519 * t554;
t623 = -rSges(5,2) * t722 + rSges(5,3) * t552;
t644 = rSges(5,2) * t723 + t554 * rSges(5,3);
t718 = t520 * t554;
t157 = t552 * (rSges(5,1) * t719 - t644) + t554 * (rSges(5,1) * t718 + t623) + t675;
t132 = t157 + t663;
t549 = t554 ^ 2;
t640 = t548 + t549;
t457 = rSges(5,1) * t519 + rSges(5,2) * t520;
t643 = rSges(5,1) * t723 + rSges(5,2) * t719;
t655 = -t549 * t457 - t552 * t643;
t236 = -t640 * t776 + t655;
t599 = t776 + t777;
t579 = t457 + t599;
t301 = t579 * t552;
t303 = t579 * t554;
t770 = rSges(5,1) * t520;
t627 = rSges(5,2) * t519 - t521 - t770;
t331 = t627 * t552;
t333 = t627 * t554;
t678 = -t301 * t331 - t303 * t333;
t75 = t132 * t236 + t678;
t888 = m(4) * t114 + m(5) * t75 + m(6) * t46;
t628 = t457 + t776;
t330 = t628 * t552;
t834 = m(5) / 0.2e1;
t577 = t446 + t598;
t280 = t577 * t552;
t282 = t577 * t554;
t851 = t552 * t280 + t282 * t554;
t597 = t628 * t554;
t856 = t554 * t597;
t871 = -m(6) / 0.2e1;
t682 = t851 * t871 + (-t552 * t330 - t856) * t834;
t284 = t552 * t598 + t393;
t285 = -t451 - t396;
t320 = pkin(3) * t715 + t643;
t833 = m(6) / 0.2e1;
t683 = (t284 * t552 - t285 * t554) * t833 + (t320 * t552 + t856) * t834;
t65 = t683 - t682;
t887 = qJD(1) * t65;
t273 = -t448 * t552 + t393;
t416 = t554 * t448;
t274 = t416 - t396;
t305 = t552 * t599 + t643;
t468 = t554 * t599;
t306 = -t457 * t554 - t468;
t686 = (t273 * t552 - t274 * t554) * t833 + (t305 * t552 - t306 * t554) * t834;
t852 = t552 * t267 + t269 * t554;
t687 = t852 * t871 + (-t552 * t301 - t303 * t554) * t834;
t55 = t687 - t686;
t886 = t55 * qJD(1);
t884 = Icges(4,5) * t526 + Icges(5,5) * t519 + Icges(4,6) * t527 + Icges(5,6) * t520;
t883 = -t554 * t890 - t889;
t873 = t364 * t718 + t405 * t710 + t552 * t890;
t487 = Icges(5,4) * t723;
t363 = Icges(5,1) * t719 - Icges(5,5) * t554 - t487;
t504 = Icges(4,4) * t715;
t404 = Icges(4,1) * t711 - Icges(4,5) * t554 - t504;
t882 = -t363 * t718 - t404 * t710 - t552 * t885;
t881 = t446 * t640;
t507 = Icges(5,4) * t520;
t454 = -Icges(5,2) * t519 + t507;
t362 = Icges(5,6) * t552 + t454 * t554;
t518 = Icges(4,4) * t527;
t471 = -Icges(4,2) * t526 + t518;
t403 = Icges(4,6) * t552 + t471 * t554;
t880 = -t362 * t723 - t403 * t715 + t883;
t361 = Icges(5,4) * t719 - Icges(5,2) * t723 - Icges(5,6) * t554;
t402 = Icges(4,4) * t711 - Icges(4,2) * t715 - Icges(4,6) * t554;
t879 = t361 * t722 + t402 * t714 + t882;
t878 = -t362 * t722 - t403 * t714 + t873;
t501 = Icges(6,4) * t511;
t436 = -Icges(6,2) * t510 + t501;
t437 = Icges(6,1) * t510 + t501;
t877 = t436 + t437;
t455 = Icges(5,1) * t519 + t507;
t876 = t455 + t454;
t472 = Icges(4,1) * t526 + t518;
t875 = t472 + t471;
t874 = t362 * t519 + t403 * t526 - t885;
t872 = t361 * t519 + t402 * t526;
t835 = m(4) / 0.2e1;
t809 = t552 / 0.2e1;
t807 = -t554 / 0.2e1;
t867 = t554 / 0.2e1;
t758 = Icges(3,4) * t551;
t493 = Icges(3,2) * t553 + t758;
t496 = Icges(3,1) * t553 - t758;
t866 = (t496 / 0.2e1 - t493 / 0.2e1) * t551;
t865 = t884 * t552;
t864 = t884 * t554;
t453 = Icges(5,2) * t520 + t756;
t470 = Icges(4,2) * t527 + t757;
t863 = (-t470 + t473) * t527 - t875 * t526 + (-t453 + t456) * t520 - t876 * t519;
t626 = t523 + t771;
t298 = -t552 * t626 + t642 - t689;
t299 = t554 * t626 - t524 + t624;
t571 = (-t298 * t554 - t299 * t552) * t475;
t625 = t476 + t770;
t271 = -t552 * t625 + t644 - t690;
t272 = t554 * t625 - t517 + t623;
t679 = t333 * t271 + t331 * t272;
t251 = -t348 + t650;
t252 = -t508 + (t427 + t769) * t554 + t622;
t685 = t283 * t251 + t281 * t252;
t634 = (t571 + (t552 * t854 - t554 * t855) * t474) * t835 + (-t273 * t282 - t274 * t280 + t685) * t833 + (-t305 * t597 - t306 * t330 + t679) * t834;
t635 = (-t449 * t854 + t450 * t855 + t571) * t835 + (-t267 * t285 - t269 * t284 + t685) * t833 + (t301 * t597 - t303 * t320 + t679) * t834;
t862 = t634 - t635;
t810 = -t552 / 0.2e1;
t861 = t809 + t810;
t651 = -t470 * t554 + t405;
t652 = -Icges(4,2) * t711 + t404 - t504;
t653 = -t472 * t554 - t403;
t654 = t472 * t552 + t402;
t659 = -t453 * t554 + t364;
t660 = -Icges(5,2) * t719 + t363 - t487;
t661 = -t455 * t554 - t362;
t662 = t455 * t552 + t361;
t860 = (-t659 * t552 + t554 * t660) * t519 + (t661 * t552 + t554 * t662) * t520 + (-t651 * t552 + t554 * t652) * t526 + (t653 * t552 + t554 * t654) * t527;
t853 = t251 * t554 + t252 * t552;
t540 = Icges(3,4) * t553;
t494 = -Icges(3,2) * t551 + t540;
t495 = Icges(3,1) * t551 + t540;
t849 = -t850 + t548 * (t448 + t599) + t554 * (t468 + t416);
t337 = Icges(6,6) * t552 + t436 * t554;
t755 = Icges(6,4) * t510;
t438 = Icges(6,1) * t511 - t755;
t339 = Icges(6,5) * t552 + t438 * t554;
t294 = t339 * t727;
t434 = Icges(6,5) * t511 - Icges(6,6) * t510;
t738 = t434 * t554;
t335 = Icges(6,3) * t552 + t738;
t613 = t335 * t554 - t294;
t170 = -t337 * t731 - t613;
t336 = Icges(6,4) * t727 - Icges(6,2) * t731 - Icges(6,6) * t554;
t334 = Icges(6,5) * t727 - Icges(6,6) * t731 - Icges(6,3) * t554;
t479 = Icges(6,4) * t731;
t338 = Icges(6,1) * t727 - Icges(6,5) * t554 - t479;
t674 = -t552 * t334 - t338 * t726;
t171 = -t336 * t730 - t674;
t673 = t552 * t335 + t339 * t726;
t172 = -t337 * t730 + t673;
t610 = t337 * t510 - t334;
t745 = t336 * t510;
t605 = (-t171 * t554 + t172 * t552) * t867 + ((t170 - t294 + (t335 + t745) * t554 + t674) * t554 + t673 * t552) * t807 + ((t552 * t610 + t170 + t171 + t613) * t552 + (t552 * (-t338 * t511 + t745) + t172 - t673 + (t334 + t610) * t554) * t554) * t809;
t435 = Icges(6,2) * t511 + t755;
t604 = t877 * t511 / 0.2e1 + (-t435 / 0.2e1 + t438 / 0.2e1) * t510;
t433 = Icges(3,5) * t552 + t496 * t554;
t646 = -t493 * t554 + t433;
t709 = t551 * t552;
t514 = Icges(3,4) * t709;
t701 = t552 * t553;
t432 = Icges(3,1) * t701 - Icges(3,5) * t554 - t514;
t647 = -Icges(3,2) * t701 + t432 - t514;
t431 = Icges(3,6) * t552 + t494 * t554;
t648 = -t495 * t554 - t431;
t430 = Icges(3,4) * t701 - Icges(3,2) * t709 - Icges(3,6) * t554;
t649 = t495 * t552 + t430;
t846 = (-t646 * t552 + t554 * t647) * t551 + (t648 * t552 + t554 * t649) * t553;
t665 = -t435 * t554 + t339;
t666 = -Icges(6,2) * t727 + t338 - t479;
t667 = -t437 * t554 - t337;
t668 = t437 * t552 + t336;
t843 = (-t665 * t552 + t554 * t666) * t510 + (t667 * t552 + t554 * t668) * t511;
t842 = (t552 * t878 + t554 * t879) * t867 + (t873 * t552 + ((t872 + t890) * t554 + t880 + t882 + t889) * t554) * t807 + ((t874 * t552 - t879 + t880 - t883) * t552 + ((t874 + t885) * t554 + (-t363 * t520 - t404 * t527 + t872) * t552 - t873 + t878) * t554) * t809;
t841 = t876 * t520 / 0.2e1 + t875 * t527 / 0.2e1 + (t473 / 0.2e1 - t470 / 0.2e1) * t526 + (t456 / 0.2e1 - t453 / 0.2e1) * t519;
t840 = 0.4e1 * qJD(1);
t839 = 2 * qJD(2);
t837 = 2 * qJD(3);
t836 = 4 * qJD(3);
t708 = t551 * t554;
t656 = t548 * (-t599 + t777) + t554 * (pkin(2) * t708 - t468);
t193 = t655 + t656;
t677 = -t330 * t331 - t333 * t597;
t826 = m(5) * (t157 * t193 + t677);
t81 = t89 * t850;
t87 = t107 * t850;
t824 = m(6) * (-t81 - t87 + ((t269 + t282) * t554 + (t267 + t280) * t552) * t447);
t130 = t244 * t159;
t572 = (-t281 * t552 - t283 * t554) * t446;
t64 = t447 * t852 - t81;
t822 = m(6) * (t130 + t572 + t64);
t142 = t656 + t849;
t71 = t447 * t851 - t87;
t558 = t572 + t71;
t820 = m(6) * (t142 * t244 + t558);
t680 = -t280 * t281 - t282 * t283;
t816 = m(6) * (t107 * t142 + t680);
t573 = t853 * t447;
t813 = m(6) * (t267 * t396 - t269 * t393 - t573);
t812 = m(6) * (-t573 + (-t273 * t554 - t274 * t552) * t446);
t811 = m(6) * (t280 * t396 - t282 * t393 - t573);
t772 = rSges(3,1) * t553;
t630 = pkin(1) + t772;
t641 = rSges(3,2) * t709 + t554 * rSges(3,3);
t349 = -t552 * t630 + t546 + t641;
t516 = rSges(3,2) * t708;
t350 = -t516 + t630 * t554 + (rSges(3,3) + pkin(6)) * t552;
t497 = rSges(3,1) * t551 + rSges(3,2) * t553;
t465 = t497 * t552;
t466 = t497 * t554;
t806 = m(3) * (t349 * t465 - t350 * t466);
t149 = t474 * t475 * t640 + t262 * t279;
t147 = m(4) * t149;
t802 = m(4) * (t298 * t855 - t299 * t854);
t801 = m(4) * (t298 * t449 - t299 * t450);
t798 = m(5) * (t271 * t305 + t272 * t306);
t797 = m(5) * (t271 * t320 - t272 * t597);
t796 = m(5) * (t271 * t554 + t272 * t552);
t789 = m(6) * (-t573 + (-t284 * t554 - t285 * t552) * t446);
t788 = m(6) * (t251 * t273 + t252 * t274);
t787 = m(6) * (t251 * t284 + t252 * t285);
t127 = -t244 * t850 + t447 * t881;
t786 = m(6) * t127;
t785 = m(6) * (t251 * t393 - t252 * t396);
t784 = m(6) * t853;
t588 = Icges(6,5) * t510 + Icges(6,6) * t511;
t387 = t588 * t552;
t388 = t554 * t588;
t773 = (-t549 * t387 + (t554 * t388 + t843) * t552) * t807 + (-t548 * t388 + (t552 * t387 + t843) * t554) * t809;
t739 = t430 * t551;
t700 = t553 * t554;
t428 = Icges(3,5) * t701 - Icges(3,6) * t709 - Icges(3,3) * t554;
t658 = -t552 * t428 - t432 * t700;
t592 = Icges(3,5) * t553 - Icges(3,6) * t551;
t429 = Icges(3,3) * t552 + t554 * t592;
t657 = t552 * t429 + t433 * t700;
t575 = t850 * t833;
t603 = m(6) * t881;
t176 = t575 + t603 / 0.2e1;
t639 = t176 * qJD(1);
t638 = t820 / 0.2e1 + t773;
t636 = -t786 + t773;
t633 = t107 * t159 + t680;
t632 = t157 * t236 + t677;
t629 = -t475 - t545;
t370 = t433 * t701;
t608 = t429 * t554 - t370;
t606 = t431 * t551 - t428;
t602 = t640 * t777;
t591 = -Icges(3,5) * t551 - Icges(3,6) * t553;
t580 = t773 + ((t552 * t865 + t860) * t554 - t864 * t548) * t809 + ((t554 * t864 + t860) * t552 - t865 * t549) * t807;
t578 = t627 - t545;
t574 = t147 + t580;
t570 = -t602 + t656;
t564 = (-t435 + t438) * t511 - t877 * t510;
t568 = -t605 + (t434 * t552 + t510 * t667 + t511 * t665 + t554 * t564) * t809 + (-t510 * t668 + t511 * t666 + t552 * t564 - t738) * t807;
t566 = t576 - t545;
t565 = -t604 + t861 * (t336 * t511 + t338 * t510);
t561 = t605 + t842;
t560 = t604 + t841;
t557 = -t842 + t568 + (t519 * t661 + t520 * t659 + t526 * t653 + t527 * t651 + t893 * t552 + t863 * t554) * t809 + (-t519 * t662 + t520 * t660 - t526 * t654 + t527 * t652 + t552 * t863 - t891) * t807;
t556 = t565 - t841 + (t361 * t520 + t363 * t519 + t402 * t527 + t404 * t526) * t861;
t499 = -rSges(3,2) * t551 + t772;
t460 = t591 * t554;
t459 = t591 * t552;
t398 = t629 * t554;
t395 = t629 * t552;
t304 = t578 * t554;
t302 = t578 * t552;
t270 = t566 * t554;
t268 = t566 * t552;
t255 = -t602 + t279;
t229 = -t431 * t708 + t657;
t228 = -t430 * t708 - t658;
t227 = -t431 * t709 - t608;
t175 = t575 - t603 / 0.2e1;
t167 = t570 + t655;
t151 = -t228 * t554 + t229 * t552;
t150 = -(-t552 * (-t432 * t553 + t739) - t428 * t554) * t554 + t227 * t552;
t146 = m(6) * (-t281 * t554 + t283 * t552) + m(5) * (-t331 * t554 + t333 * t552);
t145 = t146 * qJD(3);
t129 = t570 + t849;
t109 = t784 + t796;
t105 = t604 + t785;
t97 = t789 / 0.2e1;
t94 = t811 / 0.2e1;
t93 = t812 / 0.2e1;
t88 = t813 / 0.2e1;
t69 = (t227 - t370 + (t429 + t739) * t554 + t658) * t554 + t657 * t552;
t68 = (t554 * t606 + t229 - t657) * t554 + (t552 * t606 + t228 + t608) * t552;
t67 = t682 + t683;
t56 = t686 + t687;
t38 = t560 + t787 + t797 + t801;
t31 = t822 / 0.2e1;
t27 = t773 + t786;
t26 = t27 * qJD(5);
t24 = t560 + t866 + t806 + t802 + t798 + t788 + (t495 / 0.2e1 + t494 / 0.2e1) * t553;
t22 = t824 / 0.2e1;
t21 = m(6) * t71 + t773;
t20 = m(6) * t64 + t773;
t16 = t574 + t816 + t826;
t15 = t94 - t789 / 0.2e1 + t605;
t14 = t97 - t811 / 0.2e1 + t605;
t13 = t580 + t888;
t12 = t22 - t822 / 0.2e1 + t638;
t11 = t22 + t31 - t820 / 0.2e1 + t773;
t10 = t31 - t824 / 0.2e1 + t638;
t9 = t88 - t812 / 0.2e1 + t605;
t8 = t93 - t813 / 0.2e1 + t605;
t7 = t94 + t97 + t568;
t6 = t88 + t93 + t568;
t4 = t561 + (-t69 / 0.2e1 + t151 / 0.2e1) * t554 + (t150 / 0.2e1 + t68 / 0.2e1) * t552;
t3 = t561 + t862;
t2 = t561 - t862;
t1 = t557 + t634 + t635;
t5 = [qJD(2) * t24 + qJD(3) * t38 + qJD(4) * t109 + qJD(5) * t105, t24 * qJD(1) + t1 * qJD(3) + t56 * qJD(4) + t6 * qJD(5) + (m(3) * ((-t349 * t554 - t350 * t552) * t499 + (-t465 * t554 + t466 * t552) * t497) / 0.2e1 + (t298 * t398 + t299 * t395) * t835 + (t271 * t304 + t272 * t302 - t301 * t306 - t303 * t305) * t834 + (t251 * t270 + t252 * t268 - t267 * t274 - t269 * t273) * t833) * t839 + (t557 + (t551 * t648 + t553 * t646) * t809 + t69 * t867 + (t150 + t68) * t810 + (-t551 * t649 + t553 * t647 + t151) * t807 + (t549 / 0.2e1 + t548 / 0.2e1) * t592) * qJD(2), t38 * qJD(1) + t1 * qJD(2) + t557 * qJD(3) + t67 * qJD(4) + t7 * qJD(5) + ((t571 + (-t449 * t554 + t450 * t552) * t474) * t835 + (t679 + (-t320 + t330) * t597) * t834 + (-t280 * t285 - t282 * t284 + t685) * t833) * t837, qJD(1) * t109 + qJD(2) * t56 + qJD(3) * t67 + qJD(5) * t175, t105 * qJD(1) + t6 * qJD(2) + t7 * qJD(3) + t175 * qJD(4) + ((-t573 + (-t393 * t554 + t396 * t552) * t446) * m(6) + t568) * qJD(5); (t556 - (t495 + t494) * t553 / 0.2e1 - t866) * qJD(1) + t4 * qJD(2) + t2 * qJD(3) + t55 * qJD(4) + t9 * qJD(5) + (-t806 / 0.4e1 - t802 / 0.4e1 - t798 / 0.4e1 - t788 / 0.4e1) * t840, t4 * qJD(1) + (m(6) * (t129 * t89 - t267 * t268 - t269 * t270) + m(5) * (t132 * t167 - t301 * t302 - t303 * t304) + m(4) * (t181 * t255 - t395 * t855 - t398 * t854) + m(3) * ((t552 * (rSges(3,1) * t701 - t641) + t554 * (rSges(3,1) * t700 + rSges(3,3) * t552 - t516)) * (-t465 * t552 - t466 * t554) + t640 * t499 * t497) + (t549 * t459 + (-t554 * t460 + t846) * t552) * t807 + (t548 * t460 + (-t552 * t459 + t846) * t554) * t809 + t580) * qJD(2) + t13 * qJD(3) + t20 * qJD(5), t2 * qJD(1) + t13 * qJD(2) + t580 * qJD(3) + t11 * qJD(5) + (-t816 / 0.4e1 - t826 / 0.4e1 - t147 / 0.4e1) * t836 + ((t633 + t46) * t833 + (t632 + t75) * t834 + (t114 + t149) * t835) * t837, t886, t9 * qJD(1) + t20 * qJD(2) + t11 * qJD(3) + (m(6) * (t64 + t127) + t636) * qJD(5); t556 * qJD(1) + t3 * qJD(2) + t561 * qJD(3) - t65 * qJD(4) + t15 * qJD(5) + (-t797 / 0.4e1 - t801 / 0.4e1 - t787 / 0.4e1) * t840, t3 * qJD(1) + t16 * qJD(3) + t12 * qJD(5) + ((t107 * t129 + t142 * t89 - t268 * t280 - t270 * t282 + t681) * t833 + (t132 * t193 + t157 * t167 - t302 * t330 - t304 * t597 + t678) * t834 + (t255 * t262 + (-t395 * t552 - t398 * t554) * t474 + t114) * t835) * t839 + (t580 - t888) * qJD(2), t561 * qJD(1) + t16 * qJD(2) + t574 * qJD(3) + t21 * qJD(5) + (m(6) * t633 / 0.4e1 + m(5) * t632 / 0.4e1) * t836, -t887, t15 * qJD(1) + t12 * qJD(2) + t21 * qJD(3) + (m(6) * (t71 + t127) + t636) * qJD(5); -t55 * qJD(2) + t65 * qJD(3) + t176 * qJD(5) + (-t796 / 0.4e1 - t784 / 0.4e1) * t840, -t886 + ((-t268 * t554 + t270 * t552) * t833 + (-t302 * t554 + t304 * t552) * t834) * t839 + t145, qJD(2) * t146 + t145 + t887, 0, t639; (t565 - t785) * qJD(1) + t8 * qJD(2) + t14 * qJD(3) - t176 * qJD(4) + t605 * qJD(5), t8 * qJD(1) + ((t129 * t244 + (-t268 * t552 - t270 * t554) * t446) * m(6) + t773) * qJD(2) + t10 * qJD(3) + t26, t14 * qJD(1) + t10 * qJD(2) + ((t130 - t71 + t558) * m(6) + t773) * qJD(3) + t26, -t639, qJD(1) * t605 + t26 + (qJD(2) + qJD(3)) * t27;];
Cq = t5;
