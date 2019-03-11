% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:38:01
% EndTime: 2019-03-09 02:38:38
% DurationCPUTime: 28.04s
% Computational Cost: add. (76998->827), mult. (54187->1158), div. (0->0), fcn. (57474->12), ass. (0->479)
t523 = qJ(3) + pkin(10);
t516 = sin(t523);
t524 = qJ(1) + pkin(9);
t517 = sin(t524);
t705 = t516 * t517;
t821 = m(7) / 0.2e1;
t823 = m(6) / 0.2e1;
t633 = t821 + t823;
t513 = t517 ^ 2;
t520 = cos(t524);
t514 = t520 ^ 2;
t644 = t513 + t514;
t445 = t644 * t516;
t874 = m(6) + m(7);
t655 = t874 * t445 / 0.2e1;
t305 = -t516 * t633 + t655;
t519 = cos(t523);
t522 = pkin(11) + qJ(6);
t515 = sin(t522);
t518 = cos(t522);
t752 = rSges(7,1) * t518;
t592 = -rSges(7,2) * t515 + t752;
t387 = -rSges(7,3) * t519 + t516 * t592;
t721 = qJ(5) * t519;
t760 = pkin(4) * t516;
t528 = sin(qJ(3));
t762 = pkin(3) * t528;
t616 = -t721 + t760 + t762;
t527 = -pkin(8) - qJ(5);
t683 = qJ(5) + t527;
t526 = cos(pkin(11));
t511 = pkin(5) * t526 + pkin(4);
t758 = -pkin(4) + t511;
t567 = t516 * t758 + t519 * t683 + t387 + t616;
t244 = t567 * t517;
t246 = t567 * t520;
t525 = sin(pkin(11));
t753 = rSges(6,1) * t526;
t595 = -rSges(6,2) * t525 + t753;
t554 = t595 * t516;
t862 = rSges(6,3) * t519 - t554;
t600 = t616 - t862;
t318 = t600 * t517;
t320 = t600 * t520;
t530 = cos(qJ(3));
t761 = pkin(3) * t530;
t512 = pkin(2) + t761;
t687 = t520 * t525;
t700 = t517 * t526;
t443 = -t519 * t687 + t700;
t686 = t520 * t526;
t701 = t517 * t525;
t444 = t519 * t686 + t701;
t596 = t444 * rSges(6,1) + t443 * rSges(6,2);
t757 = -qJ(4) - pkin(7);
t497 = t517 * t757;
t763 = cos(qJ(1)) * pkin(1);
t617 = -t497 + t763;
t745 = rSges(6,3) + qJ(5);
t759 = pkin(4) * t519;
t832 = -t516 * t745 - t759;
t243 = (t512 - t832) * t520 + t596 + t617;
t692 = t519 * t520;
t702 = t517 * t519;
t646 = -t517 * t512 - t520 * t757;
t764 = sin(qJ(1)) * pkin(1);
t601 = t646 - t764;
t441 = t519 * t701 + t686;
t442 = t519 * t700 - t687;
t843 = -t442 * rSges(6,1) + t441 * rSges(6,2);
t864 = t517 * t832 + t601 + t843;
t678 = t243 * t702 + t692 * t864;
t688 = t520 * t515;
t415 = t517 * t518 - t519 * t688;
t416 = t515 * t517 + t518 * t692;
t593 = t416 * rSges(7,1) + t415 * rSges(7,2);
t703 = t516 * t520;
t602 = pkin(5) * t701 - t527 * t703;
t708 = t511 * t519;
t750 = rSges(7,3) * t516;
t211 = (t512 + t708 + t750) * t520 + t593 + t602 + t617;
t603 = -pkin(5) * t687 + t511 * t702;
t744 = -rSges(7,3) + t527;
t695 = t518 * t520;
t413 = t515 * t702 + t695;
t414 = t518 * t702 - t688;
t844 = -t414 * rSges(7,1) + t413 * rSges(7,2);
t865 = t705 * t744 + t601 - t603 + t844;
t681 = t211 * t702 + t692 * t865;
t741 = (-t318 * t703 + t320 * t705 + t678) * t823 + (-t244 * t703 + t246 * t705 + t681) * t821;
t269 = (t762 + t744 * t519 + (t511 + t592) * t516) * t517;
t649 = t516 * rSges(7,2) * t688 + rSges(7,3) * t692;
t691 = t519 * t527;
t270 = (-t762 - t691 + (-t511 - t752) * t516) * t520 + t649;
t479 = pkin(4) * t705;
t279 = t479 + (-t519 * t745 + t554 + t762) * t517;
t471 = qJ(5) * t692;
t647 = t516 * rSges(6,2) * t687 + rSges(6,3) * t692;
t280 = t471 + (-t762 + (-pkin(4) - t753) * t516) * t520 + t647;
t742 = ((t279 * t520 + t280 * t517) * t516 + t678) * t823 + ((t269 * t520 + t270 * t517) * t516 + t681) * t821;
t7 = t742 - t741;
t888 = -t7 * qJD(1) + t305 * qJD(2);
t580 = Icges(7,5) * t518 - Icges(7,6) * t515;
t375 = -Icges(7,3) * t519 + t516 * t580;
t730 = Icges(7,4) * t518;
t586 = -Icges(7,2) * t515 + t730;
t379 = -Icges(7,6) * t519 + t516 * t586;
t731 = Icges(7,4) * t515;
t588 = Icges(7,1) * t518 - t731;
t383 = -Icges(7,5) * t519 + t516 * t588;
t192 = t375 * t705 - t379 * t413 + t383 * t414;
t290 = Icges(7,5) * t414 - Icges(7,6) * t413 + Icges(7,3) * t705;
t401 = Icges(7,4) * t414;
t293 = -Icges(7,2) * t413 + Icges(7,6) * t705 + t401;
t400 = Icges(7,4) * t413;
t297 = -Icges(7,1) * t414 - Icges(7,5) * t705 + t400;
t140 = t290 * t705 - t293 * t413 - t297 * t414;
t292 = Icges(7,5) * t416 + Icges(7,6) * t415 + Icges(7,3) * t703;
t732 = Icges(7,4) * t416;
t295 = Icges(7,2) * t415 + Icges(7,6) * t703 + t732;
t402 = Icges(7,4) * t415;
t298 = Icges(7,1) * t416 + Icges(7,5) * t703 + t402;
t141 = t292 * t705 - t413 * t295 + t414 * t298;
t579 = t140 * t517 + t141 * t520;
t17 = t519 * t192 - t516 * t579;
t299 = rSges(7,3) * t705 - t844;
t389 = t519 * t592 + t750;
t176 = (t389 * t517 - t299) * t516;
t301 = rSges(7,3) * t703 + t593;
t344 = -rSges(7,1) * t516 * t695 + t649;
t177 = (-t387 * t520 - t344) * t519 + (-t389 * t520 + t301) * t516;
t109 = t176 * t517 - t177 * t520;
t824 = m(5) / 0.2e1;
t467 = rSges(5,1) * t516 + rSges(5,2) * t519;
t550 = t467 + t762;
t849 = t550 * t520;
t850 = t550 * t517;
t859 = t517 * t850 + t520 * t849;
t625 = (t269 * t517 - t270 * t520) * t821 + (t279 * t517 - t280 * t520) * t823 + t859 * t824;
t845 = t517 * t244 + t246 * t520;
t626 = -t845 * m(7) / 0.2e1 + (-t517 * t318 - t320 * t520) * t823 - m(5) * t859 / 0.2e1;
t32 = t626 - t625;
t756 = m(7) * qJD(6);
t887 = t32 * qJD(1) - t109 * t756 / 0.2e1;
t331 = -rSges(7,1) * t413 - rSges(7,2) * t414;
t332 = rSges(7,1) * t415 - rSges(7,2) * t416;
t886 = m(7) * (t211 * t332 - t331 * t865);
t194 = t375 * t703 + t379 * t415 + t383 * t416;
t142 = t290 * t703 + t415 * t293 - t297 * t416;
t143 = t292 * t703 + t415 * t295 + t416 * t298;
t578 = t142 * t517 + t143 * t520;
t884 = -t519 * t194 + t516 * t578;
t309 = Icges(6,5) * t442 - Icges(6,6) * t441 + Icges(6,3) * t705;
t312 = Icges(6,4) * t442 - Icges(6,2) * t441 + Icges(6,6) * t705;
t315 = Icges(6,1) * t442 - Icges(6,4) * t441 + Icges(6,5) * t705;
t860 = t443 * t312 + t444 * t315;
t162 = t309 * t703 + t860;
t311 = Icges(6,5) * t444 + Icges(6,6) * t443 + Icges(6,3) * t703;
t314 = Icges(6,4) * t444 + Icges(6,2) * t443 + Icges(6,6) * t703;
t317 = Icges(6,1) * t444 + Icges(6,4) * t443 + Icges(6,5) * t703;
t880 = -t142 * t520 + t143 * t517;
t883 = (t311 * t703 + t443 * t314 + t444 * t317) * t517 - t162 * t520 + t880;
t77 = -t140 * t520 + t141 * t517;
t234 = t301 * t519 + t387 * t703;
t879 = t299 * t519 + t387 * t705;
t846 = t234 * t517 - t520 * t879;
t716 = t290 * t519;
t876 = t293 * t515 + t297 * t518;
t167 = t516 * t876 + t716;
t733 = Icges(5,4) * t516;
t465 = Icges(5,1) * t519 - t733;
t386 = Icges(5,5) * t517 + t465 * t520;
t462 = Icges(5,2) * t519 + t733;
t723 = Icges(6,3) * t519;
t581 = Icges(6,5) * t526 - Icges(6,6) * t525;
t852 = t516 * t581;
t541 = t723 - t852;
t878 = t314 * t525 - t317 * t526 - t386 + (t462 + t541) * t520;
t877 = -Icges(4,5) * t528 - Icges(5,5) * t516 - Icges(4,6) * t530 - Icges(5,6) * t519;
t847 = -t211 * t517 - t520 * t865;
t799 = -t519 / 0.2e1;
t734 = Icges(4,4) * t528;
t483 = Icges(4,2) * t530 + t734;
t486 = Icges(4,1) * t530 - t734;
t866 = (t486 / 0.2e1 - t483 / 0.2e1) * t528;
t573 = -t312 * t441 + t315 * t442;
t858 = -0.2e1 * t445;
t853 = t387 * t517;
t851 = t516 * t683;
t205 = (t331 * t520 - t332 * t517) * t516;
t634 = m(6) / 0.4e1 + m(7) / 0.4e1;
t704 = t516 * t519;
t648 = t644 * t704;
t848 = t634 * (t648 - t704);
t842 = t877 * t517;
t841 = t877 * t520;
t521 = Icges(4,4) * t530;
t484 = -Icges(4,2) * t528 + t521;
t485 = Icges(4,1) * t528 + t521;
t474 = Icges(5,4) * t705;
t385 = Icges(5,1) * t702 - Icges(5,5) * t520 - t474;
t840 = -Icges(5,2) * t702 - t312 * t525 + t315 * t526 - t541 * t517 + t385 - t474;
t798 = -t520 / 0.2e1;
t838 = qJD(3) * t798;
t837 = t311 * t705 - t441 * t314 + t442 * t317;
t587 = Icges(6,4) * t526 - Icges(6,2) * t525;
t834 = -Icges(6,6) * t519 + t516 * t587;
t589 = Icges(6,1) * t526 - Icges(6,4) * t525;
t833 = -Icges(6,5) * t519 + t516 * t589;
t699 = t517 * t528;
t493 = Icges(4,4) * t699;
t698 = t517 * t530;
t425 = Icges(4,1) * t698 - Icges(4,5) * t520 - t493;
t652 = -Icges(4,2) * t698 + t425 - t493;
t423 = Icges(4,4) * t698 - Icges(4,2) * t699 - Icges(4,6) * t520;
t654 = t485 * t517 + t423;
t381 = Icges(5,4) * t702 - Icges(5,2) * t705 - Icges(5,6) * t520;
t507 = Icges(5,4) * t519;
t735 = Icges(5,1) * t516;
t590 = -t507 - t735;
t661 = -t590 * t517 + t381;
t831 = t519 * t661 + t528 * t652 + t530 * t654;
t433 = (-Icges(7,2) * t518 - t731) * t516;
t436 = (-Icges(7,1) * t515 - t730) * t516;
t830 = -t515 * (t383 / 0.2e1 + t433 / 0.2e1) + t518 * (t436 / 0.2e1 - t379 / 0.2e1);
t726 = Icges(5,2) * t516;
t382 = Icges(5,6) * t517 + (t507 - t726) * t520;
t660 = t590 * t520 - t382;
t829 = (-t309 * t520 + t311 * t517) * t519 + t660 * t702 + (t517 * t878 + t520 * t840) * t516;
t827 = 0.4e1 * qJD(1);
t826 = 0.2e1 * qJD(3);
t825 = m(4) / 0.2e1;
t819 = -t884 / 0.2e1;
t818 = t77 / 0.2e1;
t817 = t880 / 0.2e1;
t574 = t299 * t520 - t301 * t517;
t146 = t574 * t519 + (-t344 * t517 - t520 * t853) * t516;
t190 = t574 * t516;
t679 = -t234 * t702 + t692 * t879;
t813 = m(7) * (-t146 * t519 + (t176 * t520 + t177 * t517 + t190) * t516 + t679);
t812 = m(7) * (t176 * t865 + t177 * t211 - t234 * t270 + t269 * t879);
t439 = (-rSges(7,1) * t515 - rSges(7,2) * t518) * t516;
t809 = m(7) * (-t244 * t332 + t246 * t331 + t439 * t847);
t325 = -Icges(7,5) * t413 - Icges(7,6) * t414;
t671 = -Icges(7,2) * t414 - t297 - t400;
t673 = -Icges(7,1) * t413 - t293 - t401;
t129 = -t325 * t519 + (-t515 * t671 + t518 * t673) * t516;
t805 = t129 / 0.2e1;
t803 = t516 / 0.2e1;
t801 = t517 / 0.2e1;
t800 = t517 / 0.4e1;
t797 = -t520 / 0.4e1;
t796 = t520 / 0.2e1;
t795 = -t525 / 0.2e1;
t794 = t526 / 0.2e1;
t510 = t520 * pkin(7);
t755 = rSges(4,1) * t530;
t618 = pkin(2) + t755;
t645 = rSges(4,2) * t699 + t520 * rSges(4,3);
t345 = -t517 * t618 + t510 + t645 - t764;
t685 = t520 * t528;
t495 = rSges(4,2) * t685;
t346 = t763 - t495 + t618 * t520 + (rSges(4,3) + pkin(7)) * t517;
t487 = rSges(4,1) * t528 + rSges(4,2) * t530;
t452 = t487 * t517;
t453 = t487 * t520;
t792 = m(4) * (t345 * t452 - t346 * t453);
t565 = rSges(5,1) * t702 - rSges(5,2) * t705 - t520 * rSges(5,3);
t333 = -t565 + t601;
t613 = -rSges(5,2) * t703 + t517 * rSges(5,3);
t754 = rSges(5,1) * t519;
t334 = (t512 + t754) * t520 + t613 + t617;
t790 = m(5) * (t333 * t850 - t334 * t849);
t789 = m(5) * (t333 * t520 + t334 * t517);
t782 = m(6) * (t243 * t280 + t279 * t864);
t231 = t243 * t703;
t781 = m(6) * (-t705 * t864 + t231);
t780 = m(6) * (t243 * t517 + t520 * t864);
t776 = m(7) * (t190 * t445 + t679);
t775 = m(7) * (t211 * t270 + t269 * t865);
t199 = t211 * t703;
t774 = m(7) * (-t705 * t865 + t199);
t773 = m(7) * t847;
t772 = m(7) * (-t234 * t703 - t705 * t879);
t771 = m(7) * t846;
t215 = t331 * t517 + t332 * t520;
t768 = m(7) * (-t215 * t519 - t439 * t445);
t766 = m(7) * t205;
t765 = m(7) * t215;
t748 = t517 * t17;
t746 = t520 * t884;
t339 = t379 * t517;
t341 = t383 * t517;
t560 = -t375 * t517 + t876;
t121 = -t560 * t519 + (t339 * t515 - t341 * t518 + t290) * t516;
t720 = t121 * t520;
t340 = t379 * t520;
t342 = t383 * t520;
t575 = -t295 * t515 + t298 * t518;
t559 = -t375 * t520 - t575;
t122 = -t559 * t519 + (t340 * t515 - t342 * t518 + t292) * t516;
t719 = t122 * t517;
t715 = t292 * t519;
t714 = t375 * t519;
t712 = t381 * t516;
t710 = t423 * t528;
t707 = t515 * t379;
t380 = Icges(7,6) * t516 + t519 * t586;
t706 = t515 * t380;
t697 = t518 * t383;
t384 = Icges(7,5) * t516 + t519 * t588;
t696 = t518 * t384;
t430 = (-Icges(7,5) * t515 - Icges(7,6) * t518) * t516;
t693 = t519 * t430;
t684 = t520 * t530;
t676 = -t244 * t702 - t246 * t692;
t674 = -t318 * t702 - t320 * t692;
t672 = Icges(7,1) * t415 - t295 - t732;
t670 = -Icges(7,2) * t416 + t298 + t402;
t377 = Icges(5,5) * t702 - Icges(5,6) * t705 - Icges(5,3) * t520;
t668 = -t517 * t377 - t385 * t692;
t583 = Icges(5,5) * t519 - Icges(5,6) * t516;
t378 = Icges(5,3) * t517 + t520 * t583;
t667 = t517 * t378 + t386 * t692;
t665 = -t517 * (pkin(2) * t517 - t510 + t646) + t520 * (-t517 * pkin(7) - t497 + (-pkin(2) + t512) * t520);
t421 = Icges(4,5) * t698 - Icges(4,6) * t699 - Icges(4,3) * t520;
t664 = -t517 * t421 - t425 * t684;
t585 = Icges(4,5) * t530 - Icges(4,6) * t528;
t422 = Icges(4,3) * t517 + t520 * t585;
t426 = Icges(4,5) * t517 + t486 * t520;
t663 = t517 * t422 + t426 * t684;
t662 = -t379 + t436;
t659 = t383 + t433;
t424 = Icges(4,6) * t517 + t484 * t520;
t653 = -t485 * t520 - t424;
t651 = -t483 * t520 + t426;
t643 = qJD(1) * t516;
t642 = qJD(1) * t519;
t641 = qJD(6) * t516;
t640 = qJD(6) * t519;
t639 = t109 * qJD(4);
t118 = 0.2e1 * (t146 / 0.4e1 - t215 / 0.4e1) * m(7);
t638 = t118 * qJD(2);
t249 = t633 * t858;
t636 = t249 * qJD(1);
t629 = t819 + t884 / 0.2e1;
t624 = t705 / 0.4e1;
t468 = qJ(5) * t516 + t759;
t615 = -t468 - t761;
t614 = rSges(5,2) * t516 - t754 - t761;
t354 = t386 * t702;
t610 = t378 * t520 - t354;
t365 = t426 * t698;
t609 = t422 * t520 - t365;
t608 = t382 * t516 - t377;
t607 = t424 * t528 - t421;
t605 = t468 * t644 + t665;
t604 = t644 * t762;
t599 = -rSges(6,3) * t516 - t519 * t595 + t615;
t111 = (t299 - (t759 + t851) * t517 + t603) * t517 + (t301 + (-t468 + t708) * t520 + t602) * t520 + t605;
t144 = t517 * (rSges(6,3) * t705 - t843) + t520 * (rSges(6,3) * t703 + t596) + t605;
t38 = m(6) * (t144 * t445 + t674) + m(7) * (t111 * t445 + t676);
t326 = Icges(7,5) * t415 - Icges(7,6) * t416;
t130 = -t326 * t519 + (-t515 * t670 + t518 * t672) * t516;
t150 = -t413 * t659 + t414 * t662 + t430 * t705;
t151 = t415 * t659 + t416 * t662 + t430 * t703;
t598 = t809 / 0.2e1 + (t130 + t151) * t800 + (t129 + t150) * t797;
t168 = t516 * t575 - t715;
t577 = -t167 * t517 + t168 * t520;
t571 = t697 - t707;
t568 = -t511 * t516 - t691;
t566 = -t519 * t758 - t389 + t615 + t851;
t104 = t325 * t705 - t413 * t671 + t414 * t673;
t105 = t326 * t705 - t413 * t670 + t414 * t672;
t49 = -t104 * t520 + t105 * t517;
t106 = t325 * t703 + t415 * t671 + t416 * t673;
t107 = t326 * t703 + t415 * t670 + t416 * t672;
t50 = -t106 * t520 + t107 * t517;
t563 = t49 * t798 + t50 * t801;
t376 = Icges(7,3) * t516 + t519 * t580;
t556 = t376 - t571;
t551 = t517 * (qJ(5) * t702 - t479) + t520 * (-pkin(4) * t703 + t471) - t604;
t547 = -t528 * t651 + t530 * t653;
t546 = t17 * t800 + t884 * t797 - t748 / 0.4e1 + t746 / 0.4e1 + (t624 - t705 / 0.4e1) * t880;
t539 = t516 * t560 + t716;
t538 = t516 * t559 + t715;
t537 = t516 * t556 + t714;
t134 = -t380 * t413 + t384 * t414 + t517 * t537;
t135 = t380 * t415 + t384 * t416 + t520 * t537;
t164 = -t556 * t519 + (t375 + t696 - t706) * t516;
t224 = t516 * t571 - t714;
t536 = t164 * t799 + t224 * t803 + t812 / 0.2e1 + (t121 + t134) * t624 + (t122 + t135) * t703 / 0.4e1 + (-t167 + t192) * t702 / 0.4e1 + (t168 + t194) * t692 / 0.4e1;
t533 = t833 * t794 + t834 * t795 + t697 / 0.2e1 - t707 / 0.2e1 + t507 + t735 / 0.2e1 - t726 / 0.2e1 - Icges(6,3) * t516 / 0.2e1 + t581 * t799 - t376 / 0.2e1;
t408 = Icges(6,6) * t516 + t519 * t587;
t410 = Icges(6,5) * t516 + t519 * t589;
t532 = t410 * t794 + t408 * t795 + t696 / 0.2e1 - t706 / 0.2e1 + t465 / 0.2e1 - t462 / 0.2e1 + t852 / 0.2e1 - t723 / 0.2e1 + t375 / 0.2e1;
t489 = -rSges(4,2) * t528 + t755;
t393 = t614 * t520;
t391 = t614 * t517;
t353 = t833 * t520;
t352 = t833 * t517;
t351 = t834 * t520;
t350 = t834 * t517;
t336 = -t452 * t517 - t453 * t520;
t321 = t599 * t520;
t319 = t599 * t517;
t304 = t803 * t874 + t655;
t275 = -t467 * t644 - t604;
t262 = 0.4e1 * t848;
t255 = -t332 * t519 - t439 * t703;
t254 = t331 * t519 + t439 * t705;
t253 = -t424 * t685 + t663;
t252 = -t423 * t685 - t664;
t251 = -t424 * t699 - t609;
t248 = t634 * t858 + t655;
t247 = t566 * t520;
t245 = t566 * t517;
t228 = -t382 * t703 + t667;
t227 = -t381 * t703 - t668;
t226 = -t382 * t705 - t610;
t213 = -t765 / 0.2e1;
t200 = -t766 / 0.2e1;
t184 = t520 * (-rSges(6,1) * t516 * t686 + t647) + t862 * t513 + t551;
t181 = -t693 + (-t515 * t659 + t518 * t662) * t516;
t178 = t768 / 0.2e1;
t175 = -t252 * t520 + t253 * t517;
t174 = -(-t517 * (-t425 * t530 + t710) - t421 * t520) * t520 + t251 * t517;
t155 = -t771 / 0.2e1;
t154 = -t227 * t520 + t228 * t517;
t153 = -(-t517 * (-t385 * t519 + t712) - t377 * t520) * t520 + t226 * t517;
t149 = t772 / 0.2e1;
t133 = (-t471 + t344 + (t568 + t760) * t520) * t520 + (t479 - t853 + (t568 - t721) * t517) * t517 + t551;
t119 = (t146 + t215) * t821;
t103 = t109 * qJD(3) * t821;
t99 = t776 / 0.2e1;
t98 = -t340 * t415 - t342 * t416 + t520 * t538;
t97 = -t339 * t415 - t341 * t416 + t520 * t539;
t96 = t340 * t413 - t342 * t414 + t517 * t538;
t95 = t339 * t413 - t341 * t414 + t517 * t539;
t86 = -(t309 * t705 + t573) * t520 + t517 * t837;
t83 = t774 + t781;
t76 = t886 - t693 / 0.2e1 + t830 * t516;
t73 = -t224 * t519 + t516 * t577;
t72 = (t251 - t365 + (t422 + t710) * t520 + t664) * t520 + t663 * t517;
t71 = (t520 * t607 + t253 - t663) * t520 + (t517 * t607 + t252 + t609) * t517;
t70 = -t773 + t780 + t789;
t58 = (t226 - t354 + (t378 + t712) * t520 + t668) * t520 + t667 * t517;
t57 = (t520 * t608 + t228 - t667) * t520 + (t517 * t608 + t227 + t610) * t517;
t54 = t111 * t215 + t439 * t845;
t53 = t155 + t765 / 0.2e1;
t52 = t213 + t155;
t51 = t213 + t771 / 0.2e1;
t48 = t517 * t98 - t520 * t97;
t47 = t517 * t96 - t520 * t95;
t44 = t149 + t766 / 0.2e1;
t43 = t200 + t149;
t42 = t200 - t772 / 0.2e1;
t41 = t146 * t190 + t176 * t879 - t177 * t234;
t36 = t813 / 0.2e1;
t35 = -t151 * t519 + (t106 * t517 + t107 * t520) * t516;
t34 = -t150 * t519 + (t104 * t517 + t105 * t520) * t516;
t31 = t625 + t626;
t28 = t573 * t520 + (t162 - t837 - t860) * t517;
t23 = (-t164 + t577) * t519 + (t121 * t517 + t122 * t520 + t224) * t516;
t22 = (t485 / 0.2e1 + t484 / 0.2e1) * t530 + t866 + t792 + t790 + t782 + t775 + t533 * t519 + t532 * t516;
t21 = t99 + t36 - t768 / 0.2e1;
t20 = t178 + t99 - t813 / 0.2e1;
t19 = t178 + t36 - t776 / 0.2e1;
t14 = (-t135 + t578) * t519 + (t517 * t97 + t520 * t98 + t194) * t516;
t13 = (-t134 + t579) * t519 + (t517 * t95 + t520 * t96 + t192) * t516;
t10 = m(7) * t54 + t563;
t8 = t741 + t742;
t6 = t629 * t705;
t5 = m(7) * t41 + (t746 / 0.2e1 - t748 / 0.2e1 - t23 / 0.2e1) * t519 + (t14 * t796 + t13 * t801 + t73 / 0.2e1) * t516;
t4 = (t817 - t880 / 0.2e1 + t175 / 0.2e1 + t154 / 0.2e1 - t58 / 0.2e1 - t72 / 0.2e1) * t520 + (-t77 / 0.2e1 + t818 + t71 / 0.2e1 + t57 / 0.2e1 + t153 / 0.2e1 + t28 / 0.2e1 + t174 / 0.2e1 + t86 / 0.2e1) * t517;
t3 = t536 + t598;
t2 = t546 + (-t151 / 0.4e1 - t130 / 0.4e1) * t517 + (t150 / 0.4e1 + t129 / 0.4e1) * t520 + t536 - t809 / 0.2e1;
t1 = t546 + (-t224 / 0.2e1 + (-t135 / 0.4e1 - t122 / 0.4e1) * t520 + (-t134 / 0.4e1 - t121 / 0.4e1) * t517) * t516 + (t164 / 0.2e1 + (-t194 / 0.4e1 - t168 / 0.4e1) * t520 + (-t192 / 0.4e1 + t167 / 0.4e1) * t517) * t519 - t812 / 0.2e1 + t598;
t9 = [t22 * qJD(3) + t70 * qJD(4) + t83 * qJD(5) + t76 * qJD(6), 0, t22 * qJD(1) + t31 * qJD(4) + t8 * qJD(5) + t3 * qJD(6) + ((t333 * t393 + t334 * t391) * t824 + (t243 * t319 - t279 * t320 - t280 * t318 + t321 * t864) * t823 + (t211 * t245 - t244 * t270 - t246 * t269 + t247 * t865) * t821 + ((-t345 * t520 - t346 * t517) * t489 + (-t452 * t520 + t453 * t517) * t487) * t825) * t826 + (-t408 * t441 + t410 * t442 - t528 * t654 + t530 * t652 + t134 + t154 + t175 + t840 * t519 + (t350 * t525 - t352 * t526 + t309 - t661) * t516 + t883) * t838 + ((t585 + t583) * (t514 / 0.2e1 + t513 / 0.2e1) + t719 / 0.2e1 - t720 / 0.2e1 + (t72 + t58 + t883) * t796 + (t408 * t443 + t410 * t444 + t528 * t653 + t530 * t651 + t135 - t878 * t519 + (t351 * t525 - t353 * t526 + t311 + t660) * t516) * t801 - (t174 + t153 + t86 + t71 + t57 + t28) * t517 / 0.2e1) * qJD(3), qJD(1) * t70 + qJD(3) * t31 + qJD(5) * t248 + qJD(6) * t52, qJD(1) * t83 + qJD(3) * t8 + qJD(4) * t248 + qJD(6) * t43, t76 * qJD(1) + t3 * qJD(3) + t52 * qJD(4) + t43 * qJD(5) - t181 * t640 + (t211 * t255 - t234 * t332 + t254 * t865 - t331 * t879) * t756 + ((t130 / 0.2e1 + t151 / 0.2e1) * t520 + (t805 + t150 / 0.2e1 - t629) * t517) * t641; 0, 0, t304 * qJD(5) + t119 * qJD(6) + (t133 * t821 + t184 * t823 + t275 * t824 + t336 * t825) * t826, 0, t304 * qJD(3), t119 * qJD(3) + t205 * t756; t4 * qJD(3) + t32 * qJD(4) - t7 * qJD(5) + t1 * qJD(6) + (-t792 / 0.4e1 - t790 / 0.4e1 - t782 / 0.4e1 - t775 / 0.4e1) * t827 - t533 * t642 - t532 * t643 + (-(t485 + t484) * t530 / 0.2e1 - t866) * qJD(1), qJD(5) * t305 - qJD(6) * t118, t4 * qJD(1) + t38 * qJD(5) + t10 * qJD(6) + (t47 - (t350 * t441 - t352 * t442) * t520 + t842 * t514 + (t441 * t351 - t442 * t353 + t547 * t517 + (t831 - t841) * t520 + t829) * t517) * t838 + (m(7) * (t111 * t133 - t244 * t245 - t246 * t247) + m(6) * (t144 * t184 - t318 * t319 - t320 * t321) + m(5) * ((t517 * t565 + t520 * (rSges(5,1) * t692 + t613) + t665) * t275 - t850 * t391 - t849 * t393) + m(4) * (t487 * t489 * t644 + (t517 * (rSges(4,1) * t698 - t645) + t520 * (rSges(4,1) * t684 + t517 * rSges(4,3) - t495)) * t336) + (t48 + (-t351 * t443 - t353 * t444) * t517 + t841 * t513 + (t350 * t443 + t352 * t444 + t831 * t520 + (t547 - t842) * t517 + t829) * t520) * t801) * qJD(3), t887, t38 * qJD(3) + (-0.4e1 * t848 + 0.2e1 * t633 * (-t445 * t519 + t648)) * qJD(5) + t20 * qJD(6) + t888, t1 * qJD(1) - t638 + t10 * qJD(3) + t20 * qJD(5) + (t34 * t798 + t35 * t801) * qJD(6) + (t23 / 0.2e1 + (t805 + t819) * t520 + (-t130 / 0.2e1 + t17 / 0.2e1) * t517) * t640 + (-t639 / 0.2e1 + (t111 * t205 + t190 * t215 - t244 * t255 - t246 * t254 + t439 * t846 - t41) * qJD(6)) * m(7) + (-t73 / 0.2e1 + (t50 / 0.2e1 - t14 / 0.2e1) * t520 + (t49 / 0.2e1 - t13 / 0.2e1) * t517) * t641; -t32 * qJD(3) + t249 * qJD(5) + t51 * qJD(6) + (t773 / 0.4e1 - t780 / 0.4e1 - t789 / 0.4e1) * t827, 0 ((-t245 * t520 + t247 * t517) * t821 + (-t319 * t520 + t321 * t517) * t823 + (-t391 * t520 + t393 * t517) * t824) * t826 - t887, 0, t636, t51 * qJD(1) + t103 + (t254 * t517 - t255 * t520) * t756; t7 * qJD(3) - t249 * qJD(4) + t42 * qJD(6) + (-t774 / 0.4e1 - t781 / 0.4e1) * t827 + 0.2e1 * (t199 * t821 + t231 * t823 + (-t211 * t821 - t243 * t823) * t703) * qJD(1), -t305 * qJD(3) (m(7) * (-t133 * t519 + t676) + m(6) * (-t184 * t519 + t674) + 0.2e1 * ((t245 * t517 + t247 * t520 + t111) * t821 + (t319 * t517 + t321 * t520 + t144) * t823) * t516 - t38) * qJD(3) + t262 * qJD(5) + t19 * qJD(6) - t888, -t636, t262 * qJD(3), t42 * qJD(1) + t19 * qJD(3) + (-t205 * t519 + (t254 * t520 + t255 * t517) * t516) * t756; t430 * t642 / 0.2e1 + t2 * qJD(3) + t53 * qJD(4) + t44 * qJD(5) + t6 * qJD(6) - qJD(1) * t886 - t830 * t643, qJD(3) * t118, t2 * qJD(1) + t638 + (t14 * t801 + (t167 * t520 + t168 * t517) * t803 + (t719 - t720) * t799 + t702 * t818 + t47 * t705 / 0.2e1 + t692 * t817 + t48 * t703 / 0.2e1 + t13 * t798 - t563) * qJD(3) + t21 * qJD(5) + t5 * qJD(6) + ((t111 * t146 + t133 * t190 - t176 * t246 - t177 * t244 - t234 * t245 + t247 * t879 - t54) * qJD(3) + t639 / 0.2e1) * m(7), qJD(1) * t53 + t103, qJD(1) * t44 + qJD(3) * t21, t6 * qJD(1) + t5 * qJD(3) + (t519 ^ 2 * t181 / 0.2e1 + m(7) * (t190 * t205 - t234 * t255 + t254 * t879) + (t35 * t796 + t34 * t801 + (t129 * t517 + t130 * t520) * t799) * t516) * qJD(6);];
Cq  = t9;
