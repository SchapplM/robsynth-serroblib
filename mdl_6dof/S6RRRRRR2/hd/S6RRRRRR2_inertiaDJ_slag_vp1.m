% Calculate time derivative of joint inertia matrix for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:18
% EndTime: 2019-03-10 03:34:06
% DurationCPUTime: 27.95s
% Computational Cost: add. (116145->1243), mult. (92044->1650), div. (0->0), fcn. (87670->12), ass. (0->647)
t519 = qJ(2) + qJ(3);
t507 = qJ(4) + t519;
t496 = cos(t507);
t518 = qJ(5) + qJ(6);
t505 = cos(t518);
t525 = cos(qJ(1));
t756 = t505 * t525;
t503 = sin(t518);
t522 = sin(qJ(1));
t760 = t503 * t522;
t423 = -t496 * t760 - t756;
t757 = t505 * t522;
t759 = t503 * t525;
t424 = t496 * t757 - t759;
t495 = sin(t507);
t767 = t495 * t522;
t295 = Icges(7,5) * t424 + Icges(7,6) * t423 + Icges(7,3) * t767;
t297 = Icges(7,4) * t424 + Icges(7,2) * t423 + Icges(7,6) * t767;
t299 = Icges(7,1) * t424 + Icges(7,4) * t423 + Icges(7,5) * t767;
t601 = -t297 * t503 + t299 * t505;
t140 = -t295 * t496 + t495 * t601;
t613 = Icges(7,5) * t505 - Icges(7,6) * t503;
t372 = -Icges(7,3) * t496 + t495 * t613;
t795 = Icges(7,4) * t505;
t618 = -Icges(7,2) * t503 + t795;
t373 = -Icges(7,6) * t496 + t495 * t618;
t796 = Icges(7,4) * t503;
t624 = Icges(7,1) * t505 - t796;
t374 = -Icges(7,5) * t496 + t495 * t624;
t182 = t372 * t767 + t373 * t423 + t374 * t424;
t866 = -t182 - t140;
t425 = -t496 * t759 + t757;
t426 = t496 * t756 + t760;
t766 = t495 * t525;
t296 = Icges(7,5) * t426 + Icges(7,6) * t425 + Icges(7,3) * t766;
t298 = Icges(7,4) * t426 + Icges(7,2) * t425 + Icges(7,6) * t766;
t300 = Icges(7,1) * t426 + Icges(7,4) * t425 + Icges(7,5) * t766;
t600 = -t298 * t503 + t300 * t505;
t141 = -t296 * t496 + t495 * t600;
t183 = t372 * t766 + t373 * t425 + t374 * t426;
t865 = -t183 - t141;
t509 = t522 * rSges(4,3);
t504 = sin(t519);
t506 = cos(t519);
t815 = rSges(4,1) * t506;
t637 = -rSges(4,2) * t504 + t815;
t422 = t525 * t637 + t509;
t720 = qJD(1) * t522;
t681 = t495 * t720;
t513 = qJD(5) + qJD(6);
t655 = -t496 * t513 + qJD(1);
t585 = t525 * t655;
t721 = qJD(1) * t496;
t654 = -t513 + t721;
t514 = qJD(2) + qJD(3);
t502 = qJD(4) + t514;
t761 = t502 * t525;
t700 = t495 * t761;
t841 = t522 * t654 + t700;
t260 = t503 * t841 + t505 * t585;
t261 = t503 * t585 - t505 * t841;
t698 = t496 * t761;
t693 = rSges(7,1) * t261 + rSges(7,2) * t260 + rSges(7,3) * t698;
t170 = -rSges(7,3) * t681 + t693;
t460 = pkin(10) * t698;
t526 = -pkin(11) - pkin(10);
t520 = sin(qJ(5));
t715 = qJD(5) * t520;
t707 = pkin(5) * t715;
t572 = -t502 * t526 - t707;
t523 = cos(qJ(5));
t714 = qJD(5) * t523;
t706 = pkin(5) * t714;
t822 = pkin(5) * t520;
t708 = qJD(1) * t822;
t682 = t522 * t706 + t525 * t708 + t526 * t681;
t497 = pkin(5) * t523 + pkin(4);
t819 = pkin(4) - t497;
t863 = t170 - t460 + (pkin(10) * t720 + t761 * t819) * t495 + (t525 * t572 + t720 * t819) * t496 + t682;
t631 = -rSges(7,1) * t424 - rSges(7,2) * t423;
t303 = rSges(7,3) * t767 - t631;
t817 = pkin(10) + t526;
t549 = -t495 * t817 - t496 * t819;
t751 = t520 * t525;
t320 = -pkin(5) * t751 + t522 * t549;
t862 = t303 + t320;
t765 = t496 * t502;
t674 = t765 / 0.2e1;
t861 = -t525 * t674 + t681 / 0.2e1;
t719 = qJD(1) * t525;
t670 = t719 / 0.2e1;
t860 = -t495 * t670 - t522 * t674;
t763 = t502 * t522;
t699 = t496 * t763;
t553 = t495 * t719 + t699;
t241 = t618 * t765 + (Icges(7,6) * t502 + (-Icges(7,2) * t505 - t796) * t513) * t495;
t859 = -t373 * t505 * t513 + (-t374 * t513 - t241) * t503;
t527 = -pkin(8) - pkin(7);
t521 = sin(qJ(2));
t717 = qJD(2) * t521;
t710 = pkin(2) * t717;
t858 = qJD(1) * t527 + t710;
t454 = rSges(5,1) * t495 + rSges(5,2) * t496;
t570 = t454 * t502;
t524 = cos(qJ(2));
t483 = rSges(3,1) * t521 + rSges(3,2) * t524;
t568 = qJD(2) * t483;
t857 = t522 * t568;
t803 = Icges(3,4) * t524;
t623 = -Icges(3,2) * t521 + t803;
t440 = Icges(3,6) * t522 + t525 * t623;
t804 = Icges(3,4) * t521;
t629 = Icges(3,1) * t524 - t804;
t442 = Icges(3,5) * t522 + t525 * t629;
t589 = t440 * t521 - t442 * t524;
t856 = t522 * t589;
t801 = Icges(4,4) * t506;
t621 = -Icges(4,2) * t504 + t801;
t418 = Icges(4,6) * t522 + t525 * t621;
t802 = Icges(4,4) * t504;
t627 = Icges(4,1) * t506 - t802;
t420 = Icges(4,5) * t522 + t525 * t627;
t591 = t418 * t504 - t420 * t506;
t855 = t522 * t591;
t799 = Icges(5,4) * t496;
t620 = -Icges(5,2) * t495 + t799;
t396 = Icges(5,6) * t522 + t525 * t620;
t800 = Icges(5,4) * t495;
t626 = Icges(5,1) * t496 - t800;
t398 = Icges(5,5) * t522 + t525 * t626;
t593 = t396 * t495 - t398 * t496;
t854 = t522 * t593;
t498 = pkin(2) * t524 + pkin(1);
t820 = pkin(1) - t498;
t853 = t522 * t820;
t439 = -Icges(3,6) * t525 + t522 * t623;
t441 = -Icges(3,5) * t525 + t522 * t629;
t590 = t439 * t521 - t441 * t524;
t852 = t525 * t590;
t417 = -Icges(4,6) * t525 + t522 * t621;
t419 = -Icges(4,5) * t525 + t522 * t627;
t592 = t417 * t504 - t419 * t506;
t851 = t525 * t592;
t395 = -Icges(5,6) * t525 + t522 * t620;
t397 = -Icges(5,5) * t525 + t522 * t626;
t594 = t395 * t495 - t397 * t496;
t850 = t525 * t594;
t240 = t613 * t765 + (Icges(7,3) * t502 + (-Icges(7,5) * t503 - Icges(7,6) * t505) * t513) * t495;
t784 = t373 * t503;
t849 = -t502 * t784 - t240;
t614 = Icges(6,5) * t523 - Icges(6,6) * t520;
t269 = t614 * t765 + (Icges(6,3) * t502 + (-Icges(6,5) * t520 - Icges(6,6) * t523) * qJD(5)) * t495;
t797 = Icges(6,4) * t523;
t619 = -Icges(6,2) * t520 + t797;
t378 = -Icges(6,6) * t496 + t495 * t619;
t781 = t378 * t520;
t848 = -t502 * t781 - t269;
t304 = rSges(7,1) * t426 + rSges(7,2) * t425 + rSges(7,3) * t766;
t847 = -t303 * t522 - t304 * t525;
t749 = t523 * t525;
t752 = t520 * t522;
t446 = -t496 * t752 - t749;
t750 = t522 * t523;
t447 = t496 * t750 - t751;
t634 = -rSges(6,1) * t447 - rSges(6,2) * t446;
t339 = rSges(6,3) * t767 - t634;
t448 = -t496 * t751 + t750;
t449 = t496 * t749 + t752;
t340 = rSges(6,1) * t449 + rSges(6,2) * t448 + rSges(6,3) * t766;
t846 = -t339 * t522 - t340 * t525;
t615 = Icges(5,5) * t496 - Icges(5,6) * t495;
t393 = -Icges(5,3) * t525 + t522 * t615;
t845 = qJD(1) * t393;
t616 = Icges(4,5) * t506 - Icges(4,6) * t504;
t415 = -Icges(4,3) * t525 + t522 * t616;
t844 = qJD(1) * t415;
t617 = Icges(3,5) * t524 - Icges(3,6) * t521;
t437 = -Icges(3,3) * t525 + t522 * t617;
t651 = -qJD(5) + t721;
t842 = t522 * t651 + t700;
t701 = t495 * t763;
t840 = t525 * t654 - t701;
t451 = Icges(5,2) * t496 + t800;
t452 = Icges(5,1) * t495 + t799;
t588 = t451 * t495 - t452 * t496;
t839 = qJD(1) * t588 + t502 * t615;
t463 = Icges(4,2) * t506 + t802;
t464 = Icges(4,1) * t504 + t801;
t587 = t463 * t504 - t464 * t506;
t838 = qJD(1) * t587 + t514 * t616;
t515 = -pkin(9) + t527;
t726 = t515 - t527;
t471 = pkin(3) * t506 + t498;
t731 = t471 - t498;
t365 = t522 * t731 + t525 * t726;
t837 = 2 * m(3);
t836 = 2 * m(4);
t835 = 2 * m(5);
t834 = 2 * m(6);
t833 = 2 * m(7);
t516 = t522 ^ 2;
t517 = t525 ^ 2;
t832 = -t496 / 0.2e1;
t831 = t522 / 0.2e1;
t830 = -t525 / 0.2e1;
t829 = -rSges(6,3) - pkin(10);
t828 = m(3) * t483;
t811 = rSges(4,2) * t506;
t465 = rSges(4,1) * t504 + t811;
t827 = m(4) * t465;
t826 = m(5) * t454;
t825 = pkin(2) * t521;
t824 = pkin(3) * t504;
t823 = pkin(4) * t496;
t821 = t522 * pkin(7);
t512 = t525 * pkin(7);
t818 = -pkin(7) - t527;
t816 = rSges(3,1) * t524;
t814 = rSges(5,1) * t496;
t813 = rSges(3,2) * t521;
t810 = rSges(3,3) * t525;
t586 = t522 * t655;
t262 = -t503 * t840 + t505 * t586;
t263 = t503 * t586 + t505 * t840;
t165 = Icges(7,5) * t263 + Icges(7,6) * t262 + Icges(7,3) * t553;
t167 = Icges(7,4) * t263 + Icges(7,2) * t262 + Icges(7,6) * t553;
t169 = Icges(7,1) * t263 + Icges(7,4) * t262 + Icges(7,5) * t553;
t39 = (t502 * t601 - t165) * t496 + (t295 * t502 + (-t297 * t513 + t169) * t505 + (-t299 * t513 - t167) * t503) * t495;
t809 = t39 * t525;
t552 = -t681 + t698;
t164 = Icges(7,5) * t261 + Icges(7,6) * t260 + Icges(7,3) * t552;
t166 = Icges(7,4) * t261 + Icges(7,2) * t260 + Icges(7,6) * t552;
t168 = Icges(7,1) * t261 + Icges(7,4) * t260 + Icges(7,5) * t552;
t40 = (t502 * t600 - t164) * t496 + (t296 * t502 + (-t298 * t513 + t168) * t505 + (-t300 * t513 - t166) * t503) * t495;
t808 = t40 * t522;
t652 = -qJD(5) * t496 + qJD(1);
t582 = t523 * t652;
t308 = t522 * t582 + (-t525 * t651 + t701) * t520;
t581 = t652 * t520;
t762 = t502 * t523;
t309 = t651 * t749 + (-t495 * t762 + t581) * t522;
t193 = Icges(6,5) * t309 + Icges(6,6) * t308 + Icges(6,3) * t553;
t195 = Icges(6,4) * t309 + Icges(6,2) * t308 + Icges(6,6) * t553;
t197 = Icges(6,1) * t309 + Icges(6,4) * t308 + Icges(6,5) * t553;
t329 = Icges(6,5) * t447 + Icges(6,6) * t446 + Icges(6,3) * t767;
t331 = Icges(6,4) * t447 + Icges(6,2) * t446 + Icges(6,6) * t767;
t333 = Icges(6,1) * t447 + Icges(6,4) * t446 + Icges(6,5) * t767;
t599 = -t331 * t520 + t333 * t523;
t52 = (t502 * t599 - t193) * t496 + (-t195 * t520 + t197 * t523 + t329 * t502 + (-t331 * t523 - t333 * t520) * qJD(5)) * t495;
t807 = t52 * t525;
t510 = t522 * rSges(3,3);
t508 = t522 * rSges(5,3);
t306 = t520 * t842 + t525 * t582;
t307 = -t523 * t842 + t525 * t581;
t192 = Icges(6,5) * t307 + Icges(6,6) * t306 + Icges(6,3) * t552;
t194 = Icges(6,4) * t307 + Icges(6,2) * t306 + Icges(6,6) * t552;
t196 = Icges(6,1) * t307 + Icges(6,4) * t306 + Icges(6,5) * t552;
t330 = Icges(6,5) * t449 + Icges(6,6) * t448 + Icges(6,3) * t766;
t332 = Icges(6,4) * t449 + Icges(6,2) * t448 + Icges(6,6) * t766;
t334 = Icges(6,1) * t449 + Icges(6,4) * t448 + Icges(6,5) * t766;
t598 = -t332 * t520 + t334 * t523;
t53 = (t502 * t598 - t192) * t496 + (-t194 * t520 + t196 * t523 + t330 * t502 + (-t332 * t523 - t334 * t520) * qJD(5)) * t495;
t806 = t53 * t522;
t805 = -rSges(7,3) + t526;
t798 = Icges(6,4) * t520;
t270 = t619 * t765 + (Icges(6,6) * t502 + (-Icges(6,2) * t523 - t798) * qJD(5)) * t495;
t785 = t270 * t520;
t783 = t374 * t505;
t636 = -rSges(5,2) * t495 + t814;
t402 = t636 * t502;
t780 = t402 * t522;
t779 = t402 * t525;
t433 = t637 * t514;
t778 = t433 * t522;
t777 = t439 * t524;
t776 = t440 * t524;
t775 = t441 * t521;
t774 = t442 * t521;
t773 = t451 * t502;
t772 = t452 * t502;
t771 = t463 * t514;
t770 = t464 * t514;
t769 = t495 * t502;
t764 = t496 * t525;
t758 = t504 * t514;
t755 = t506 * t514;
t754 = t514 * t522;
t753 = t514 * t525;
t748 = t525 * t515;
t632 = rSges(7,1) * t263 + rSges(7,2) * t262;
t171 = rSges(7,3) * t553 + t632;
t747 = t171 * t766 + t303 * t698;
t630 = rSges(7,1) * t505 - rSges(7,2) * t503;
t243 = t630 * t765 + (rSges(7,3) * t502 + (-rSges(7,1) * t503 - rSges(7,2) * t505) * t513) * t495;
t677 = t495 * t715;
t268 = -pkin(5) * t677 + t502 * t549;
t745 = -t243 - t268;
t375 = -rSges(7,3) * t496 + t495 * t630;
t367 = t375 * t720;
t744 = t304 * t769 + t367 * t495;
t633 = rSges(6,1) * t523 - rSges(6,2) * t520;
t276 = t633 * t765 + (rSges(6,3) * t502 + (-rSges(6,1) * t520 - rSges(6,2) * t523) * qJD(5)) * t495;
t640 = pkin(10) * t495 + t823;
t406 = t640 * t502;
t743 = -t276 - t406;
t479 = pkin(4) * t764;
t428 = pkin(10) * t766 + t479;
t494 = pkin(5) * t752;
t577 = t497 * t764 - t526 * t766 + t494;
t321 = t577 - t428;
t741 = -t304 - t321;
t740 = -t340 - t428;
t216 = t303 * t496 + t375 * t767;
t456 = t525 * t471;
t484 = t525 * t498;
t366 = -t522 * t726 + t456 - t484;
t739 = t365 * t522 + t366 * t525;
t404 = rSges(5,1) * t764 - rSges(5,2) * t766 + t508;
t738 = -t366 - t404;
t368 = -t495 * t819 + t496 * t817;
t737 = -t368 - t375;
t383 = -rSges(6,3) * t496 + t495 * t633;
t371 = t383 * t720;
t455 = pkin(4) * t495 - pkin(10) * t496;
t434 = t455 * t720;
t736 = t371 + t434;
t403 = -rSges(5,3) * t525 + t522 * t636;
t310 = t403 * t522 + t404 * t525;
t735 = -t383 - t455;
t409 = t525 * t527 + t512 - t853;
t410 = -pkin(1) * t525 + t522 * t818 + t484;
t734 = t409 * t522 + t410 * t525;
t421 = -rSges(4,3) * t525 + t522 * t637;
t322 = t421 * t522 + t422 * t525;
t427 = t640 * t522;
t733 = t427 * t522 + t428 * t525;
t679 = t504 * t720;
t474 = pkin(3) * t679;
t732 = t454 * t720 + t474;
t730 = rSges(5,2) * t681 + rSges(5,3) * t719;
t729 = rSges(4,2) * t679 + rSges(4,3) * t719;
t728 = t858 * t522;
t727 = t525 * t816 + t510;
t725 = t516 + t517;
t394 = Icges(5,3) * t522 + t525 * t615;
t724 = qJD(1) * t394;
t416 = Icges(4,3) * t522 + t525 * t616;
t723 = qJD(1) * t416;
t438 = Icges(3,3) * t522 + t525 * t617;
t722 = qJD(1) * t438;
t716 = qJD(2) * t524;
t32 = t165 * t767 + t167 * t423 + t169 * t424 + t262 * t297 + t263 * t299 + t295 * t553;
t33 = t164 * t767 + t166 * t423 + t168 * t424 + t262 * t298 + t263 * t300 + t296 * t553;
t124 = t295 * t767 + t297 * t423 + t299 * t424;
t125 = t296 * t767 + t298 * t423 + t300 * t424;
t611 = t124 * t522 + t125 * t525;
t16 = qJD(1) * t611 - t32 * t525 + t33 * t522;
t45 = t193 * t767 + t195 * t446 + t197 * t447 + t308 * t331 + t309 * t333 + t329 * t553;
t46 = t192 * t767 + t194 * t446 + t196 * t447 + t308 * t332 + t309 * t334 + t330 * t553;
t142 = t329 * t767 + t331 * t446 + t333 * t447;
t143 = t330 * t767 + t332 * t446 + t334 * t447;
t605 = t142 * t522 + t143 * t525;
t22 = qJD(1) * t605 - t45 * t525 + t46 * t522;
t222 = -t393 * t525 - t522 * t594;
t223 = -t394 * t525 - t854;
t450 = Icges(5,5) * t495 + Icges(5,6) * t496;
t564 = t502 * t450;
t285 = -t525 * t564 - t845;
t286 = -t522 * t564 + t724;
t287 = -qJD(1) * t395 - t525 * t773;
t289 = -qJD(1) * t397 - t525 * t772;
t290 = qJD(1) * t398 - t522 * t772;
t662 = t395 * t502 - t290;
t288 = qJD(1) * t396 - t522 * t773;
t664 = t397 * t502 + t288;
t29 = (t525 * t286 + (t223 + t850) * qJD(1)) * t525 + (t222 * qJD(1) + (-t287 * t495 + t289 * t496 - t396 * t765 - t398 * t769 + t724) * t522 + (-t285 + t662 * t496 + t664 * t495 + (-t393 - t593) * qJD(1)) * t525) * t522;
t713 = -t16 - t22 - t29;
t712 = pkin(3) * t755;
t711 = t525 * t813;
t709 = pkin(2) * t716;
t702 = t497 * t769;
t242 = t624 * t765 + (Icges(7,5) * t502 + (-Icges(7,1) * t503 - t795) * t513) * t495;
t697 = t242 * t495 * t505 + t372 * t769 + t765 * t783;
t696 = -t406 + t745;
t625 = Icges(6,1) * t523 - t798;
t271 = t625 * t765 + (Icges(6,5) * t502 + (-Icges(6,1) * t520 - t797) * qJD(5)) * t495;
t377 = -Icges(6,3) * t496 + t495 * t614;
t379 = -Icges(6,5) * t496 + t495 * t625;
t695 = t271 * t495 * t523 + t379 * t496 * t762 + t377 * t769;
t453 = -pkin(3) * t758 - t710;
t443 = t525 * t453;
t491 = t515 * t720;
t653 = t525 * t710;
t694 = t522 * (t453 * t522 + t719 * t731 - t491 + t728) + t525 * (-qJD(1) * t365 + t443 + t653) + t365 * t719;
t554 = -t496 * t720 - t700;
t692 = t522 * (-t522 * t570 + (t525 * t636 + t508) * qJD(1)) + t525 * (rSges(5,1) * t554 - rSges(5,2) * t698 + t730) + t403 * t719;
t691 = rSges(6,1) * t307 + rSges(6,2) * t306 + rSges(6,3) * t698;
t690 = -t428 + t741;
t459 = pkin(4) * t701;
t689 = t522 * (pkin(10) * t553 + qJD(1) * t479 - t459) + t525 * (pkin(4) * t554 - pkin(10) * t681 + t460) + t427 * t719;
t569 = t465 * t514;
t688 = t522 * (qJD(1) * t422 - t522 * t569) + t525 * (-t753 * t811 + (-t504 * t753 - t506 * t720) * rSges(4,1) + t729) + t421 * t719;
t687 = -t366 + t740;
t351 = t368 * t720;
t686 = t351 + t367 + t434;
t685 = t522 * ((-t525 * t820 - t821) * qJD(1) - t728) + t525 * (-t653 + (t525 * t818 + t853) * qJD(1)) + t409 * t719;
t684 = -t455 + t737;
t683 = t474 + t736;
t678 = t521 * t720;
t676 = t767 / 0.2e1;
t675 = t766 / 0.2e1;
t160 = -t329 * t496 + t495 * t599;
t205 = t377 * t767 + t378 * t446 + t379 * t447;
t673 = t160 / 0.2e1 + t205 / 0.2e1;
t161 = -t330 * t496 + t495 * t598;
t206 = t377 * t766 + t378 * t448 + t379 * t449;
t672 = t161 / 0.2e1 + t206 / 0.2e1;
t671 = t720 / 0.2e1;
t669 = -t465 - t825;
t668 = -t454 - t824;
t667 = -t455 - t824;
t666 = t737 * t525;
t338 = t735 * t525;
t665 = t398 * t502 + t287;
t663 = -t396 * t502 + t289;
t325 = -qJD(1) * t417 - t525 * t771;
t661 = t420 * t514 + t325;
t326 = qJD(1) * t418 - t522 * t771;
t660 = t419 * t514 + t326;
t327 = -qJD(1) * t419 - t525 * t770;
t659 = -t418 * t514 + t327;
t328 = qJD(1) * t420 - t522 * t770;
t658 = t417 * t514 - t328;
t657 = -t515 * t522 + t456;
t656 = -t496 * t497 - t471;
t650 = t496 * t171 + t243 * t767 + t375 * t553;
t649 = -t366 + t690;
t648 = t474 + t686;
t202 = t310 + t739;
t181 = t733 - t846;
t643 = -t406 - t712;
t642 = -t383 + t667;
t229 = t684 * t525;
t641 = -t824 - t825;
t207 = -t372 * t496 + (t783 - t784) * t495;
t30 = t165 * t766 + t167 * t425 + t169 * t426 + t260 * t297 + t261 * t299 + t295 * t552;
t31 = t164 * t766 + t166 * t425 + t168 * t426 + t260 * t298 + t261 * t300 + t296 * t552;
t126 = t295 * t766 + t297 * t425 + t299 * t426;
t127 = t296 * t766 + t298 * t425 + t300 * t426;
t609 = t126 * t522 + t127 * t525;
t610 = t126 * t525 - t127 * t522;
t66 = t240 * t766 + t241 * t425 + t242 * t426 + t260 * t373 + t261 * t374 + t372 * t552;
t5 = (t502 * t609 - t66) * t496 + (qJD(1) * t610 + t183 * t502 + t30 * t522 + t31 * t525) * t495;
t612 = t124 * t525 - t125 * t522;
t67 = t240 * t767 + t241 * t423 + t242 * t424 + t262 * t373 + t263 * t374 + t372 * t553;
t6 = (t502 * t611 - t67) * t496 + (qJD(1) * t612 + t182 * t502 + t32 * t522 + t33 * t525) * t495;
t60 = -t182 * t496 + t495 * t611;
t607 = t140 * t522 + t141 * t525;
t61 = -t183 * t496 + t495 * t609;
t639 = t5 * t766 + t6 * t767 + t61 * t698 + (-t207 * t496 + t495 * t607) * t769 + t553 * t60;
t638 = -t813 + t816;
t635 = rSges(6,1) * t309 + rSges(6,2) * t308;
t628 = Icges(3,1) * t521 + t803;
t622 = Icges(3,2) * t524 + t804;
t462 = Icges(4,5) * t504 + Icges(4,6) * t506;
t608 = t140 * t525 - t141 * t522;
t606 = t142 * t525 - t143 * t522;
t144 = t329 * t766 + t331 * t448 + t333 * t449;
t145 = t330 * t766 + t332 * t448 + t334 * t449;
t604 = t144 * t525 - t145 * t522;
t603 = t144 * t522 + t145 * t525;
t602 = t160 * t522 + t161 * t525;
t597 = t339 * t525 - t340 * t522;
t584 = -t276 + t643;
t583 = t667 + t737;
t580 = -pkin(1) - t638;
t579 = -t454 + t641;
t578 = -t455 + t641;
t292 = t642 * t525;
t576 = -t498 - t637;
t575 = -t471 - t636;
t198 = -rSges(6,3) * t681 + t691;
t199 = rSges(6,3) * t553 + t635;
t574 = t198 * t525 + t199 * t522 + t339 * t719 + t689;
t573 = t692 + t694;
t121 = t181 + t739;
t115 = t320 * t522 + t321 * t525 + t733 - t847;
t571 = t643 + t745;
t567 = -t383 + t578;
t561 = t514 * t462;
t560 = t320 * t525 + t522 * t741;
t559 = -t709 - t712;
t558 = qJD(2) * t628;
t557 = qJD(2) * t622;
t556 = qJD(2) * (-Icges(3,5) * t521 - Icges(3,6) * t524);
t219 = t583 * t525;
t555 = t495 * t829 - t471 - t823;
t370 = t579 * t525;
t551 = t495 * t805 + t656;
t550 = t578 + t737;
t548 = -t402 + t559;
t547 = -t406 + t559;
t106 = t115 + t739;
t273 = t567 * t525;
t15 = qJD(1) * t609 - t30 * t525 + t31 * t522;
t43 = t193 * t766 + t195 * t448 + t197 * t449 + t306 * t331 + t307 * t333 + t329 * t552;
t44 = t192 * t766 + t194 * t448 + t196 * t449 + t306 * t332 + t307 * t334 + t330 * t552;
t21 = qJD(1) * t603 - t43 * t525 + t44 * t522;
t224 = t393 * t522 - t850;
t225 = t394 * t522 - t525 * t593;
t546 = (-t222 * t525 - t606 - t612) * t720 + (-t224 * t525 - t604 - t610) * t719 + (t15 + t21 + t223 * t720 + t225 * t719 + (t225 * qJD(1) + (t288 * t495 - t290 * t496 + t395 * t765 + t397 * t769 - t845) * t525) * t525 + ((t224 + t854) * qJD(1) + (-t286 + t663 * t496 - t665 * t495 + (t394 - t594) * qJD(1)) * t525 + t522 * t285) * t522) * t522;
t545 = -t276 + t547;
t191 = -t525 * t706 + t459 + (-t702 + (-t502 * t817 - t707) * t496) * t522 + (t525 * t549 + t494) * qJD(1);
t544 = t689 + t862 * t719 + t863 * t525 + (t171 + t191) * t522;
t543 = t574 + t694;
t214 = t550 * t525;
t542 = t547 + t745;
t541 = t15 * t675 + t16 * t676 + t5 * t831 + t6 * t830 + (qJD(1) * t607 + t808 - t809) * t832 + t60 * t671 + t61 * t670 - t608 * t769 / 0.2e1 + t861 * t610 + t860 * t612;
t540 = rSges(3,2) * t678 + rSges(3,3) * t719 - t525 * t568;
t232 = -t415 * t525 - t522 * t592;
t233 = -t416 * t525 - t855;
t234 = t415 * t522 - t851;
t235 = t416 * t522 - t525 * t591;
t323 = -t525 * t561 - t844;
t324 = -t522 * t561 + t723;
t539 = (-t232 * t525 + t233 * t522) * t720 + (-t234 * t525 + t235 * t522) * t719 + t522 * ((t522 * t323 + (t234 + t855) * qJD(1)) * t522 + (t235 * qJD(1) + (t326 * t504 - t328 * t506 + t417 * t755 + t419 * t758 - t844) * t525 + (-t324 + t659 * t506 - t661 * t504 + (t416 - t592) * qJD(1)) * t522) * t525) + t546;
t538 = t522 * t555 - t748;
t203 = t207 * t769;
t80 = t495 * t859 + t496 * t849 + t697;
t7 = t203 + (t502 * t607 - t80) * t496 + (qJD(1) * t608 + t39 * t522 + t40 * t525) * t495;
t537 = -t496 * t7 - t61 * t681 + t639;
t536 = t544 + t694;
t535 = t203 + (t39 + t67) * t676 + (t40 + t66) * t675 + t865 * t861 + t866 * t860;
t534 = t525 * t713 + t546;
t400 = t620 * t502;
t401 = t626 * t502;
t533 = qJD(1) * t450 + (t401 - t773) * t496 + (-t400 - t772) * t495;
t430 = t621 * t514;
t431 = t627 * t514;
t532 = qJD(1) * t462 + (t431 - t771) * t506 + (-t430 - t770) * t504;
t42 = (t525 * t324 + (t233 + t851) * qJD(1)) * t525 + (t232 * qJD(1) + (-t325 * t504 + t327 * t506 - t418 * t755 - t420 * t758 + t723) * t522 + (-t323 + t658 * t506 + t660 * t504 + (-t415 - t591) * qJD(1)) * t525) * t522;
t531 = t539 + (-t42 + t713) * t525;
t78 = t269 * t766 + t270 * t448 + t271 * t449 + t306 * t378 + t307 * t379 + t377 * t552;
t10 = (t502 * t603 - t78) * t496 + (qJD(1) * t604 + t206 * t502 + t43 * t522 + t44 * t525) * t495;
t79 = t269 * t767 + t270 * t446 + t271 * t447 + t308 * t378 + t309 * t379 + t377 * t553;
t11 = (t502 * t605 - t79) * t496 + (qJD(1) * t606 + t205 * t502 + t45 * t522 + t46 * t525) * t495;
t73 = -t205 * t496 + t495 * t605;
t74 = -t206 * t496 + t495 * t603;
t530 = t10 * t831 + t11 * t830 + t21 * t675 + t22 * t676 + (-t160 * t525 + t161 * t522) * t769 / 0.2e1 + (qJD(1) * t602 + t806 - t807) * t832 + t541 + t73 * t671 + t74 * t670 + t861 * t604 + t860 * t606;
t529 = -t809 / 0.2e1 + t808 / 0.2e1 - t807 / 0.2e1 + t806 / 0.2e1 + (t495 * t663 + t496 * t665 + t522 * t839 + t525 * t533 + t66 + t78) * t831 + (-t495 * t662 + t496 * t664 + t522 * t533 - t525 * t839 + t67 + t79) * t830 + (t395 * t496 + t397 * t495 - t450 * t525 - t522 * t588 + t160 + t205 - t866) * t671 + (t396 * t496 + t398 * t495 + t450 * t522 - t525 * t588 + t161 + t206 - t865) * t670;
t528 = t529 + (t504 * t659 + t506 * t661 + t522 * t838 + t525 * t532) * t831 + (-t504 * t658 + t506 * t660 + t522 * t532 - t525 * t838) * t830 + (t417 * t506 + t419 * t504 - t462 * t525 - t522 * t587) * t671 + (t418 * t506 + t420 * t504 + t462 * t522 - t525 * t587) * t670;
t489 = pkin(2) * t678;
t470 = t638 * qJD(2);
t445 = -t711 + t727;
t444 = t522 * t638 - t810;
t408 = t669 * t525;
t407 = t669 * t522;
t388 = t821 + (pkin(1) - t813) * t525 + t727;
t387 = t522 * t580 + t512 + t810;
t385 = t668 * t525;
t384 = t668 * t522;
t369 = t579 * t522;
t364 = -t522 * t527 + t422 + t484;
t363 = (rSges(4,3) - t527) * t525 + t576 * t522;
t355 = t522 * t556 + t722;
t354 = -qJD(1) * t437 + t525 * t556;
t343 = t404 + t657;
t342 = (rSges(5,3) - t515) * t525 + t575 * t522;
t337 = t735 * t522;
t336 = t857 + ((-rSges(3,3) - pkin(7)) * t522 + t580 * t525) * qJD(1);
t335 = (t512 + (-pkin(1) - t816) * t522) * qJD(1) + t540;
t302 = -t465 * t719 - t778 + (-t521 * t719 - t522 * t716) * pkin(2);
t301 = t465 * t720 + t489 + (-t433 - t709) * t525;
t291 = t642 * t522;
t272 = t567 * t522;
t267 = t303 * t766;
t259 = -t454 * t719 - t780 + (-t504 * t719 - t506 * t754) * pkin(3);
t258 = (-t402 - t712) * t525 + t732;
t252 = t438 * t522 - t525 * t589;
t251 = t437 * t522 - t852;
t250 = -t438 * t525 - t856;
t249 = -t437 * t525 - t522 * t590;
t245 = t465 * t754 + (t525 * t576 - t509) * qJD(1) + t728;
t244 = (-t498 - t815) * t720 + (-t569 - t858) * t525 + t729;
t238 = qJD(1) * t370 + t522 * t548;
t237 = t525 * t548 + t489 + t732;
t231 = t657 - t740;
t230 = t538 + t634;
t228 = t684 * t522;
t227 = -t340 * t496 - t383 * t766;
t226 = t339 * t496 + t383 * t767;
t221 = t491 + (-t453 + t570) * t522 + (t525 * t575 - t508) * qJD(1);
t220 = t443 - t525 * t570 + (-t748 + (-t471 - t814) * t522) * qJD(1) + t730;
t218 = t583 * t522;
t217 = -t304 * t496 - t375 * t766;
t215 = t322 + t734;
t213 = t550 * t522;
t212 = -t377 * t496 + (t379 * t523 - t781) * t495;
t211 = t577 + t657 + t304;
t210 = (-t515 + t822) * t525 + t551 * t522 + t631;
t209 = t597 * t495;
t208 = t212 * t769;
t204 = -t304 * t767 + t267;
t187 = qJD(1) * t338 + t522 * t743;
t186 = t525 * t743 + t736;
t176 = -t422 * t720 + t688;
t173 = qJD(1) * t292 + t522 * t584;
t172 = t525 * t584 + t683;
t153 = -t404 * t720 + t692;
t152 = qJD(1) * t273 + t522 * t545;
t151 = t525 * t545 + t489 + t683;
t148 = t202 + t734;
t135 = t495 * t666 + t496 * t741;
t134 = t320 * t496 + t368 * t767 + t216;
t120 = t459 + t491 + (t765 * t829 - t453) * t522 + t555 * t719 - t635;
t119 = -pkin(4) * t700 + qJD(1) * t538 + t443 + t460 + t691;
t118 = t495 * t560 + t267;
t117 = qJD(1) * t229 + t522 * t696;
t116 = t525 * t696 + t686;
t114 = t121 + t734;
t113 = qJD(1) * t219 + t522 * t571;
t112 = t525 * t571 + t648;
t111 = (-t410 - t422) * t720 + t685 + t688;
t110 = t491 + (qJD(1) * t551 + t706) * t525 + (-t708 + t702 - t453 + (t502 * t805 + t707) * t496) * t522 - t632;
t109 = t443 + (t496 * t572 - t702) * t525 + (-t748 + (-rSges(7,3) * t495 + t656) * t522) * qJD(1) + t682 + t693;
t108 = qJD(1) * t214 + t522 * t542;
t107 = t525 * t542 + t489 + t648;
t105 = (t383 * t763 + t199) * t496 + (t276 * t522 - t339 * t502 + t383 * t719) * t495;
t104 = (-t383 * t761 - t198) * t496 + (-t276 * t525 + t340 * t502 + t371) * t495;
t98 = t720 * t738 + t573;
t97 = t106 + t734;
t92 = -t303 * t769 + t650;
t91 = -t243 * t766 + (-t375 * t761 - t170) * t496 + t744;
t90 = t848 * t496 + (-t785 + (-t378 * t523 - t379 * t520) * qJD(5)) * t495 + t695;
t75 = (-t410 + t738) * t720 + t573 + t685;
t72 = t597 * t765 + (qJD(1) * t846 - t198 * t522 + t199 * t525) * t495;
t69 = t720 * t740 + t574;
t57 = -t304 * t699 + (qJD(1) * t847 - t170 * t522) * t495 + t747;
t49 = t687 * t720 + t543;
t48 = (t368 * t763 + t191) * t496 + (t268 * t522 + t368 * t719 - t502 * t862) * t495 + t650;
t47 = (t502 * t666 - t863) * t496 + (t321 * t502 + t525 * t745 + t351) * t495 + t744;
t34 = (-t410 + t687) * t720 + t543 + t685;
t27 = t690 * t720 + t544;
t26 = t560 * t765 + (t191 * t525 - t863 * t522 + (-t522 * t862 + t525 * t741) * qJD(1)) * t495 + t747;
t25 = t649 * t720 + t536;
t24 = (-t410 + t649) * t720 + t536 + t685;
t1 = [t504 * t431 + t506 * t430 + (t109 * t211 + t110 * t210) * t833 + (t119 * t231 + t120 * t230) * t834 + (t220 * t343 + t221 * t342) * t835 + (t244 * t364 + t245 * t363) * t836 + (t335 * t388 + t336 * t387) * t837 - t463 * t758 + t452 * t765 - t379 * t677 + t464 * t755 - t451 * t769 + t695 + t697 + (-t622 + t629) * t717 + (t623 + t628) * t716 + (t400 + t848 + t849) * t496 + (-t378 * t714 + t401 - t785 + t859) * t495; (t517 / 0.2e1 + t516 / 0.2e1) * t617 * qJD(2) + ((t776 / 0.2e1 + t774 / 0.2e1 - t388 * t828) * t525 + (t387 * t828 + t777 / 0.2e1 + t775 / 0.2e1) * t522) * qJD(1) + t528 + m(5) * (t220 * t369 + t221 * t370 + t237 * t342 + t238 * t343) + m(6) * (t119 * t272 + t120 * t273 + t151 * t230 + t152 * t231) + m(7) * (t107 * t210 + t108 * t211 + t109 * t213 + t110 * t214) + m(3) * ((-t335 * t522 - t336 * t525) * t483 + (-t387 * t525 - t388 * t522) * t470) + (-qJD(2) * t589 + (-qJD(1) * t439 - t525 * t557) * t524 + (-qJD(1) * t441 - t525 * t558) * t521) * t831 + (-qJD(2) * t590 + (qJD(1) * t440 - t522 * t557) * t524 + (qJD(1) * t442 - t522 * t558) * t521) * t830 + m(4) * (t244 * t407 + t245 * t408 + t301 * t363 + t302 * t364); -t525 * t16 - t525 * t22 - t525 * t29 - t525 * t42 + t539 + ((t444 * t522 + t445 * t525) * ((qJD(1) * t444 + t540) * t525 + (-t857 + (-t445 - t711 + t510) * qJD(1)) * t522) + t725 * t483 * t470) * t837 + (t111 * t215 + t301 * t408 + t302 * t407) * t836 + (t148 * t75 + t237 * t370 + t238 * t369) * t835 + (t114 * t34 + t151 * t273 + t152 * t272) * t834 + (t107 * t214 + t108 * t213 + t24 * t97) * t833 + (-t249 * t525 + t250 * t522) * t720 + (-t251 * t525 + t252 * t522) * t719 + t522 * ((t522 * t354 + (t251 + t856) * qJD(1)) * t522 + (t252 * qJD(1) + (t439 * t716 + t441 * t717) * t525 + (-t355 + (-t774 - t776) * qJD(2) + (t438 - t590) * qJD(1)) * t522) * t525) - t525 * ((t525 * t355 + (t250 + t852) * qJD(1)) * t525 + (t249 * qJD(1) + (-t440 * t716 - t442 * t717 + t722) * t522 + (-t354 + (t775 + t777) * qJD(2) - t589 * qJD(1)) * t525) * t522); m(5) * (t220 * t384 + t221 * t385 + t258 * t342 + t259 * t343) + m(6) * (t119 * t291 + t120 * t292 + t172 * t230 + t173 * t231) + m(7) * (t109 * t218 + t110 * t219 + t112 * t210 + t113 * t211) + t528 + (-t244 * t522 - t245 * t525 + (t363 * t522 - t364 * t525) * qJD(1)) * t827 + m(4) * (-t363 * t525 - t364 * t522) * t433; t531 + m(4) * (-t408 * t433 * t525 + t111 * t322 + t176 * t215 - t407 * t778) + m(7) * (t106 * t24 + t107 * t219 + t108 * t218 + t112 * t214 + t113 * t213 + t25 * t97) + m(6) * (t114 * t49 + t121 * t34 + t151 * t292 + t152 * t291 + t172 * t273 + t173 * t272) + m(5) * (t148 * t98 + t202 * t75 + t237 * t385 + t238 * t384 + t258 * t370 + t259 * t369) + (-t301 * t525 - t302 * t522 + (-t407 * t525 + t408 * t522) * qJD(1)) * t827; (t433 * t465 * t725 + t176 * t322) * t836 + t531 + (t106 * t25 + t112 * t219 + t113 * t218) * t833 + (t121 * t49 + t172 * t292 + t173 * t291) * t834 + (t202 * t98 + t258 * t385 + t259 * t384) * t835; t529 + m(6) * (t119 * t337 + t120 * t338 + t186 * t230 + t187 * t231) + m(7) * (t109 * t228 + t110 * t229 + t116 * t210 + t117 * t211) + m(5) * (-t342 * t525 - t343 * t522) * t402 + (-t220 * t522 - t221 * t525 + (t342 * t522 - t343 * t525) * qJD(1)) * t826; m(5) * (t148 * t153 + t310 * t75 - t369 * t780 - t370 * t779) + m(7) * (t107 * t229 + t108 * t228 + t115 * t24 + t116 * t214 + t117 * t213 + t27 * t97) + m(6) * (t114 * t69 + t151 * t338 + t152 * t337 + t181 * t34 + t186 * t273 + t187 * t272) + t534 + (-t237 * t525 - t238 * t522 + (-t369 * t525 + t370 * t522) * qJD(1)) * t826; t534 + (-t258 * t525 - t259 * t522 + (-t384 * t525 + t385 * t522) * qJD(1)) * t826 + m(7) * (t106 * t27 + t112 * t229 + t113 * t228 + t115 * t25 + t116 * t219 + t117 * t218) + m(6) * (t121 * t69 + t172 * t338 + t173 * t337 + t181 * t49 + t186 * t292 + t187 * t291) + m(5) * (t153 * t202 + t310 * t98 - t384 * t780 - t385 * t779); (t115 * t27 + t116 * t229 + t117 * t228) * t833 + (t181 * t69 + t186 * t338 + t187 * t337) * t834 + (t402 * t454 * t725 + t153 * t310) * t835 + t534; t535 + t208 + ((t53 / 0.2e1 + t78 / 0.2e1) * t525 + (t52 / 0.2e1 + t79 / 0.2e1) * t522 + (-t522 * t672 + t525 * t673) * qJD(1)) * t495 + (-t80 - t90 + (t522 * t673 + t525 * t672) * t502) * t496 + m(6) * (t104 * t231 + t105 * t230 + t119 * t227 + t120 * t226) + m(7) * (t109 * t135 + t110 * t134 + t210 * t48 + t211 * t47); t530 + m(6) * (t104 * t272 + t105 * t273 + t114 * t72 + t151 * t226 + t152 * t227 + t209 * t34) + m(7) * (t107 * t134 + t108 * t135 + t118 * t24 + t213 * t47 + t214 * t48 + t26 * t97); t530 + m(6) * (t104 * t291 + t105 * t292 + t121 * t72 + t172 * t226 + t173 * t227 + t209 * t49) + m(7) * (t106 * t26 + t112 * t134 + t113 * t135 + t118 * t25 + t218 * t47 + t219 * t48); t530 + m(6) * (t104 * t337 + t105 * t338 + t181 * t72 + t186 * t226 + t187 * t227 + t209 * t69) + m(7) * (t115 * t26 + t116 * t134 + t117 * t135 + t118 * t27 + t228 * t47 + t229 * t48); (t118 * t26 + t134 * t48 + t135 * t47) * t833 + (t104 * t227 + t105 * t226 + t209 * t72) * t834 + (t90 * t496 - t208 - t7 + (-t496 * t602 + t522 * t73 + t525 * t74) * t502) * t496 + (t525 * t10 + t522 * t11 + t602 * t769 + (-t212 * t502 - t52 * t522 - t53 * t525) * t496 + ((-t160 * t496 + t73) * t525 + (t161 * t496 - t61 - t74) * t522) * qJD(1)) * t495 + t639; -t80 * t496 + m(7) * (t109 * t217 + t110 * t216 + t210 * t92 + t211 * t91) + t535; t541 + m(7) * (t107 * t216 + t108 * t217 + t204 * t24 + t213 * t91 + t214 * t92 + t57 * t97); t541 + m(7) * (t106 * t57 + t112 * t216 + t113 * t217 + t204 * t25 + t218 * t91 + t219 * t92); t541 + m(7) * (t115 * t57 + t116 * t216 + t117 * t217 + t204 * t27 + t228 * t91 + t229 * t92); m(7) * (t118 * t57 + t134 * t92 + t135 * t91 + t204 * t26 + t216 * t48 + t217 * t47) + t537; (t204 * t57 + t216 * t92 + t217 * t91) * t833 + t537;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
