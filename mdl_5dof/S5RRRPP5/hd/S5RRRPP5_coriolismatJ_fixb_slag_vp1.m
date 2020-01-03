% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:37
% EndTime: 2019-12-31 20:58:00
% DurationCPUTime: 15.44s
% Computational Cost: add. (34554->595), mult. (42021->761), div. (0->0), fcn. (39046->6), ass. (0->363)
t876 = Icges(5,1) + Icges(6,1);
t547 = qJ(2) + qJ(3);
t528 = cos(t547);
t527 = sin(t547);
t724 = Icges(4,4) * t527;
t482 = Icges(4,1) * t528 - t724;
t549 = sin(qJ(1));
t551 = cos(qJ(1));
t382 = Icges(4,5) * t549 + t482 * t551;
t521 = Icges(5,5) * t527;
t828 = Icges(5,1) * t528 + t521;
t875 = Icges(5,4) * t549 + t551 * t828 + t382;
t522 = Icges(6,4) * t527;
t874 = t521 + t522 + (-Icges(6,2) - Icges(5,3)) * t528;
t523 = Icges(4,4) * t528;
t481 = Icges(4,1) * t527 + t523;
t722 = Icges(5,5) * t528;
t723 = Icges(6,4) * t528;
t873 = t527 * t876 + t481 - t722 - t723;
t861 = Icges(4,2) * t528 + t724 - t874;
t483 = pkin(3) * t527 - qJ(4) * t528;
t748 = pkin(4) * t527;
t590 = rSges(6,1) * t527 - rSges(6,2) * t528 + t483 + t748;
t548 = sin(qJ(2));
t750 = pkin(2) * t548;
t563 = t590 + t750;
t275 = t563 * t549;
t277 = t563 * t551;
t630 = rSges(5,1) * t527 - rSges(5,3) * t528 + t483;
t592 = t630 + t750;
t302 = t592 * t549;
t304 = t592 * t551;
t543 = t551 * rSges(5,2);
t550 = cos(qJ(2));
t749 = pkin(2) * t550;
t524 = pkin(1) + t749;
t810 = -pkin(7) - pkin(6);
t625 = -t524 * t549 - t551 * t810;
t734 = rSges(5,3) + qJ(4);
t780 = rSges(5,1) + pkin(3);
t823 = -t527 * t734 - t528 * t780;
t264 = t549 * t823 + t543 + t625;
t525 = t549 * t810;
t742 = t549 * rSges(5,2);
t265 = t742 - t525 + (t524 - t823) * t551;
t695 = t528 * t551;
t696 = t528 * t549;
t662 = t264 * t695 + t265 * t696;
t619 = rSges(6,1) + pkin(3) + pkin(4);
t735 = rSges(6,2) + qJ(4);
t824 = -t527 * t735 - t528 * t619;
t856 = rSges(6,3) + qJ(5);
t238 = t549 * t824 - t551 * t856 + t625;
t239 = -t525 + (t524 - t824) * t551 - t856 * t549;
t663 = t238 * t695 + t239 * t696;
t703 = t527 * t551;
t704 = t527 * t549;
t811 = m(6) / 0.2e1;
t812 = m(5) / 0.2e1;
t731 = (-t302 * t703 + t304 * t704 + t662) * t812 + (-t275 * t703 + t277 * t704 + t663) * t811;
t564 = -t528 * t735 + t748;
t510 = rSges(6,1) * t704;
t515 = pkin(3) * t704;
t624 = t510 + t515;
t268 = (t564 + t750) * t549 + t624;
t594 = t619 * t527;
t502 = qJ(4) * t695;
t514 = rSges(6,2) * t695;
t626 = t502 + t514;
t269 = (-t594 - t750) * t551 + t626;
t608 = t734 * t528;
t511 = rSges(5,1) * t704;
t623 = t511 + t515;
t289 = (-t608 + t750) * t549 + t623;
t611 = t780 * t527;
t513 = rSges(5,3) * t695;
t627 = t502 + t513;
t290 = (-t611 - t750) * t551 + t627;
t732 = ((t289 * t551 + t290 * t549) * t527 + t662) * t812 + ((t268 * t551 + t269 * t549) * t527 + t663) * t811;
t12 = t732 - t731;
t872 = t12 * qJD(1);
t470 = Icges(5,3) * t527 + t722;
t367 = -Icges(5,6) * t551 + t470 * t549;
t473 = Icges(6,2) * t527 + t723;
t371 = Icges(6,6) * t551 + t473 * t549;
t871 = t367 + t371;
t827 = Icges(6,1) * t528 + t522;
t377 = Icges(6,5) * t551 + t549 * t827;
t379 = -Icges(5,4) * t551 + t549 * t828;
t870 = -t377 - t379;
t291 = t590 * t549;
t293 = t590 * t551;
t311 = t630 * t549;
t313 = t630 * t551;
t667 = (-t291 * t703 + t293 * t704 + t663) * t811 + (-t311 * t703 + t313 * t704 + t662) * t812;
t286 = t549 * t564 + t624;
t287 = -t551 * t594 + t626;
t308 = -t549 * t608 + t623;
t309 = -t551 * t611 + t627;
t730 = ((t308 * t551 + t309 * t549) * t527 + t662) * t812 + ((t286 * t551 + t287 * t549) * t527 + t663) * t811;
t15 = t730 - t667;
t869 = t15 * qJD(1);
t505 = Icges(5,5) * t695;
t368 = Icges(5,6) * t549 + Icges(5,3) * t703 + t505;
t471 = Icges(4,5) * t528 - Icges(4,6) * t527;
t714 = t471 * t551;
t370 = Icges(4,3) * t549 + t714;
t474 = Icges(5,4) * t528 + Icges(5,6) * t527;
t713 = t474 * t551;
t868 = t368 * t703 + t875 * t695 + (Icges(5,2) * t549 + t370 + t713) * t549;
t867 = t470 + t473;
t369 = Icges(4,5) * t696 - Icges(4,6) * t704 - Icges(4,3) * t551;
t507 = Icges(4,4) * t704;
t381 = Icges(4,1) * t696 - Icges(4,5) * t551 - t507;
t866 = -t549 * t369 - t371 * t703 + (-t377 - t381) * t695;
t476 = -Icges(4,2) * t527 + t523;
t864 = t476 + t873;
t863 = t482 + t827 + t828;
t862 = (-Icges(4,6) + Icges(5,6) - Icges(6,6)) * t528 + (-Icges(5,4) - Icges(4,5) + Icges(6,5)) * t527;
t378 = -Icges(6,5) * t549 + t551 * t827;
t860 = -t551 * t861 + t378 + t875;
t859 = -Icges(4,2) * t696 + t549 * t874 + t381 - t507 - t870;
t506 = Icges(6,4) * t695;
t372 = Icges(6,2) * t703 - Icges(6,6) * t549 + t506;
t376 = Icges(4,6) * t549 + t476 * t551;
t858 = -t481 * t551 - t703 * t876 + t368 + t372 - t376 + t505 + t506;
t375 = Icges(4,4) * t696 - Icges(4,2) * t704 - Icges(4,6) * t551;
t857 = t549 * t873 + t375 - t871;
t468 = Icges(6,5) * t528 + Icges(6,6) * t527;
t683 = t549 * t468;
t365 = Icges(6,3) * t551 + t683;
t684 = t549 * t365;
t855 = t375 * t703 + t684 + t866;
t854 = -t376 * t703 + t868;
t545 = t549 ^ 2;
t546 = t551 ^ 2;
t621 = t545 + t546;
t570 = t371 * t527 + t377 * t528;
t832 = (t367 * t527 + t379 * t528) * t549;
t852 = -t365 * t551 - t549 * t570 - t832 - t868;
t813 = m(4) / 0.2e1;
t783 = t549 / 0.2e1;
t781 = -t551 / 0.2e1;
t847 = t551 / 0.2e1;
t682 = t549 * t474;
t340 = t549 * (-Icges(5,2) * t551 + t682);
t219 = t367 * t703 + t379 * t695 + t340;
t846 = t219 * t551;
t725 = Icges(3,4) * t548;
t494 = Icges(3,2) * t550 + t725;
t497 = Icges(3,1) * t550 - t725;
t845 = (t497 / 0.2e1 - t494 / 0.2e1) * t548;
t486 = rSges(4,1) * t527 + rSges(4,2) * t528;
t385 = rSges(4,1) * t696 - rSges(4,2) * t704 - rSges(4,3) * t551;
t306 = -t385 + t625;
t607 = -rSges(4,2) * t703 + t549 * rSges(4,3);
t743 = rSges(4,1) * t528;
t307 = -t525 + (t524 + t743) * t551 + t607;
t490 = -rSges(4,2) * t527 + t743;
t561 = (-t306 * t551 - t307 * t549) * t490;
t487 = pkin(3) * t528 + qJ(4) * t527;
t489 = rSges(5,1) * t528 + rSges(5,3) * t527;
t629 = -t487 - t489;
t312 = t629 * t549;
t314 = t629 * t551;
t664 = t264 * t314 + t265 * t312;
t488 = rSges(6,1) * t528 + rSges(6,2) * t527;
t589 = -pkin(4) * t528 - t487 - t488;
t292 = t589 * t549;
t294 = t589 * t551;
t665 = t238 * t294 + t239 * t292;
t560 = t486 + t750;
t830 = t560 * t551;
t831 = t560 * t549;
t615 = (t561 + (t549 * t830 - t551 * t831) * t486) * t813 + (-t268 * t293 - t269 * t291 + t665) * t811 + (-t289 * t313 - t290 * t311 + t664) * t812;
t451 = t486 * t549;
t453 = t486 * t551;
t616 = (-t451 * t830 + t453 * t831 + t561) * t813 + (-t275 * t287 - t277 * t286 + t665) * t811 + (-t302 * t309 - t304 * t308 + t664) * t812;
t843 = t615 - t616;
t636 = t621 * t487;
t165 = t636 + (pkin(4) * t696 + t488 * t549) * t549 + (pkin(4) * t695 + t488 * t551) * t551;
t202 = t549 * (t489 * t549 - t543) + t551 * (t489 * t551 + t742) + t636;
t637 = t549 * (qJ(4) * t696 - t515) + t551 * (-pkin(3) * t703 + t502);
t223 = t549 * (rSges(5,3) * t696 - t511) + t551 * (-rSges(5,1) * t703 + t513) + t637;
t587 = -t223 * t528 + t312 * t704 + t314 * t703;
t194 = -t621 * t748 + t549 * (rSges(6,2) * t696 - t510) + t551 * (-rSges(6,1) * t703 + t514) + t637;
t588 = -t194 * t528 + t292 * t704 + t294 * t703;
t658 = -t311 * t696 - t313 * t695;
t685 = t549 * t291;
t660 = -t293 * t695 - t528 * t685;
t566 = (t165 * t527 + t588 + t660) * t811 + (t202 * t527 + t587 + t658) * t812;
t458 = t621 * t527;
t105 = t165 * t458 + t660;
t124 = t202 * t458 + t658;
t544 = t551 * pkin(6);
t652 = -t549 * (pkin(1) * t549 - t544 + t625) + t551 * (-t549 * pkin(6) - t525 + (-pkin(1) + t524) * t551);
t149 = t165 + t652;
t686 = t549 * t275;
t661 = -t277 * t695 - t528 * t686;
t79 = t149 * t458 + t661;
t161 = t202 + t652;
t659 = -t302 * t696 - t304 * t695;
t98 = t161 * t458 + t659;
t747 = (t105 + t79) * t811 + (t124 + t98) * t812;
t841 = t747 - t566;
t840 = t862 * t549;
t839 = t862 * t551;
t838 = (-t861 + t863) * t528 + (-t864 + t867) * t527;
t836 = (t549 * t858 + t551 * t857) * t528 + (-t549 * t860 + t551 * t859) * t527;
t705 = t527 * t528;
t628 = t621 * t705;
t829 = (m(5) / 0.4e1 + m(6) / 0.4e1) * (t628 - t705);
t539 = Icges(3,4) * t550;
t495 = -Icges(3,2) * t548 + t539;
t496 = Icges(3,1) * t548 + t539;
t612 = t161 * t527 + t659;
t613 = t149 * t527 + t661;
t746 = (t588 + t613) * t811 + (t587 + t612) * t812;
t431 = Icges(3,5) * t549 + t497 * t551;
t631 = -t494 * t551 + t431;
t694 = t548 * t549;
t518 = Icges(3,4) * t694;
t681 = t549 * t550;
t430 = Icges(3,1) * t681 - Icges(3,5) * t551 - t518;
t632 = -Icges(3,2) * t681 + t430 - t518;
t429 = Icges(3,6) * t549 + t495 * t551;
t633 = -t496 * t551 - t429;
t428 = Icges(3,4) * t681 - Icges(3,2) * t694 - Icges(3,6) * t551;
t634 = t496 * t549 + t428;
t822 = (-t549 * t631 + t551 * t632) * t548 + (t549 * t633 + t551 * t634) * t550;
t319 = t382 * t696;
t598 = t370 * t551 - t319;
t216 = -t376 * t704 - t598;
t715 = t468 * t551;
t366 = -Icges(6,3) * t549 + t715;
t596 = t376 * t527 - t369;
t656 = t372 * t703 + t378 * t695;
t717 = t375 * t527;
t556 = (-t846 + t855 * t551 + (-t366 * t549 + t656 + t854) * t549) * t847 + (-t846 + (t216 - t319 + (t370 + t717) * t551 + t866) * t551 + (t656 - t832 + (-t366 - t570) * t549 - t852) * t549) * t781 + (t365 * t546 + ((t369 + t596) * t551 + t852 + t854) * t551 + (-(t381 * t528 - t717) * t551 + t216 + t684 + t219 - t340 + t598 + t549 * t596 - t855) * t549) * t783;
t555 = (-t867 / 0.2e1 + t864 / 0.2e1) * t528 + (-t861 / 0.2e1 + t863 / 0.2e1) * t527;
t818 = 0.4e1 * qJD(1);
t817 = 2 * qJD(2);
t816 = 4 * qJD(2);
t815 = 2 * qJD(3);
t814 = 4 * qJD(3);
t67 = t161 * t223 - t302 * t312 - t304 * t314;
t803 = m(5) * t67;
t78 = t202 * t223 - t311 * t312 - t313 * t314;
t799 = m(5) * t78;
t798 = m(5) * t98;
t48 = t149 * t194 - t275 * t292 - t277 * t294;
t793 = m(6) * t48;
t64 = t165 * t194 - t291 * t292 - t293 * t294;
t789 = m(6) * t64;
t788 = m(6) * t79;
t784 = -t549 / 0.2e1;
t744 = rSges(3,1) * t550;
t610 = pkin(1) + t744;
t622 = rSges(3,2) * t694 + rSges(3,3) * t551;
t332 = -t549 * t610 + t544 + t622;
t693 = t548 * t551;
t520 = rSges(3,2) * t693;
t333 = -t520 + t610 * t551 + (rSges(3,3) + pkin(6)) * t549;
t498 = rSges(3,1) * t548 + rSges(3,2) * t550;
t465 = t498 * t549;
t466 = t498 * t551;
t779 = m(3) * (t332 * t465 - t333 * t466);
t266 = t549 * t385 + t551 * (rSges(4,1) * t695 + t607);
t190 = t266 + t652;
t288 = -t451 * t549 - t453 * t551;
t111 = t190 * t288 + (t549 * t831 + t551 * t830) * t490;
t778 = m(4) * t111;
t155 = t486 * t490 * t621 + t266 * t288;
t775 = m(4) * t155;
t774 = m(4) * (t306 * t831 - t307 * t830);
t773 = m(4) * (t306 * t451 - t307 * t453);
t767 = m(5) * (t264 * t289 + t265 * t290);
t766 = m(5) * t124;
t764 = m(5) * (t264 * t308 + t265 * t309);
t763 = m(5) * (-t264 * t704 + t265 * t703);
t760 = m(6) * t105;
t759 = m(6) * (t238 * t268 + t239 * t269);
t757 = m(6) * (t238 * t286 + t239 * t287);
t756 = m(6) * (-t238 * t704 + t239 * t703);
t755 = m(6) * (-t238 * t551 - t239 * t549);
t754 = m(6) * (-t268 * t549 + t269 * t551);
t753 = m(6) * (t277 * t551 + t686);
t752 = m(6) * (-t286 * t549 + t287 * t551);
t751 = m(6) * (t293 * t551 + t685);
t745 = m(6) * qJD(2);
t716 = t428 * t548;
t680 = t550 * t551;
t426 = Icges(3,5) * t681 - Icges(3,6) * t694 - Icges(3,3) * t551;
t651 = -t549 * t426 - t430 * t680;
t578 = Icges(3,5) * t550 - Icges(3,6) * t548;
t427 = Icges(3,3) * t549 + t551 * t578;
t650 = t427 * t549 + t431 * t680;
t295 = m(6) * t458;
t620 = t295 * qJD(1);
t614 = qJD(1) * t755;
t609 = -t490 - t749;
t355 = t431 * t681;
t597 = t551 * t427 - t355;
t595 = t429 * t548 - t426;
t593 = t621 * t750;
t591 = t629 - t749;
t585 = ((-t549 * t840 + t836) * t551 + t839 * t545) * t783 + ((-t551 * t839 + t836) * t549 + t840 * t546) * t781;
t577 = -Icges(3,5) * t548 - Icges(3,6) * t550;
t562 = t589 - t749;
t553 = -t556 + (t549 * t471 + t527 * t858 + t528 * t860 + t838 * t551 + t682 - t683) * t783 + (-t527 * t857 + t528 * t859 + t838 * t549 - t713 - t714 + t715) * t781;
t552 = -t555 + ((t375 + t871) * t528 + (t381 + t870) * t527) * (t783 + t784);
t500 = -rSges(3,2) * t548 + t744;
t460 = t577 * t551;
t459 = t577 * t549;
t363 = t609 * t551;
t361 = t609 * t549;
t305 = t591 * t551;
t303 = t591 * t549;
t278 = t562 * t551;
t276 = t562 * t549;
t260 = -t593 + t288;
t237 = -t429 * t693 + t650;
t236 = -t428 * t693 - t651;
t235 = -t429 * t694 - t597;
t228 = 0.4e1 * t829;
t227 = t228 * qJD(4);
t207 = t292 * t551 - t294 * t549;
t201 = t751 / 0.2e1;
t197 = m(6) * t207 * qJD(3);
t196 = -t593 + t223;
t195 = t752 / 0.2e1;
t177 = t753 / 0.2e1;
t170 = t754 / 0.2e1;
t169 = -t593 + t194;
t157 = -t236 * t551 + t237 * t549;
t156 = -(-t549 * (-t430 * t550 + t716) - t551 * t426) * t551 + t235 * t549;
t153 = (-0.4e1 * t829 + 0.2e1 * (t811 + t812) * (-t458 * t528 + t628)) * qJD(4);
t87 = t756 + t763;
t85 = t201 - t752 / 0.2e1;
t84 = t201 + t195;
t83 = t195 - t751 / 0.2e1;
t74 = t177 - t754 / 0.2e1;
t73 = t177 + t170;
t72 = t170 - t753 / 0.2e1;
t54 = t760 + t766;
t53 = (t235 - t355 + (t427 + t716) * t551 + t651) * t551 + t650 * t549;
t52 = (t551 * t595 + t237 - t650) * t551 + (t549 * t595 + t236 + t597) * t549;
t24 = t788 + t798;
t23 = t555 + t757 + t764 + t773;
t20 = t779 + t774 + t767 + t759 + (t496 / 0.2e1 + t495 / 0.2e1) * t550 + t555 + t845;
t16 = t667 + t730;
t13 = t731 + t732;
t11 = t585 + t775 + t789 + t799;
t10 = t11 * qJD(3);
t9 = t566 + t747 - t746;
t8 = t746 + t841;
t7 = t746 - t841;
t6 = t585 + t778 + t793 + t803;
t4 = (-t53 / 0.2e1 + t157 / 0.2e1) * t551 + (t156 / 0.2e1 + t52 / 0.2e1) * t549 + t556;
t3 = t556 - t843;
t2 = t556 + t843;
t1 = t553 + t615 + t616;
t5 = [qJD(2) * t20 + qJD(3) * t23 + qJD(4) * t87 + qJD(5) * t755, t20 * qJD(1) + t1 * qJD(3) + t13 * qJD(4) + t73 * qJD(5) + ((t306 * t363 + t307 * t361) * t813 + (t264 * t305 + t265 * t303 - t289 * t304 - t290 * t302) * t812 + (t238 * t278 + t239 * t276 - t268 * t277 - t269 * t275) * t811 + m(3) * ((-t332 * t551 - t333 * t549) * t500 + (-t465 * t551 + t466 * t549) * t498) / 0.2e1) * t817 + (t553 + t53 * t847 + (t548 * t633 + t550 * t631) * t783 + (t156 + t52) * t784 + (-t548 * t634 + t550 * t632 + t157) * t781 + (t546 / 0.2e1 + t545 / 0.2e1) * t578) * qJD(2), t23 * qJD(1) + t1 * qJD(2) + t553 * qJD(3) + t16 * qJD(4) + t84 * qJD(5) + ((-t308 * t313 - t309 * t311 + t664) * t812 + (-t286 * t293 - t287 * t291 + t665) * t811 + (t561 + (-t451 * t551 + t453 * t549) * t486) * t813) * t815, qJD(1) * t87 + qJD(2) * t13 + qJD(3) * t16, qJD(2) * t73 + qJD(3) * t84 + t614; (t552 - (t496 + t495) * t550 / 0.2e1 - t845) * qJD(1) + t4 * qJD(2) + t3 * qJD(3) - t12 * qJD(4) + t74 * qJD(5) + (-t779 / 0.4e1 - t774 / 0.4e1 - t767 / 0.4e1 - t759 / 0.4e1) * t818, t4 * qJD(1) + (m(3) * ((t549 * (rSges(3,1) * t681 - t622) + t551 * (rSges(3,1) * t680 + t549 * rSges(3,3) - t520)) * (-t465 * t549 - t466 * t551) + t621 * t500 * t498) + m(6) * (t149 * t169 - t275 * t276 - t277 * t278) + m(5) * (t161 * t196 - t302 * t303 - t304 * t305) + m(4) * (t190 * t260 - t361 * t831 - t363 * t830) + (t546 * t459 + (-t551 * t460 + t822) * t549) * t781 + (t545 * t460 + (-t549 * t459 + t822) * t551) * t783 + t585) * qJD(2) + t6 * qJD(3) + t24 * qJD(4), t3 * qJD(1) + t6 * qJD(2) + t585 * qJD(3) + t8 * qJD(4) + (-t789 / 0.4e1 - t799 / 0.4e1 - t775 / 0.4e1) * t814 + ((t64 + t48) * t811 + (t78 + t67) * t812 + (t111 + t155) * t813) * t815, qJD(2) * t24 + qJD(3) * t8 + t153 - t872, t74 * qJD(1); t552 * qJD(1) + t2 * qJD(2) + t556 * qJD(3) - t15 * qJD(4) + t85 * qJD(5) + (-t773 / 0.4e1 - t764 / 0.4e1 - t757 / 0.4e1) * t818, t2 * qJD(1) + t585 * qJD(2) + t10 + t9 * qJD(4) + (-t793 / 0.4e1 - t803 / 0.4e1 - t778 / 0.4e1) * t816 + ((t165 * t169 - t276 * t291 - t278 * t293 + t48) * t811 + (t196 * t202 - t303 * t311 - t305 * t313 + t67) * t812 + (t266 * t260 + (-t361 * t549 - t363 * t551) * t486 + t111) * t813) * t817, qJD(1) * t556 + qJD(2) * t11 + qJD(4) * t54 + t10, qJD(2) * t9 + qJD(3) * t54 + t153 - t869, t85 * qJD(1); t12 * qJD(2) + t15 * qJD(3) - t295 * qJD(5) + (-t756 / 0.4e1 - t763 / 0.4e1) * t818, t872 + t7 * qJD(3) + t227 + (-t788 / 0.4e1 - t798 / 0.4e1) * t816 + ((-t528 * t169 + (t276 * t549 + t278 * t551) * t527 + t613) * t811 + (-t528 * t196 + (t303 * t549 + t305 * t551) * t527 + t612) * t812) * t817, t869 + t7 * qJD(2) + t227 + (-t760 / 0.4e1 - t766 / 0.4e1) * t814 + t566 * t815, (qJD(2) + qJD(3)) * t228, -t620; qJD(2) * t72 + qJD(3) * t83 + qJD(4) * t295 - t614, t72 * qJD(1) + (t276 * t551 - t278 * t549) * t745 + t197, qJD(1) * t83 + t207 * t745 + t197, t620, 0;];
Cq = t5;
