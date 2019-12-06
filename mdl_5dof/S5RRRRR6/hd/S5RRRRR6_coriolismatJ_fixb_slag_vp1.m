% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR6
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:04
% EndTime: 2019-12-05 19:00:25
% DurationCPUTime: 13.45s
% Computational Cost: add. (107289->687), mult. (67925->817), div. (0->0), fcn. (62840->10), ass. (0->434)
t587 = qJ(1) + qJ(2);
t579 = sin(t587);
t586 = qJ(3) + qJ(4);
t583 = qJ(5) + t586;
t574 = cos(t583);
t581 = cos(t587);
t728 = t574 * t581;
t573 = sin(t583);
t732 = t573 * t581;
t418 = Icges(6,5) * t728 - Icges(6,6) * t732 + Icges(6,3) * t579;
t566 = Icges(6,4) * t574;
t880 = Icges(6,2) * t573 - t566;
t419 = Icges(6,6) * t581 + t880 * t579;
t733 = t573 * t579;
t502 = Icges(6,5) * t574 - Icges(6,6) * t573;
t746 = t502 * t579;
t677 = t581 * (Icges(6,3) * t581 - t746) + t419 * t733;
t420 = Icges(6,4) * t728 - Icges(6,2) * t732 + Icges(6,6) * t579;
t529 = Icges(6,4) * t732;
t422 = Icges(6,1) * t728 + Icges(6,5) * t579 - t529;
t884 = (t420 * t573 - t422 * t574) * t581;
t902 = t579 * t418 + t677 - t884;
t580 = cos(t586);
t716 = t580 * t581;
t578 = sin(t586);
t724 = t578 * t581;
t435 = Icges(5,5) * t716 - Icges(5,6) * t724 + Icges(5,3) * t579;
t567 = Icges(5,4) * t580;
t879 = Icges(5,2) * t578 - t567;
t436 = Icges(5,6) * t581 + t879 * t579;
t725 = t578 * t579;
t517 = Icges(5,5) * t580 - Icges(5,6) * t578;
t743 = t517 * t579;
t675 = t581 * (Icges(5,3) * t581 - t743) + t436 * t725;
t437 = Icges(5,4) * t716 - Icges(5,2) * t724 + Icges(5,6) * t579;
t541 = Icges(5,4) * t724;
t439 = Icges(5,1) * t716 + Icges(5,5) * t579 - t541;
t885 = (t437 * t578 - t439 * t580) * t581;
t901 = t579 * t435 + t675 - t885;
t590 = cos(qJ(3));
t584 = t590 * pkin(3);
t575 = t584 + pkin(2);
t784 = pkin(4) * t580;
t536 = t575 + t784;
t592 = -pkin(8) - pkin(7);
t585 = -pkin(9) + t592;
t654 = -rSges(6,2) * t733 - t581 * rSges(6,3);
t778 = rSges(6,1) * t574;
t351 = t581 * t585 + (t536 + t778) * t579 + t654;
t783 = sin(qJ(1)) * pkin(1);
t338 = t351 + t783;
t900 = -t338 + t351;
t501 = t581 * t536;
t632 = -rSges(6,1) * t728 + rSges(6,2) * t732;
t352 = -t501 + (-rSges(6,3) + t585) * t579 + t632;
t787 = cos(qJ(1)) * pkin(1);
t339 = t352 - t787;
t899 = t339 - t352;
t779 = rSges(5,1) * t580;
t640 = t575 + t779;
t653 = -rSges(5,2) * t725 - t581 * rSges(5,3);
t708 = t581 * t592;
t369 = t579 * t640 + t653 + t708;
t362 = t369 + t783;
t898 = -t362 + t369;
t897 = t436 * t724;
t896 = t419 * t732;
t855 = m(4) / 0.2e1;
t854 = m(5) / 0.2e1;
t853 = m(6) / 0.2e1;
t895 = -t579 / 0.2e1;
t838 = t579 / 0.2e1;
t836 = t581 / 0.2e1;
t835 = m(3) * (-t787 * (rSges(3,1) * t579 + rSges(3,2) * t581) - (-rSges(3,1) * t581 + t579 * rSges(3,2)) * t783);
t576 = t579 ^ 2;
t577 = t581 ^ 2;
t649 = t576 + t577;
t565 = t579 * t592;
t639 = rSges(5,2) * t724 - rSges(5,3) * t579;
t370 = -t581 * t640 + t565 + t639;
t524 = -rSges(5,2) * t578 + t779;
t740 = t524 * t579;
t330 = t370 * t740;
t522 = rSges(5,1) * t578 + rSges(5,2) * t580;
t588 = sin(qJ(3));
t786 = pkin(3) * t588;
t625 = (t522 + t786) * t581;
t741 = t522 * t579;
t368 = t625 * t741;
t729 = t574 * t579;
t463 = rSges(6,1) * t733 + rSges(6,2) * t729;
t785 = pkin(4) * t578;
t537 = -t785 - t786;
t384 = -t537 * t579 + t463;
t464 = rSges(6,1) * t732 + rSges(6,2) * t728;
t385 = -t537 * t581 + t464;
t507 = rSges(6,1) * t573 + rSges(6,2) * t574;
t481 = t579 * t507;
t546 = pkin(4) * t725;
t409 = t546 + t481;
t411 = (t507 + t785) * t581;
t685 = -t411 * t384 + t409 * t385;
t508 = -rSges(6,2) * t573 + t778;
t482 = t579 * t508;
t721 = t579 * t580;
t410 = pkin(4) * t721 + t482;
t641 = -t508 - t784;
t412 = t641 * t581;
t691 = -t412 * t351 + t410 * t352;
t479 = rSges(5,1) * t725 + rSges(5,2) * t721;
t720 = t579 * t588;
t563 = pkin(3) * t720;
t425 = t563 + t479;
t751 = t425 * t522;
t760 = t369 * t524;
t699 = (t685 + t691) * t853 + (t330 + t368 + (-t751 + t760) * t581) * t854;
t442 = t563 + t741;
t480 = t522 * t581;
t679 = t442 * t480 - t479 * t625;
t380 = t563 + t409;
t382 = (t507 - t537) * t581;
t399 = t546 + t463;
t400 = pkin(4) * t724 + t464;
t688 = t380 * t400 - t382 * t399;
t739 = t524 * t581;
t700 = (t688 + t691) * t853 + (t369 * t739 + t330 + t679) * t854;
t894 = t699 - t700;
t363 = t370 - t787;
t321 = t363 * t740;
t694 = -t412 * t338 + t410 * t339;
t761 = t362 * t524;
t701 = (t685 + t694) * t853 + (t321 + t368 + (-t751 + t761) * t581) * t854;
t702 = (t688 + t694) * t853 + (t362 * t739 + t321 + t679) * t854;
t893 = t701 - t702;
t837 = -t581 / 0.2e1;
t892 = t836 + t837;
t667 = -Icges(6,2) * t728 + t422 - t529;
t528 = Icges(6,4) * t733;
t421 = -Icges(6,1) * t729 + Icges(6,5) * t581 + t528;
t668 = Icges(6,2) * t729 + t421 + t528;
t883 = Icges(6,1) * t573 + t566;
t669 = t581 * t883 + t420;
t670 = -t579 * t883 + t419;
t891 = -(t579 * t667 + t581 * t668) * t573 - (t579 * t669 + t581 * t670) * t574;
t663 = -Icges(5,2) * t716 + t439 - t541;
t540 = Icges(5,4) * t725;
t438 = -Icges(5,1) * t721 + Icges(5,5) * t581 + t540;
t664 = Icges(5,2) * t721 + t438 + t540;
t882 = Icges(5,1) * t578 + t567;
t665 = t581 * t882 + t437;
t666 = -t579 * t882 + t436;
t890 = -(t579 * t663 + t581 * t664) * t578 - (t579 * t665 + t581 * t666) * t580;
t202 = -t351 * t463 + t352 * t464;
t889 = m(6) * t202;
t780 = rSges(4,1) * t590;
t555 = -rSges(4,2) * t588 + t780;
t887 = t555 * t855;
t709 = t581 * t590;
t710 = t581 * t588;
t453 = Icges(4,4) * t709 - Icges(4,2) * t710 + Icges(4,6) * t579;
t560 = Icges(4,4) * t710;
t455 = Icges(4,1) * t709 + Icges(4,5) * t579 - t560;
t886 = (t453 * t588 - t455 * t590) * t581;
t291 = t625 * t363;
t304 = t625 * t370;
t172 = -t352 * t338 + t339 * t351;
t190 = -t370 * t362 + t363 * t369;
t572 = t581 * pkin(7);
t642 = pkin(2) + t780;
t650 = -rSges(4,2) * t720 - t581 * rSges(4,3);
t390 = t579 * t642 - t572 + t650;
t378 = t390 + t783;
t562 = rSges(4,2) * t710;
t391 = t562 - t642 * t581 + (-rSges(4,3) - pkin(7)) * t579;
t379 = t391 - t787;
t205 = -t391 * t378 + t379 * t390;
t582 = Icges(4,4) * t590;
t881 = Icges(4,1) * t588 + t582;
t878 = Icges(4,2) * t588 - t582;
t782 = pkin(2) - t575;
t393 = t581 * (pkin(7) * t579 + t782 * t581 + t565);
t404 = t581 * (rSges(5,1) * t716 - t639);
t415 = t782 * t579 - t572 - t708;
t440 = -rSges(5,1) * t721 - t653;
t213 = -t393 + t404 + (-t415 - t440) * t579;
t356 = -t479 * t579 - t581 * t480;
t135 = t213 * t356 + t442 * t740 + t625 * t739;
t398 = t581 * (rSges(6,3) * t579 - t632);
t673 = t398 - t581 * (t575 * t581 + t579 * t585 - t501 - t565);
t426 = -rSges(6,1) * t729 - t654;
t678 = -(-t585 + t592) * t581 - (-t536 + t575) * t579 - t426;
t161 = -t393 + (-t415 + t678) * t579 + t673;
t433 = t581 * t464;
t342 = -t463 * t579 - t433;
t287 = -t649 * t785 + t342;
t686 = t380 * t410 - t382 * t412;
t88 = t161 * t287 + t686;
t877 = m(5) * t135 + m(6) * t88;
t648 = qJD(1) + qJD(2);
t703 = (t900 * t409 + t899 * t411) * t853 + ((t363 - t370) * t480 + t898 * t741) * t854;
t193 = -t338 * t399 + t339 * t400;
t197 = -t351 * t399 + t352 * t400;
t209 = -t362 * t479 + t363 * t480;
t214 = -t369 * t479 + t370 * t480;
t876 = (t197 + t193) * t853 + (t214 + t209) * t854;
t553 = rSges(4,1) * t588 + rSges(4,2) * t590;
t498 = t553 * t581;
t497 = t553 * t579;
t756 = t390 * t497;
t758 = t378 * t497;
t645 = (t898 * t442 + t291 - t304) * t854 + (-t758 + t756 + (t379 - t391) * t498) * t855 + (t900 * t380 + t899 * t382) * t853;
t189 = -t338 * t384 + t339 * t385;
t192 = -t351 * t384 + t352 * t385;
t201 = -t362 * t425 + t291;
t204 = -t369 * t425 + t304;
t229 = t379 * t498 - t758;
t236 = t391 * t498 - t756;
t875 = (t204 + t201) * t854 + (t236 + t229) * t855 + (t192 + t189) * t853;
t768 = Icges(6,4) * t573;
t503 = Icges(6,2) * t574 + t768;
t506 = Icges(6,1) * t574 - t768;
t874 = (-t880 + t883) * t573 + (t503 - t506) * t574;
t769 = Icges(5,4) * t578;
t518 = Icges(5,2) * t580 + t769;
t521 = Icges(5,1) * t580 - t769;
t873 = (-t879 + t882) * t578 + (t518 - t521) * t580;
t770 = Icges(4,4) * t588;
t549 = Icges(4,2) * t590 + t770;
t552 = Icges(4,1) * t590 - t770;
t868 = (-t878 + t881) * t588 + (t549 - t552) * t590;
t659 = -Icges(4,2) * t709 + t455 - t560;
t661 = t581 * t881 + t453;
t867 = t588 * t659 + t590 * t661;
t559 = Icges(4,4) * t720;
t719 = t579 * t590;
t454 = -Icges(4,1) * t719 + Icges(4,5) * t581 + t559;
t660 = Icges(4,2) * t719 + t454 + t559;
t452 = Icges(4,6) * t581 + t878 * t579;
t662 = -t579 * t881 + t452;
t866 = -t588 * t660 - t590 * t662;
t865 = (t881 / 0.2e1 - t878 / 0.2e1) * t590 + (t552 / 0.2e1 - t549 / 0.2e1) * t588;
t237 = -t438 * t721 + t675;
t238 = t581 * t435 + t437 * t725 - t439 * t721;
t637 = -t438 * t580 - t435;
t750 = t436 * t578;
t864 = (t237 * t581 + t238 * t579) * t895 + ((t885 + t901) * t581 + ((t637 - t750) * t581 + t238 + t897) * t579) * t838 + ((t238 + (-t435 + t750) * t581 - t897) * t581 + (t579 * t637 - t237 + t901) * t579) * t836;
t863 = (t882 / 0.2e1 - t879 / 0.2e1) * t580 + (t521 / 0.2e1 - t518 / 0.2e1) * t578;
t224 = -t421 * t729 + t677;
t225 = t581 * t418 + t420 * t733 - t422 * t729;
t638 = -t421 * t574 - t418;
t752 = t419 * t573;
t32 = (t224 * t581 + t225 * t579) * t895 + ((t884 + t902) * t581 + ((t638 - t752) * t581 + t225 + t896) * t579) * t838 + ((t225 + (-t418 + t752) * t581 - t896) * t581 + (t579 * t638 - t224 + t902) * t579) * t836;
t634 = (-t880 / 0.2e1 + t883 / 0.2e1) * t574 + (-t503 / 0.2e1 + t506 / 0.2e1) * t573;
t861 = 4 * qJD(1);
t859 = 4 * qJD(2);
t858 = 2 * qJD(3);
t857 = 2 * qJD(4);
t199 = t579 * t678 + t673;
t697 = t199 * t342 + t409 * t482;
t698 = t161 * t342 + t380 * t482;
t744 = t508 * t581;
t849 = m(6) * ((t382 + t411) * t744 + t697 + t698);
t315 = -t426 * t579 + t398;
t207 = t315 * t287;
t359 = t410 * t481;
t753 = t412 * t507;
t847 = m(6) * (t207 + t359 + (t382 * t508 - t753) * t581 + t698);
t228 = t577 * (t537 + t786) - t433 + (t563 - t384) * t579;
t601 = t359 + (t411 * t508 - t753) * t581 + t697;
t845 = m(6) * (t228 * t315 + t601);
t682 = t409 * t410 - t411 * t412;
t841 = m(6) * (t199 * t228 + t682);
t831 = m(4) * t205;
t829 = m(4) * t229;
t828 = m(4) * t236;
t818 = m(5) * t190;
t323 = -t440 * t579 + t404;
t196 = t649 * t522 * t524 + t323 * t356;
t194 = m(5) * t196;
t816 = m(5) * t201;
t815 = m(5) * t204;
t814 = m(5) * t209;
t813 = m(5) * t214;
t200 = -t338 * t463 + t339 * t464;
t809 = m(6) * (t202 + t200);
t805 = m(6) * (t899 * t581 * t507 + t900 * t481);
t289 = t339 * t482;
t629 = t338 * t744 + t289;
t683 = t380 * t464 - t382 * t463;
t804 = m(6) * (t629 + t683);
t297 = t352 * t482;
t628 = t351 * t744 + t297;
t803 = m(6) * (t628 + t683);
t336 = t385 * t481;
t757 = t384 * t507;
t764 = t338 * t508;
t802 = m(6) * (t289 + t336 + (-t757 + t764) * t581);
t680 = t409 * t464 - t411 * t463;
t801 = m(6) * (t629 + t680);
t763 = t351 * t508;
t800 = m(6) * (t297 + t336 + (-t757 + t763) * t581);
t353 = t400 * t481;
t754 = t399 * t507;
t799 = m(6) * (t289 + t353 + (-t754 + t764) * t581);
t798 = m(6) * (t628 + t680);
t797 = m(6) * (t297 + t353 + (-t754 + t763) * t581);
t796 = m(6) * t172;
t186 = t649 * t507 * t508 + t315 * t342;
t794 = m(6) * t186;
t793 = m(6) * t189;
t792 = m(6) * t192;
t791 = m(6) * t193;
t790 = m(6) * t197;
t789 = m(6) * t200;
t618 = Icges(6,5) * t573 + Icges(6,6) * t574;
t456 = t579 * t618;
t457 = t618 * t581;
t781 = (-t576 * t457 + (t579 * t456 + t891) * t581) * t838 + (t577 * t456 + (-t581 * t457 - t891) * t579) * t836;
t749 = t452 * t588;
t748 = t463 * t507;
t747 = t479 * t522;
t548 = Icges(4,5) * t590 - Icges(4,6) * t588;
t736 = t548 * t579;
t690 = t380 * t385 - t382 * t384;
t684 = -t411 * t399 + t409 * t400;
t681 = (-t425 + t442) * t625;
t450 = Icges(4,3) * t581 - t736;
t672 = t581 * t450 + t452 * t720;
t671 = t579 * t450 + t454 * t709;
t647 = t845 / 0.2e1 + t781;
t646 = -t794 + t781;
t643 = t199 * t287 + t682;
t451 = Icges(4,5) * t709 - Icges(4,6) * t710 + Icges(4,3) * t579;
t636 = -t454 * t590 - t451;
t619 = Icges(5,5) * t578 + Icges(5,6) * t580;
t473 = t579 * t619;
t474 = t619 * t581;
t635 = (-t576 * t474 + (t579 * t473 + t890) * t581) * t838 + (t577 * t473 + (-t581 * t474 - t890) * t579) * t836 + t781;
t633 = t649 * t786;
t626 = t194 + t635;
t624 = t809 / 0.2e1 + t634;
t620 = Icges(4,5) * t588 + Icges(4,6) * t590;
t89 = t382 * t744 + t698;
t112 = t411 * t744 + t697;
t256 = t581 * t451 + t453 * t720 - t455 * t719;
t491 = t579 * t620;
t610 = t32 + t864;
t609 = -t32 + (-t573 * t669 + t574 * t667 - t874 * t581 + t746) * t838 + (t502 * t581 - t573 * t670 + t574 * t668 + t874 * t579) * t836;
t608 = -t634 + t892 * (t420 * t574 + t422 * t573);
t607 = t634 + t863;
t602 = t607 + t876;
t600 = t607 + t865;
t599 = t600 + t875;
t598 = t609 - t864 + (-t578 * t665 + t580 * t663 - t873 * t581 + t743) * t838 + (t517 * t581 - t578 * t666 + t873 * t579 + t580 * t664) * t836;
t597 = t608 - t863 + t892 * (t437 * t580 + t439 * t578);
t595 = t598 * qJD(4);
t594 = t597 - t865 + t892 * (t590 * t453 + t588 * t455);
t255 = -t454 * t719 + t672;
t187 = t255 * t581 + t256 * t579;
t257 = -t452 * t710 + t671;
t258 = t579 * t451 - t886;
t188 = t257 * t581 + t258 * t579;
t78 = (t258 + t672 + t886) * t581 + (-t257 + (t636 - t749) * t581 + t256 + t671) * t579;
t79 = (t256 + (-t451 + t749) * t581 - t671) * t581 + (t579 * t636 - t255 + t672) * t579;
t593 = (t598 + t78 * t895 + (t188 + t79) * t837 + (t548 * t581 + t868 * t579 - t588 * t662 + t590 * t660) * t836 + (-t868 * t581 - t588 * t661 + t590 * t659 + t187 + t736) * t838) * qJD(3);
t564 = pkin(3) * t719;
t492 = t620 * t581;
t445 = (-t524 - t584) * t581;
t443 = t564 + t740;
t386 = t480 * t741;
t383 = (t641 - t584) * t581;
t381 = t564 + t410;
t372 = t464 * t481;
t313 = t356 - t633;
t215 = -t633 + t228;
t156 = t634 + t889;
t150 = t634 + t789;
t147 = t797 / 0.2e1;
t146 = t798 / 0.2e1;
t141 = t799 / 0.2e1;
t140 = t800 / 0.2e1;
t139 = t801 / 0.2e1;
t134 = t802 / 0.2e1;
t133 = t803 / 0.2e1;
t130 = t804 / 0.2e1;
t123 = t805 / 0.2e1;
t85 = t607 + t790 + t813;
t81 = t607 + t791 + t814;
t67 = t796 + t818 + t831 + t835;
t63 = t847 / 0.2e1;
t58 = t600 + t792 + t815 + t828;
t57 = t600 + t793 + t816 + t829;
t46 = t849 / 0.2e1;
t45 = -t805 / 0.2e1 + t624;
t44 = t123 + t624;
t43 = t781 + t794;
t42 = t43 * qJD(5);
t41 = m(6) * t112 + t781;
t37 = m(6) * t89 + t781;
t36 = t123 - t809 / 0.2e1 + t608;
t35 = t602 + t703;
t34 = t602 - t703;
t33 = t626 + t841;
t31 = t32 * qJD(5);
t30 = t635 + t877;
t29 = t597 + t703 - t876;
t28 = t63 - t849 / 0.2e1 + t647;
t27 = t46 - t847 / 0.2e1 + t647;
t26 = t46 + t63 - t845 / 0.2e1 + t781;
t25 = t599 + t645;
t24 = t599 - t645;
t23 = t146 - t797 / 0.2e1 + t32;
t22 = t147 - t798 / 0.2e1 + t32;
t21 = t139 - t799 / 0.2e1 + t32;
t20 = t141 - t801 / 0.2e1 + t32;
t19 = t133 - t800 / 0.2e1 + t32;
t18 = t140 - t803 / 0.2e1 + t32;
t17 = t130 - t802 / 0.2e1 + t32;
t16 = t134 - t804 / 0.2e1 + t32;
t15 = t146 + t147 + t609;
t14 = t139 + t141 + t609;
t13 = t594 + t645 - t875;
t12 = t133 + t140 + t609;
t11 = t130 + t134 + t609;
t9 = t610 * qJD(4);
t8 = t610 + t894;
t7 = t610 - t894;
t6 = (t79 / 0.2e1 + t188 / 0.2e1) * t581 + (-t187 / 0.2e1 + t78 / 0.2e1) * t579 + t610;
t5 = t6 * qJD(3);
t4 = t610 + t893;
t3 = t610 - t893;
t2 = t598 + t699 + t700;
t1 = t598 + t701 + t702;
t10 = [qJD(2) * t67 + qJD(3) * t57 + qJD(4) * t81 + qJD(5) * t150, t67 * qJD(1) + t25 * qJD(3) + t35 * qJD(4) + t44 * qJD(5) + 0.2e1 * (t835 / 0.2e1 + t172 * t853 + t190 * t854 + t205 * t855) * qJD(2), t57 * qJD(1) + t25 * qJD(2) + t593 + t1 * qJD(4) + t11 * qJD(5) + ((t378 * t581 + t379 * t579) * t887 + (-t362 * t445 + t363 * t443 + t681) * t854 + (-t338 * t383 + t339 * t381 + t690) * t853) * t858, t81 * qJD(1) + t35 * qJD(2) + t1 * qJD(3) + t595 + t14 * qJD(5) + ((t684 + t694) * t853 + (t321 + t386 + (-t747 + t761) * t581) * t854) * t857, t150 * qJD(1) + t44 * qJD(2) + t11 * qJD(3) + t14 * qJD(4) + ((t289 + t372 + (-t748 + t764) * t581) * m(6) + t609) * qJD(5); t24 * qJD(3) + t34 * qJD(4) + t45 * qJD(5) + (-t796 / 0.4e1 - t818 / 0.4e1 - t835 / 0.4e1 - t831 / 0.4e1) * t861, qJD(3) * t58 + qJD(4) * t85 + qJD(5) * t156, t24 * qJD(1) + t58 * qJD(2) + t593 + t2 * qJD(4) + t12 * qJD(5) + ((t390 * t581 + t391 * t579) * t887 + (-t369 * t445 + t370 * t443 + t681) * t854 + (-t351 * t383 + t352 * t381 + t690) * t853) * t858, t34 * qJD(1) + t85 * qJD(2) + t2 * qJD(3) + t595 + t15 * qJD(5) + ((t684 + t691) * t853 + (t330 + t386 + (-t747 + t760) * t581) * t854) * t857, t45 * qJD(1) + t156 * qJD(2) + t12 * qJD(3) + t15 * qJD(4) + ((t297 + t372 + (-t748 + t763) * t581) * m(6) + t609) * qJD(5); t594 * qJD(1) + t13 * qJD(2) + t5 + t3 * qJD(4) + t17 * qJD(5) + (-t793 / 0.4e1 - t816 / 0.4e1 - t829 / 0.4e1) * t861, t13 * qJD(1) + t594 * qJD(2) + t5 + t7 * qJD(4) + t19 * qJD(5) + (-t792 / 0.4e1 - t815 / 0.4e1 - t828 / 0.4e1) * t859, (m(6) * (t161 * t215 + t380 * t381 - t382 * t383) + m(5) * (t213 * t313 + t442 * t443 - t445 * t625) + (-t576 * t492 + (t866 * t581 + (t491 - t867) * t579) * t581) * t838 + (t577 * t491 + (t867 * t579 + (-t492 - t866) * t581) * t579) * t836 + m(4) * ((t581 * (rSges(4,1) * t709 + rSges(4,3) * t579 - t562) - t579 * (-rSges(4,1) * t719 - t650)) * (-t497 * t579 - t498 * t581) + t649 * t555 * t553) + t635) * qJD(3) + t30 * qJD(4) + t37 * qJD(5) + t648 * t6, t3 * qJD(1) + t7 * qJD(2) + t30 * qJD(3) + t26 * qJD(5) + ((t643 + t88) * t853 + (t135 + t196) * t854) * t857 + (t635 - t194 - t841) * qJD(4), t17 * qJD(1) + t19 * qJD(2) + t37 * qJD(3) + t26 * qJD(4) + (m(6) * (t89 + t186) + t646) * qJD(5); t29 * qJD(2) + t4 * qJD(3) + t9 + t21 * qJD(5) + (-t814 / 0.4e1 - t791 / 0.4e1) * t861 + t597 * qJD(1), t29 * qJD(1) + t597 * qJD(2) + t8 * qJD(3) + t9 + t23 * qJD(5) + (-t790 / 0.4e1 - t813 / 0.4e1) * t859, t4 * qJD(1) + t8 * qJD(2) + t33 * qJD(4) + t27 * qJD(5) + ((t161 * t228 + t199 * t215 + t381 * t409 - t383 * t411 + t686) * t853 + (t313 * t323 + (t443 * t579 - t445 * t581) * t522 + t135) * t854) * t858 + (t635 - t877) * qJD(3), t33 * qJD(3) + (m(6) * t643 + t626) * qJD(4) + t41 * qJD(5) + t648 * t610, t21 * qJD(1) + t23 * qJD(2) + t27 * qJD(3) + t41 * qJD(4) + (m(6) * (t112 + t186) + t646) * qJD(5); (t608 - t789) * qJD(1) + t36 * qJD(2) + t16 * qJD(3) + t20 * qJD(4) + t31, t36 * qJD(1) + (t608 - t889) * qJD(2) + t18 * qJD(3) + t22 * qJD(4) + t31, t16 * qJD(1) + t18 * qJD(2) + ((t215 * t315 + (t381 * t579 - t383 * t581) * t507) * m(6) + t781) * qJD(3) + t28 * qJD(4) + t42, t20 * qJD(1) + t22 * qJD(2) + t28 * qJD(3) + ((t207 - t112 + t601) * m(6) + t781) * qJD(4) + t42, t42 + (qJD(3) + qJD(4)) * t43 + t648 * t32;];
Cq = t10;
