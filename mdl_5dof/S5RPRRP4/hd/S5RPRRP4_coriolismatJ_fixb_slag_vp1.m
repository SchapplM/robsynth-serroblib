% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:52
% EndTime: 2022-01-23 09:32:26
% DurationCPUTime: 23.92s
% Computational Cost: add. (60581->603), mult. (78489->819), div. (0->0), fcn. (84655->8), ass. (0->379)
t899 = Icges(5,4) + Icges(6,4);
t893 = Icges(5,1) + Icges(6,1);
t883 = Icges(5,5) + Icges(6,5);
t892 = Icges(5,2) + Icges(6,2);
t881 = Icges(5,6) + Icges(6,6);
t880 = Icges(5,3) + Icges(6,3);
t553 = qJ(3) + qJ(4);
t546 = sin(t553);
t555 = cos(pkin(8));
t557 = sin(qJ(1));
t685 = t555 * t557;
t547 = cos(t553);
t559 = cos(qJ(1));
t693 = t547 * t559;
t497 = t546 * t685 + t693;
t897 = t899 * t497;
t681 = t559 * t546;
t498 = t547 * t685 - t681;
t896 = t899 * t498;
t895 = t899 * t547;
t894 = t899 * t546;
t554 = sin(pkin(8));
t692 = t554 * t557;
t891 = -t881 * t497 + t883 * t498 + t880 * t692;
t890 = -t892 * t497 + t881 * t692 + t896;
t889 = -t893 * t498 - t883 * t692 + t897;
t887 = t883 * t555 + (-t893 * t547 + t894) * t554;
t861 = -t887 + (-t892 * t547 - t894) * t554;
t888 = t861 * t546;
t821 = -t881 * t555 + (-t892 * t546 + t895) * t554;
t816 = -t547 * t557 + t555 * t681;
t886 = t899 * t816;
t819 = (-t893 * t546 - t895) * t554;
t500 = t546 * t557 + t555 * t693;
t885 = t899 * t500;
t691 = t554 * t559;
t372 = t500 * rSges(5,1) - rSges(5,2) * t816 + rSges(5,3) * t691;
t462 = -rSges(5,3) * t555 + (rSges(5,1) * t547 - rSges(5,2) * t546) * t554;
t619 = t462 * t691;
t278 = t372 * t555 + t619;
t558 = cos(qJ(3));
t550 = t558 * pkin(3);
t560 = pkin(7) + pkin(6);
t503 = t554 * t550 + (pkin(6) - t560) * t555;
t450 = t503 * t691;
t723 = pkin(2) * t555 + pkin(1);
t527 = pkin(6) * t554 + t723;
t709 = qJ(2) * t557;
t496 = t550 * t555 + t554 * t560 + t723;
t556 = sin(qJ(3));
t725 = pkin(3) * t556;
t543 = qJ(2) + t725;
t820 = t496 * t559 + t543 * t557;
t792 = -t527 * t559 - t709 + t820;
t664 = t555 * t792 + t450;
t209 = t278 + t664;
t662 = -t792 - t372;
t210 = t555 * t662 - t450 - t619;
t879 = t209 + t210;
t878 = -t881 * t691 + t892 * t816 - t885;
t877 = t893 * t500 + t883 * t691 - t886;
t860 = -t819 + t821;
t876 = t891 * t691;
t875 = t889 * t500 + t890 * t816;
t716 = Icges(4,4) * t558;
t487 = -Icges(4,6) * t555 + (-Icges(4,2) * t556 + t716) * t554;
t717 = Icges(4,4) * t556;
t488 = -Icges(4,5) * t555 + (Icges(4,1) * t558 - t717) * t554;
t515 = (-Icges(4,2) * t558 - t717) * t554;
t516 = (-Icges(4,1) * t556 - t716) * t554;
t514 = (-Icges(4,5) * t556 - Icges(4,6) * t558) * t554;
t686 = t555 * t514;
t874 = t686 / 0.2e1 - (-(t488 / 0.2e1 + t515 / 0.2e1) * t556 + (t516 / 0.2e1 - t487 / 0.2e1) * t558) * t554;
t552 = t554 ^ 2;
t548 = t559 * qJ(2);
t724 = pkin(4) * t546;
t525 = t724 + t725;
t523 = t559 * t525;
t551 = -qJ(5) - t560;
t795 = t551 * t692 + t523;
t608 = t548 + t795;
t794 = -t496 * t557 + t543 * t559;
t526 = pkin(4) * t547 + t550;
t524 = pkin(2) + t526;
t595 = t524 * t555 + pkin(1);
t797 = -t498 * rSges(6,1) + t497 * rSges(6,2);
t853 = rSges(6,3) * t692 + t557 * t595 - t797;
t670 = -t608 + t794 + t853;
t589 = -t548 + t794;
t823 = t589 - t795 + t853;
t873 = t557 * (t670 - t823);
t345 = t527 * t557 + t589;
t449 = t503 * t692;
t822 = (t551 + t560 - rSges(6,3)) * t555 + (rSges(6,1) * t547 - rSges(6,2) * t546 - pkin(2) + t524 - t550) * t554;
t651 = t822 * t692;
t134 = t449 + (-t345 + t670) * t555 + t651;
t161 = t823 * t555 + t651;
t339 = t555 * t345;
t135 = t449 - t339 + t161;
t678 = t134 - t135;
t160 = t555 * t670 + t651;
t677 = t160 - t161;
t872 = -t875 + t876;
t600 = t691 / 0.2e1;
t602 = t692 / 0.2e1;
t603 = -t692 / 0.2e1;
t859 = t872 - t876;
t807 = (t859 + t875) * t691;
t490 = (-Icges(6,5) * t546 - Icges(6,6) * t547) * t554;
t491 = (-Icges(5,5) * t546 - Icges(5,6) * t547) * t554;
t862 = -t490 - t491;
t826 = t500 * t860 + t691 * t862 + t816 * t861;
t827 = t497 * t861 + t498 * t860 + t692 * t862;
t812 = -t816 * t893 + t878 - t885;
t814 = -t500 * t892 + t877 - t886;
t865 = -t881 * t500 - t883 * t816;
t828 = -t865 * t555 + (-t546 * t814 + t547 * t812) * t554;
t811 = -t497 * t893 - t890 - t896;
t813 = -t498 * t892 - t889 - t897;
t864 = -t883 * t497 - t881 * t498;
t829 = -t864 * t555 + (-t546 * t813 + t547 * t811) * t554;
t851 = (t821 * t816 + t887 * t500 + (t880 * t555 + (t881 * t546 - t883 * t547) * t554) * t691) * t555;
t868 = t878 * t816 + t877 * t500 + (t883 * t500 + t880 * t691 - t881 * t816) * t691;
t839 = t851 + (t872 * t557 + t868 * t559) * t554;
t840 = t851 + ((t891 * t692 + t868) * t559 + t859 * t557) * t554;
t562 = t840 * t603 - t807 * t691 / 0.2e1 + (-t826 + t828) * t600 + (-t827 + t829 + t839) * t602;
t687 = t555 * t491;
t688 = t555 * t490;
t674 = t687 + t688 + (t860 * t547 + t888) * t554;
t809 = t674 * t555;
t871 = t562 + t809;
t579 = rSges(6,3) * t554 + t595;
t583 = (qJ(2) + t525) * t557 - t551 * t691;
t584 = rSges(6,1) * t500 - rSges(6,2) * t816;
t260 = t559 * t579 + t583 + t584;
t818 = -t557 * t579 + t608 + t797;
t870 = m(6) * (t260 * t691 - t692 * t818);
t371 = rSges(6,3) * t691 + t584;
t610 = t555 * t371 + t822 * t691;
t793 = t559 * t595 + t583 - t820;
t162 = t555 * t793 + t610;
t869 = t162 * t557;
t368 = t498 * rSges(5,1) - t497 * rSges(5,2) + rSges(5,3) * t692;
t433 = t462 * t692;
t643 = t433 + t449;
t207 = (-t345 + t368) * t555 + t643;
t342 = t555 * t368;
t208 = -t339 + t342 + t643;
t867 = t207 - t208;
t863 = t433 + t342;
t682 = t558 * t559;
t684 = t556 * t557;
t517 = t555 * t684 + t682;
t680 = t559 * t556;
t683 = t557 * t558;
t518 = t555 * t683 - t680;
t796 = -t518 * rSges(4,1) + t517 * rSges(4,2);
t409 = rSges(4,3) * t692 - t796;
t489 = -rSges(4,3) * t555 + (rSges(4,1) * t558 - rSges(4,2) * t556) * t554;
t858 = t409 * t555 + t489 * t692;
t857 = -t554 / 0.2e1;
t855 = qJD(1) * t870;
t852 = -t368 + t794;
t507 = Icges(4,4) * t518;
t399 = -Icges(4,2) * t517 + Icges(4,6) * t692 + t507;
t506 = Icges(4,4) * t517;
t403 = -Icges(4,1) * t518 - Icges(4,5) * t692 + t506;
t519 = -t555 * t680 + t683;
t520 = t555 * t682 + t684;
t843 = -t519 * t399 + t403 * t520;
t776 = m(4) / 0.2e1;
t775 = m(5) / 0.2e1;
t774 = m(6) / 0.2e1;
t396 = Icges(4,5) * t518 - Icges(4,6) * t517 + Icges(4,3) * t692;
t831 = t396 * t691;
t407 = -rSges(6,1) * t816 - t500 * rSges(6,2);
t645 = t816 * pkin(4) - t407;
t830 = t557 * t645;
t594 = rSges(4,3) * t554 + t527;
t817 = -t557 * t594 + t548 + t796;
t122 = t208 * t559 + t210 * t557;
t587 = rSges(4,1) * t520 + rSges(4,2) * t519;
t411 = rSges(4,3) * t691 + t587;
t295 = t411 * t555 + t489 * t691;
t669 = -t793 - t792;
t137 = t555 * t669 - t450 - t610;
t76 = t135 * t559 + t137 * t557;
t623 = t122 * t775 + (-t295 * t557 + t559 * t858) * t776 + t76 * t774;
t136 = t162 + t664;
t628 = (-t134 * t559 + t136 * t557 + t76) * t774 + (-t207 * t559 + t209 * t557 + t122) * t775;
t815 = t623 - t628;
t105 = t161 * t559 - t869;
t679 = t105 * t774 + (-t278 * t557 + t559 * t863) * t775;
t719 = (-t160 * t559 + t105 + t869) * t774;
t810 = t679 - t719;
t571 = -t525 * t685 - t526 * t559;
t633 = t517 * pkin(3);
t378 = t571 + t633;
t466 = (-t525 + t725) * t554;
t501 = (-rSges(6,1) * t546 - rSges(6,2) * t547) * t554;
t585 = rSges(6,1) * t497 + rSges(6,2) * t498;
t650 = t501 * t692 - t555 * t585;
t222 = t555 * t378 + t466 * t692 + t650;
t590 = (-t466 - t501) * t691;
t292 = -t555 * t523 + t557 * t526 + t407;
t510 = t519 * pkin(3);
t652 = t510 - t292;
t223 = t555 * t652 + t590;
t291 = -t571 + t585;
t406 = -rSges(5,1) * t497 - rSges(5,2) * t498;
t320 = -t406 + t633;
t408 = -rSges(5,1) * t816 - t500 * rSges(5,2);
t644 = -t408 - t510;
t286 = t372 + t820;
t502 = (-rSges(5,1) * t546 - rSges(5,2) * t547) * t554;
t296 = t555 * t406 + t502 * t692;
t618 = t502 * t691;
t297 = -t408 * t555 - t618;
t675 = t297 * t286 + t296 * t852;
t721 = (t161 * t291 - t162 * t292 + t222 * t818 + t223 * t260) * t774 + (t278 * t644 + t320 * t863 + t675) * t775;
t473 = t497 * pkin(4);
t313 = t585 + t473;
t629 = t552 * t724;
t251 = -t473 * t555 - t557 * t629 + t650;
t252 = (-t501 * t554 + t629) * t559 + t645 * t555;
t676 = t251 * t818 + t252 * t260;
t722 = (t135 * t313 - t137 * t645 + t676) * t774 + (-t208 * t406 + t210 * t408 + t675) * t775;
t808 = t721 - t722;
t806 = -t687 / 0.2e1 - t688 / 0.2e1 + t857 * t888;
t805 = t691 * t865 + t692 * t864;
t804 = t554 / 0.2e1;
t752 = -t555 / 0.2e1;
t803 = m(6) * t554;
t202 = t831 - t843;
t798 = t202 - t831;
t673 = (t313 * t557 + t559 * t645) * t774 + (-t406 * t557 - t408 * t559) * t775;
t328 = t345 * t691;
t668 = -t793 - t371;
t671 = t823 * t691;
t108 = -t328 + (-t792 + t668) * t692 + t671;
t138 = t668 * t692 + t671;
t269 = t793 * t691;
t67 = t269 + t792 * t691 + (t669 * t559 + t873) * t554;
t98 = t269 + (-t559 * t793 + t873) * t554;
t786 = (-t278 * t867 + t863 * t879) * t775 + (t108 * t98 + t136 * t161 + t137 * t160 + t138 * t67 - t678 * t162) * t774;
t421 = -rSges(4,1) * t517 - rSges(4,2) * t518;
t422 = rSges(4,1) * t519 - rSges(4,2) * t520;
t611 = (t291 * t557 - t292 * t559) * t774 + (t320 * t557 + t559 * t644) * t775 + (-t421 * t557 - t422 * t559) * t776;
t782 = 0.2e1 * qJD(1);
t781 = 0.4e1 * qJD(1);
t780 = 2 * qJD(3);
t779 = 4 * qJD(3);
t778 = 2 * qJD(4);
t777 = 4 * qJD(4);
t770 = m(5) * (t207 * t210 + t208 * t209);
t337 = t368 * t691;
t153 = t662 * t692 - t328 + t337;
t377 = t406 * t691;
t264 = -t408 * t692 + t377;
t63 = t153 * t264 + t208 * t296 + t210 * t297;
t769 = m(5) * t63;
t250 = -t372 * t692 + t337;
t96 = t250 * t264 - t278 * t297 + t296 * t863;
t95 = m(5) * t96;
t764 = m(6) * (t108 * t67 + t134 * t137 + t135 * t136);
t130 = t137 * t691;
t763 = m(6) * (-t555 * t67 + t130 + (t136 * t559 + t557 * t678) * t554);
t762 = m(6) * (t138 * t98 - t677 * t162);
t154 = t162 * t691;
t761 = m(6) * (-t555 * t98 - t154 + (t162 * t559 + t557 * t677) * t554);
t376 = t585 * t691;
t224 = -t376 + (-t473 * t559 + t830) * t554;
t30 = t108 * t224 + t135 * t251 + t137 * t252;
t760 = m(6) * t30;
t665 = t378 * t691 - t376;
t169 = t652 * t692 + t665;
t759 = m(6) * (t138 * t169 + t161 * t222 - t162 * t223);
t754 = m(6) * (-t135 * t692 + t130);
t751 = m(3) * ((rSges(3,2) * t692 + rSges(3,3) * t559 + t548) * t559 + (-rSges(3,2) * t691 + (rSges(3,3) + qJ(2)) * t557) * t557);
t317 = t559 * t594 + t587 + t709;
t749 = m(4) * (t317 * t422 - t421 * t817);
t747 = m(4) * (t317 * t557 + t559 * t817);
t743 = m(5) * (-t286 * t644 + t320 * t852);
t742 = m(5) * (t286 * t408 - t406 * t852);
t740 = m(5) * (t286 * t557 + t559 * t852);
t737 = m(6) * (-t161 * t692 - t154);
t735 = m(6) * (-t224 * t555 + (t251 * t559 + t252 * t557) * t554);
t734 = m(6) * (t260 * t292 + t291 * t818);
t733 = m(6) * (-t260 * t645 + t313 * t818);
t155 = t251 * t557 - t252 * t559;
t732 = m(6) * t155;
t731 = m(6) * (t260 * t557 + t559 * t818);
t730 = (t291 * t559 + t292 * t557) * t803;
t728 = (t313 * t559 - t830) * t803;
t718 = Icges(4,4) * t520;
t705 = t222 * t557;
t704 = t223 * t559;
t695 = t547 * t554;
t649 = -Icges(4,1) * t517 - t399 - t507;
t401 = Icges(4,2) * t519 + Icges(4,6) * t691 + t718;
t648 = Icges(4,1) * t519 - t401 - t718;
t647 = -Icges(4,2) * t518 - t403 - t506;
t508 = Icges(4,4) * t519;
t404 = Icges(4,1) * t520 + Icges(4,5) * t691 + t508;
t646 = -Icges(4,2) * t520 + t404 + t508;
t637 = -t487 + t516;
t636 = t488 + t515;
t414 = (-t557 ^ 2 - t559 ^ 2) * t803;
t631 = t414 * qJD(1);
t630 = t552 * t725;
t624 = t735 / 0.2e1;
t614 = (t798 + t843) * t804 * t559 + (t752 + t555 / 0.2e1) * ((-Icges(4,3) * t555 + (Icges(4,5) * t558 - Icges(4,6) * t556) * t554) * t692 - t487 * t517 + t488 * t518);
t203 = (Icges(4,5) * t520 + Icges(4,6) * t519 + Icges(4,3) * t691) * t691 + t519 * t401 + t520 * t404;
t613 = ((t396 * t692 + t203) * t559 + t798 * t557) * t804 + (t202 * t557 + t203 * t559) * t857;
t612 = t138 * t224 + t161 * t251 - t162 * t252;
t609 = -t510 + t652;
t605 = -t695 / 0.2e1;
t604 = t695 / 0.2e1;
t582 = ((t829 * t557 + t828 * t559) * t554 + t809) * t752 + ((-t814 * t497 + t812 * t498) * t691 + t827 * t555 + (-t813 * t497 + t811 * t498 + t805) * t692) * t602 + ((t811 * t500 - t813 * t816) * t692 + t826 * t555 + (t812 * t500 - t814 * t816 + t805) * t691) * t600;
t564 = m(6) * (-t555 * t169 + (t222 * t559 + t223 * t557) * t554);
t55 = t624 - t564 / 0.2e1;
t81 = 0.2e1 * (t155 / 0.4e1 - t705 / 0.4e1 + t704 / 0.4e1) * m(6);
t581 = t81 * qJD(2) + t55 * qJD(5);
t580 = -t555 * t633 - t557 * t630;
t578 = t95 + t582;
t568 = t807 * t600 + t840 * t602 + t839 * t603;
t567 = t819 * t604 + t821 * t605 + t806;
t566 = t568 + t786;
t561 = t821 * t604 + t819 * t605 - t806;
t531 = t559 * t630;
t521 = (-rSges(4,1) * t556 - rSges(4,2) * t558) * t554;
t463 = t633 * t691;
t416 = Icges(4,5) * t519 - Icges(4,6) * t520;
t415 = -Icges(4,5) * t517 - Icges(4,6) * t518;
t316 = -t422 * t555 - t521 * t691;
t315 = t421 * t555 + t521 * t692;
t257 = t555 * t644 + t531 - t618;
t256 = t580 + t296;
t248 = -t686 + (-t556 * t636 + t558 * t637) * t554;
t229 = t728 / 0.2e1;
t226 = t644 * t692 + t377 - t463;
t218 = m(5) * (t296 * t557 - t297 * t559);
t214 = t514 * t691 + t519 * t636 + t520 * t637;
t213 = t514 * t692 - t517 * t636 + t518 * t637;
t206 = t730 / 0.2e1;
t188 = t555 * t609 + t531 + t590;
t187 = t580 + t222;
t159 = -t416 * t555 + (-t556 * t646 + t558 * t648) * t554;
t158 = -t415 * t555 + (-t556 * t647 + t558 * t649) * t554;
t148 = t609 * t692 - t463 + t665;
t99 = t737 / 0.2e1;
t78 = t731 + t740 + t747 + t751;
t72 = t754 / 0.2e1;
t66 = t218 + t732 / 0.2e1 + (-t704 + t705) * t774;
t58 = t567 + t733 + t742;
t54 = t624 + t564 / 0.2e1;
t43 = t567 + t749 + t743 + t734 - t874;
t28 = t761 / 0.2e1;
t23 = t763 / 0.2e1;
t22 = t99 + t229 - t761 / 0.2e1;
t21 = t99 + t28 - t728 / 0.2e1;
t20 = t229 + t28 - t737 / 0.2e1;
t19 = t679 + t719 - t673;
t18 = t673 - t810;
t17 = t673 + t810;
t15 = t72 + t206 - t763 / 0.2e1;
t14 = t72 + t23 - t730 / 0.2e1;
t13 = t206 + t23 - t754 / 0.2e1;
t12 = t623 + t628 - t611;
t11 = t611 - t815;
t10 = t611 + t815;
t7 = t578 + t759;
t6 = t582 + t760 + t769;
t5 = t568 + t762;
t4 = t770 + t764 + (t557 * t613 + t559 * t614) * t554 + t568;
t3 = t566 + t808;
t2 = t566 - t808;
t1 = t721 + t722 - t786 + t871;
t8 = [t78 * qJD(2) + t43 * qJD(3) + t58 * qJD(4) + qJD(5) * t870, qJD(1) * t78 + qJD(3) * t10 + qJD(4) * t17, t43 * qJD(1) + t10 * qJD(2) + t1 * qJD(4) + t15 * qJD(5) + (-t770 / 0.4e1 - t764 / 0.4e1) * t779 + ((t135 * t291 + t137 * t292 + t187 * t818 + t188 * t260) * t774 + (t208 * t320 - t210 * t644 + t256 * t852 + t257 * t286) * t775 + (-t295 * t422 + t315 * t817 + t316 * t317 - t421 * t858) * t776) * t780 + (t562 + (-t248 + t674) * t555 + ((t159 / 0.2e1 + t214 / 0.2e1 - t614) * t559 + (t158 / 0.2e1 + t213 / 0.2e1 - t613) * t557) * t554) * qJD(3), t58 * qJD(1) + t17 * qJD(2) + t1 * qJD(3) + t871 * qJD(4) + t22 * qJD(5) - t762 * t777 / 0.4e1 + ((t161 * t313 + t162 * t645 + t676) * t774 + (-t278 * t408 - t406 * t863 + t675) * t775) * t778, t15 * qJD(3) + t22 * qJD(4) + t855; t11 * qJD(3) + t18 * qJD(4) + t414 * qJD(5) + (-t731 / 0.4e1 - t740 / 0.4e1 - t747 / 0.4e1 - t751 / 0.4e1) * t781, 0, t11 * qJD(1) + t66 * qJD(4) + ((t315 * t557 - t316 * t559) * t776 + (t256 * t557 - t257 * t559) * t775 + (t187 * t557 - t188 * t559) * t774) * t780, t18 * qJD(1) + t66 * qJD(3) + (t218 + t732) * qJD(4), t631; t12 * qJD(2) + t4 * qJD(3) + t2 * qJD(4) + t14 * qJD(5) + (-t749 / 0.4e1 - t743 / 0.4e1 - t734 / 0.4e1) * t781 + ((t286 * t867 + t852 * t879) * t775 + (t260 * t678 + (t136 + t137) * t818) * t774) * t782 + (t561 + t874) * qJD(1), qJD(1) * t12 + qJD(4) * t81, t4 * qJD(1) + (((t416 * t692 - t517 * t646 + t518 * t648) * t691 + (t415 * t692 - t517 * t647 + t518 * t649) * t692 - t213 * t555) * t602 + ((t416 * t691 + t519 * t646 + t520 * t648) * t691 + (t415 * t691 + t519 * t647 + t520 * t649) * t692 - t214 * t555) * t600 + m(6) * (t108 * t148 + t135 * t187 + t137 * t188) + m(5) * (t153 * t226 + t208 * t256 + t210 * t257) + (-t248 * t555 + (t158 * t557 + t159 * t559) * t554) * t752 + m(4) * (t858 * t315 - t295 * t316 + (t409 * t559 - t411 * t557) * t552 * (t421 * t559 - t422 * t557)) + t582) * qJD(3) + t6 * qJD(4), t2 * qJD(1) + t6 * qJD(3) + t582 * qJD(4) + (-t95 / 0.4e1 - t759 / 0.4e1) * t777 + ((t612 + t30) * t774 + (t96 + t63) * t775) * t778 + t581, qJD(1) * t14 + qJD(4) * t55; t561 * qJD(1) + t19 * qJD(2) + t3 * qJD(3) + t5 * qJD(4) + t21 * qJD(5) + (-t733 / 0.4e1 - t742 / 0.4e1) * t781 + t260 * t677 * t774 * t782, qJD(1) * t19 - qJD(3) * t81, t3 * qJD(1) + t582 * qJD(3) + t7 * qJD(4) + (-t760 / 0.4e1 - t769 / 0.4e1) * t779 + ((t108 * t169 + t135 * t222 + t137 * t223 + t138 * t148 + t161 * t187 - t162 * t188) * t774 + (t226 * t250 + t256 * t863 - t257 * t278 + t63) * t775) * t780 - t581, t5 * qJD(1) + t7 * qJD(3) + (m(6) * t612 + t578) * qJD(4), qJD(1) * t21 - qJD(3) * t55; -t414 * qJD(2) + t13 * qJD(3) + t20 * qJD(4) - t855, -t631, t13 * qJD(1) + m(6) * (-t148 * t555 + (t187 * t559 + t188 * t557) * t554) * qJD(3) + t54 * qJD(4), t20 * qJD(1) + t54 * qJD(3) + qJD(4) * t735, 0;];
Cq = t8;
