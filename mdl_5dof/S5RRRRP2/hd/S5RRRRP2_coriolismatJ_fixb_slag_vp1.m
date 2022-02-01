% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:46
% EndTime: 2022-01-20 11:49:07
% DurationCPUTime: 15.18s
% Computational Cost: add. (62785->477), mult. (46522->615), div. (0->0), fcn. (42674->8), ass. (0->319)
t880 = Icges(5,4) + Icges(6,4);
t878 = Icges(5,1) + Icges(6,1);
t877 = Icges(5,5) + Icges(6,5);
t872 = Icges(5,6) + Icges(6,6);
t541 = qJ(3) + qJ(4);
t534 = sin(t541);
t881 = t880 * t534;
t876 = Icges(5,2) + Icges(6,2);
t536 = cos(t541);
t879 = -t872 * t534 + t877 * t536;
t856 = t878 * t536 - t881;
t875 = Icges(5,3) + Icges(6,3);
t874 = t880 * t536;
t542 = qJ(1) + qJ(2);
t535 = sin(t542);
t537 = cos(t542);
t866 = t535 * t877 + t537 * t856;
t871 = -t534 * t876 + t874;
t683 = t534 * t535;
t870 = t880 * t683;
t869 = t879 * t537;
t545 = cos(qJ(3));
t539 = t545 * pkin(3);
t531 = t539 + pkin(2);
t730 = pkin(4) * t536;
t479 = t531 + t730;
t669 = t536 * t537;
t682 = t534 * t537;
t793 = -pkin(8) - pkin(7);
t864 = rSges(6,3) + qJ(5) - t793;
t312 = rSges(6,1) * t669 - rSges(6,2) * t682 + t537 * t479 + t864 * t535;
t723 = rSges(6,1) * t536;
t824 = rSges(6,2) * t683 + t864 * t537;
t812 = (-t479 - t723) * t535 + t824;
t177 = t312 * t535 + t537 * t812;
t868 = t177 * m(6) * qJD(2);
t676 = t535 * t536;
t858 = -t875 * t537 + t877 * t676 - t872 * t683;
t867 = t875 * t535 + t869;
t842 = -t877 * t537 + t878 * t676 - t870;
t857 = t876 * t536 + t881;
t865 = t878 * t534 + t874;
t733 = cos(qJ(1)) * pkin(1);
t302 = t733 + t312;
t863 = -t302 + t312;
t862 = t866 * t676;
t861 = -t872 * t537 + t676 * t880 - t683 * t876;
t860 = t872 * t535 + t871 * t537;
t734 = sin(qJ(1)) * pkin(1);
t817 = t812 - t734;
t859 = -t812 + t817;
t855 = t877 * t534 + t872 * t536;
t854 = -t867 * t537 + t862;
t853 = -t857 * t537 + t866;
t852 = -t676 * t876 + t842 - t870;
t851 = -t865 * t537 - t860;
t850 = t865 * t535 + t861;
t849 = -t858 * t535 - t842 * t669;
t837 = t867 * t535 + t866 * t669;
t848 = -t865 - t871;
t847 = m(6) * qJD(5);
t174 = t302 * t535 + t537 * t817;
t846 = t174 * m(6) * qJD(1);
t845 = -t860 * t683 + t854;
t844 = t861 * t682 + t849;
t843 = -t860 * t682 + t837;
t538 = Icges(4,4) * t545;
t543 = sin(qJ(3));
t497 = -Icges(4,2) * t543 + t538;
t498 = Icges(4,1) * t543 + t538;
t839 = t497 + t498;
t838 = t860 * t534 - t858;
t836 = t861 * t534;
t141 = t302 * t812 - t312 * t817;
t796 = m(4) / 0.2e1;
t795 = m(5) / 0.2e1;
t794 = m(6) / 0.2e1;
t771 = t535 / 0.2e1;
t770 = -t537 / 0.2e1;
t832 = t537 / 0.2e1;
t768 = m(3) * (t733 * (-rSges(3,1) * t535 - rSges(3,2) * t537) + (t537 * rSges(3,1) - t535 * rSges(3,2)) * t734);
t826 = t855 * t535;
t825 = t855 * t537;
t532 = t535 ^ 2;
t533 = t537 ^ 2;
t615 = t532 + t533;
t823 = (t856 - t857) * t536 + t848 * t534;
t474 = rSges(5,1) * t534 + rSges(5,2) * t536;
t732 = pkin(3) * t543;
t565 = t474 + t732;
t813 = t565 * t537;
t814 = t565 * t535;
t568 = (t535 * t813 - t537 * t814) * t474;
t395 = rSges(5,1) * t676 - rSges(5,2) * t683 - t537 * rSges(5,3);
t619 = t535 * t531 + t537 * t793;
t324 = -t395 - t619;
t512 = t535 * t793;
t603 = -rSges(5,2) * t682 + t535 * rSges(5,3);
t724 = rSges(5,1) * t536;
t325 = -t512 + (t531 + t724) * t537 + t603;
t477 = -rSges(5,2) * t534 + t724;
t569 = (-t324 * t537 - t325 * t535) * t477;
t731 = pkin(4) * t534;
t480 = -t731 - t732;
t617 = rSges(6,1) * t683 + rSges(6,2) * t676;
t341 = -t480 * t535 + t617;
t454 = t537 * t480;
t722 = rSges(6,2) * t536;
t473 = rSges(6,1) * t534 + t722;
t342 = -t473 * t537 + t454;
t606 = t473 + t731;
t378 = t606 * t535;
t380 = t606 * t537;
t645 = -t380 * t341 - t378 * t342;
t605 = rSges(6,2) * t534 - t723 - t730;
t379 = t605 * t535;
t381 = t605 * t537;
t652 = t379 * t312 + t381 * t812;
t713 = (t569 + t568) * t795 + (t645 + t652) * t794;
t425 = t474 * t535;
t426 = t474 * t537;
t641 = -t425 * t813 + t426 * t814;
t573 = t473 - t480;
t337 = t573 * t535;
t339 = t573 * t537;
t365 = pkin(4) * t683 + t617;
t366 = (-t722 + (-rSges(6,1) - pkin(4)) * t534) * t537;
t647 = -t337 * t366 - t339 * t365;
t714 = (t569 + t641) * t795 + (t647 + t652) * t794;
t821 = t713 - t714;
t314 = t324 - t734;
t315 = t325 + t733;
t570 = (-t314 * t537 - t315 * t535) * t477;
t653 = t379 * t302 + t381 * t817;
t715 = (t570 + t568) * t795 + (t645 + t653) * t794;
t716 = (t570 + t641) * t795 + (t647 + t653) * t794;
t820 = t715 - t716;
t772 = -t535 / 0.2e1;
t819 = t771 + t772;
t818 = (t851 * t535 + t850 * t537) * t536 + (-t853 * t535 + t852 * t537) * t534;
t159 = t314 * t814 - t315 * t813;
t160 = t324 * t814 - t325 * t813;
t149 = -t325 * t314 + t315 * t324;
t530 = t537 * pkin(7);
t725 = rSges(4,1) * t545;
t608 = pkin(2) + t725;
t675 = t535 * t543;
t616 = rSges(4,2) * t675 + t537 * rSges(4,3);
t343 = -t535 * t608 + t530 + t616;
t329 = t343 - t734;
t663 = t537 * t543;
t508 = rSges(4,2) * t663;
t344 = -t508 + t608 * t537 + (rSges(4,3) + pkin(7)) * t535;
t330 = t344 + t733;
t161 = -t344 * t329 + t330 * t343;
t158 = (-t537 * t531 + t312 + t512) * t537 + (rSges(6,1) * t676 + t479 * t535 - t619 - t824) * t535;
t635 = -t535 * (pkin(2) * t535 - t530 - t619) + t537 * (-t535 * pkin(7) - t512 + (-pkin(2) + t531) * t537);
t120 = t158 + t635;
t624 = -t533 * t473 - t535 * t617;
t247 = -t615 * t731 + t624;
t646 = -t337 * t379 - t339 * t381;
t55 = t120 * t247 + t646;
t278 = t535 * t395 + t537 * (rSges(5,1) * t669 + t603);
t164 = t278 + t635;
t308 = -t535 * t425 - t537 * t426;
t92 = t164 * t308 + (t535 * t814 + t537 * t813) * t477;
t811 = -m(5) * t92 - m(6) * t55;
t614 = qJD(1) + qJD(2);
t717 = ((-t315 + t325) * t537 + (t314 - t324) * t535) * t474 * t795 + (t859 * t378 + t863 * t380) * t794;
t156 = t302 * t366 + t365 * t817;
t157 = t312 * t366 + t365 * t812;
t163 = t314 * t425 - t315 * t426;
t166 = t324 * t425 - t325 * t426;
t809 = (t157 + t156) * t794 + (t166 + t163) * t795;
t500 = rSges(4,1) * t543 + rSges(4,2) * t545;
t612 = ((-t330 + t344) * t537 + (t329 - t343) * t535) * t500 * t796 + (t859 * t337 + t863 * t339) * t794 + (t159 - t160) * t795;
t151 = t302 * t342 + t341 * t817;
t155 = t312 * t342 + t341 * t812;
t449 = t500 * t535;
t450 = t500 * t537;
t176 = t329 * t449 - t330 * t450;
t185 = t343 * t449 - t344 * t450;
t807 = (t185 + t176) * t796 + (t155 + t151) * t794 + (t160 + t159) * t795;
t706 = Icges(4,4) * t543;
t496 = Icges(4,2) * t545 + t706;
t499 = Icges(4,1) * t545 - t706;
t806 = t839 * t545 / 0.2e1 + (t499 / 0.2e1 - t496 / 0.2e1) * t543;
t564 = (t843 * t535 + t844 * t537) * t832 + (t837 * t535 + ((t836 + t867) * t537 + t845 + t849 - t862) * t537) * t770 + ((t838 * t535 - t844 + t845 - t854) * t535 + ((t838 + t858) * t537 + (-t842 * t536 + t836) * t535 - t837 + t843) * t537) * t771;
t562 = -t848 * t536 / 0.2e1 + (-t857 / 0.2e1 + t856 / 0.2e1) * t534;
t802 = 0.4e1 * qJD(1);
t800 = 0.4e1 * qJD(2);
t799 = 2 * qJD(3);
t798 = 2 * qJD(4);
t182 = t532 * (t480 + t732) + t537 * (pkin(3) * t663 + t454) + t624;
t643 = -t378 * t379 - t380 * t381;
t782 = m(6) * (t158 * t182 + t643);
t409 = Icges(4,6) * t535 + t497 * t537;
t411 = Icges(4,5) * t535 + t499 * t537;
t674 = t535 * t545;
t362 = t411 * t674;
t495 = Icges(4,5) * t545 - Icges(4,6) * t543;
t688 = t495 * t537;
t407 = Icges(4,3) * t535 + t688;
t595 = t407 * t537 - t362;
t207 = -t409 * t675 - t595;
t406 = Icges(4,5) * t674 - Icges(4,6) * t675 - Icges(4,3) * t537;
t506 = Icges(4,4) * t675;
t410 = Icges(4,1) * t674 - Icges(4,5) * t537 - t506;
t408 = Icges(4,4) * t674 - Icges(4,2) * t675 - Icges(4,6) * t537;
t693 = t408 * t543;
t147 = -(-t535 * (-t410 * t545 + t693) - t406 * t537) * t537 + t207 * t535;
t662 = t537 * t545;
t634 = -t535 * t406 - t410 * t662;
t208 = -t408 * t663 - t634;
t633 = t535 * t407 + t411 * t662;
t209 = -t409 * t663 + t633;
t148 = -t208 * t537 + t209 * t535;
t592 = t409 * t543 - t406;
t46 = (t537 * t592 + t209 - t633) * t537 + (t535 * t592 + t208 + t595) * t535;
t47 = (t207 - t362 + (t407 + t693) * t537 + t634) * t537 + t633 * t535;
t8 = (-t47 / 0.2e1 + t148 / 0.2e1) * t537 + (t147 / 0.2e1 + t46 / 0.2e1) * t535 + t564;
t738 = m(6) * (-t535 * t337 - t339 * t537);
t218 = t738 / 0.2e1;
t739 = m(6) * (t341 * t535 - t342 * t537);
t99 = t218 - t739 / 0.2e1;
t769 = t8 * qJD(3) + t99 * qJD(5);
t764 = m(4) * t161;
t762 = m(4) * t176;
t761 = m(4) * t185;
t755 = m(5) * t149;
t154 = t615 * t474 * t477 + t278 * t308;
t152 = m(5) * t154;
t753 = m(5) * t159;
t752 = m(5) * t160;
t751 = m(5) * t163;
t750 = m(5) * t166;
t749 = m(6) * (-t177 + t174);
t748 = m(6) * (t177 + t174);
t746 = m(6) * t141;
t744 = m(6) * t151;
t743 = m(6) * t155;
t742 = m(6) * t156;
t741 = m(6) * t157;
t737 = m(6) * (t365 * t535 - t366 * t537);
t735 = m(6) * (-t535 * t378 - t380 * t537);
t257 = t735 / 0.2e1;
t139 = t257 - t737 / 0.2e1;
t729 = qJD(4) * t564 + t139 * qJD(5);
t726 = m(6) * qJD(3);
t238 = t737 / 0.2e1;
t137 = t238 - t735 / 0.2e1;
t217 = t739 / 0.2e1;
t98 = t217 - t738 / 0.2e1;
t712 = t98 * qJD(3) + t137 * qJD(4);
t100 = t218 + t217;
t138 = t257 + t238;
t711 = t100 * qJD(3) + t138 * qJD(4);
t649 = -t337 * t342 - t339 * t341;
t644 = -t380 * t365 - t378 * t366;
t623 = t498 * t535 + t408;
t622 = -t498 * t537 - t409;
t621 = -Icges(4,2) * t674 + t410 - t506;
t620 = -t496 * t537 + t411;
t610 = t158 * t247 + t643;
t607 = -t477 - t539;
t604 = ((t826 * t535 + t818) * t537 - t825 * t532) * t771 + ((t825 * t537 + t818) * t535 - t826 * t533) * t770;
t590 = t615 * t732;
t589 = t152 + t604;
t583 = Icges(4,5) * t543 + Icges(4,6) * t545;
t572 = t605 - t539;
t567 = (-t425 * t537 + t426 * t535) * t474;
t566 = (-t449 * t537 + t450 * t535) * t500;
t561 = t543 * t621 + t545 * t623;
t560 = -t543 * t620 + t545 * t622;
t557 = (-t496 + t499) * t545 - t839 * t543;
t555 = t562 + t809;
t554 = t562 + t806;
t553 = t554 + t807;
t552 = -t564 + (t851 * t534 + t535 * t879 + t853 * t536 + t823 * t537) * t771 + (-t850 * t534 + t823 * t535 + t852 * t536 - t869) * t770;
t551 = -t562 + (t842 * t534 + t861 * t536) * t819;
t549 = t552 * qJD(4) + t138 * qJD(5);
t548 = t551 - t806 + t819 * (t408 * t545 + t410 * t543);
t547 = t100 * qJD(5) + (t47 * t832 + t552 + (t535 * t557 - t543 * t623 + t545 * t621 + t148 - t688) * t770 + (t46 + t147) * t772 + (t495 * t535 + t537 * t557 + t543 * t622 + t545 * t620) * t771) * qJD(3);
t502 = -rSges(4,2) * t543 + t725;
t444 = t537 * t583;
t443 = t583 * t535;
t399 = t607 * t537;
t397 = t607 * t535;
t340 = t572 * t537;
t338 = t572 * t535;
t270 = -t379 * t537 + t381 * t535;
t262 = -t590 + t308;
t249 = m(6) * t270 * qJD(4);
t168 = -t590 + t182;
t105 = t748 / 0.2e1;
t104 = t749 / 0.2e1;
t54 = t562 + t741 + t750;
t53 = t562 + t742 + t751;
t41 = t746 + t755 + t764 + t768;
t40 = t554 + t743 + t752 + t761;
t27 = t554 + t744 + t753 + t762;
t24 = t105 - t749 / 0.2e1;
t23 = t105 + t104;
t22 = t104 - t748 / 0.2e1;
t18 = t555 - t717;
t17 = t555 + t717;
t16 = t589 + t782;
t15 = t604 - t811;
t14 = t551 + t717 - t809;
t13 = t553 - t612;
t12 = t553 + t612;
t11 = t548 + t612 - t807;
t6 = t564 + t821;
t5 = t564 - t821;
t4 = t564 + t820;
t3 = t564 - t820;
t2 = t552 + t713 + t714;
t1 = t552 + t715 + t716;
t7 = [t41 * qJD(2) + t27 * qJD(3) + t53 * qJD(4) + t174 * t847, t41 * qJD(1) + t12 * qJD(3) + t17 * qJD(4) + t23 * qJD(5) + 0.2e1 * (t141 * t794 + t149 * t795 + t161 * t796 + t768 / 0.2e1) * qJD(2), t27 * qJD(1) + t12 * qJD(2) + t1 * qJD(4) + (((-t329 * t537 - t330 * t535) * t502 + t566) * t796 + (t302 * t338 + t340 * t817 + t649) * t794 + (t314 * t399 + t315 * t397) * t795) * t799 + t547, t53 * qJD(1) + t17 * qJD(2) + t1 * qJD(3) + ((t644 + t653) * t794 + (t570 + t567) * t795) * t798 + t549, t23 * qJD(2) + t711 + t846; t13 * qJD(3) + t18 * qJD(4) + t24 * qJD(5) + (-t746 / 0.4e1 - t755 / 0.4e1 - t764 / 0.4e1 - t768 / 0.4e1) * t802, t40 * qJD(3) + t54 * qJD(4) + t177 * t847, t13 * qJD(1) + t40 * qJD(2) + t2 * qJD(4) + (((-t343 * t537 - t344 * t535) * t502 + t566) * t796 + (t312 * t338 + t340 * t812 + t649) * t794 + (t324 * t399 + t325 * t397) * t795) * t799 + t547, t18 * qJD(1) + t54 * qJD(2) + t2 * qJD(3) + ((t644 + t652) * t794 + (t569 + t567) * t795) * t798 + t549, t24 * qJD(1) + t711 + t868; t548 * qJD(1) + t11 * qJD(2) + t3 * qJD(4) + (-t744 / 0.4e1 - t753 / 0.4e1 - t762 / 0.4e1) * t802 + t769, t11 * qJD(1) + t548 * qJD(2) + t5 * qJD(4) + (-t743 / 0.4e1 - t752 / 0.4e1 - t761 / 0.4e1) * t800 + t769, (m(6) * (t120 * t168 - t337 * t338 - t339 * t340) + m(5) * (t164 * t262 - t397 * t814 - t399 * t813) + m(4) * ((t535 * (rSges(4,1) * t674 - t616) + t537 * (rSges(4,1) * t662 + t535 * rSges(4,3) - t508)) * (-t449 * t535 - t450 * t537) + t615 * t502 * t500) + (-t443 * t533 + (t560 * t535 + (t444 + t561) * t537) * t535) * t770 + (-t444 * t532 + (t561 * t537 + (t443 + t560) * t535) * t537) * t771 + t604) * qJD(3) + t15 * qJD(4) + t614 * t8, t3 * qJD(1) + t5 * qJD(2) + t15 * qJD(3) + ((t610 + t55) * t794 + (t92 + t154) * t795) * t798 + (t604 - t152 - t782) * qJD(4), t614 * t99; t551 * qJD(1) + t14 * qJD(2) + t4 * qJD(3) + (-t742 / 0.4e1 - t751 / 0.4e1) * t802 + t729, t14 * qJD(1) + t551 * qJD(2) + t6 * qJD(3) + (-t741 / 0.4e1 - t750 / 0.4e1) * t800 + t729, t4 * qJD(1) + t6 * qJD(2) + t16 * qJD(4) + ((t120 * t182 + t158 * t168 - t338 * t378 - t340 * t380 + t646) * t794 + (t262 * t278 + (-t397 * t535 - t399 * t537) * t474 + t92) * t795) * t799 + (t604 + t811) * qJD(3), t16 * qJD(3) + (m(6) * t610 + t589) * qJD(4) + t614 * t564, t614 * t139; t22 * qJD(2) + t712 - t846, t22 * qJD(1) + t712 - t868, (-t338 * t537 + t340 * t535) * t726 + t249 + t614 * t98, t137 * t614 + t270 * t726 + t249, 0;];
Cq = t7;
