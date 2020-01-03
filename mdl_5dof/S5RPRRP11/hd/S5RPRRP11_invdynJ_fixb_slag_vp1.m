% Calculate vector of inverse dynamics joint torques for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP11_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:27
% EndTime: 2019-12-31 18:54:31
% DurationCPUTime: 54.39s
% Computational Cost: add. (25986->1011), mult. (37358->1291), div. (0->0), fcn. (35554->8), ass. (0->477)
t845 = Icges(5,1) + Icges(6,1);
t859 = Icges(6,4) + Icges(5,5);
t858 = Icges(5,6) - Icges(6,6);
t456 = sin(qJ(4));
t884 = (Icges(5,4) - Icges(6,5)) * t456;
t844 = Icges(5,2) + Icges(6,3);
t883 = Icges(6,2) + Icges(5,3);
t458 = cos(qJ(4));
t882 = -t858 * t456 + t859 * t458;
t881 = t845 * t458 - t884;
t452 = pkin(8) + qJ(3);
t436 = sin(t452);
t437 = cos(t452);
t713 = Icges(5,4) * t458;
t544 = -Icges(5,2) * t456 + t713;
t880 = t436 * t544 - t858 * t437;
t688 = t436 * t458;
t390 = Icges(6,5) * t688;
t690 = t436 * t456;
t857 = Icges(6,3) * t690 + t390 - t880;
t707 = Icges(6,5) * t458;
t541 = Icges(6,3) * t456 + t707;
t810 = (t541 - t544) * t437 - t858 * t436;
t879 = t436 * t883 + t437 * t882;
t848 = t881 * t436 - t859 * t437;
t809 = t859 * t436 + t881 * t437;
t878 = (t458 * t844 + t884) * t436;
t877 = (-t859 * t456 - t858 * t458) * t436;
t459 = cos(qJ(1));
t627 = qJD(4) * t459;
t457 = sin(qJ(1));
t634 = qJD(3) * t457;
t353 = t436 * t627 + t634;
t629 = qJD(4) * t457;
t632 = qJD(3) * t459;
t354 = -t436 * t629 + t632;
t630 = qJD(4) * t437;
t407 = qJD(1) - t630;
t679 = t459 * t456;
t681 = t457 * t458;
t345 = t437 * t679 - t681;
t680 = t458 * t459;
t615 = t437 * t680;
t682 = t456 * t457;
t346 = t615 + t682;
t687 = t436 * t459;
t850 = t436 * t882 - t437 * t883;
t798 = t345 * t857 + t346 * t848 + t687 * t850;
t166 = Icges(5,5) * t346 - Icges(5,6) * t345 + Icges(5,3) * t687;
t169 = Icges(6,4) * t346 + Icges(6,2) * t687 + Icges(6,6) * t345;
t851 = t166 + t169;
t709 = Icges(6,5) * t345;
t175 = Icges(6,1) * t346 + Icges(6,4) * t687 + t709;
t311 = Icges(5,4) * t345;
t178 = Icges(5,1) * t346 + Icges(5,5) * t687 - t311;
t861 = t175 + t178;
t308 = Icges(6,5) * t346;
t163 = Icges(6,6) * t687 + Icges(6,3) * t345 + t308;
t715 = Icges(5,4) * t346;
t172 = -Icges(5,2) * t345 + Icges(5,6) * t687 + t715;
t863 = t163 - t172;
t820 = t863 * t345 + t861 * t346 + t851 * t687;
t343 = t437 * t682 + t680;
t344 = t437 * t681 - t679;
t689 = t436 * t457;
t164 = Icges(5,5) * t344 - Icges(5,6) * t343 + Icges(5,3) * t689;
t167 = Icges(6,4) * t344 + Icges(6,2) * t689 + Icges(6,6) * t343;
t817 = t164 + t167;
t307 = Icges(6,5) * t344;
t162 = -Icges(6,6) * t689 - Icges(6,3) * t343 - t307;
t310 = Icges(5,4) * t344;
t170 = -Icges(5,2) * t343 + Icges(5,6) * t689 + t310;
t818 = t162 + t170;
t306 = Icges(6,5) * t343;
t173 = Icges(6,1) * t344 + Icges(6,4) * t689 + t306;
t309 = Icges(5,4) * t343;
t177 = -Icges(5,1) * t344 - Icges(5,5) * t689 + t309;
t862 = t173 - t177;
t821 = -t345 * t818 + t346 * t862 + t687 * t817;
t827 = t353 * t820 - t354 * t821 + t407 * t798;
t799 = t343 * t857 + t344 * t848 + t689 * t850;
t822 = t863 * t343 + t861 * t344 + t851 * t689;
t823 = -t343 * t818 + t344 * t862 + t689 * t817;
t828 = t353 * t822 - t354 * t823 + t799 * t407;
t603 = t436 * t632;
t147 = qJD(1) * t343 - qJD(4) * t615 + (t603 - t629) * t456;
t524 = t407 * t456;
t580 = qJD(1) * t437 - qJD(4);
t148 = t459 * t524 + (-t457 * t580 - t603) * t458;
t602 = t437 * t632;
t637 = qJD(1) * t457;
t605 = t436 * t637;
t497 = t602 - t605;
t73 = Icges(6,5) * t148 + Icges(6,6) * t497 - Icges(6,3) * t147;
t79 = Icges(5,4) * t148 + Icges(5,2) * t147 + Icges(5,6) * t497;
t872 = t73 - t79;
t604 = t436 * t634;
t628 = qJD(4) * t458;
t636 = qJD(1) * t459;
t685 = t437 * t457;
t149 = t628 * t685 - t458 * t637 + (t437 * t636 - t604 - t627) * t456;
t633 = qJD(3) * t458;
t150 = t580 * t680 + (-t436 * t633 + t524) * t457;
t498 = t436 * t636 + t437 * t634;
t74 = Icges(6,5) * t150 + Icges(6,6) * t498 + Icges(6,3) * t149;
t80 = Icges(5,4) * t150 - Icges(5,2) * t149 + Icges(5,6) * t498;
t871 = t74 - t80;
t75 = Icges(5,5) * t148 + Icges(5,6) * t147 + Icges(5,3) * t497;
t77 = Icges(6,4) * t148 + Icges(6,2) * t497 - Icges(6,6) * t147;
t870 = t75 + t77;
t76 = Icges(5,5) * t150 - Icges(5,6) * t149 + Icges(5,3) * t498;
t78 = Icges(6,4) * t150 + Icges(6,2) * t498 + Icges(6,6) * t149;
t869 = t76 + t78;
t81 = Icges(6,1) * t148 + Icges(6,4) * t497 - Icges(6,5) * t147;
t83 = Icges(5,1) * t148 + Icges(5,4) * t147 + Icges(5,5) * t497;
t868 = t81 + t83;
t82 = Icges(6,1) * t150 + Icges(6,4) * t498 + Icges(6,5) * t149;
t84 = Icges(5,1) * t150 - Icges(5,4) * t149 + Icges(5,5) * t498;
t867 = t82 + t84;
t866 = t810 * qJD(3) + t878 * qJD(4);
t865 = t879 * qJD(3) + t877 * qJD(4);
t301 = (-Icges(5,1) * t456 - t713) * t436;
t631 = qJD(4) * t436;
t864 = (-Icges(6,1) * t456 + t707) * t631 + qJD(4) * t301 + t809 * qJD(3);
t860 = t857 * t456 + t848 * t458;
t568 = t344 * rSges(5,1) - t343 * rSges(5,2);
t182 = -rSges(5,3) * t689 - t568;
t736 = rSges(5,1) * t458;
t567 = -rSges(5,2) * t456 + t736;
t263 = -rSges(5,3) * t437 + t436 * t567;
t856 = t182 * t407 - t263 * t354;
t305 = qJD(5) * t345;
t819 = rSges(6,3) + qJ(5);
t835 = rSges(6,1) + pkin(4);
t815 = -t343 * t819 - t344 * t835;
t676 = rSges(6,2) * t689 - t815;
t855 = t407 * t676 - t305;
t623 = qJD(1) * qJD(3);
t372 = -qJDD(3) * t459 + t457 * t623;
t622 = qJDD(4) * t436;
t191 = qJD(4) * t498 + t457 * t622 + t372;
t342 = qJD(3) * t631 - qJDD(4) * t437 + qJDD(1);
t378 = pkin(3) * t604;
t683 = t437 * t459;
t412 = pkin(3) * t683;
t195 = pkin(7) * t498 + qJD(1) * t412 - t378;
t424 = t437 * pkin(3);
t786 = t436 * pkin(7) + t424;
t352 = t786 * qJD(3);
t363 = pkin(3) * t436 - pkin(7) * t437;
t337 = t786 * t457;
t443 = t459 * qJ(2);
t382 = pkin(1) * t457 - t443;
t454 = cos(pkin(8));
t425 = pkin(2) * t454 + pkin(1);
t455 = -pkin(6) - qJ(2);
t430 = t459 * t455;
t645 = -t457 * t425 - t430;
t266 = t382 + t645;
t663 = t266 - t382;
t608 = -t337 + t663;
t624 = qJD(1) * qJD(2);
t641 = qJDD(2) * t457 + t459 * t624;
t442 = t457 * qJ(2);
t384 = t459 * pkin(1) + t442;
t441 = qJD(2) * t459;
t350 = qJD(1) * t384 - t441;
t414 = t455 * t637;
t743 = pkin(1) - t425;
t667 = t414 - (-t459 * t743 - t442) * qJD(1) - t350;
t465 = t372 * t363 + (-t195 + t667) * qJD(1) + t608 * qJDD(1) - t352 * t632 + t641;
t564 = rSges(6,1) * t458 + rSges(6,3) * t456;
t800 = (-pkin(4) * t458 - qJ(5) * t456 - t564) * t436;
t666 = -rSges(6,2) * t437 - t800;
t626 = qJD(5) * t456;
t380 = t436 * t626;
t422 = t436 * rSges(6,2);
t635 = qJD(3) * t437;
t496 = t436 * t628 + t456 * t635;
t678 = t380 + t496 * qJ(5) + (t437 * t633 - t456 * t631) * pkin(4) + (-rSges(6,1) * t456 + rSges(6,3) * t458) * t631 + (t437 * t564 + t422) * qJD(3);
t304 = qJD(5) * t343;
t796 = -t819 * t149 - t150 * t835 - t304;
t739 = rSges(6,2) * t498 - t796;
t6 = -qJD(5) * t147 + qJDD(5) * t345 + t191 * t666 - t342 * t676 - t354 * t678 - t407 * t739 + t465;
t854 = -g(1) + t6;
t834 = t818 * t147 + t862 * t148 + t871 * t345 + t867 * t346 + t817 * t497 + t869 * t687;
t833 = -t863 * t147 + t861 * t148 + t872 * t345 + t868 * t346 + t851 * t497 + t870 * t687;
t832 = -t818 * t149 + t862 * t150 + t871 * t343 + t867 * t344 + t817 * t498 + t869 * t689;
t831 = t863 * t149 + t861 * t150 + t872 * t343 + t868 * t344 + t851 * t498 + t870 * t689;
t825 = -t857 * t147 + t848 * t148 + t866 * t345 + t864 * t346 + t850 * t497 + t865 * t687;
t824 = t857 * t149 + t848 * t150 + t866 * t343 + t864 * t344 + t850 * t498 + t865 * t689;
t536 = -t162 * t456 + t173 * t458;
t67 = -t167 * t437 + t436 * t536;
t534 = -t170 * t456 - t177 * t458;
t69 = -t164 * t437 + t436 * t534;
t853 = t67 + t69;
t535 = t163 * t456 + t175 * t458;
t68 = -t169 * t437 + t436 * t535;
t533 = -t172 * t456 + t178 * t458;
t70 = -t166 * t437 + t436 * t533;
t852 = t68 + t70;
t797 = t860 * t436 - t850 * t437;
t849 = -t436 * t541 + t880;
t847 = (-t860 + t879) * t407 + (t850 * t457 + t534 + t536) * t354 + (-t850 * t459 - t533 - t535) * t353;
t655 = -t819 * t688 + t690 * t835;
t846 = (t860 * qJD(3) - t865) * t437 + (t864 * t458 + t866 * t456 + (-t456 * t848 + t458 * t857) * qJD(4) + t850 * qJD(3)) * t436;
t843 = -t877 * t407 + (-t343 * t859 - t344 * t858) * t354 + (t345 * t859 + t346 * t858) * t353;
t371 = qJDD(3) * t457 + t459 * t623;
t190 = qJD(4) * t497 + t459 * t622 + t371;
t379 = pkin(7) * t602;
t499 = -t437 * t637 - t603;
t194 = pkin(3) * t499 - pkin(7) * t605 + t379;
t339 = pkin(7) * t687 + t412;
t581 = t459 * t425 - t455 * t457;
t267 = t581 - t384;
t431 = qJ(2) * t636;
t440 = qJD(2) * t457;
t640 = t431 + t440;
t520 = -qJDD(2) * t459 + qJD(1) * (-pkin(1) * t637 + t640) + qJDD(1) * t384 + t457 * t624;
t489 = qJD(1) * (-t431 + (t457 * t743 - t430) * qJD(1)) + qJDD(1) * t267 + t520;
t466 = qJD(1) * t194 + qJDD(1) * t339 - t352 * t634 - t363 * t371 + t489;
t674 = rSges(6,2) * t687 + t345 * t819 + t346 * t835;
t793 = rSges(6,2) * t602 - t819 * t147 + t148 * t835 + t305;
t740 = -rSges(6,2) * t605 + t793;
t7 = qJD(5) * t149 + qJDD(5) * t343 - t190 * t666 + t342 * t674 - t353 * t678 + t407 * t740 + t466;
t838 = -g(2) + t7;
t837 = t190 * t820 + t191 * t821 + t342 * t798 + t353 * t833 - t354 * t834 + t407 * t825;
t836 = t190 * t822 + t191 * t823 + t342 * t799 + t353 * t831 - t354 * t832 + t407 * t824;
t16 = (qJD(3) * t536 - t78) * t437 + (qJD(3) * t167 + t456 * t74 + t458 * t82 + (-t162 * t458 - t173 * t456) * qJD(4)) * t436;
t18 = (qJD(3) * t534 - t76) * t437 + (qJD(3) * t164 - t456 * t80 + t458 * t84 + (-t170 * t458 + t177 * t456) * qJD(4)) * t436;
t830 = t16 + t18;
t17 = (qJD(3) * t535 - t77) * t437 + (qJD(3) * t169 + t456 * t73 + t458 * t81 + (t163 * t458 - t175 * t456) * qJD(4)) * t436;
t19 = (qJD(3) * t533 - t75) * t437 + (qJD(3) * t166 - t456 * t79 + t458 * t83 + (-t172 * t458 - t178 * t456) * qJD(4)) * t436;
t829 = t17 + t19;
t826 = t353 * t852 - t354 * t853 + t407 * t797;
t814 = t849 * t457;
t813 = t849 * t459;
t812 = t848 * t457;
t811 = t848 * t459;
t808 = t380 - t352 - t678;
t807 = t847 * t436;
t806 = t353 * t851 - t354 * t817 + t407 * t850;
t782 = t457 * t853 + t459 * t852;
t805 = t457 * t852 - t459 * t853;
t781 = t457 * t821 + t459 * t820;
t804 = t457 * t820 - t459 * t821;
t780 = t457 * t823 + t459 * t822;
t803 = t457 * t822 - t459 * t823;
t662 = t267 + t384;
t500 = (t339 + t662) * qJD(1) - t363 * t634 - t441;
t45 = -t353 * t666 + t674 * t407 + t304 + t500;
t584 = t45 * t666;
t658 = t337 * t634 + t339 * t632;
t41 = t353 * t676 + t354 * t674 + t380 + t658;
t586 = t41 * t676;
t802 = t584 - t586;
t408 = pkin(7) * t685;
t336 = -pkin(3) * t689 + t408;
t411 = pkin(7) * t683;
t338 = -pkin(3) * t687 + t411;
t801 = t459 * t194 + t457 * t195 - t336 * t634 + t337 * t636 - t338 * t632;
t374 = qJD(1) * t382;
t788 = qJD(1) * t266 - t374;
t420 = Icges(4,4) * t437;
t545 = -Icges(4,2) * t436 + t420;
t359 = Icges(4,1) * t436 + t420;
t785 = (-t848 + t878) * t407 + (-t344 * t844 + t306 - t309 + t862) * t354 + (t346 * t844 + t311 - t709 - t861) * t353;
t784 = (-Icges(6,1) * t690 + t301 + t390 + t857) * t407 + (t343 * t845 - t307 + t310 + t818) * t354 + (-t345 * t845 + t308 - t715 + t863) * t353;
t783 = t843 * t436;
t314 = (-rSges(5,1) * t456 - rSges(5,2) * t458) * t436;
t421 = t436 * rSges(5,3);
t152 = qJD(4) * t314 + (t437 * t567 + t421) * qJD(3);
t569 = rSges(5,1) * t150 - rSges(5,2) * t149;
t88 = rSges(5,3) * t498 + t569;
t21 = -t152 * t354 + t182 * t342 + t191 * t263 - t407 * t88 + t465;
t184 = t346 * rSges(5,1) - t345 * rSges(5,2) + rSges(5,3) * t687;
t614 = t148 * rSges(5,1) + t147 * rSges(5,2) + rSges(5,3) * t602;
t86 = -rSges(5,3) * t605 + t614;
t22 = -t152 * t353 + t184 * t342 - t190 * t263 + t407 * t86 + t466;
t775 = t21 * t457 - t22 * t459;
t774 = t342 * t797 + t407 * t846;
t453 = sin(pkin(8));
t735 = rSges(3,2) * t453;
t738 = rSges(3,1) * t454;
t291 = t457 * rSges(3,3) + (-t735 + t738) * t459;
t746 = g(1) * t457;
t773 = -g(2) * t459 + t746;
t772 = t436 * (g(1) * t459 + g(2) * t457);
t701 = Icges(4,3) * t459;
t270 = Icges(4,5) * t685 - Icges(4,6) * t689 - t701;
t391 = Icges(4,4) * t689;
t711 = Icges(4,5) * t459;
t274 = Icges(4,1) * t685 - t391 - t711;
t704 = Icges(4,6) * t459;
t272 = Icges(4,4) * t685 - Icges(4,2) * t689 - t704;
t694 = t272 * t436;
t529 = -t274 * t437 + t694;
t102 = -t270 * t459 - t457 * t529;
t356 = Icges(4,5) * t437 - Icges(4,6) * t436;
t355 = Icges(4,5) * t436 + Icges(4,6) * t437;
t501 = qJD(3) * t355;
t716 = Icges(4,4) * t436;
t360 = Icges(4,1) * t437 - t716;
t275 = Icges(4,5) * t457 + t360 * t459;
t273 = Icges(4,6) * t457 + t459 * t545;
t693 = t273 * t436;
t528 = -t275 * t437 + t693;
t770 = -t459 * t501 + (-t356 * t457 + t528 + t701) * qJD(1);
t271 = Icges(4,3) * t457 + t356 * t459;
t639 = qJD(1) * t271;
t769 = qJD(1) * t529 - t457 * t501 + t639;
t357 = Icges(4,2) * t437 + t716;
t525 = t357 * t436 - t359 * t437;
t768 = t525 * qJD(1) + t356 * qJD(3);
t522 = t194 * t632 + t195 * t634 + t371 * t337 - t339 * t372;
t5 = qJD(5) * t496 + qJDD(5) * t690 + t190 * t676 - t191 * t674 + t353 * t739 + t354 * t740 + t522;
t767 = t41 * t740 + t5 * t674;
t766 = t457 * (-t357 * t459 + t275) - t459 * (-Icges(4,2) * t685 + t274 - t391);
t763 = t190 / 0.2e1;
t762 = t191 / 0.2e1;
t761 = t342 / 0.2e1;
t760 = -t353 / 0.2e1;
t759 = t353 / 0.2e1;
t758 = -t354 / 0.2e1;
t757 = t354 / 0.2e1;
t756 = t371 / 0.2e1;
t755 = t372 / 0.2e1;
t754 = -t407 / 0.2e1;
t753 = t407 / 0.2e1;
t751 = t457 / 0.2e1;
t750 = -t459 / 0.2e1;
t748 = -rSges(6,2) - pkin(7);
t747 = -rSges(5,3) - pkin(7);
t737 = rSges(4,1) * t437;
t734 = t16 * t354;
t733 = t17 * t353;
t732 = t18 * t354;
t731 = t19 * t353;
t447 = t457 * rSges(4,3);
t63 = t184 * t407 - t263 * t353 + t500;
t722 = t457 * t63;
t721 = t67 * t191;
t720 = t68 * t190;
t719 = t69 * t191;
t718 = t70 * t190;
t361 = rSges(4,1) * t436 + rSges(4,2) * t437;
t571 = -t361 * t632 + t440;
t276 = rSges(4,1) * t685 - rSges(4,2) * t689 - t459 * rSges(4,3);
t609 = -t276 + t663;
t108 = qJD(1) * t609 + t571;
t698 = t108 * t457;
t277 = rSges(4,1) * t683 - rSges(4,2) * t687 + t447;
t109 = -t361 * t634 - t441 + (t277 + t662) * qJD(1);
t318 = t361 * t459;
t697 = t109 * t318;
t692 = t355 * t457;
t691 = t355 * t459;
t686 = t437 * t456;
t684 = t437 * t458;
t677 = -t152 - t352;
t673 = -t343 * t835 + t344 * t819;
t672 = -t345 * t835 + t346 * t819;
t397 = rSges(6,2) * t685;
t671 = t457 * t800 + t397;
t404 = rSges(6,2) * t683;
t670 = t459 * t800 + t404;
t669 = -t457 * t270 - t274 * t683;
t668 = t457 * t271 + t275 * t683;
t665 = -t263 - t363;
t664 = -t684 * t835 - t686 * t819 - t422;
t656 = t457 * t337 + t459 * t339;
t244 = t291 + t384;
t651 = -t357 + t360;
t650 = t359 + t545;
t618 = rSges(5,2) * t690;
t649 = rSges(5,3) * t685 + t457 * t618;
t648 = rSges(5,3) * t683 + t459 * t618;
t647 = rSges(4,2) * t605 + rSges(4,3) * t636;
t621 = t457 * t738;
t415 = t457 * t735;
t642 = t459 * rSges(3,3) + t415;
t290 = t621 - t642;
t646 = -t382 - t290;
t644 = rSges(3,3) * t636 + qJD(1) * t415;
t643 = t414 + t441;
t638 = qJD(1) * t356;
t124 = -t457 * t525 - t691;
t625 = t124 * qJD(1);
t620 = rSges(5,1) * t688;
t610 = -t363 - t666;
t606 = t378 + t643;
t597 = -pkin(1) - t738;
t594 = t636 / 0.2e1;
t592 = -t634 / 0.2e1;
t591 = t634 / 0.2e1;
t590 = -t632 / 0.2e1;
t589 = t632 / 0.2e1;
t570 = -t363 * t632 + t440;
t475 = qJD(1) * t608 + t570;
t44 = -t354 * t666 + t475 - t855;
t585 = t44 * t666;
t238 = t275 * t685;
t583 = t271 * t459 - t238;
t582 = -t270 + t693;
t572 = qJD(1) * t338 - t634 * t786;
t385 = rSges(2,1) * t459 - rSges(2,2) * t457;
t383 = rSges(2,1) * t457 + rSges(2,2) * t459;
t362 = -rSges(4,2) * t436 + t737;
t127 = t273 * t437 + t275 * t436;
t502 = qJD(3) * t357;
t155 = -t459 * t502 + (-t457 * t545 + t704) * qJD(1);
t503 = qJD(3) * t359;
t157 = -t459 * t503 + (-t360 * t457 + t711) * qJD(1);
t468 = -qJD(3) * t127 - t155 * t436 + t157 * t437 + t639;
t126 = t272 * t437 + t274 * t436;
t156 = qJD(1) * t273 - t457 * t502;
t158 = qJD(1) * t275 - t457 * t503;
t469 = qJD(1) * t270 - qJD(3) * t126 - t156 * t436 + t158 * t437;
t562 = -(t457 * t769 + t469 * t459) * t459 + (t457 * t770 + t468 * t459) * t457;
t561 = -(t469 * t457 - t459 * t769) * t459 + (t468 * t457 - t459 * t770) * t457;
t103 = -t273 * t689 - t583;
t540 = -t102 * t459 + t103 * t457;
t104 = -t272 * t687 - t669;
t105 = -t273 * t687 + t668;
t539 = -t104 * t459 + t105 * t457;
t538 = -t108 * t459 - t109 * t457;
t159 = rSges(4,1) * t499 - rSges(4,2) * t602 + t647;
t316 = t361 * t457;
t160 = -qJD(3) * t316 + (t362 * t459 + t447) * qJD(1);
t537 = t159 * t459 + t160 * t457;
t532 = -t182 * t459 - t184 * t457;
t527 = t276 * t457 + t277 * t459;
t526 = t357 * t437 + t359 * t436;
t265 = rSges(5,1) * t684 - rSges(5,2) * t686 + t421;
t523 = -t425 - t786;
t521 = t581 + t339;
t519 = t436 * t748 - t424;
t518 = t436 * t747 - t424;
t505 = -qJD(1) * t336 - t632 * t786;
t504 = -pkin(3) * t603 + t379 + t440;
t495 = t41 * t739 + t5 * t676;
t490 = -qJD(1) * t337 + t570 + t788;
t488 = -t44 * t676 + t45 * t674;
t487 = -t41 * t674 + t585;
t486 = t272 * t459 - t273 * t457;
t476 = (-t436 * t650 + t437 * t651) * qJD(1);
t62 = t475 + t856;
t64 = -t182 * t353 + t184 * t354 + t658;
t470 = t64 * t532 + (t457 * t62 - t459 * t63) * t263;
t348 = t545 * qJD(3);
t349 = t360 * qJD(3);
t467 = qJD(1) * t355 - qJD(3) * t526 - t348 * t436 + t349 * t437;
t464 = t487 * t457 - t459 * t802;
t463 = -t436 * t766 + t486 * t437;
t351 = t362 * qJD(3);
t341 = t363 * t637;
t236 = -t459 * t620 + t648;
t234 = -t457 * t620 + t649;
t220 = qJD(1) * t244 - t441;
t219 = qJD(1) * t646 + t440;
t217 = -rSges(5,1) * t345 - rSges(5,2) * t346;
t212 = -rSges(5,1) * t343 - rSges(5,2) * t344;
t128 = t527 * qJD(3);
t125 = -t459 * t525 + t692;
t123 = t125 * qJD(1);
t107 = qJDD(1) * t291 + qJD(1) * (-qJD(1) * t621 + t644) + t520;
t106 = t646 * qJDD(1) + (-qJD(1) * t291 - t350) * qJD(1) + t641;
t66 = -qJD(3) * t528 + t155 * t437 + t157 * t436;
t65 = -t529 * qJD(3) + t156 * t437 + t158 * t436;
t61 = t467 * t457 - t459 * t768;
t60 = t457 * t768 + t467 * t459;
t49 = qJD(1) * t159 + qJDD(1) * t277 - t351 * t634 - t361 * t371 + t489;
t48 = -t351 * t632 + t361 * t372 + t609 * qJDD(1) + (-t160 + t667) * qJD(1) + t641;
t47 = qJD(3) * t539 + t123;
t46 = qJD(3) * t540 + t625;
t20 = -t182 * t190 - t184 * t191 + t353 * t88 + t354 * t86 + t522;
t1 = [(-t625 + ((t459 * t582 + t105 - t668) * t459 + (t457 * t582 + t104 + t583) * t457) * qJD(3) + t46) * t592 + (t125 + t127) * t756 + (t66 + t60) * t591 + (t584 * t354 + t838 * (t521 + t674) + (t6 * t457 - t746) * t519 + (t606 + t796 + t635 * t748 * t457 + (t523 - t422) * t636 + t45) * t44 + (-t490 + t504 + t793 + ((-t425 + t519) * t457 - t430) * qJD(1) + t855) * t45 + t854 * (t645 + t815)) * m(6) + t827 * t757 + (-(-qJD(1) * t276 - t108 + t571 + t788) * t109 + t108 * t643 + t109 * (t440 + t647) + (t361 * t698 - t697) * qJD(3) + ((-t108 * rSges(4,3) + t109 * (-t425 - t737)) * t457 + (t108 * (-t362 - t425) - t109 * t455) * t459) * qJD(1) + (-g(2) + t49) * (t277 + t581) + (-g(1) + t48) * (-t276 + t645)) * m(4) + (t219 * t441 + t220 * (t640 + t644) + (t219 * (t597 + t735) * t459 + (t219 * (-rSges(3,3) - qJ(2)) + t220 * t597) * t457) * qJD(1) - (-qJD(1) * t290 - t219 - t374 + t440) * t220 + (-g(2) + t107) * t244 + (-g(1) + t106) * (t597 * t457 + t443 + t642)) * m(3) + (t817 * t437 + (t456 * t818 - t458 * t862) * t436 + t853) * t353 * t754 + (t124 + t126) * t755 + t825 * t759 + (t824 + t827) * t758 + t774 + (t123 + ((t103 - t238 + (t271 + t694) * t459 + t669) * t459 + t668 * t457) * qJD(3)) * t589 + t721 / 0.2e1 + t733 / 0.2e1 - t734 / 0.2e1 + t731 / 0.2e1 - t732 / 0.2e1 + (Icges(3,2) * t454 ^ 2 + (Icges(3,1) * t453 + 0.2e1 * Icges(3,4) * t454) * t453 + m(2) * (t383 ^ 2 + t385 ^ 2) + Icges(2,3) + t526) * qJDD(1) + (-(t490 - t62 + t856) * t63 - t518 * t746 + t62 * (-t569 + t606) + t63 * (t504 + t614) + (t62 * t635 * t747 + t21 * t518) * t457 + ((-t425 + t518) * t722 + (t62 * (t523 - t421) - t63 * t455) * t459) * qJD(1) + (-g(2) + t22) * (t521 + t184) + (-g(1) + t21) * (-t568 + t645)) * m(5) + t798 * t763 + t799 * t762 + (t65 + t61 + t47) * t590 - m(2) * (-g(1) * t383 + g(2) * t385) + (-qJD(3) * t525 + t348 * t437 + t349 * t436) * qJD(1) + t718 / 0.2e1 + t719 / 0.2e1 + t720 / 0.2e1; (-m(3) - m(4) - m(6)) * t773 + 0.2e1 * (t6 * t751 + t7 * t750) * m(6) + 0.2e1 * (t48 * t751 + t49 * t750) * m(4) + 0.2e1 * (t106 * t751 + t107 * t750) * m(3) + (-t773 + t775) * m(5); (qJD(1) * t61 + qJD(3) * t561 + qJDD(1) * t124 + t102 * t372 + t103 * t371 + t836) * t750 + (qJD(1) * t60 + qJD(3) * t562 + qJDD(1) * t125 + t104 * t372 + t105 * t371 + t837) * t751 + ((-t632 * t692 - t638) * t459 + (t476 + (t459 * t691 + t463) * qJD(3)) * t457) * t589 + ((-t634 * t691 + t638) * t457 + (t476 + (t457 * t692 + t463) * qJD(3)) * t459) * t592 - qJD(1) * ((t436 * t651 + t437 * t650) * qJD(1) + (t486 * t436 + t437 * t766) * qJD(3)) / 0.2e1 + t803 * t762 + t804 * t763 + t805 * t761 - (t457 * t828 + t459 * t827) * t630 / 0.2e1 + (-g(1) * (t404 + t411) - g(2) * (t397 + t408) - g(3) * (t786 - t664) - (-t456 * t819 - t458 * t835 - pkin(3)) * t772 - (t436 * t488 + t437 * t464) * qJD(4) + t5 * t656 + (qJD(1) * t586 + t6 * t610 + t767) * t459 + (qJD(1) * t585 + t7 * t610 + t495) * t457 + (-t664 * t353 - t670 * t407 + t457 * t808 + t610 * t636 - t572) * t45 + (-t437 * t626 - t670 * t354 - t671 * t353 + (-t339 - t674) * t637 + t801) * t41 + (-t664 * t354 + t671 * t407 + t459 * t808 + t341 - t505) * t44) * m(6) - t826 * t631 / 0.2e1 + (t47 + t827) * t594 + (t46 + t828) * t637 / 0.2e1 + (qJD(1) * t782 + t457 * t829 - t459 * t830) * t753 + (qJD(1) * t780 + t457 * t831 - t459 * t832) * t758 + (-g(1) * (t411 + t648) - g(2) * (t408 + t649) - g(3) * (t265 + t786) - (-pkin(3) - t736) * t772 - t62 * (-t234 * t407 - t265 * t354 + t505) - t63 * (t236 * t407 - t265 * t353 + t572) - ((t182 * t62 + t184 * t63) * t436 + t470 * t437) * qJD(4) + t62 * t341 + t20 * t656 + (t20 * t184 + t62 * t677 + (qJD(1) * t63 + t21) * t665) * t459 + (qJD(1) * t263 * t62 - t182 * t20 + t22 * t665 + t63 * t677) * t457 + (-t234 * t353 - t236 * t354 + (-qJD(1) * t182 + t86) * t459 + (t88 + (-t184 - t339) * qJD(1)) * t457 + t801) * t64) * m(5) + ((qJD(4) * t782 - t847) * t437 + ((t456 * t810 + t458 * t809 + t850) * t407 + (-t456 * t814 + t458 * t812 - t817) * t354 + (t456 * t813 - t458 * t811 + t851) * t353 + t797 * qJD(4)) * t436) * t754 + ((t436 * t799 + t683 * t822) * qJD(4) + ((qJD(4) * t823 + t806) * t437 + t807) * t457 + (t343 * t810 + t344 * t809) * t407 + (-t343 * t814 + t344 * t812) * t354 + (t343 * t813 - t344 * t811) * t353) * t757 + ((t436 * t798 + t685 * t821) * qJD(4) + ((qJD(4) * t820 + t806) * t437 + t807) * t459 + (t345 * t810 + t346 * t809) * t407 + (-t345 * t814 + t346 * t812) * t354 + (t345 * t813 - t346 * t811) * t353) * t760 + (g(1) * t318 + g(2) * t316 - g(3) * t362 - (t108 * t316 - t697) * qJD(1) - (t128 * (-t316 * t457 - t318 * t459) + t538 * t362) * qJD(3) + (qJD(3) * t537 + t276 * t371 - t277 * t372) * t527 + t128 * ((t276 * t459 - t277 * t457) * qJD(1) + t537) + t538 * t351 + (-t49 * t457 - t48 * t459 + (-t109 * t459 + t698) * qJD(1)) * t361) * m(4) + (qJD(1) * t781 + t457 * t833 - t459 * t834) * t759 + qJD(1) * (t457 * t66 - t459 * t65 + (t126 * t457 + t127 * t459) * qJD(1)) / 0.2e1 + ((t102 * t457 + t103 * t459) * qJD(1) + t561) * t590 + ((t104 * t457 + t105 * t459) * qJD(1) + t562) * t591 + t540 * t755 + t539 * t756 + qJDD(1) * (-t126 * t459 + t127 * t457) / 0.2e1; (t436 * t781 - t437 * t798) * t763 + (t780 * t436 - t437 * t799) * t762 + (t436 * t782 - t437 * t797) * t761 + (t345 * t785 + t346 * t784 - t459 * t783) * t760 + ((qJD(3) * t781 - t825) * t437 + (-qJD(1) * t804 + t798 * qJD(3) + t457 * t834 + t459 * t833) * t436) * t759 + ((qJD(3) * t780 - t824) * t437 + (-qJD(1) * t803 + t799 * qJD(3) + t457 * t832 + t459 * t831) * t436) * t758 + (t343 * t785 + t344 * t784 - t457 * t783) * t757 + (t843 * t437 + (t456 * t785 + t458 * t784) * t436) * t754 + ((qJD(3) * t782 - t846) * t437 + (-qJD(1) * t805 + t797 * qJD(3) + t457 * t830 + t459 * t829) * t436) * t753 - (t720 + t721 + t733 - t734 + t718 + t719 + t731 - t732 + t774) * t437 / 0.2e1 + t836 * t689 / 0.2e1 + t837 * t687 / 0.2e1 + t826 * qJD(3) * t436 / 0.2e1 + ((qJD(3) * t464 + t44 * t739 - t45 * t740 + t6 * t676 - t674 * t7) * t437 + (t488 * qJD(3) + (qJD(1) * t487 - t45 * t678 - t666 * t7 + t495) * t459 + (qJD(1) * t802 + t44 * t678 + t6 * t666 - t767) * t457) * t436 - g(1) * t672 - g(2) * t673 + g(3) * t655 - (t344 * t45 + t346 * t44 + t41 * t688) * qJD(5) - (-t44 * t673 + t45 * t672) * t407 - (t41 * t672 + t44 * t655) * t354 - (t41 * t673 + t45 * t655) * t353) * m(6) + (-g(1) * t217 - g(2) * t212 - g(3) * t314 + (qJD(3) * t470 - t182 * t21 - t22 * t184 + t62 * t88 - t63 * t86) * t437 + (t62 * (qJD(3) * t182 + t152 * t457) + t63 * (qJD(3) * t184 - t152 * t459) + t20 * t532 + t64 * (t182 * t637 - t184 * t636 - t457 * t86 + t459 * t88) + ((t459 * t62 + t722) * qJD(1) + t775) * t263) * t436 - t62 * (-t212 * t407 - t314 * t354) - t63 * (t217 * t407 - t314 * t353) - t64 * (t212 * t353 + t217 * t354)) * m(5) + t828 * (t436 * t594 + t437 * t591) + t827 * (t437 * t589 - t605 / 0.2e1); (-t147 * t44 + t149 * t45 + t496 * t41 + (t353 * t45 + t354 * t44 - g(3) + t5) * t690 + (-t41 * t354 - t45 * t407 + t854) * t345 + (-t353 * t41 + t407 * t44 + t838) * t343) * m(6);];
tau = t1;
