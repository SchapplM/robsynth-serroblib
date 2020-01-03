% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:41
% EndTime: 2019-12-31 21:51:01
% DurationCPUTime: 13.00s
% Computational Cost: add. (69128->567), mult. (52092->708), div. (0->0), fcn. (47910->8), ass. (0->369)
t564 = qJ(3) + qJ(4);
t561 = cos(t564);
t559 = sin(t564);
t731 = Icges(5,4) * t559;
t504 = Icges(5,1) * t561 - t731;
t565 = qJ(1) + qJ(2);
t560 = sin(t565);
t562 = cos(t565);
t417 = Icges(5,5) * t560 + t504 * t562;
t546 = Icges(6,5) * t559;
t846 = Icges(6,1) * t561 + t546;
t884 = Icges(6,4) * t560 + t562 * t846 + t417;
t550 = Icges(5,4) * t561;
t503 = Icges(5,1) * t559 + t550;
t730 = Icges(6,5) * t561;
t883 = Icges(6,1) * t559 + t503 - t730;
t601 = Icges(6,3) * t561 - t546;
t880 = Icges(5,2) * t561 + t601 + t731;
t509 = pkin(4) * t561 + qJ(5) * t559;
t510 = rSges(6,1) * t561 + rSges(6,3) * t559;
t869 = t509 + t510;
t373 = t869 * t560;
t802 = rSges(6,1) + pkin(4);
t742 = rSges(6,3) + qJ(5);
t696 = t561 * t562;
t516 = Icges(6,5) * t696;
t709 = t559 * t562;
t407 = Icges(6,6) * t560 + Icges(6,3) * t709 + t516;
t497 = Icges(5,5) * t561 - Icges(5,6) * t559;
t721 = t497 * t562;
t409 = Icges(5,3) * t560 + t721;
t498 = Icges(6,4) * t561 + Icges(6,6) * t559;
t719 = t498 * t562;
t882 = t407 * t709 + t884 * t696 + (Icges(6,2) * t560 + t409 + t719) * t560;
t500 = -Icges(5,2) * t559 + t550;
t881 = t500 + t883;
t879 = (Icges(5,6) - Icges(6,6)) * t561 + (Icges(6,4) + Icges(5,5)) * t559;
t878 = -t880 * t562 + t884;
t414 = -Icges(6,4) * t562 + t560 * t846;
t711 = t559 * t560;
t517 = Icges(5,4) * t711;
t703 = t560 * t561;
t416 = Icges(5,1) * t703 - Icges(5,5) * t562 - t517;
t877 = -Icges(5,2) * t703 - t601 * t560 + t414 + t416 - t517;
t413 = Icges(5,6) * t560 + t500 * t562;
t876 = -Icges(6,1) * t709 - t503 * t562 + t407 - t413 + t516;
t496 = Icges(6,3) * t559 + t730;
t406 = -Icges(6,6) * t562 + t496 * t560;
t412 = Icges(5,4) * t703 - Icges(5,2) * t711 - Icges(5,6) * t562;
t875 = t883 * t560 - t406 + t412;
t874 = t846 + t504;
t625 = t742 * t561;
t628 = t802 * t559;
t873 = -t413 * t709 + t882;
t554 = t562 * rSges(6,2);
t568 = cos(qJ(3));
t754 = pkin(3) * t568;
t556 = pkin(2) + t754;
t824 = -pkin(8) - pkin(7);
t638 = -t560 * t556 - t562 * t824;
t838 = -t742 * t559 - t802 * t561;
t311 = t838 * t560 + t554 + t638;
t757 = sin(qJ(1)) * pkin(1);
t303 = t311 - t757;
t872 = t303 - t311;
t538 = t560 * t824;
t747 = t560 * rSges(6,2);
t312 = t747 - t538 + (t556 - t838) * t562;
t756 = cos(qJ(1)) * pkin(1);
t304 = t312 + t756;
t871 = -t304 + t312;
t563 = Icges(4,4) * t568;
t566 = sin(qJ(3));
t526 = -Icges(4,2) * t566 + t563;
t527 = Icges(4,1) * t566 + t563;
t868 = t527 + t526;
t849 = t560 * (t406 * t559 + t414 * t561);
t867 = t849 + t882;
t827 = m(4) / 0.2e1;
t826 = m(5) / 0.2e1;
t825 = m(6) / 0.2e1;
t804 = t560 / 0.2e1;
t803 = -t562 / 0.2e1;
t863 = t562 / 0.2e1;
t801 = m(3) * (t756 * (-rSges(3,1) * t560 - rSges(3,2) * t562) + (t562 * rSges(3,1) - t560 * rSges(3,2)) * t757);
t630 = t312 * t709;
t194 = -t311 * t711 + t630;
t861 = t194 * m(6) * qJD(2);
t720 = t498 * t560;
t384 = t560 * (-Icges(6,2) * t562 + t720);
t237 = t406 * t709 + t414 * t696 + t384;
t860 = t237 * t562;
t859 = t879 * t562;
t858 = t879 * t560;
t639 = t742 * t696;
t637 = t802 * t711;
t557 = t560 ^ 2;
t558 = t562 ^ 2;
t635 = t557 + t558;
t857 = (t874 - t880) * t561 + (t496 - t881) * t559;
t507 = rSges(5,1) * t559 + rSges(5,2) * t561;
t755 = pkin(3) * t566;
t587 = t507 + t755;
t847 = t587 * t562;
t848 = t587 * t560;
t590 = (t560 * t847 - t562 * t848) * t507;
t419 = rSges(5,1) * t703 - rSges(5,2) * t711 - t562 * rSges(5,3);
t351 = -t419 + t638;
t624 = -rSges(5,2) * t709 + t560 * rSges(5,3);
t748 = rSges(5,1) * t561;
t352 = -t538 + (t556 + t748) * t562 + t624;
t511 = -rSges(5,2) * t559 + t748;
t591 = (-t351 * t562 - t352 * t560) * t511;
t339 = (-t625 + t755) * t560 + t637;
t340 = (-t628 - t755) * t562 + t639;
t641 = t628 - t625;
t372 = t641 * t560;
t374 = t641 * t562;
t672 = -t374 * t339 - t372 * t340;
t375 = t869 * t562;
t675 = -t375 * t311 - t373 * t312;
t737 = (t591 + t590) * t826 + (t672 + t675) * t825;
t455 = t507 * t560;
t457 = t507 * t562;
t662 = -t455 * t847 + t457 * t848;
t613 = t641 + t755;
t353 = t613 * t560;
t355 = t613 * t562;
t358 = -t560 * t625 + t637;
t359 = -t562 * t628 + t639;
t673 = -t353 * t359 - t355 * t358;
t738 = (t591 + t662) * t826 + (t673 + t675) * t825;
t856 = t737 - t738;
t341 = t351 - t757;
t342 = t352 + t756;
t592 = (-t341 * t562 - t342 * t560) * t511;
t677 = -t375 * t303 - t373 * t304;
t739 = (t592 + t590) * t826 + (t672 + t677) * t825;
t740 = (t592 + t662) * t826 + (t673 + t677) * t825;
t855 = t739 - t740;
t805 = -t560 / 0.2e1;
t854 = t804 + t805;
t853 = (t876 * t560 + t875 * t562) * t561 + (-t878 * t560 + t877 * t562) * t559;
t188 = t341 * t848 - t342 * t847;
t191 = t351 * t848 - t352 * t847;
t126 = -t312 * t303 + t304 * t311;
t182 = -t352 * t341 + t342 * t351;
t555 = t562 * pkin(7);
t749 = rSges(4,1) * t568;
t627 = pkin(2) + t749;
t702 = t560 * t566;
t636 = rSges(4,2) * t702 + t562 * rSges(4,3);
t376 = -t560 * t627 + t555 + t636;
t360 = t376 - t757;
t686 = t562 * t566;
t537 = rSges(4,2) * t686;
t377 = -t537 + t627 * t562 + (rSges(4,3) + pkin(7)) * t560;
t361 = t377 + t756;
t192 = -t377 * t360 + t361 * t376;
t215 = t562 * (t510 * t562 + t747) + t558 * t509 + (-t554 + t373) * t560;
t245 = (-t802 * t709 + t639) * t562 + (t742 * t703 - t637) * t560;
t100 = t215 * t245 + t372 * t373 + t374 * t375;
t313 = t560 * t419 + t562 * (rSges(5,1) * t696 + t624);
t335 = -t560 * t455 - t562 * t457;
t187 = t635 * t507 * t511 + t313 * t335;
t845 = m(5) * t187 + m(6) * t100;
t657 = -t560 * (pkin(2) * t560 - t555 + t638) + t562 * (-t560 * pkin(7) - t538 + (-pkin(2) + t556) * t562);
t200 = t313 + t657;
t122 = t200 * t335 + (t560 * t848 + t562 * t847) * t511;
t184 = t215 + t657;
t75 = t184 * t245 + t353 * t373 + t355 * t375;
t844 = m(5) * t122 + m(6) * t75;
t842 = qJD(1) + qJD(2);
t741 = ((-t342 + t352) * t562 + (t341 - t351) * t560) * t507 * t826 + (t872 * t372 + t871 * t374) * t825;
t164 = t303 * t358 + t304 * t359;
t172 = t311 * t358 + t312 * t359;
t196 = t341 * t455 - t342 * t457;
t202 = t351 * t455 - t352 * t457;
t841 = (t202 + t196) * t826 + (t172 + t164) * t825;
t529 = rSges(4,1) * t566 + rSges(4,2) * t568;
t632 = (t188 - t191) * t826 + ((-t361 + t377) * t562 + (t360 - t376) * t560) * t529 * t827 + (t872 * t353 + t871 * t355) * t825;
t142 = t303 * t339 + t304 * t340;
t145 = t311 * t339 + t312 * t340;
t481 = t529 * t560;
t482 = t529 * t562;
t219 = t360 * t481 - t361 * t482;
t232 = t376 * t481 - t377 * t482;
t839 = (t191 + t188) * t826 + (t232 + t219) * t827 + (t145 + t142) * t825;
t732 = Icges(4,4) * t566;
t525 = Icges(4,2) * t568 + t732;
t528 = Icges(4,1) * t568 - t732;
t837 = t868 * t568 / 0.2e1 + (t528 / 0.2e1 - t525 / 0.2e1) * t566;
t365 = t417 * t703;
t619 = t562 * t409 - t365;
t236 = -t413 * t711 - t619;
t408 = Icges(5,5) * t703 - Icges(5,6) * t711 - Icges(5,3) * t562;
t660 = -t560 * t408 - t416 * t696;
t239 = -t412 * t709 - t660;
t617 = t413 * t559 - t408;
t723 = t412 * t559;
t586 = (-t239 * t562 + t873 * t560 - t860) * t863 + (-t860 + (t236 - t365 + (t409 + t723) * t562 + t660) * t562 + (-t849 + t867) * t560) * t803 + (((t408 + t617) * t562 - t867 + t873) * t562 + (-(t416 * t561 - t723) * t562 + t236 + t237 - t384 + t239 + t619 + t560 * t617) * t560) * t804;
t585 = (-t496 / 0.2e1 + t881 / 0.2e1) * t561 + (-t880 / 0.2e1 + t874 / 0.2e1) * t559;
t833 = 0.4e1 * qJD(1);
t831 = 0.4e1 * qJD(2);
t830 = 2 * qJD(3);
t829 = 2 * qJD(4);
t458 = t635 * t559;
t664 = -t353 * t703 - t355 * t696;
t105 = t184 * t458 + t664;
t661 = -t372 * t703 - t374 * t696;
t156 = t215 * t458 + t661;
t819 = m(6) * (t156 + t105);
t609 = -t245 * t561 - t373 * t711 - t375 * t709;
t629 = t559 * t184 + t664;
t818 = m(6) * (t609 + t629);
t85 = t215 * t559 + t609 + t661;
t812 = m(6) * t85;
t797 = m(4) * t192;
t795 = m(4) * t219;
t794 = m(4) * t232;
t784 = m(5) * t182;
t781 = m(5) * t188;
t780 = m(5) * t191;
t779 = m(5) * t196;
t778 = m(5) * t202;
t277 = t304 * t709;
t776 = m(6) * (-t872 * t711 + t277 - t630);
t775 = m(6) * (t277 + t630 + (-t303 - t311) * t711);
t666 = t339 * t709 + t340 * t711;
t671 = t303 * t696 + t304 * t703;
t773 = m(6) * (t666 + t671);
t772 = m(6) * t126;
t668 = t311 * t696 + t312 * t703;
t771 = m(6) * (t666 + t668);
t611 = -t353 * t709 + t355 * t711;
t770 = m(6) * (t611 + t671);
t663 = t358 * t709 + t359 * t711;
t768 = m(6) * (t663 + t671);
t767 = m(6) * (t611 + t668);
t766 = m(6) * (t663 + t668);
t610 = -t372 * t709 + t374 * t711;
t765 = m(6) * (t610 + t671);
t764 = m(6) * t142;
t763 = m(6) * (t610 + t668);
t762 = m(6) * t145;
t761 = m(6) * t164;
t760 = m(6) * t172;
t752 = m(6) * qJD(3);
t751 = m(6) * qJD(4);
t750 = m(6) * qJD(5);
t701 = t560 * t568;
t438 = Icges(4,4) * t701 - Icges(4,2) * t702 - Icges(4,6) * t562;
t722 = t438 * t566;
t524 = Icges(4,5) * t568 - Icges(4,6) * t566;
t716 = t524 * t562;
t710 = t559 * t561;
t685 = t562 * t568;
t674 = -t355 * t339 - t353 * t340;
t670 = -t374 * t358 - t372 * t359;
t436 = Icges(4,5) * t701 - Icges(4,6) * t702 - Icges(4,3) * t562;
t535 = Icges(4,4) * t702;
t440 = Icges(4,1) * t701 - Icges(4,5) * t562 - t535;
t656 = -t560 * t436 - t440 * t685;
t437 = Icges(4,3) * t560 + t716;
t441 = Icges(4,5) * t560 + t528 * t562;
t655 = t560 * t437 + t441 * t685;
t646 = t527 * t560 + t438;
t439 = Icges(4,6) * t560 + t526 * t562;
t645 = -t527 * t562 - t439;
t644 = -Icges(4,2) * t701 + t440 - t535;
t643 = -t525 * t562 + t441;
t642 = t635 * t710;
t189 = -t303 * t711 + t277;
t633 = m(6) * t189 * qJD(1);
t626 = -t511 - t754;
t394 = t441 * t701;
t618 = t562 * t437 - t394;
t616 = t439 * t566 - t436;
t615 = ((t858 * t560 + t853) * t562 - t859 * t557) * t804 + ((t859 * t562 + t853) * t560 - t858 * t558) * t803;
t614 = t635 * t755;
t612 = -t869 - t754;
t603 = Icges(4,5) * t566 + Icges(4,6) * t568;
t589 = (-t455 * t562 + t457 * t560) * t507;
t588 = (-t481 * t562 + t482 * t560) * t529;
t584 = t566 * t644 + t568 * t646;
t583 = -t566 * t643 + t568 * t645;
t580 = (-t525 + t528) * t568 - t868 * t566;
t578 = t585 + t841;
t577 = t585 + t837;
t576 = t577 + t839;
t575 = -t586 + (t497 * t560 + t876 * t559 + t878 * t561 + t857 * t562 + t720) * t804 + (-t875 * t559 + t857 * t560 + t877 * t561 - t719 - t721) * t803;
t574 = -t585 + ((t406 + t412) * t561 + (-t414 + t416) * t559) * t854;
t572 = t575 * qJD(4);
t571 = t574 - t837 + t854 * (t568 * t438 + t566 * t440);
t260 = -t439 * t702 - t618;
t179 = -(-t560 * (-t440 * t568 + t722) - t562 * t436) * t562 + t260 * t560;
t261 = -t438 * t686 - t656;
t262 = -t439 * t686 + t655;
t180 = -t261 * t562 + t262 * t560;
t60 = (t562 * t616 + t262 - t655) * t562 + (t560 * t616 + t261 + t618) * t560;
t61 = (t260 - t394 + (t437 + t722) * t562 + t656) * t562 + t655 * t560;
t570 = (t61 * t863 + t575 + (t60 + t179) * t805 + (t524 * t560 + t562 * t580 + t566 * t645 + t568 * t643) * t804 + (t560 * t580 - t566 * t646 + t568 * t644 + t180 - t716) * t803) * qJD(3);
t531 = -rSges(4,2) * t566 + t749;
t476 = t562 * t603;
t475 = t603 * t560;
t423 = t626 * t562;
t421 = t626 * t560;
t362 = t642 - t710;
t357 = t362 * t750;
t356 = t612 * t562;
t354 = t612 * t560;
t302 = -t614 + t335;
t253 = (-t458 * t561 - t362 + t642) * t750;
t211 = -t614 + t245;
t141 = t763 / 0.2e1;
t138 = t765 / 0.2e1;
t137 = t766 / 0.2e1;
t133 = t767 / 0.2e1;
t130 = t768 / 0.2e1;
t127 = t770 / 0.2e1;
t125 = t771 / 0.2e1;
t120 = t773 / 0.2e1;
t108 = t775 / 0.2e1;
t107 = t776 / 0.2e1;
t83 = t812 / 0.2e1;
t68 = t818 / 0.2e1;
t67 = t585 + t760 + t778;
t66 = t585 + t761 + t779;
t64 = t819 / 0.2e1;
t43 = t772 + t784 + t797 + t801;
t42 = t577 + t762 + t780 + t794;
t41 = t141 - t766 / 0.2e1;
t40 = t141 + t137;
t39 = t137 - t763 / 0.2e1;
t38 = t577 + t764 + t781 + t795;
t37 = t138 - t768 / 0.2e1;
t36 = t138 + t130;
t35 = t130 - t765 / 0.2e1;
t34 = t108 - t776 / 0.2e1;
t33 = t108 + t107;
t32 = t107 - t775 / 0.2e1;
t31 = t133 - t771 / 0.2e1;
t30 = t133 + t125;
t29 = t125 - t767 / 0.2e1;
t28 = t127 - t773 / 0.2e1;
t27 = t127 + t120;
t26 = t120 - t770 / 0.2e1;
t22 = t64 + t83 - t818 / 0.2e1;
t21 = t64 + t68 - t812 / 0.2e1;
t20 = t68 + t83 - t819 / 0.2e1;
t19 = t615 + t845;
t18 = t19 * qJD(4);
t17 = t578 + t741;
t16 = t578 - t741;
t15 = t615 + t844;
t14 = t574 + t741 - t841;
t13 = t576 + t632;
t12 = t576 - t632;
t11 = t571 + t632 - t839;
t9 = t586 * qJD(4);
t8 = (t180 / 0.2e1 - t61 / 0.2e1) * t562 + (t179 / 0.2e1 + t60 / 0.2e1) * t560 + t586;
t7 = t8 * qJD(3);
t6 = t586 + t856;
t5 = t586 - t856;
t4 = t586 + t855;
t3 = t586 - t855;
t2 = t575 + t737 + t738;
t1 = t575 + t739 + t740;
t10 = [t43 * qJD(2) + t38 * qJD(3) + t66 * qJD(4) + t189 * t750, t43 * qJD(1) + t13 * qJD(3) + t17 * qJD(4) + t33 * qJD(5) + 0.2e1 * (t801 / 0.2e1 + t126 * t825 + t182 * t826 + t192 * t827) * qJD(2), t38 * qJD(1) + t13 * qJD(2) + t570 + t1 * qJD(4) + t27 * qJD(5) + (((-t360 * t562 - t361 * t560) * t531 + t588) * t827 + (t341 * t423 + t342 * t421) * t826 + (t303 * t356 + t304 * t354 + t674) * t825) * t830, t66 * qJD(1) + t17 * qJD(2) + t1 * qJD(3) + t572 + t36 * qJD(5) + ((t670 + t677) * t825 + (t592 + t589) * t826) * t829, t33 * qJD(2) + t27 * qJD(3) + t36 * qJD(4) + t633; t12 * qJD(3) + t16 * qJD(4) + t34 * qJD(5) + (-t772 / 0.4e1 - t784 / 0.4e1 - t797 / 0.4e1 - t801 / 0.4e1) * t833, t42 * qJD(3) + t67 * qJD(4) + t194 * t750, t12 * qJD(1) + t42 * qJD(2) + t570 + t2 * qJD(4) + t30 * qJD(5) + (((-t376 * t562 - t377 * t560) * t531 + t588) * t827 + (t351 * t423 + t352 * t421) * t826 + (t311 * t356 + t312 * t354 + t674) * t825) * t830, t16 * qJD(1) + t67 * qJD(2) + t2 * qJD(3) + t572 + t40 * qJD(5) + ((t670 + t675) * t825 + (t591 + t589) * t826) * t829, t34 * qJD(1) + t30 * qJD(3) + t40 * qJD(4) + t861; t571 * qJD(1) + t11 * qJD(2) + t7 + t3 * qJD(4) + t28 * qJD(5) + (-t795 / 0.4e1 - t781 / 0.4e1 - t764 / 0.4e1) * t833, t11 * qJD(1) + t571 * qJD(2) + t7 + t5 * qJD(4) + t31 * qJD(5) + (-t794 / 0.4e1 - t780 / 0.4e1 - t762 / 0.4e1) * t831, (m(6) * (t184 * t211 - t353 * t354 - t355 * t356) + m(5) * (t200 * t302 - t421 * t848 - t423 * t847) + (-t557 * t476 + (t584 * t562 + (t475 + t583) * t560) * t562) * t804 + (-t558 * t475 + (t583 * t560 + (t476 + t584) * t562) * t560) * t803 + m(4) * ((t560 * (rSges(4,1) * t701 - t636) + t562 * (rSges(4,1) * t685 + t560 * rSges(4,3) - t537)) * (-t481 * t560 - t482 * t562) + t635 * t531 * t529) + t615) * qJD(3) + t15 * qJD(4) + t105 * t750 + t842 * t8, t3 * qJD(1) + t5 * qJD(2) + t15 * qJD(3) + t21 * qJD(5) + ((t100 + t75) * t825 + (t122 + t187) * t826) * t829 + (t615 - t845) * qJD(4), t28 * qJD(1) + t31 * qJD(2) + t21 * qJD(4) + t105 * t752 + t253; t574 * qJD(1) + t14 * qJD(2) + t4 * qJD(3) + t9 + t37 * qJD(5) + (-t761 / 0.4e1 - t779 / 0.4e1) * t833, t14 * qJD(1) + t574 * qJD(2) + t6 * qJD(3) + t9 + t41 * qJD(5) + (-t760 / 0.4e1 - t778 / 0.4e1) * t831, t4 * qJD(1) + t6 * qJD(2) + t18 + t22 * qJD(5) + ((t211 * t215 - t354 * t372 - t356 * t374 + t75) * t825 + (t302 * t313 + (-t421 * t560 - t423 * t562) * t507 + t122) * t826) * t830 + (t615 - t844) * qJD(3), t19 * qJD(3) + t156 * t750 + t586 * t842 + t18, t37 * qJD(1) + t41 * qJD(2) + t22 * qJD(3) + t156 * t751 + t253; t32 * qJD(2) + t26 * qJD(3) + t35 * qJD(4) - t633, t32 * qJD(1) + t29 * qJD(3) + t39 * qJD(4) - t861, t26 * qJD(1) + t29 * qJD(2) + (-t211 * t561 + (t354 * t560 + t356 * t562) * t559 - t105 + t629) * t752 + t20 * qJD(4) + t357, t35 * qJD(1) + t39 * qJD(2) + t20 * qJD(3) + (t85 - t156) * t751 + t357, 0.4e1 * (qJD(3) / 0.4e1 + qJD(4) / 0.4e1) * t362 * m(6);];
Cq = t10;
