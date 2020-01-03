% Calculate vector of inverse dynamics joint torques for
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR13_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR13_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:57
% EndTime: 2019-12-31 18:32:51
% DurationCPUTime: 45.96s
% Computational Cost: add. (17792->970), mult. (25016->1236), div. (0->0), fcn. (22395->8), ass. (0->468)
t819 = Icges(5,4) - Icges(4,5);
t818 = Icges(5,5) - Icges(4,6);
t817 = Icges(5,1) + Icges(4,3);
t415 = pkin(8) + qJ(3);
t397 = sin(t415);
t398 = cos(t415);
t790 = t818 * t397 - t819 * t398;
t675 = Icges(4,4) * t397;
t303 = Icges(4,2) * t398 + t675;
t664 = Icges(5,6) * t397;
t490 = Icges(5,3) * t398 + t664;
t807 = -t303 - t490;
t382 = Icges(4,4) * t398;
t305 = Icges(4,1) * t397 + t382;
t663 = Icges(5,6) * t398;
t492 = Icges(5,2) * t397 + t663;
t816 = t305 + t492;
t424 = cos(qJ(1));
t815 = t817 * t424;
t422 = sin(qJ(1));
t648 = t398 * t422;
t650 = t397 * t422;
t793 = t819 * t648 - t818 * t650 + t815;
t801 = t817 * t422 + t790 * t424;
t665 = Icges(4,6) * t424;
t212 = Icges(4,4) * t648 - Icges(4,2) * t650 - t665;
t351 = Icges(5,6) * t650;
t673 = Icges(5,4) * t424;
t219 = Icges(5,2) * t648 - t351 + t673;
t814 = t212 * t397 - t219 * t398;
t306 = Icges(4,1) * t398 - t675;
t215 = Icges(4,5) * t422 + t306 * t424;
t491 = -Icges(5,3) * t397 + t663;
t216 = Icges(5,5) * t422 - t424 * t491;
t813 = -t215 * t648 - t216 * t650;
t668 = Icges(5,5) * t424;
t217 = Icges(5,6) * t648 - Icges(5,3) * t650 + t668;
t812 = t212 + t217;
t498 = -Icges(4,2) * t397 + t382;
t213 = Icges(4,6) * t422 + t424 * t498;
t811 = t213 - t216;
t356 = Icges(4,4) * t650;
t669 = Icges(4,5) * t424;
t214 = Icges(4,1) * t648 - t356 - t669;
t810 = t214 + t219;
t649 = t397 * t424;
t352 = Icges(5,6) * t649;
t647 = t398 * t424;
t674 = Icges(5,4) * t422;
t218 = -Icges(5,2) * t647 + t352 + t674;
t809 = t215 - t218;
t299 = Icges(4,5) * t397 + Icges(4,6) * t398;
t496 = Icges(5,4) * t397 + Icges(5,5) * t398;
t808 = t299 - t496;
t493 = Icges(5,2) * t398 - t664;
t805 = t306 + t493;
t804 = t807 * qJD(3);
t803 = t816 * qJD(3);
t786 = -t214 * t398 + t217 * t397 + t814;
t802 = t491 + t498;
t475 = t303 * t397 - t305 * t398;
t788 = -t397 * t490 + t398 * t492 - t475;
t800 = t424 * t801 + t813;
t739 = t215 * t647 + t216 * t649 + t422 * t801;
t799 = -t214 * t647 + t217 * t649 + t422 * t793;
t798 = t213 * t397 + t218 * t398;
t759 = -t786 * t422 + t793 * t424;
t758 = -t213 * t650 - t218 * t648 - t800;
t757 = -t212 * t649 + t219 * t647 - t799;
t756 = -t213 * t649 - t218 * t647 + t739;
t750 = t810 * t397 + t812 * t398;
t749 = t809 * t397 + t811 * t398;
t797 = t804 * t424 + (-t422 * t802 + t665 - t668) * qJD(1);
t796 = t811 * qJD(1) + t804 * t422;
t795 = -t803 * t424 + (-t422 * t805 + t669 - t673) * qJD(1);
t794 = t803 * t422 + (-t424 * t493 - t215 + t674) * qJD(1);
t792 = t802 * qJD(3);
t791 = t805 * qJD(3);
t789 = t808 * qJD(3);
t787 = -t397 * t816 + t398 * t807;
t785 = t215 * t398 + t216 * t397 - t798;
t737 = t808 * t422;
t651 = t299 * t424;
t96 = -t422 * t475 - t651;
t254 = t496 * t424;
t99 = t490 * t650 - t492 * t648 - t254;
t784 = t96 - t99;
t755 = t788 * t424 + t737;
t783 = t801 * qJD(1);
t423 = cos(qJ(5));
t643 = t422 * t423;
t421 = sin(qJ(5));
t644 = t421 * t424;
t282 = t397 * t643 + t644;
t642 = t423 * t424;
t645 = t421 * t422;
t283 = -t397 * t645 + t642;
t613 = t283 * rSges(6,1) - t282 * rSges(6,2);
t154 = rSges(6,3) * t648 - t613;
t520 = rSges(6,1) * t421 + rSges(6,2) * t423;
t741 = t398 * t520;
t206 = rSges(6,3) * t397 - t741;
t579 = qJD(5) * t398;
t584 = qJD(3) * t424;
t294 = -t422 * t579 + t584;
t580 = qJD(5) * t397;
t371 = qJD(1) + t580;
t307 = pkin(3) * t397 - qJ(4) * t398;
t541 = -pkin(7) * t397 - t307;
t581 = qJD(4) * t424;
t341 = t397 * t581;
t401 = qJD(2) * t422;
t602 = t341 + t401;
t782 = t154 * t371 + t294 * t206 - t541 * t584 - t602;
t781 = t424 ^ 2;
t780 = qJD(1) * t808 + t787 * qJD(3) - t792 * t397 + t791 * t398;
t779 = t793 * qJD(1) + t750 * qJD(3) + t796 * t397 + t794 * t398;
t778 = -t749 * qJD(3) - t797 * t397 + t795 * t398 + t783;
t777 = -t811 * t422 + t812 * t424;
t776 = t756 * t422 - t757 * t424;
t775 = t758 * t422 - t759 * t424;
t774 = t802 + t816;
t773 = t805 + t807;
t772 = (t351 + t356 + (Icges(4,2) + Icges(5,3)) * t648 - t810) * t424 + (-Icges(5,3) * t647 - t303 * t424 - t352 + t809) * t422;
t771 = t788 * qJD(1) - t790 * qJD(3);
t770 = -t789 * t424 + (-t790 * t422 - t785 + t815) * qJD(1);
t769 = t786 * qJD(1) - t789 * t422 + t783;
t768 = t755 * qJD(1);
t384 = t398 * rSges(6,3);
t205 = t397 * t520 + t384;
t314 = pkin(4) * t424 - pkin(7) * t648;
t419 = cos(pkin(8));
t386 = pkin(2) * t419 + pkin(1);
t420 = -pkin(6) - qJ(2);
t391 = t424 * t420;
t600 = -t422 * t386 - t391;
t767 = t314 + t600 + t613;
t766 = t784 * qJD(1);
t585 = qJD(3) * t422;
t293 = t424 * t579 + t585;
t280 = t397 * t642 - t645;
t281 = t397 * t644 + t643;
t143 = Icges(6,5) * t281 + Icges(6,6) * t280 + Icges(6,3) * t647;
t672 = Icges(6,4) * t281;
t146 = Icges(6,2) * t280 + Icges(6,6) * t647 + t672;
t260 = Icges(6,4) * t280;
t149 = Icges(6,1) * t281 + Icges(6,5) * t647 + t260;
t41 = t143 * t647 + t280 * t146 + t281 * t149;
t145 = -Icges(6,5) * t283 + Icges(6,6) * t282 + Icges(6,3) * t648;
t262 = Icges(6,4) * t283;
t148 = Icges(6,2) * t282 + Icges(6,6) * t648 - t262;
t261 = Icges(6,4) * t282;
t150 = Icges(6,1) * t283 - Icges(6,5) * t648 - t261;
t42 = t145 * t647 + t280 * t148 - t150 * t281;
t494 = Icges(6,5) * t421 + Icges(6,6) * t423;
t444 = -Icges(6,3) * t397 + t398 * t494;
t671 = Icges(6,4) * t421;
t495 = Icges(6,2) * t423 + t671;
t445 = -Icges(6,6) * t397 + t398 * t495;
t670 = Icges(6,4) * t423;
t499 = Icges(6,1) * t421 + t670;
t446 = -Icges(6,5) * t397 + t398 * t499;
t68 = -t280 * t445 - t281 * t446 - t444 * t647;
t12 = t293 * t41 - t294 * t42 + t68 * t371;
t43 = t143 * t648 + t282 * t146 - t283 * t149;
t44 = t145 * t648 + t148 * t282 + t150 * t283;
t69 = -t282 * t445 + t283 * t446 - t444 * t648;
t13 = t293 * t43 - t294 * t44 + t371 * t69;
t487 = t148 * t423 - t150 * t421;
t54 = t145 * t397 - t398 * t487;
t765 = qJD(3) * t775 + t766;
t764 = qJD(3) * t776 + t768;
t763 = -t422 * t771 + t424 * t780;
t762 = t422 * t780 + t424 * t771;
t761 = t786 * qJD(3) + t794 * t397 - t796 * t398;
t760 = t785 * qJD(3) + t795 * t397 + t797 * t398;
t748 = -t397 * t772 + t398 * t777;
t747 = (-t397 * t774 + t398 * t773) * qJD(1);
t746 = t779 * t781 + (t770 * t422 + (-t769 + t778) * t424) * t422;
t745 = t769 * t781 + (t778 * t422 + (-t770 + t779) * t424) * t422;
t744 = t793 + t798;
t576 = qJD(3) * qJD(4);
t743 = qJDD(4) * t397 + t398 * t576;
t555 = t397 * t584;
t589 = qJD(1) * t422;
t558 = t398 * t589;
t742 = t555 + t558;
t740 = t790 * qJD(1);
t738 = t651 - t254;
t404 = t424 * qJ(2);
t343 = pkin(1) * t422 - t404;
t208 = t343 + t600;
t324 = qJD(1) * t343;
t735 = qJD(1) * t208 - t324;
t381 = t397 * qJ(4);
t733 = t398 * pkin(3) + t381;
t734 = -pkin(7) * t398 - t733;
t383 = t397 * rSges(5,3);
t691 = rSges(5,2) * t398;
t519 = t383 - t691;
t413 = t422 * pkin(4);
t310 = pkin(7) * t647 + t413;
t403 = t422 * qJ(2);
t345 = t424 * pkin(1) + t403;
t266 = t733 * t422;
t625 = t208 - t343;
t564 = -t266 + t625;
t526 = t314 + t564;
t39 = qJD(1) * t526 - t782;
t152 = t281 * rSges(6,1) + t280 * rSges(6,2) + rSges(6,3) * t647;
t556 = t397 * t585;
t334 = pkin(7) * t556;
t582 = qJD(4) * t422;
t551 = t397 * t582;
t348 = qJ(4) * t649;
t271 = pkin(3) * t647 + t348;
t358 = t424 * t386;
t534 = -t420 * t422 + t358;
t209 = t534 - t345;
t624 = t209 + t345;
t563 = t271 + t624;
t402 = qJD(2) * t424;
t612 = -t307 * t585 - t402;
t40 = t551 + t152 * t371 - t206 * t293 - t334 + (t310 + t563) * qJD(1) + t612;
t732 = t39 * t424 + t40 * t422;
t730 = -qJD(1) * t266 + t735;
t418 = sin(pkin(8));
t693 = rSges(3,2) * t418;
t696 = rSges(3,1) * t419;
t244 = t422 * rSges(3,3) + (-t693 + t696) * t424;
t252 = (Icges(6,2) * t421 - t670) * t398;
t439 = t293 * (-Icges(6,2) * t281 + t149 + t260) - t294 * (Icges(6,2) * t283 - t150 + t261) + t371 * (-t446 + t252);
t257 = (-Icges(6,1) * t423 + t671) * t398;
t440 = t293 * (-Icges(6,1) * t280 + t146 + t672) - t294 * (-Icges(6,1) * t282 + t148 - t262) + t371 * (-t445 - t257);
t721 = m(5) / 0.2e1;
t720 = m(6) / 0.2e1;
t719 = -m(5) - m(6);
t577 = qJD(1) * qJD(3);
t321 = qJDD(3) * t422 + t424 * t577;
t574 = qJDD(5) * t398;
t155 = -qJD(5) * t742 + t424 * t574 + t321;
t718 = t155 / 0.2e1;
t322 = -qJDD(3) * t424 + t422 * t577;
t588 = qJD(1) * t424;
t557 = t398 * t588;
t456 = -t556 + t557;
t156 = qJD(5) * t456 + t422 * t574 + t322;
t717 = t156 / 0.2e1;
t279 = qJD(3) * t579 + qJDD(5) * t397 + qJDD(1);
t716 = t279 / 0.2e1;
t715 = -t293 / 0.2e1;
t714 = t293 / 0.2e1;
t713 = -t294 / 0.2e1;
t712 = t294 / 0.2e1;
t711 = t321 / 0.2e1;
t710 = t322 / 0.2e1;
t709 = -t371 / 0.2e1;
t708 = t371 / 0.2e1;
t707 = t422 / 0.2e1;
t706 = -t424 / 0.2e1;
t705 = -rSges(6,3) - pkin(3);
t703 = g(1) * t422;
t488 = t146 * t423 + t149 * t421;
t529 = qJD(1) * t397 + qJD(5);
t553 = t398 * t584;
t447 = -t422 * t529 + t553;
t472 = t371 * t421;
t124 = t423 * t447 - t424 * t472;
t473 = t423 * t371;
t125 = t421 * t447 + t424 * t473;
t61 = Icges(6,5) * t125 + Icges(6,6) * t124 - Icges(6,3) * t742;
t63 = Icges(6,4) * t125 + Icges(6,2) * t124 - Icges(6,6) * t742;
t65 = Icges(6,1) * t125 + Icges(6,4) * t124 - Icges(6,5) * t742;
t8 = (qJD(3) * t488 + t61) * t397 + (qJD(3) * t143 - t421 * t65 - t423 * t63 + (t146 * t421 - t149 * t423) * qJD(5)) * t398;
t702 = t8 * t293;
t471 = t529 * t424;
t586 = qJD(3) * t398;
t122 = t423 * t471 + (t423 * t586 - t472) * t422;
t554 = t398 * t585;
t123 = t422 * t473 + (t471 + t554) * t421;
t60 = Icges(6,5) * t123 + Icges(6,6) * t122 + Icges(6,3) * t456;
t62 = Icges(6,4) * t123 + Icges(6,2) * t122 + Icges(6,6) * t456;
t64 = Icges(6,1) * t123 + Icges(6,4) * t122 + Icges(6,5) * t456;
t9 = (qJD(3) * t487 + t60) * t397 + (qJD(3) * t145 - t421 * t64 - t423 * t62 + (t148 * t421 + t150 * t423) * qJD(5)) * t398;
t701 = t9 * t294;
t698 = pkin(1) - t386;
t198 = Icges(6,3) * t398 + t397 * t494;
t249 = (-Icges(6,5) * t423 + Icges(6,6) * t421) * t398;
t119 = qJD(3) * t198 + qJD(5) * t249;
t200 = Icges(6,6) * t398 + t397 * t495;
t120 = qJD(3) * t200 + qJD(5) * t252;
t202 = Icges(6,5) * t398 + t397 * t499;
t121 = qJD(3) * t202 + qJD(5) * t257;
t485 = -t421 * t446 - t423 * t445;
t21 = (qJD(3) * t485 + t119) * t397 + (-qJD(3) * t444 - t120 * t423 - t121 * t421 + (-t421 * t445 + t423 * t446) * qJD(5)) * t398;
t72 = -t397 * t444 - t398 * t485;
t697 = t21 * t371 + t72 * t279;
t695 = rSges(4,1) * t398;
t692 = rSges(5,2) * t397;
t267 = (-rSges(6,1) * t423 + rSges(6,2) * t421) * t398;
t126 = qJD(3) * t205 + qJD(5) * t267;
t194 = qJD(1) * t310 - t334;
t583 = qJD(4) * t398;
t229 = qJD(3) * t733 - t583;
t241 = t397 * t588 + t554;
t335 = pkin(3) * t556;
t118 = pkin(3) * t557 + qJ(4) * t241 - t335 + t551;
t290 = qJD(1) * t345 - t402;
t377 = t420 * t589;
t628 = t377 - (-t424 * t698 - t403) * qJD(1) - t290;
t465 = -t118 - t551 + t628;
t578 = qJD(1) * qJD(2);
t596 = qJDD(2) * t422 + t424 * t578;
t501 = t322 * t307 + t424 * t743 + t596;
t646 = t398 * qJD(3) ^ 2;
t522 = rSges(6,1) * t123 + rSges(6,2) * t122;
t66 = rSges(6,3) * t456 + t522;
t10 = -t229 * t584 - t294 * t126 - t279 * t154 + t156 * t206 - t371 * t66 + (t322 * t397 - t424 * t646) * pkin(7) + t526 * qJDD(1) + (-t194 + t465) * qJD(1) + t501;
t689 = t10 * t422;
t396 = pkin(4) * t588;
t195 = -pkin(7) * t742 + t396;
t328 = qJ(4) * t553;
t559 = t397 * t589;
t117 = -pkin(3) * t742 - qJ(4) * t559 + t328 + t341;
t392 = qJ(2) * t588;
t595 = t392 + t401;
t468 = -qJDD(2) * t424 + qJD(1) * (-pkin(1) * t589 + t595) + qJDD(1) * t345 + t422 * t578;
t452 = qJD(1) * (-t392 + (t422 * t698 - t391) * qJD(1)) + qJDD(1) * t209 + t468;
t438 = qJDD(1) * t271 + t452 + t743 * t422 + (t117 + t341) * qJD(1);
t639 = t125 * rSges(6,1) + t124 * rSges(6,2);
t67 = -rSges(6,3) * t742 + t639;
t11 = t438 - t229 * t585 + (-t321 * t397 - t422 * t646) * pkin(7) + qJD(1) * t195 - t155 * t206 + t279 * t152 - t293 * t126 + qJDD(1) * t310 - t321 * t307 + t371 * t67;
t688 = t11 * t424;
t309 = rSges(4,1) * t397 + rSges(4,2) * t398;
t270 = t309 * t424;
t408 = t422 * rSges(4,3);
t231 = rSges(4,1) * t647 - rSges(4,2) * t649 + t408;
t87 = -t309 * t585 - t402 + (t231 + t624) * qJD(1);
t687 = t270 * t87;
t410 = t422 * rSges(5,1);
t524 = -t309 * t584 + t401;
t230 = rSges(4,1) * t648 - rSges(4,2) * t650 - t424 * rSges(4,3);
t565 = -t230 + t625;
t86 = qJD(1) * t565 + t524;
t681 = t422 * t86;
t53 = t143 * t397 - t398 * t488;
t680 = t53 * t155;
t679 = t54 * t156;
t523 = t266 * t585 + t271 * t584 - t583;
t38 = t152 * t294 + t154 * t293 + (t310 * t424 - t314 * t422) * qJD(3) + t523;
t659 = qJD(3) * t38;
t634 = t152 + t310;
t633 = t154 - t314;
t618 = -t519 * qJD(3) - t229;
t232 = -rSges(5,2) * t647 + rSges(5,3) * t649 + t410;
t617 = -t232 - t271;
t616 = t422 * t266 + t424 * t271;
t349 = qJ(4) * t647;
t268 = -pkin(3) * t649 + t349;
t615 = qJD(1) * t268 + t398 * t582;
t189 = t244 + t345;
t614 = -t271 - t310;
t518 = rSges(5,3) * t398 + t692;
t607 = -t307 + t518;
t606 = -t733 - t519;
t605 = t741 * t422;
t604 = t741 * t424;
t603 = rSges(4,2) * t559 + rSges(4,3) * t588;
t572 = t422 * t696;
t378 = t422 * t693;
t597 = t424 * rSges(3,3) + t378;
t243 = t572 - t597;
t601 = -t343 - t243;
t264 = rSges(5,2) * t650 + rSges(5,3) * t648;
t269 = rSges(5,2) * t649 + rSges(5,3) * t647;
t599 = rSges(3,3) * t588 + qJD(1) * t378;
t598 = t377 + t402;
t594 = t422 ^ 2 + t781;
t587 = qJD(3) * t397;
t380 = qJD(4) * t397;
t573 = -pkin(7) + t705;
t567 = t40 * t588;
t566 = t424 * t117 + t422 * t118 + t266 * t588;
t347 = qJ(4) * t648;
t263 = -pkin(3) * t650 + t347;
t562 = t263 * t585 + t268 * t584 + t380;
t561 = t328 + t602;
t560 = t335 + t598;
t550 = -pkin(1) - t696;
t547 = t588 / 0.2e1;
t546 = -t585 / 0.2e1;
t545 = t585 / 0.2e1;
t544 = -t584 / 0.2e1;
t543 = t584 / 0.2e1;
t539 = rSges(5,1) * t424 - rSges(5,3) * t650;
t533 = qJD(3) * t618;
t532 = -qJD(1) * t263 + t398 * t581;
t531 = t422 * t573;
t530 = t573 * t424;
t233 = rSges(5,2) * t648 + t539;
t528 = t233 + t564;
t527 = rSges(5,1) * t588 + rSges(5,2) * t742 + rSges(5,3) * t553;
t525 = -t206 + t541;
t346 = rSges(2,1) * t424 - rSges(2,2) * t422;
t344 = rSges(2,1) * t422 + rSges(2,2) * t424;
t313 = -rSges(4,2) * t397 + t695;
t512 = t41 * t424 + t42 * t422;
t511 = t41 * t422 - t42 * t424;
t510 = t422 * t44 + t424 * t43;
t509 = t422 * t43 - t424 * t44;
t508 = t422 * t54 + t424 * t53;
t507 = t422 * t53 - t424 * t54;
t502 = -t422 * t87 - t424 * t86;
t139 = -rSges(4,1) * t742 - rSges(4,2) * t553 + t603;
t265 = t309 * t422;
t140 = -qJD(3) * t265 + (t313 * t424 + t408) * qJD(1);
t489 = t139 * t424 + t140 * t422;
t486 = t152 * t422 - t154 * t424;
t479 = t230 * t422 + t231 * t424;
t470 = -pkin(7) * t586 - t126 - t229;
t469 = t534 + t271;
t458 = -qJDD(4) * t398 + t117 * t584 + t118 * t585 + t321 * t266 + t397 * t576;
t457 = t584 * t607 + t602;
t454 = -t143 * t293 + t145 * t294 + t371 * t444;
t453 = (Icges(6,5) * t280 - Icges(6,6) * t281) * t293 - (Icges(6,5) * t282 + Icges(6,6) * t283) * t294 + t249 * t371;
t449 = -t386 - t733 - t383;
t448 = t398 * t453;
t437 = t38 * t486 + (-t39 * t422 + t40 * t424) * t206;
t428 = (t444 * t424 + t488) * t293 - (t444 * t422 + t487) * t294 + (t198 + t485) * t371;
t427 = t428 * t398;
t292 = t313 * qJD(3);
t278 = t307 * t589;
t242 = t553 - t559;
t240 = t594 * t587;
t178 = -rSges(6,3) * t649 + t604;
t177 = -rSges(6,3) * t650 + t605;
t176 = t446 * t424;
t175 = t446 * t422;
t174 = t445 * t424;
t173 = t445 * t422;
t170 = qJD(1) * t189 - t402;
t169 = qJD(1) * t601 + t401;
t168 = rSges(6,1) * t282 + rSges(6,2) * t283;
t167 = rSges(6,1) * t280 - rSges(6,2) * t281;
t142 = -rSges(5,3) * t559 + t527;
t141 = t518 * t585 + (t424 * t519 + t410) * qJD(1);
t108 = t479 * qJD(3);
t85 = qJDD(1) * t244 + qJD(1) * (-qJD(1) * t572 + t599) + t468;
t84 = t601 * qJDD(1) + (-qJD(1) * t244 - t290) * qJD(1) + t596;
t71 = (t232 * t424 - t233 * t422) * qJD(3) + t523;
t58 = (qJD(3) * t518 + t380) * t422 + (t232 + t563) * qJD(1) + t612;
t57 = qJD(1) * t528 + t457;
t35 = qJD(1) * t139 + qJDD(1) * t231 - t292 * t585 - t309 * t321 + t452;
t34 = -t292 * t584 + t309 * t322 + t565 * qJDD(1) + (-t140 + t628) * qJD(1) + t596;
t19 = -t233 * t321 + t617 * t322 + (t141 * t422 + t142 * t424) * qJD(3) + t458;
t18 = qJD(1) * t142 + qJDD(1) * t232 + t321 * t607 + t422 * t533 + t438;
t17 = -t518 * t322 + t424 * t533 + t528 * qJDD(1) + (-t141 + t465) * qJD(1) + t501;
t16 = t119 * t647 + t120 * t280 + t121 * t281 - t124 * t445 - t125 * t446 + t444 * t742;
t15 = t119 * t648 + t120 * t282 - t121 * t283 - t122 * t445 - t123 * t446 - t444 * t456;
t14 = t293 * t53 - t294 * t54 + t371 * t72;
t7 = t124 * t148 - t125 * t150 - t145 * t742 + t280 * t62 + t281 * t64 + t60 * t647;
t6 = t124 * t146 + t125 * t149 - t143 * t742 + t280 * t63 + t281 * t65 + t61 * t647;
t5 = t122 * t148 - t123 * t150 + t145 * t456 + t282 * t62 - t283 * t64 + t60 * t648;
t4 = t122 * t146 + t123 * t149 + t143 * t456 + t282 * t63 - t283 * t65 + t61 * t648;
t3 = -t152 * t156 + t154 * t155 + t293 * t66 + t294 * t67 - t314 * t321 + t614 * t322 + (t194 * t422 + t195 * t424) * qJD(3) + t458;
t2 = t155 * t41 + t156 * t42 + t16 * t371 + t279 * t68 + t293 * t6 - t294 * t7;
t1 = t15 * t371 + t155 * t43 + t156 * t44 + t279 * t69 + t293 * t4 - t294 * t5;
t20 = [(-(qJD(1) * t233 + t457 - t57 + t730) * t58 + t57 * t560 + t58 * (-pkin(3) * t555 + t527 + t561) + t57 * (-t380 + (-t692 + (-rSges(5,3) - qJ(4)) * t398) * qJD(3)) * t422 + ((-t57 * rSges(5,1) + t449 * t58) * t422 + (t57 * (t449 + t691) - t58 * t420) * t424) * qJD(1) + (t18 - g(2)) * (t232 + t469) + (t17 - g(1)) * ((-t381 + (rSges(5,2) - pkin(3)) * t398) * t422 + t539 + t600)) * m(5) + (-(-qJD(1) * t230 + t524 + t735 - t86) * t87 + t86 * t598 + t87 * (t401 + t603) + (t309 * t681 - t687) * qJD(3) + ((-t86 * rSges(4,3) + t87 * (-t386 - t695)) * t422 + (t86 * (-t313 - t386) - t87 * t420) * t424) * qJD(1) + (t35 - g(2)) * (t231 + t534) + (t34 - g(1)) * (-t230 + t600)) * m(4) + (-(-qJD(1) * t243 - t169 - t324 + t401) * t170 + t169 * t402 + t170 * (t595 + t599) + (t169 * (t550 + t693) * t424 + (t169 * (-rSges(3,3) - qJ(2)) + t170 * t550) * t422) * qJD(1) + (t85 - g(2)) * t189 + (t84 - g(1)) * (t422 * t550 + t404 + t597)) * m(3) + (t788 * qJD(3) + t791 * t397 + t792 * t398) * qJD(1) + t16 * t714 + t69 * t717 + t68 * t718 + (t749 + t755) * t711 + (-g(1) * t767 - (t398 * t705 - t381) * t703 + ((-t733 - t384) * t422 + t767) * t10 + (t334 - t522 + t560 + (rSges(6,3) * t587 - qJ(4) * t586 - t380) * t422 + (t398 * t530 - t348 - t358 - t413) * qJD(1)) * t39 + (-g(2) + t11) * (t469 + t634) + (t530 * t587 + t39 + t396 + t561 + t639 - t730 + (-qJ(4) * t650 + t398 * t531 - t314 + t600) * qJD(1) + t782) * t40) * m(6) + (t15 + t12) * t713 + t679 / 0.2e1 + t680 / 0.2e1 + t12 * t712 + (Icges(3,2) * t419 ^ 2 + (Icges(3,1) * t418 + 0.2e1 * Icges(3,4) * t419) * t418 + m(2) * (t344 ^ 2 + t346 ^ 2) + Icges(2,3) - t787) * qJDD(1) + (t96 + t750) * t710 + (((t424 * t744 - t739 + t756) * t424 + (t422 * t744 + t757 + t800) * t422) * qJD(3) + t765 - t766) * t546 + t702 / 0.2e1 - t701 / 0.2e1 + t697 + ((t739 * t422 + ((t801 + t814) * t424 + t758 + t799 + t813) * t424) * qJD(3) + t768) * t543 + (t760 + t763) * t545 + (-t761 + t762 + t764) * t544 - t322 * t99 / 0.2e1 - m(2) * (-g(1) * t344 + g(2) * t346); (-m(3) - m(4) + t719) * (-g(2) * t424 + t703) + 0.2e1 * (t689 / 0.2e1 - t688 / 0.2e1) * m(6) + 0.2e1 * (t17 * t707 + t18 * t706) * m(5) + 0.2e1 * (t34 * t707 + t35 * t706) * m(4) + 0.2e1 * (t706 * t85 + t707 * t84) * m(3); (t12 * t424 + t13 * t422) * t580 / 0.2e1 + ((t174 * t282 - t176 * t283) * t293 - (t173 * t282 - t175 * t283) * t294 + (t200 * t282 - t202 * t283) * t371 + (t398 * t69 - t43 * t649) * qJD(5) + ((-qJD(5) * t44 + t454) * t397 + t427) * t422) * t712 + (qJD(1) * t510 + t4 * t422 - t424 * t5) * t713 + (qJD(1) * t512 + t422 * t6 - t424 * t7) * t714 + ((t174 * t280 + t176 * t281) * t293 - (t173 * t280 + t175 * t281) * t294 + (t200 * t280 + t202 * t281) * t371 + (t398 * t68 - t42 * t650) * qJD(5) + ((-qJD(5) * t41 + t454) * t397 + t427) * t424) * t715 + t507 * t716 + t509 * t717 + t511 * t718 + (-t57 * (-qJD(1) * t264 + t532) - t58 * (qJD(1) * t269 + t615) - t71 * t562 - ((t71 * t269 + t57 * t606) * t424 + (t71 * t264 + t58 * t606) * t422) * qJD(3) - g(1) * (t268 + t269) - g(2) * (t263 + t264) + g(3) * t606 + t57 * t278 + t19 * t616 + t71 * t566 + (t17 * t607 + t57 * t618 + t19 * t232 + t71 * t142 + (-t71 * t233 + t58 * t607) * qJD(1)) * t424 + (t18 * t607 + t58 * t618 - t19 * t233 + t71 * t141 + (-t518 * t57 + t617 * t71) * qJD(1)) * t422) * m(5) + ((-t584 * t737 - t740) * t424 + ((t424 * t738 + t748) * qJD(3) + t747) * t422) * t543 + ((-t585 * t738 + t740) * t422 + ((t422 * t737 + t748) * qJD(3) + t747) * t424) * t546 + (((-t174 * t423 - t176 * t421 + t143) * t293 - (-t173 * t423 - t175 * t421 + t145) * t294 + (-t200 * t423 - t202 * t421 - t444) * t371 + t72 * qJD(5)) * t398 + (-qJD(5) * t508 + t428) * t397) * t709 + (t762 * qJD(1) + t745 * qJD(3) + qJDD(1) * t784 + t758 * t321 + t759 * t322 + t1) * t706 - t14 * t579 / 0.2e1 + (t422 * t749 - t424 * t750) * qJDD(1) / 0.2e1 + (t13 + t765) * t589 / 0.2e1 + ((qJD(3) * t489 + t230 * t321 - t231 * t322) * t479 + t108 * ((t230 * t424 - t231 * t422) * qJD(1) + t489) + t502 * t292 + (-t34 * t424 - t35 * t422 + (-t424 * t87 + t681) * qJD(1)) * t309 - (t265 * t86 - t687) * qJD(1) - (t108 * (-t265 * t422 - t270 * t424) + t502 * t313) * qJD(3) + g(1) * t270 + g(2) * t265 - g(3) * t313) * m(4) + (qJD(1) * t508 + t422 * t8 - t424 * t9) * t708 + (t3 * t616 + (t11 * t525 + t3 * t633 + t40 * t470) * t422 + (t3 * t634 + (qJD(1) * t40 + t10) * t525) * t424 - g(1) * (t349 + t604) - g(2) * (t347 + t605) - g(3) * (t205 - t734) - t40 * (t152 * t579 + t178 * t371 - t205 * t293 + t615) - t732 * qJD(3) * t734 + (-g(1) * t530 - g(2) * t531 - (-t594 * t659 - t567) * pkin(7) - t437 * qJD(5)) * t397 + (t154 * t579 + t177 * t371 + t205 * t294 + t206 * t589 + t424 * t470 + t278 - t532) * t39 + (t566 + (t194 + t66 + (-t152 + t614) * qJD(1)) * t422 + (qJD(1) * t633 + t195 + t67) * t424 - t177 * t293 - t178 * t294 - t562) * t38) * m(6) + t775 * t710 + t776 * t711 - ((t397 * t777 + t398 * t772) * qJD(3) + (t773 * t397 + t774 * t398) * qJD(1)) * qJD(1) / 0.2e1 + ((t422 * t757 + t424 * t756) * qJD(1) + t746) * t545 + ((t422 * t759 + t424 * t758) * qJD(1) + t745) * t544 + (t761 * t424 + t760 * t422 + (t422 * t750 + t424 * t749) * qJD(1)) * qJD(1) / 0.2e1 + (qJD(1) * t763 + qJD(3) * t746 + qJDD(1) * t755 + t321 * t756 + t322 * t757 + t2) * t707 + (t12 + t764) * t547; t719 * (-g(3) * t398 + (g(1) * t424 + g(2) * t422) * t397) - m(5) * (t240 * t71 + t241 * t58 + t242 * t57) - m(6) * (t240 * t38 + t241 * t40 + t242 * t39) + 0.2e1 * ((t57 * t584 + t58 * t585 - t19) * t721 + (t39 * t584 + t40 * t585 - t3) * t720) * t398 + 0.2e1 * ((qJD(3) * t71 + t17 * t424 + t18 * t422 - t57 * t589 + t58 * t588) * t721 + (t10 * t424 + t11 * t422 - t39 * t589 + t567 + t659) * t720) * t397; t2 * t647 / 0.2e1 + (t397 * t68 + t398 * t512) * t718 + ((-qJD(3) * t512 + t16) * t397 + (-qJD(1) * t511 + qJD(3) * t68 + t422 * t7 + t424 * t6) * t398) * t714 + t1 * t648 / 0.2e1 + (t397 * t69 + t398 * t510) * t717 + ((-qJD(3) * t510 + t15) * t397 + (-qJD(1) * t509 + qJD(3) * t69 + t4 * t424 + t422 * t5) * t398) * t713 + t14 * t586 / 0.2e1 + t397 * (t679 + t680 + t697 - t701 + t702) / 0.2e1 + (t397 * t72 + t398 * t508) * t716 + ((-qJD(3) * t508 + t21) * t397 + (-qJD(1) * t507 + qJD(3) * t72 + t422 * t9 + t424 * t8) * t398) * t708 + (t439 * t280 - t281 * t440 + t424 * t448) * t715 + (t282 * t439 + t283 * t440 + t422 * t448) * t712 + (t453 * t397 + (t440 * t421 - t423 * t439) * t398) * t709 + (t397 * t546 + t398 * t547) * t13 + (-t558 / 0.2e1 + t397 * t544) * t12 + ((qJD(3) * t437 - t10 * t154 + t11 * t152 - t39 * t66 + t40 * t67) * t397 + (t39 * (-qJD(3) * t154 + t126 * t422) + t40 * (qJD(3) * t152 - t126 * t424) - t3 * t486 + t38 * (-t152 * t588 - t154 * t589 - t422 * t67 + t424 * t66) + (qJD(1) * t732 - t688 + t689) * t206) * t398 - t39 * (-t168 * t371 - t267 * t294) - t40 * (t167 * t371 - t267 * t293) - t38 * (t167 * t294 + t168 * t293) - g(1) * t167 - g(2) * t168 - g(3) * t267) * m(6);];
tau = t20;
