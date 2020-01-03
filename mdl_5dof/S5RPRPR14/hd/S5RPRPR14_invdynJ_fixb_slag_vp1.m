% Calculate vector of inverse dynamics joint torques for
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR14_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:31
% EndTime: 2019-12-31 18:35:15
% DurationCPUTime: 36.91s
% Computational Cost: add. (15189->956), mult. (23516->1213), div. (0->0), fcn. (21116->8), ass. (0->466)
t762 = Icges(4,3) + Icges(5,3);
t385 = qJ(3) + pkin(8);
t361 = sin(t385);
t362 = cos(t385);
t390 = sin(qJ(3));
t393 = cos(qJ(3));
t761 = Icges(4,5) * t390 + Icges(5,5) * t361 + Icges(4,6) * t393 + Icges(5,6) * t362;
t391 = sin(qJ(1));
t394 = cos(qJ(1));
t759 = t761 * t391 + t762 * t394;
t617 = Icges(5,4) * t361;
t462 = Icges(5,2) * t362 + t617;
t201 = -Icges(5,6) * t391 + t394 * t462;
t616 = Icges(5,4) * t362;
t465 = Icges(5,1) * t361 + t616;
t203 = -Icges(5,5) * t391 + t394 * t465;
t619 = Icges(4,4) * t390;
t463 = Icges(4,2) * t393 + t619;
t216 = -Icges(4,6) * t391 + t394 * t463;
t618 = Icges(4,4) * t393;
t466 = Icges(4,1) * t390 + t618;
t218 = -Icges(4,5) * t391 + t394 * t466;
t738 = t201 * t362 + t203 * t361 + t216 * t393 + t218 * t390;
t760 = t738 * t394;
t694 = -t762 * t391 + t761 * t394;
t758 = Icges(4,5) * t393 + Icges(5,5) * t362 - Icges(4,6) * t390 - Icges(5,6) * t361;
t200 = Icges(5,6) * t394 + t391 * t462;
t601 = t362 * t391;
t319 = Icges(5,4) * t601;
t603 = t361 * t391;
t611 = Icges(5,5) * t394;
t202 = Icges(5,1) * t603 + t319 + t611;
t215 = Icges(4,6) * t394 + t391 * t463;
t595 = t391 * t393;
t345 = Icges(4,4) * t595;
t598 = t390 * t391;
t612 = Icges(4,5) * t394;
t217 = Icges(4,1) * t598 + t345 + t612;
t692 = t200 * t362 + t202 * t361 + t215 * t393 + t217 * t390;
t757 = pkin(1) * t391;
t756 = t394 * pkin(1);
t389 = sin(qJ(5));
t392 = cos(qJ(5));
t232 = (-rSges(6,1) * t389 - rSges(6,2) * t392) * t362;
t350 = t362 * rSges(6,3);
t633 = rSges(6,2) * t389;
t636 = rSges(6,1) * t392;
t484 = -t633 + t636;
t110 = qJD(5) * t232 + (-t361 * t484 + t350) * qJD(3);
t596 = t391 * t392;
t599 = t389 * t394;
t255 = t361 * t599 + t596;
t594 = t392 * t394;
t597 = t391 * t389;
t256 = t361 * t594 - t597;
t571 = t256 * rSges(6,1) - t255 * rSges(6,2);
t600 = t362 * t394;
t130 = rSges(6,3) * t600 - t571;
t536 = qJD(1) * qJD(3);
t294 = qJDD(3) * t391 + t394 * t536;
t544 = qJD(3) * t394;
t515 = t361 * t544;
t551 = qJD(1) * t391;
t420 = -t362 * t551 - t515;
t535 = qJDD(5) * t362;
t139 = qJD(5) * t420 + t394 * t535 + t294;
t509 = t362 * t544;
t141 = t420 * pkin(7) + (t361 * t551 - t509) * pkin(4);
t631 = rSges(6,3) * t361;
t193 = t362 * t484 + t631;
t541 = qJD(5) * t362;
t251 = qJD(3) * t541 + qJDD(5) * t361 + qJDD(1);
t539 = qJD(5) * t394;
t510 = t362 * t539;
t545 = qJD(3) * t391;
t270 = t510 + t545;
t645 = t362 * pkin(4);
t284 = pkin(7) * t361 + t645;
t542 = qJD(5) * t361;
t328 = qJD(1) + t542;
t537 = qJD(1) * qJD(2);
t561 = qJDD(2) * t391 + t394 * t537;
t648 = pkin(6) * qJD(1) ^ 2;
t440 = -t394 * t648 + t561;
t652 = pkin(3) * t393;
t422 = qJDD(4) * t394 + t294 * t652 + t440;
t332 = pkin(7) * t600;
t602 = t361 * t394;
t243 = pkin(4) * t602 - t332;
t388 = -qJ(4) - pkin(6);
t640 = pkin(6) + t388;
t653 = pkin(3) * t390;
t249 = t391 * t640 + t394 * t653;
t372 = t394 * qJ(2);
t312 = -t372 + t757;
t650 = pkin(6) * t391;
t498 = -t312 - t650;
t492 = t249 + t498;
t438 = t243 + t492;
t550 = qJD(1) * t394;
t344 = t388 * t550;
t348 = pkin(3) * t598;
t511 = t393 * t544;
t543 = qJD(4) * t391;
t496 = -pkin(3) * t511 + t543;
t649 = pkin(6) * t394;
t163 = -t344 + (t348 - t649) * qJD(1) + t496;
t316 = t391 * qJ(2) + t756;
t370 = qJD(2) * t394;
t260 = qJD(1) * t316 - t370;
t490 = -t163 - t260 - t543;
t533 = qJD(3) ^ 2 * t653;
t351 = t362 * pkin(7);
t651 = pkin(4) * t361;
t283 = t351 - t651;
t262 = t283 * qJD(3);
t548 = qJD(3) * t262;
t437 = t394 * t328;
t495 = qJD(1) * t361 + qJD(5);
t685 = t391 * t495 - t509;
t106 = -t389 * t685 + t392 * t437;
t107 = t389 * t437 + t392 * t685;
t486 = rSges(6,1) * t107 + rSges(6,2) * t106;
t62 = rSges(6,3) * t420 + t486;
t10 = t270 * t110 - t251 * t130 + t139 * t193 + t294 * t284 - t328 * t62 + (-t533 + t548) * t391 + t438 * qJDD(1) + (-t141 + t490) * qJD(1) + t422;
t755 = t10 * t394;
t695 = t759 * t391;
t277 = -Icges(5,2) * t361 + t616;
t279 = Icges(5,1) * t362 - t617;
t308 = -Icges(4,2) * t390 + t618;
t310 = Icges(4,1) * t393 - t619;
t751 = t277 * t362 + t279 * t361 + t308 * t393 + t310 * t390;
t754 = t130 * t328 - t193 * t270 - t284 * t545;
t713 = t200 * t601 + t202 * t603 + t215 * t595 + t217 * t598 + t759 * t394;
t712 = -t201 * t601 - t203 * t603 - t216 * t595 - t218 * t598 - t694 * t394;
t711 = -t692 * t394 + t695;
t710 = -t694 * t391 + t760;
t705 = t201 * t361 - t203 * t362 + t216 * t390 - t218 * t393;
t727 = t758 * t391;
t752 = t277 * t361 - t279 * t362 + t308 * t390 - t310 * t393;
t706 = t200 * t361 - t202 * t362 + t215 * t390 - t217 * t393;
t723 = t758 * t394;
t704 = t751 * t391 + t723;
t703 = t394 * t751 - t727;
t253 = -t361 * t597 + t594;
t254 = t361 * t596 + t599;
t119 = Icges(6,5) * t254 + Icges(6,6) * t253 - Icges(6,3) * t601;
t121 = -Icges(6,5) * t256 + Icges(6,6) * t255 + Icges(6,3) * t600;
t230 = Icges(6,4) * t256;
t124 = Icges(6,2) * t255 + Icges(6,6) * t600 - t230;
t229 = Icges(6,4) * t255;
t126 = Icges(6,1) * t256 - Icges(6,5) * t600 - t229;
t456 = -t255 * t124 - t256 * t126;
t615 = Icges(6,4) * t254;
t122 = Icges(6,2) * t253 - Icges(6,6) * t601 + t615;
t228 = Icges(6,4) * t253;
t125 = Icges(6,1) * t254 - Icges(6,5) * t601 + t228;
t638 = t253 * t122 + t254 * t125;
t748 = t456 + t638 + (-t119 * t391 - t121 * t394) * t362;
t747 = t694 * qJD(1);
t746 = t759 * qJD(1);
t540 = qJD(5) * t391;
t271 = -t362 * t540 + t544;
t41 = t119 * t600 + t255 * t122 - t256 * t125;
t458 = Icges(6,5) * t392 - Icges(6,6) * t389;
t187 = Icges(6,3) * t361 + t362 * t458;
t613 = Icges(6,4) * t392;
t461 = -Icges(6,2) * t389 + t613;
t189 = Icges(6,6) * t361 + t362 * t461;
t614 = Icges(6,4) * t389;
t464 = Icges(6,1) * t392 - t614;
t191 = Icges(6,5) * t361 + t362 * t464;
t65 = t187 * t600 + t189 * t255 - t191 * t256;
t745 = -t271 * t41 - t328 * t65;
t744 = t394 ^ 2;
t742 = qJD(1) * t751 - qJD(3) * t761;
t741 = qJD(1) * t692 + qJD(3) * t727 + t747;
t740 = -qJD(1) * t738 - qJD(3) * t723 + t746;
t737 = t710 * t391 + t711 * t394;
t736 = t712 * t391 + t713 * t394;
t114 = qJD(1) * t201 + t277 * t545;
t226 = t279 * t391;
t116 = qJD(1) * t203 + qJD(3) * t226;
t146 = qJD(1) * t216 + t308 * t545;
t268 = t310 * t391;
t148 = qJD(1) * t218 + qJD(3) * t268;
t735 = qJD(3) * t706 - t114 * t362 - t116 * t361 - t146 * t393 - t148 * t390 + t746;
t258 = t462 * qJD(3);
t259 = t465 * qJD(3);
t287 = t463 * qJD(3);
t288 = t466 * qJD(3);
t734 = t758 * qJD(1) + t752 * qJD(3) + t258 * t362 + t259 * t361 + t287 * t393 + t288 * t390;
t224 = t277 * t394;
t113 = qJD(1) * t200 - qJD(3) * t224;
t227 = t279 * t394;
t115 = -qJD(3) * t227 + (t391 * t465 + t611) * qJD(1);
t267 = t308 * t394;
t145 = qJD(1) * t215 - qJD(3) * t267;
t269 = t310 * t394;
t147 = -qJD(3) * t269 + (t391 * t466 + t612) * qJD(1);
t733 = qJD(3) * t705 + t113 * t362 + t115 * t361 + t145 * t393 + t147 * t390 + t747;
t566 = t308 + t466;
t567 = -t463 + t310;
t568 = t277 + t465;
t569 = -t462 + t279;
t732 = (t361 * t568 - t362 * t569 + t390 * t566 - t393 * t567) * qJD(1);
t455 = t124 * t389 + t126 * t392;
t48 = t121 * t361 - t362 * t455;
t40 = -t121 * t601 + t253 * t124 - t126 * t254;
t731 = rSges(4,2) * t390;
t730 = t704 * qJD(1);
t487 = rSges(5,1) * t361 + rSges(5,2) * t362;
t424 = t487 + t653;
t508 = t362 * t545;
t726 = t361 * t550 + t508;
t725 = t703 * qJD(1);
t724 = t761 * qJD(1);
t408 = t391 * (t218 + t267) - t394 * (-Icges(4,2) * t598 + t217 + t345);
t409 = t391 * (t216 - t269) - t394 * (t215 - t268);
t410 = t391 * (t203 + t224) - t394 * (-Icges(5,2) * t603 + t202 + t319);
t411 = t391 * (t201 - t227) - t394 * (t200 - t226);
t721 = -t411 * t361 + t410 * t362 - t409 * t390 + t408 * t393;
t719 = qJD(3) * t736 + t730;
t718 = qJD(3) * t737 - t725;
t717 = t391 * t742 + t394 * t734;
t716 = -t391 * t734 + t394 * t742;
t715 = -qJD(3) * t692 - t114 * t361 + t116 * t362 - t146 * t390 + t148 * t393;
t714 = qJD(3) * t738 - t113 * t361 + t115 * t362 - t145 * t390 + t147 * t393;
t709 = rSges(4,2) * t393;
t369 = qJD(2) * t391;
t368 = qJD(4) * t394;
t512 = t393 * t545;
t563 = pkin(3) * t512 + t368;
t523 = t369 + t563;
t37 = qJD(1) * t438 + t523 - t754;
t708 = t37 * t391;
t382 = t394 * rSges(5,3);
t204 = rSges(5,1) * t603 + rSges(5,2) * t601 + t382;
t489 = -t370 + t496;
t250 = -t394 * t640 + t348;
t497 = t316 + t649;
t491 = t250 + t497;
t634 = rSges(5,2) * t361;
t637 = rSges(5,1) * t362;
t282 = -t634 + t637;
t517 = t282 * t544;
t73 = -t517 + (t204 + t491) * qJD(1) + t489;
t707 = t394 * t73;
t64 = -t187 * t601 + t189 * t253 + t191 * t254;
t702 = t270 * t40 + t64 * t328;
t654 = rSges(6,3) + pkin(7);
t701 = t361 * t654;
t186 = Icges(6,3) * t362 - t361 * t458;
t452 = t189 * t389 - t191 * t392;
t457 = t122 * t389 - t125 * t392;
t397 = t270 * (-t187 * t394 + t455) + t271 * (t187 * t391 + t457) + t328 * (t186 + t452);
t700 = t397 * t362;
t383 = t394 * rSges(4,3);
t238 = rSges(4,1) * t598 + rSges(4,2) * t595 + t383;
t693 = t238 + t497;
t317 = -rSges(3,2) * t394 + t391 * rSges(3,3);
t691 = t741 * t744 + (t733 * t391 + (-t735 + t740) * t394) * t391;
t690 = t735 * t744 + (t740 * t391 + (-t733 + t741) * t394) * t391;
t315 = rSges(4,1) * t393 - t731;
t273 = t315 * t394;
t488 = rSges(4,1) * t390 + t709;
t149 = -qJD(3) * t273 + (t391 * t488 + t383) * qJD(1);
t289 = t488 * qJD(3);
t239 = -t391 * rSges(4,3) + t394 * t488;
t493 = t239 + t498;
t49 = -t289 * t545 + t294 * t315 + (-t149 - t260) * qJD(1) + t493 * qJDD(1) + t440;
t519 = t390 * t550;
t524 = t550 * t709 + (t512 + t519) * rSges(4,1);
t546 = qJD(3) * t390;
t150 = (-rSges(4,2) * t546 - rSges(4,3) * qJD(1)) * t391 + t524;
t366 = qJDD(3) * t394;
t295 = -t391 * t536 + t366;
t560 = qJ(2) * t550 + t369;
t527 = qJD(1) * (-pkin(1) * t551 + t560) + qJDD(1) * t316 + t391 * t537;
t423 = qJDD(1) * t649 - t391 * t648 + t527;
t50 = qJD(1) * t150 + qJDD(1) * t238 - t295 * t315 + (qJD(3) * t289 - qJDD(2)) * t394 + t423;
t689 = t49 * t391 - t50 * t394;
t582 = t193 + t284;
t688 = qJD(1) * t582;
t646 = g(2) * t394;
t687 = g(1) * t391 - t646;
t117 = -t517 + (t391 * t487 + t382) * qJD(1);
t520 = t362 * t550;
t526 = rSges(5,1) * t726 + rSges(5,2) * t520;
t118 = (-rSges(5,3) * qJD(1) - qJD(3) * t634) * t391 + t526;
t494 = pkin(3) * t519 + t388 * t551 + t563;
t162 = pkin(6) * t551 + t494;
t679 = t117 * t394 + (-t118 - t162) * t391;
t222 = (-Icges(6,2) * t392 - t614) * t362;
t399 = t270 * (Icges(6,2) * t256 - t126 + t229) + t271 * (-Icges(6,2) * t254 + t125 + t228) + t328 * (t191 + t222);
t225 = (-Icges(6,1) * t389 - t613) * t362;
t676 = t270 * (-Icges(6,1) * t255 + t124 - t230) + t271 * (-Icges(6,1) * t253 + t122 + t615) + t328 * (t189 - t225);
t386 = t391 ^ 2;
t39 = -t119 * t601 + t638;
t12 = t271 * t39 + t702;
t675 = -t12 / 0.2e1;
t674 = -m(5) - m(6);
t673 = -pkin(1) - pkin(6);
t672 = t139 / 0.2e1;
t140 = -qJD(1) * t510 + t366 + (-t535 + (-qJD(1) + t542) * qJD(3)) * t391;
t671 = t140 / 0.2e1;
t670 = t251 / 0.2e1;
t669 = -t270 / 0.2e1;
t668 = t270 / 0.2e1;
t667 = -t271 / 0.2e1;
t666 = t271 / 0.2e1;
t664 = t294 / 0.2e1;
t663 = t295 / 0.2e1;
t662 = -t328 / 0.2e1;
t661 = t328 / 0.2e1;
t660 = t361 / 0.2e1;
t659 = t391 / 0.2e1;
t658 = -t394 / 0.2e1;
t656 = rSges(3,2) - pkin(1);
t655 = -rSges(5,3) - pkin(1);
t647 = g(1) * t394;
t108 = -t328 * t596 + (-t394 * t495 - t508) * t389;
t547 = qJD(3) * t362;
t109 = t495 * t594 + (-t328 * t389 + t392 * t547) * t391;
t516 = t361 * t545;
t421 = t516 - t520;
t57 = Icges(6,5) * t109 + Icges(6,6) * t108 + Icges(6,3) * t421;
t59 = Icges(6,4) * t109 + Icges(6,2) * t108 + Icges(6,6) * t421;
t61 = Icges(6,1) * t109 + Icges(6,4) * t108 + Icges(6,5) * t421;
t8 = (qJD(3) * t457 + t57) * t361 + (qJD(3) * t119 - t389 * t59 + t392 * t61 + (-t122 * t392 - t125 * t389) * qJD(5)) * t362;
t644 = t8 * t271;
t56 = Icges(6,5) * t107 + Icges(6,6) * t106 + Icges(6,3) * t420;
t58 = Icges(6,4) * t107 + Icges(6,2) * t106 + Icges(6,6) * t420;
t60 = Icges(6,1) * t107 + Icges(6,4) * t106 + Icges(6,5) * t420;
t9 = (qJD(3) * t455 + t56) * t361 + (qJD(3) * t121 - t389 * t58 + t392 * t60 + (-t124 * t392 + t126 * t389) * qJD(5)) * t362;
t643 = t9 * t270;
t219 = (-Icges(6,5) * t389 - Icges(6,6) * t392) * t362;
t103 = qJD(3) * t186 + qJD(5) * t219;
t188 = Icges(6,6) * t362 - t361 * t461;
t104 = qJD(3) * t188 + qJD(5) * t222;
t190 = Icges(6,5) * t362 - t361 * t464;
t105 = qJD(3) * t190 + qJD(5) * t225;
t18 = (qJD(3) * t452 + t103) * t361 + (qJD(3) * t187 - t104 * t389 + t105 * t392 + (-t189 * t392 - t191 * t389) * qJD(5)) * t362;
t71 = t187 * t361 - t362 * t452;
t639 = t18 * t328 + t71 * t251;
t632 = rSges(3,3) * t394;
t401 = qJDD(1) * t250 + qJDD(4) * t391 + t394 * t533 + t423 + (t162 + t368) * qJD(1);
t500 = -t282 - t652;
t261 = t487 * qJD(3);
t549 = qJD(3) * t261;
t24 = (-qJDD(2) + t549) * t394 + qJD(1) * t118 + qJDD(1) * t204 + t500 * t295 + t401;
t630 = t24 * t394;
t572 = t254 * rSges(6,1) + t253 * rSges(6,2);
t128 = -rSges(6,3) * t601 + t572;
t330 = pkin(4) * t603;
t241 = -pkin(7) * t601 + t330;
t38 = -t284 * t544 + t128 * t328 - t193 * t271 + (t241 + t491) * qJD(1) + t489;
t627 = t38 * t394;
t245 = t282 * t545;
t379 = t391 * rSges(5,3);
t205 = t394 * t487 - t379;
t439 = t205 + t492;
t72 = qJD(1) * t439 + t245 + t523;
t626 = t394 * t72;
t280 = t315 * t545;
t98 = qJD(1) * t493 + t280 + t369;
t625 = t394 * t98;
t47 = t119 * t361 - t362 * t457;
t624 = t47 * t140;
t623 = t48 * t139;
t592 = -t110 - t262;
t587 = t130 - t243;
t525 = pkin(4) * t726 + pkin(7) * t516;
t142 = -pkin(7) * t520 + t525;
t586 = -t142 - t162;
t585 = t163 * t544 - t295 * t249;
t577 = -t204 - t250;
t570 = -t241 - t250;
t313 = rSges(3,2) * t391 + t632;
t565 = -t312 + t313;
t248 = t316 + t317;
t564 = t361 * t633 + t350;
t242 = pkin(4) * t601 + pkin(7) * t603;
t559 = rSges(3,2) * t551 + rSges(3,3) * t550;
t297 = qJD(1) * t312;
t558 = t369 - t297;
t534 = -rSges(4,3) + t673;
t349 = pkin(3) * t595;
t532 = t362 * t636;
t531 = t362 * t633;
t529 = t109 * rSges(6,1) + t108 * rSges(6,2) + rSges(6,3) * t516;
t528 = -t128 + t570;
t522 = t362 * t654;
t514 = t390 * t545;
t507 = -pkin(4) - t636;
t506 = -t551 / 0.2e1;
t504 = -t545 / 0.2e1;
t503 = t545 / 0.2e1;
t502 = -t544 / 0.2e1;
t501 = t544 / 0.2e1;
t499 = -t284 - t652;
t231 = rSges(5,1) * t601 - rSges(5,2) * t603;
t318 = rSges(2,1) * t394 - rSges(2,2) * t391;
t314 = rSges(2,1) * t391 + rSges(2,2) * t394;
t63 = -rSges(6,3) * t520 + t529;
t11 = -t140 * t193 + qJD(1) * t142 + t401 + (-qJDD(2) - t548) * t394 + t499 * t295 + qJDD(1) * t241 + t251 * t128 - t271 * t110 + t328 * t63;
t483 = t11 * t391 + t755;
t478 = t39 * t394 + t391 * t40;
t477 = t39 * t391 - t394 * t40;
t42 = t121 * t600 - t456;
t476 = t391 * t42 + t394 * t41;
t475 = t391 * t41 - t394 * t42;
t474 = t391 * t48 + t394 * t47;
t473 = t391 * t47 - t394 * t48;
t99 = qJD(1) * t693 - t315 * t544 - t370;
t468 = t391 * t98 - t394 * t99;
t467 = qJD(1) * t249 - t297 + t523;
t454 = t128 * t394 + t130 * t391;
t453 = t149 * t394 - t150 * t391;
t445 = -t238 * t391 - t239 * t394;
t436 = -t388 * t394 + t316 + t348;
t435 = t344 - t489;
t434 = t494 + t560;
t170 = rSges(6,3) * t603 + (-t531 + t532) * t391;
t428 = (-t386 - t744) * t652;
t427 = -t205 * t394 + t391 * t577;
t425 = -t350 + t651 + t653;
t419 = t119 * t271 + t121 * t270 + t187 * t328;
t418 = (Icges(6,5) * t253 - Icges(6,6) * t254) * t271 + (Icges(6,5) * t255 + Icges(6,6) * t256) * t270 + t219 * t328;
t209 = t249 * t544;
t31 = -t128 * t270 + t130 * t271 - t209 + (-t243 * t394 + t391 * t570) * qJD(3);
t398 = t31 * t454 + (-t37 * t394 - t38 * t391) * t193;
t325 = rSges(5,2) * t602;
t293 = t394 * t531;
t272 = t315 * t391;
t244 = t284 * t394;
t233 = -rSges(5,1) * t600 + t325;
t210 = t394 * t249;
t192 = -t361 * t636 + t564;
t182 = qJD(1) * t248 - t370;
t181 = qJD(1) * t565 + t369;
t171 = t293 + (-t532 - t631) * t394;
t169 = t191 * t394;
t168 = t191 * t391;
t167 = t189 * t394;
t166 = t189 * t391;
t160 = t394 * t163;
t159 = rSges(6,1) * t255 + rSges(6,2) * t256;
t158 = rSges(6,1) * t253 - rSges(6,2) * t254;
t131 = t445 * qJD(3);
t84 = qJD(1) * t559 + qJDD(1) * t317 - qJDD(2) * t394 + t527;
t83 = t565 * qJDD(1) + (-qJD(1) * t317 - t260) * qJD(1) + t561;
t70 = qJD(3) * t427 - t209;
t23 = t282 * t294 + (-t533 - t549) * t391 + t439 * qJDD(1) + (-t117 + t490) * qJD(1) + t422;
t16 = -t103 * t601 + t104 * t253 + t105 * t254 + t108 * t189 + t109 * t191 + t187 * t421;
t15 = t103 * t600 + t104 * t255 - t105 * t256 + t106 * t189 + t107 * t191 + t187 * t420;
t14 = t270 * t48 + t271 * t47 + t328 * t71;
t13 = t270 * t42 - t745;
t7 = t108 * t124 - t109 * t126 + t121 * t421 + t253 * t58 + t254 * t60 - t56 * t601;
t6 = t108 * t122 + t109 * t125 + t119 * t421 + t253 * t59 + t254 * t61 - t57 * t601;
t5 = t106 * t124 - t107 * t126 + t121 * t420 + t255 * t58 - t256 * t60 + t56 * t600;
t4 = t106 * t122 + t107 * t125 + t119 * t420 + t255 * t59 - t256 * t61 + t57 * t600;
t3 = -t128 * t139 + t130 * t140 - t243 * t295 - t270 * t63 + t271 * t62 + t570 * t294 + (t141 * t394 + t391 * t586) * qJD(3) + t585;
t2 = t139 * t40 + t140 * t39 + t16 * t328 + t251 * t64 + t270 * t7 + t271 * t6;
t1 = t139 * t42 + t140 * t41 + t15 * t328 + t251 * t65 + t270 * t5 + t271 * t4;
t17 = [((t694 * t386 + ((t692 + t694) * t394 - t695 + t712) * t394) * qJD(3) + t718 + t725) * t502 + (((t695 - t711 + t712) * t391 + (-t760 + (-t692 + t694) * t391 + t710 + t713) * t394) * qJD(3) + t730) * t504 + t643 / 0.2e1 + (-(t245 - t72 + (t205 - t650) * qJD(1) + t467) * t73 + t72 * (rSges(5,1) * t509 - rSges(5,2) * t515 + t435) + t73 * (-rSges(5,2) * t516 + t434 + t526) + (t655 * t626 + (t72 * (-qJ(2) - t424) + t73 * t655) * t391) * qJD(1) + (-g(2) + t24) * (t436 + t204) + (-g(1) + t23) * (t372 - t379 + (-pkin(1) + t388) * t391 + t424 * t394)) * m(5) + (-(qJD(1) * t313 - t181 + t558) * t182 + t181 * t370 + t182 * (t559 + t560) + (t181 * t656 * t394 + (t181 * (-rSges(3,3) - qJ(2)) - t182 * pkin(1)) * t391) * qJD(1) + (t84 - g(2)) * t248 + (t83 - g(1)) * (t391 * t656 + t372 + t632)) * m(3) + ((t11 - g(2)) * (-t391 * t522 + t330 + t436 + t572) + (t10 - g(1)) * (t388 * t391 - t312 - t332 + t571) + (-t647 + t755) * t425 + (t435 - t486 + (t645 + t701) * t544 + (-t756 + (-qJ(2) - t425 + t351) * t391) * qJD(1)) * t37 + (t434 + t525 + t529 + t37 - t467 + (-t522 * t394 - t243 + t650 - t757) * qJD(1) + t754) * t38) * m(6) + t64 * t671 + t65 * t672 + (t15 + t12) * t668 + (-(t280 - t98 + (t239 - t650) * qJD(1) + t558) * t99 + t98 * (rSges(4,1) * t511 - t544 * t731 + t370) + t99 * (-rSges(4,2) * t514 + t524 + t560) + (t534 * t625 + (t98 * (-qJ(2) - t488) + t99 * t534) * t391) * qJD(1) + (-g(2) + t50) * t693 + (-g(1) + t49) * (t391 * t673 + t239 + t372)) * m(4) + (t715 + t716) * t501 + (t714 + t717 + t719) * t503 + t16 * t666 + t705 * t664 + (t704 - t706) * t663 - t703 * t294 / 0.2e1 + t644 / 0.2e1 + t639 + ((t42 + t748) * t271 + t702) * t669 + ((-t39 + t748) * t270 + t13 + t745) * t667 + (-qJD(3) * t751 + t258 * t361 - t259 * t362 + t287 * t390 - t288 * t393) * qJD(1) + (m(2) * (t314 ^ 2 + t318 ^ 2) + Icges(2,3) + Icges(3,1) - t752) * qJDD(1) - m(2) * (-g(1) * t314 + g(2) * t318) + t623 / 0.2e1 + t624 / 0.2e1; (-m(3) + t674) * t687 + 0.2e1 * (t10 * t659 + t11 * t658) * m(6) + 0.2e1 * (t23 * t659 - t630 / 0.2e1) * m(5) + 0.2e1 * (t658 * t84 + t659 * t83) * m(3) + (-t687 + t689) * m(4); -((t361 * t410 + t362 * t411 + t390 * t408 + t393 * t409) * qJD(3) + (-t361 * t569 - t362 * t568 - t390 * t567 - t393 * t566) * qJD(1)) * qJD(1) / 0.2e1 + t478 * t671 + t476 * t672 + t736 * t663 + t737 * t664 + ((-t391 * t713 + t394 * t712) * qJD(1) + t691) * t501 + (t715 * t394 + t714 * t391 + (t391 * t706 + t394 * t705) * qJD(1)) * qJD(1) / 0.2e1 + (qJD(1) * t716 + qJD(3) * t691 + qJDD(1) * t704 + t294 * t712 + t295 * t713 + t2) * t394 / 0.2e1 + (qJD(1) * t717 + qJD(3) * t690 - qJDD(1) * t703 + t294 * t710 + t295 * t711 + t1) * t659 + (t13 + t718) * t550 / 0.2e1 + (t12 + t719) * t506 + ((-t391 * t711 + t394 * t710) * qJD(1) + t690) * t503 + (-qJD(1) * t477 + t391 * t7 + t394 * t6) * t666 + (-qJD(1) * t475 + t391 * t5 + t394 * t4) * t668 + t474 * t670 + t13 * t539 * t660 - t14 * t541 / 0.2e1 + (-g(1) * (t170 + t349 + t242) - g(2) * t293 - g(3) * (t361 * t507 + t351 + t564 - t653) - (t362 * t507 - t652 - t701) * t646 + t10 * t349 - t3 * t210 + (t11 * (-t193 + t499) + t38 * t592 + t3 * t587 + t37 * t688) * t394 + (t10 * t582 + t37 * (-pkin(3) * t546 - t592) + t3 * t528 + t38 * t688) * t391 - t37 * (qJD(1) * t244 - t171 * t328 + t192 * t270) - t38 * (qJD(1) * t242 + t170 * t328 - t192 * t271) - ((t128 * t38 - t130 * t37) * t362 + t398 * t361) * qJD(5) - (-t283 * t627 + (t283 - t653) * t708) * qJD(3) + (t160 + (qJD(1) * t528 + t141 + t62) * t394 + (-t63 + t586 + (t249 - t587) * qJD(1)) * t391 + t170 * t270 - t171 * t271 - (-t242 * t391 - t244 * t394 + t428) * qJD(3)) * t31) * m(6) + (t391 * t705 - t394 * t706) * qJDD(1) / 0.2e1 + (-g(1) * (t231 + t349) - g(2) * (t325 + (-t637 - t652) * t394) + g(3) * t424 + t23 * (t282 * t391 + t349) + t72 * (-pkin(3) * t514 - t261 * t391) + t500 * t630 + t261 * t707 + (qJD(3) * t679 - t205 * t295 + t577 * t294 + t585) * (-t210 + t427) + t70 * (t160 + t679) - (t487 * t707 + t70 * (t233 * t394 + t428) + (-t70 * t231 - t424 * t72) * t391) * qJD(3) + (t70 * (t577 * t394 + (t205 + t249) * t391) + (t391 * t73 + t626) * t282 + t72 * t233 - t73 * t231) * qJD(1)) * m(5) + t361 * t540 * t675 + ((qJD(3) * t453 - t238 * t294 - t239 * t295) * t445 + t131 * ((-t238 * t394 + t239 * t391) * qJD(1) + t453) - t468 * t289 + ((t391 * t99 + t625) * qJD(1) + t689) * t315 - (t272 * t99 + t273 * t98) * qJD(1) - (t131 * (-t272 * t391 - t273 * t394) - t468 * t488) * qJD(3) - g(1) * t272 + g(2) * t273 + g(3) * t488) * m(4) + ((t166 * t253 + t168 * t254) * t271 + (-t167 * t253 - t169 * t254) * t270 + (t188 * t253 + t190 * t254) * t328 + (t362 * t64 - t40 * t602) * qJD(5) + ((qJD(5) * t39 + t419) * t361 - t700) * t391) * t667 + ((t166 * t255 - t168 * t256) * t271 + (-t167 * t255 + t169 * t256) * t270 + (t188 * t255 - t190 * t256) * t328 + (t362 * t65 + t41 * t603) * qJD(5) + ((-qJD(5) * t42 - t419) * t361 + t700) * t394) * t669 + (((-t166 * t389 + t168 * t392 + t119) * t271 + (t167 * t389 - t169 * t392 + t121) * t270 + (-t188 * t389 + t190 * t392 + t187) * t328 + t71 * qJD(5)) * t362 + (qJD(5) * t473 + t397) * t361) * t662 + ((t544 * t727 - t724) * t394 + ((-t394 * t723 - t721) * qJD(3) - t732) * t391) * t502 + ((-t545 * t723 - t724) * t391 + ((t391 * t727 + t721) * qJD(3) + t732) * t394) * t504 + (-qJD(1) * t473 + t391 * t9 + t394 * t8) * t661; t674 * (g(2) * t391 + t647) + m(5) * (t23 * t394 + t24 * t391) + m(6) * t483; t520 * t675 + t361 * t12 * t503 - t2 * t601 / 0.2e1 + (t361 * t64 - t362 * t477) * t671 + ((qJD(3) * t477 + t16) * t361 + (-qJD(1) * t478 + qJD(3) * t64 - t391 * t6 + t394 * t7) * t362) * t666 + t1 * t600 / 0.2e1 + (t361 * t65 - t362 * t475) * t672 + ((qJD(3) * t475 + t15) * t361 + (-qJD(1) * t476 + qJD(3) * t65 - t391 * t4 + t394 * t5) * t362) * t668 + t14 * t547 / 0.2e1 + (t623 + t624 + t639 + t643 + t644) * t660 + (t361 * t71 - t362 * t473) * t670 + ((qJD(3) * t473 + t18) * t361 + (-qJD(1) * t474 + qJD(3) * t71 - t391 * t8 + t394 * t9) * t362) * t661 + (t253 * t399 - t254 * t676 - t418 * t601) * t667 + (t399 * t255 + t256 * t676 + t418 * t600) * t669 + (t418 * t361 + (-t389 * t399 - t392 * t676) * t362) * t662 + (t361 * t502 + t362 * t506) * t13 + ((qJD(3) * t398 - t10 * t130 + t11 * t128 - t37 * t62 + t38 * t63) * t361 + (t37 * (-qJD(3) * t130 + t110 * t394) + t38 * (qJD(3) * t128 + t110 * t391) - t3 * t454 + t31 * (t128 * t551 - t130 * t550 - t391 * t62 - t394 * t63) + ((t627 - t708) * qJD(1) + t483) * t193) * t362 - t37 * (-t159 * t328 + t232 * t270) - t38 * (t158 * t328 - t232 * t271) - t31 * (-t158 * t270 + t159 * t271) - g(1) * t158 - g(2) * t159 - g(3) * t232) * m(6);];
tau = t17;
