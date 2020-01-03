% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:20
% EndTime: 2019-12-31 18:12:50
% DurationCPUTime: 23.34s
% Computational Cost: add. (9709->602), mult. (13312->706), div. (0->0), fcn. (10166->6), ass. (0->348)
t727 = Icges(5,4) - Icges(4,5);
t726 = Icges(5,5) - Icges(4,6);
t725 = Icges(5,1) + Icges(4,3);
t345 = pkin(7) + qJ(3);
t328 = sin(t345);
t329 = cos(t345);
t724 = t726 * t328 - t727 * t329;
t244 = Icges(6,4) * t328 + Icges(6,5) * t329;
t351 = sin(qJ(1));
t352 = cos(qJ(1));
t158 = Icges(6,1) * t351 + t244 * t352;
t689 = t725 * t351 + t724 * t352;
t723 = t158 + t689;
t316 = Icges(6,6) * t328;
t239 = -Icges(6,2) * t329 + t316;
t573 = Icges(5,6) * t328;
t423 = Icges(5,3) * t329 + t573;
t722 = t239 - t423;
t571 = Icges(6,6) * t329;
t572 = Icges(5,6) * t329;
t712 = t571 - t572 + (-Icges(5,2) - Icges(6,3)) * t328;
t721 = t725 * t352;
t559 = t328 * t351;
t298 = Icges(4,4) * t559;
t557 = t329 * t351;
t579 = Icges(4,5) * t352;
t148 = Icges(4,1) * t557 - t298 - t579;
t287 = Icges(6,6) * t559;
t577 = Icges(6,5) * t352;
t151 = -Icges(6,3) * t557 - t287 + t577;
t720 = -t148 + t151;
t584 = Icges(4,4) * t328;
t251 = Icges(4,1) * t329 - t584;
t149 = Icges(4,5) * t351 + t251 * t352;
t421 = Icges(6,3) * t329 + t316;
t150 = Icges(6,5) * t351 + t352 * t421;
t710 = t149 + t150;
t424 = -Icges(5,3) * t328 + t572;
t152 = Icges(5,5) * t351 - t352 * t424;
t556 = t329 * t352;
t289 = Icges(6,6) * t556;
t558 = t328 * t352;
t581 = Icges(6,4) * t351;
t154 = Icges(6,2) * t558 + t289 + t581;
t709 = t152 + t154;
t578 = Icges(5,5) * t352;
t153 = Icges(5,6) * t557 - Icges(5,3) * t559 + t578;
t288 = Icges(6,6) * t557;
t580 = Icges(6,4) * t352;
t155 = -Icges(6,2) * t559 - t288 + t580;
t719 = t153 + t155;
t684 = t727 * t557 - t726 * t559 + t721;
t317 = Icges(4,4) * t329;
t250 = Icges(4,1) * t328 + t317;
t696 = -t250 + t712;
t717 = (-Icges(6,4) - Icges(5,5)) * t329 + (-Icges(5,4) + Icges(6,5)) * t328;
t574 = Icges(4,6) * t352;
t146 = Icges(4,4) * t557 - Icges(4,2) * t559 - t574;
t290 = Icges(5,6) * t559;
t582 = Icges(5,4) * t352;
t157 = Icges(5,2) * t557 - t290 + t582;
t716 = t146 * t328 - t157 * t329;
t238 = Icges(6,2) * t328 + t571;
t431 = -Icges(4,2) * t328 + t317;
t697 = -t238 + t424 + t431;
t427 = Icges(5,2) * t329 - t573;
t715 = t251 + t421 + t427;
t242 = Icges(4,5) * t328 + Icges(4,6) * t329;
t691 = t242 + t717;
t714 = -t149 * t557 - t152 * t559;
t713 = t148 * t329 - t153 * t328 - t716;
t248 = Icges(4,2) * t329 + t584;
t688 = -t239 + t248;
t711 = -t423 - t688;
t410 = t248 * t328 - t250 * t329;
t692 = t722 * t328 - t712 * t329 - t410;
t440 = -t150 * t557 - t154 * t559 + t158 * t352;
t147 = Icges(4,6) * t351 + t352 * t431;
t291 = Icges(5,6) * t558;
t583 = Icges(5,4) * t351;
t156 = -Icges(5,2) * t556 + t291 + t583;
t707 = -t689 * t352 - t714;
t660 = -t147 * t559 - t156 * t557 + t707;
t708 = -t440 + t660;
t706 = t717 * t352;
t705 = t723 * t351 + t710 * t556 + t709 * t558;
t585 = Icges(6,1) * t352;
t159 = -Icges(6,4) * t559 - Icges(6,5) * t557 + t585;
t137 = t351 * t159;
t704 = t684 * t351 + t720 * t556 + t719 * t558 + t137;
t703 = t146 + t719;
t702 = t147 - t709;
t701 = t157 - t720;
t700 = -t156 + t710;
t699 = t697 * qJD(3);
t698 = t715 * qJD(3);
t677 = t244 + t724;
t694 = (t248 - t722) * qJD(3);
t693 = t696 * qJD(3);
t621 = t691 * t351;
t690 = t712 * t557 - t559 * t722 + t706;
t687 = t147 * t328 + t156 * t329;
t554 = t352 * t159;
t417 = t151 * t329 + t155 * t328;
t633 = t351 * t417;
t60 = t554 - t633;
t686 = t351 * t713 + t352 * t684 + t60;
t645 = -t146 * t558 + t157 * t556 - t704;
t644 = -t147 * t558 - t156 * t556 + t705;
t685 = t352 * t692 + t621;
t601 = rSges(6,1) + pkin(4);
t593 = rSges(6,2) * t328;
t253 = rSges(6,3) * t329 + t593;
t683 = -qJ(5) * t329 - t253;
t568 = qJ(4) * t328;
t256 = pkin(3) * t329 + t568;
t209 = t256 * t351;
t181 = qJD(1) * t209;
t496 = qJD(3) * t352;
t467 = t329 * t496;
t266 = qJ(4) * t467;
t333 = t352 * qJ(2);
t278 = pkin(1) * t351 - t333;
t349 = cos(pkin(7));
t318 = pkin(2) * t349 + pkin(1);
t350 = -pkin(6) - qJ(2);
t321 = t352 * t350;
t513 = -t351 * t318 - t321;
t142 = t278 + t513;
t261 = qJD(1) * t278;
t631 = qJD(1) * t142 - t261;
t682 = t266 + t181 - t631;
t681 = t694 * t351 + (t238 * t352 - t147 + t152 + t581) * qJD(1);
t680 = t694 * t352 + (t697 * t351 - t574 + t578 + t580) * qJD(1);
t679 = -t693 * t351 + (-t352 * t427 + t583 - t710) * qJD(1);
t678 = t693 * t352 + (-t351 * t715 + t577 + t579 - t582) * qJD(1);
t643 = t701 * t328 + t703 * t329;
t642 = t700 * t328 + t702 * t329;
t676 = t698 * t329 - t699 * t328 + (t696 * t328 + t711 * t329) * qJD(3) + t691 * qJD(1);
t675 = t691 * qJD(3);
t674 = t417 - t713;
t673 = t709 * t328 + t710 * t329 - t687;
t672 = t692 * qJD(1) - t677 * qJD(3);
t671 = (Icges(5,3) * t556 + t688 * t352 + t291 - t700) * t351 + (t287 - t290 - t298 + (-Icges(4,2) - Icges(6,2) - Icges(5,3)) * t557 + t701) * t352;
t659 = rSges(6,3) + qJ(5);
t670 = t685 * qJD(1);
t669 = (Icges(6,3) * t559 - t288 + t703) * t352 + (-Icges(6,3) * t558 + t289 - t702) * t351;
t668 = (t708 * t351 - t686 * t352) * qJD(3);
t667 = (t351 * t644 - t352 * t645) * qJD(3);
t666 = t723 * qJD(1);
t665 = -t696 + t697;
t664 = t715 + t711;
t560 = t242 * t352;
t68 = -t351 * t410 - t560;
t661 = (-t68 + t690) * qJD(1);
t529 = -rSges(6,2) * t559 + t352 * t601 - t659 * t557;
t589 = t328 * rSges(6,3);
t434 = rSges(6,2) * t329 - t589;
t210 = t434 * t351;
t654 = qJ(5) * t556 + t601 * t351;
t653 = t684 + t687;
t468 = t328 * t496;
t500 = qJD(1) * t351;
t652 = t329 * t500 + t468;
t651 = -t661 + t668;
t650 = t667 + t670;
t649 = -t351 * t672 + t352 * t676;
t648 = t351 * t676 + t352 * t672;
t647 = qJD(3) * t673 + t328 * t678 - t329 * t680;
t646 = qJD(3) * t674 + t328 * t679 + t329 * t681;
t641 = t679 * t329 - t681 * t328 + t643 * qJD(3) + (t159 + t684) * qJD(1);
t640 = -qJD(3) * t642 + t328 * t680 + t329 * t678 + t666;
t639 = -t554 + t705;
t638 = -t675 * t352 + (-t351 * t677 + t585 - t673 + t721) * qJD(1);
t637 = qJD(1) * t674 - t351 * t675 + t666;
t495 = qJD(4) * t329;
t169 = qJD(3) * t256 - t495;
t636 = qJD(3) * t683 - t169;
t635 = 0.2e1 * qJD(3);
t309 = pkin(3) * t556;
t214 = qJ(4) * t558 + t309;
t497 = qJD(3) * t351;
t438 = t209 * t497 + t214 * t496 - t495;
t492 = qJD(5) * t328;
t531 = rSges(6,2) * t558 + rSges(6,3) * t556 + t654;
t41 = t492 + (-t351 * t529 + t352 * t531) * qJD(3) + t438;
t634 = qJD(3) * t41;
t315 = qJD(4) * t328;
t590 = t328 * rSges(5,2);
t435 = rSges(5,3) * t329 + t590;
t632 = (-qJD(3) * t435 - t315) * t351;
t337 = t351 * rSges(4,3);
t171 = rSges(4,1) * t556 - rSges(4,2) * t558 + t337;
t300 = t352 * t318;
t449 = -t350 * t351 + t300;
t630 = t171 + t449;
t332 = t351 * qJ(2);
t280 = t352 * pkin(1) + t332;
t595 = rSges(3,2) * sin(pkin(7));
t597 = rSges(3,1) * t349;
t407 = t351 * rSges(3,3) + (-t595 + t597) * t352;
t629 = t280 + t407;
t499 = qJD(1) * t352;
t627 = rSges(6,2) * t467 + t601 * t499;
t507 = t351 ^ 2 + t352 ^ 2;
t252 = pkin(3) * t328 - qJ(4) * t329;
t567 = qJ(5) * t328;
t454 = -t434 + t567;
t442 = -t252 - t454;
t403 = t442 * qJD(3);
t489 = qJD(5) * t352;
t274 = t329 * t489;
t493 = qJD(4) * t352;
t276 = t328 * t493;
t330 = qJD(2) * t351;
t514 = t276 + t330;
t475 = t274 + t514;
t625 = t352 * t403 + t475;
t486 = -pkin(3) - t659;
t624 = t486 * t329 - t318 - t568 - t593;
t623 = qJD(1) * t529;
t620 = t560 + t706;
t619 = t328 * t671 + t669 * t329;
t618 = (-t665 * t328 + t664 * t329) * qJD(1);
t617 = t677 * qJD(1);
t607 = -qJD(3) * (-t256 + t683) - t492 + t636;
t606 = 0.2e1 * t329;
t605 = m(5) / 0.2e1;
t604 = m(6) / 0.2e1;
t603 = t351 / 0.2e1;
t602 = -t352 / 0.2e1;
t599 = qJD(1) / 0.2e1;
t598 = pkin(1) - t318;
t596 = rSges(4,1) * t329;
t594 = rSges(5,2) * t329;
t592 = rSges(5,3) * t328;
t255 = rSges(4,1) * t328 + rSges(4,2) * t329;
t213 = t255 * t352;
t331 = qJD(2) * t352;
t471 = t255 * t497;
t64 = qJD(1) * t630 - t331 - t471;
t591 = t213 * t64;
t340 = t351 * rSges(5,1);
t170 = rSges(4,1) * t557 - rSges(4,2) * t559 - t352 * rSges(4,3);
t470 = t255 * t496;
t439 = t330 - t470;
t542 = t142 - t278;
t63 = (-t170 + t542) * qJD(1) + t439;
t588 = t351 * t63;
t469 = t328 * t497;
t265 = qJ(5) * t469;
t490 = qJD(5) * t351;
t549 = qJD(3) * t210 + t329 * t490 - t265 + (t253 * t352 + t654) * qJD(1);
t473 = t328 * t500;
t548 = -rSges(6,2) * t473 - t652 * t659 + t274 + t627;
t228 = qJD(1) * t280 - t331;
t312 = t350 * t500;
t543 = t312 - (-t352 * t598 - t332) * qJD(1) - t228;
t258 = t592 - t594;
t532 = -t258 * qJD(3) - t169;
t173 = -rSges(5,2) * t556 + rSges(5,3) * t558 + t340;
t530 = -t173 - t214;
t528 = t351 * t209 + t352 * t214;
t211 = t252 * t352;
t494 = qJD(4) * t351;
t527 = -qJD(1) * t211 + t329 * t494;
t488 = qJD(1) * qJD(2);
t322 = qJ(2) * t499;
t508 = t322 + t330;
t526 = qJD(1) * (-pkin(1) * t500 + t508) + t351 * t488;
t217 = t252 * t497;
t525 = -t217 - t331;
t518 = -t252 + t435;
t517 = -t256 - t258;
t515 = rSges(4,2) * t473 + rSges(4,3) * t499;
t313 = t351 * t595;
t511 = rSges(3,3) * t499 + qJD(1) * t313;
t510 = t312 + t331;
t509 = t352 * rSges(3,3) + t313;
t498 = qJD(3) * t328;
t491 = qJD(5) * t329;
t487 = qJD(3) * qJD(4);
t86 = -pkin(3) * t652 - qJ(4) * t473 + t266 + t276;
t185 = t328 * t499 + t329 * t497;
t273 = pkin(3) * t469;
t87 = qJ(4) * t185 + qJD(1) * t309 + t328 * t494 - t273;
t485 = t209 * t499 + t351 * t87 + t352 * t86;
t484 = t351 * t597;
t482 = -t87 + t543;
t481 = qJD(1) * (-t322 + (t351 * t598 - t321) * qJD(1)) + t526;
t480 = -t209 + t542;
t206 = t252 * t351;
t479 = -t206 * t497 - t211 * t496 + t315;
t478 = -t214 - t531;
t320 = t352 * t488;
t464 = t329 * t487;
t477 = qJD(1) * t217 + t352 * t464 + t320;
t476 = t273 + t510;
t474 = t300 + t214;
t465 = -pkin(1) - t597;
t461 = -t497 / 0.2e1;
t458 = t496 / 0.2e1;
t456 = rSges(5,1) * t352 - rSges(5,3) * t559;
t455 = t518 * t352;
t448 = t328 * t487 + t87 * t497 + (t181 + t86) * t496;
t447 = rSges(5,1) * t499 + t652 * rSges(5,2) + rSges(5,3) * t467;
t436 = -rSges(4,2) * t328 + t596;
t433 = -t351 * t64 - t352 * t63;
t408 = -t491 - t315;
t405 = t351 * t464 + t481 + (t276 + t86) * qJD(1);
t208 = t255 * t351;
t207 = t435 * t351;
t82 = (t170 * t351 + t171 * t352) * qJD(3);
t374 = (-0.2e1 * t492 + t636) * qJD(3);
t2 = t374 * t352 + ((qJD(3) * t454 + t408) * t351 + t482 - t549) * qJD(1) + t477;
t3 = t374 * t351 + ((t403 + t491) * t352 + t548) * qJD(1) + t405;
t386 = t2 * t352 + t3 * t351 + t634;
t384 = qJD(3) * t455 + t514;
t379 = -t256 - t318 - t592;
t277 = t329 * t493;
t231 = t436 * qJD(3);
t218 = t252 * t500;
t215 = t434 * t352;
t212 = t435 * t352;
t187 = t484 - t509;
t186 = t467 - t473;
t184 = t507 * t498;
t175 = rSges(5,2) * t557 + t456;
t117 = qJD(1) * t629 - t331;
t116 = t330 + (-t187 - t278) * qJD(1);
t113 = -rSges(5,3) * t473 + t447;
t111 = qJD(3) * t207 + (t258 * t352 + t340) * qJD(1);
t109 = -qJD(3) * t208 + (t352 * t436 + t337) * qJD(1);
t108 = -rSges(4,1) * t652 - rSges(4,2) * t467 + t515;
t89 = t320 + (-qJD(1) * t407 - t228) * qJD(1);
t88 = qJD(1) * (-qJD(1) * t484 + t511) + t526;
t48 = (t173 * t352 - t175 * t351) * qJD(3) + t438;
t47 = -t632 + (-t530 + t449) * qJD(1) + t525;
t46 = (t175 + t480) * qJD(1) + t384;
t43 = -t231 * t496 + t320 + (-t109 + t471 + t543) * qJD(1);
t42 = -t231 * t497 + (t108 - t470) * qJD(1) + t481;
t34 = -t265 + (qJD(3) * t434 - t408) * t351 + (-t478 + t449) * qJD(1) + t525;
t33 = (t480 + t529) * qJD(1) + t625;
t18 = t532 * t496 + (-t111 + t482 + t632) * qJD(1) + t477;
t17 = qJD(1) * t113 + (qJD(1) * t455 + t351 * t532) * qJD(3) + t405;
t4 = (t111 * t351 + t113 * t352 + (-t175 * t352 + t351 * t530) * qJD(1)) * qJD(3) + t448;
t1 = (t491 + t548 * t352 + t549 * t351 + (t351 * t478 - t352 * t529) * qJD(1)) * qJD(3) + t448;
t5 = [(t33 * (t265 + t476) + t3 * (t474 + t531) + (-t3 * t350 + (t408 + (t589 + (-rSges(6,2) - qJ(4)) * t329) * qJD(3)) * t33) * t351 + (-t209 + t513 + t529) * t2 + (t486 * t498 * t352 + t33 + t475 - t623 - t625 + t627 + t682) * t34) * m(6) + (t18 * (t456 + t513) + t46 * t476 + t17 * (t173 + t474) + (t18 * (-t256 + t594) - t17 * t350 + (-t315 + (-t590 + (-rSges(5,3) - qJ(4)) * t329) * qJD(3)) * t46) * t351 + (-pkin(3) * t468 - t384 + t447 + t46 + t514 + t682) * t47) * m(5) + (t43 * (-t170 + t513) + t63 * t510 + t42 * t630 + (t255 * t588 - t591) * qJD(3) + (t330 + t515 - t439 + t63 - t631) * t64) * m(4) + (t89 * (t351 * t465 + t333 + t509) + t116 * t331 + t88 * t629 + (t508 + t511 + t116 + t261 - t330) * t117) * m(3) + (((t60 + t633 + t639) * t351 + ((t689 + t716) * t352 + t660 + t704 + t714) * t352) * qJD(3) + t670) * t458 + (((t653 * t352 - t639 + t644) * t352 + (t653 * t351 + t137 + t440 + t645 - t707) * t351) * qJD(3) + t651 + t661) * t461 + (t647 + t649) * t497 / 0.2e1 - (-t646 + t648 + t650) * t496 / 0.2e1 + ((t68 + t643) * t351 + (t642 + t685) * t352) * qJD(3) * t599 + (m(6) * ((t33 * t624 - t34 * t350) * t352 + (-t33 * t601 + t34 * t624) * t351) + t699 * t329 + t698 * t328 + t692 * qJD(3) + (-t175 * t47 + (-t46 * rSges(5,1) + t379 * t47) * t351 + (t46 * (t379 + t594) - t47 * t350) * t352) * m(5) + (t170 * t64 + (-t63 * rSges(4,3) + t64 * (-t318 - t596)) * t351 + (t63 * (-t318 - t436) - t64 * t350) * t352) * m(4) + (t116 * (t465 + t595) * t352 + (t116 * (-rSges(3,3) - qJ(2)) + t117 * t465) * t351 + t187 * t117) * m(3) + t690 * t461) * qJD(1); 0.2e1 * (t2 * t603 + t3 * t602) * m(6) + 0.2e1 * (t17 * t602 + t18 * t603) * m(5) + 0.2e1 * (t42 * t602 + t43 * t603) * m(4) + 0.2e1 * (t602 * t88 + t603 * t89) * m(3); (t1 * t528 + (-t1 * t529 + t3 * t442) * t351 + (t1 * t531 + t2 * t442) * t352 + (t328 * t490 - t527 + t607 * t351 + (qJ(5) * t558 + t352 * t442 - t215) * qJD(1)) * t34 + (-t206 * qJD(1) + t328 * t489 + t607 * t352 + t218 - t277) * t33 + (t485 + (qJD(1) * t478 + t549) * t351 + (t548 - t623) * t352 - t479 - t491 - (t210 * t351 + t215 * t352 - t507 * t567) * qJD(3)) * t41) * m(6) + (t46 * t218 + t4 * t528 + t48 * t485 + (t18 * t518 + t46 * t532 + t4 * t173 + t48 * t113 + (-t48 * t175 + t47 * t518) * qJD(1)) * t352 + (t17 * t518 + t47 * t532 - t4 * t175 + t48 * t111 + (-t435 * t46 + t48 * t530) * qJD(1)) * t351 - t46 * (t277 + (t206 - t207) * qJD(1)) - t47 * (qJD(1) * t212 + t527) - t48 * t479 - ((t48 * t212 + t46 * t517) * t352 + (t48 * t207 + t47 * t517) * t351) * qJD(3)) * m(5) + (0.2e1 * t82 * (t108 * t352 + t109 * t351 + (t170 * t352 - t171 * t351) * qJD(1)) + t433 * t231 + (-t42 * t351 - t43 * t352 + (-t352 * t64 + t588) * qJD(1)) * t255 - (t208 * t63 - t591) * qJD(1) - (t82 * (-t208 * t351 - t213 * t352) + t433 * t436) * qJD(3)) * m(4) - ((t669 * t328 - t329 * t671) * qJD(3) + (t664 * t328 + t665 * t329) * qJD(1)) * qJD(1) / 0.2e1 + (t646 * t352 + t647 * t351 + (t643 * t351 + t642 * t352) * qJD(1)) * t599 + ((-t497 * t620 + t617) * t351 + ((t351 * t621 + t619) * qJD(3) + t618) * t352) * t461 + ((-t496 * t621 - t617) * t352 + ((t352 * t620 + t619) * qJD(3) + t618) * t351) * t458 + (t649 * qJD(1) + ((t644 * qJD(1) + t641 * t352) * t352 + (t638 * t351 + t645 * qJD(1) + (-t637 + t640) * t352) * t351) * t635) * t603 + (t648 * qJD(1) + ((t708 * qJD(1) + t637 * t352) * t352 + (t640 * t351 + t686 * qJD(1) + (-t638 + t641) * t352) * t351) * t635) * t602 + (t651 + t668) * t500 / 0.2e1 + (t650 + t667) * t499 / 0.2e1; -m(5) * (t184 * t48 + t185 * t47 + t186 * t46) - m(6) * (t184 * t41 + t185 * t34 + t186 * t33) + ((t46 * t496 + t47 * t497 - t4) * t605 + (t33 * t496 + t34 * t497 - t1) * t604) * t606 + 0.2e1 * ((qJD(3) * t48 + t17 * t351 + t18 * t352 - t46 * t500 + t47 * t499) * t605 + (-t33 * t500 + t34 * t499 + t386) * t604) * t328; m(6) * t1 * t328 + (t386 * t604 - m(6) * t507 * t634 / 0.2e1) * t606;];
tauc = t5(:);
