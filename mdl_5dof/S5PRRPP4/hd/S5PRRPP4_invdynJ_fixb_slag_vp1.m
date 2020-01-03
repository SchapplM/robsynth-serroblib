% Calculate vector of inverse dynamics joint torques for
% S5PRRPP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:47
% EndTime: 2019-12-31 17:41:11
% DurationCPUTime: 21.66s
% Computational Cost: add. (10061->599), mult. (12276->712), div. (0->0), fcn. (9541->4), ass. (0->327)
t684 = Icges(5,6) - Icges(6,6);
t683 = Icges(5,1) + Icges(6,1);
t682 = Icges(5,4) - Icges(6,5);
t681 = Icges(6,2) + Icges(5,3);
t327 = sin(qJ(3));
t319 = Icges(6,4) * t327;
t328 = cos(qJ(3));
t388 = Icges(6,1) * t328 + t319;
t318 = Icges(5,5) * t327;
t389 = Icges(5,1) * t328 + t318;
t680 = t388 + t389;
t326 = pkin(7) + qJ(2);
t314 = sin(t326);
t512 = t314 * t327;
t282 = Icges(4,4) * t512;
t511 = t314 * t328;
t315 = cos(t326);
t539 = Icges(4,5) * t315;
t152 = Icges(4,1) * t511 - t282 - t539;
t148 = Icges(6,5) * t315 + t314 * t388;
t150 = -Icges(5,4) * t315 + t314 * t389;
t670 = -t148 - t150;
t664 = -t152 + t670;
t542 = Icges(4,4) * t327;
t261 = Icges(4,1) * t328 - t542;
t153 = Icges(4,5) * t314 + t261 * t315;
t674 = t314 * t682 + t680 * t315;
t663 = t153 + t674;
t509 = t315 * t328;
t679 = (Icges(6,4) + Icges(5,5)) * t509;
t678 = t684 * t314;
t538 = Icges(5,5) * t328;
t247 = Icges(5,3) * t327 + t538;
t138 = -Icges(5,6) * t315 + t247 * t314;
t540 = Icges(6,4) * t328;
t251 = Icges(6,2) * t327 + t540;
t142 = Icges(6,6) * t315 + t251 * t314;
t672 = t138 + t142;
t510 = t315 * t327;
t671 = t510 * t681 + t678 + t679;
t249 = Icges(4,5) * t328 - Icges(4,6) * t327;
t141 = Icges(4,3) * t314 + t249 * t315;
t253 = Icges(5,4) * t328 + Icges(5,6) * t327;
t145 = Icges(5,2) * t314 + t253 * t315;
t677 = t141 + t145;
t645 = -t328 * t681 + t318 + t319;
t659 = -Icges(4,2) * t328 - t542 + t645;
t320 = Icges(4,4) * t328;
t387 = -Icges(4,2) * t327 + t320;
t669 = t247 + t251;
t676 = -t387 + t669;
t260 = Icges(4,1) * t327 + t320;
t643 = t327 * t683 + t260 - t538 - t540;
t662 = (-Icges(4,6) + t684) * t328 + (-Icges(4,5) - t682) * t327;
t673 = t261 + t680;
t528 = Icges(4,3) * t315;
t140 = Icges(4,5) * t511 - Icges(4,6) * t512 - t528;
t533 = Icges(4,6) * t315;
t146 = Icges(4,4) * t511 - Icges(4,2) * t512 - t533;
t520 = t146 * t327;
t375 = -t152 * t328 + t520;
t378 = t142 * t327 + t148 * t328;
t144 = -Icges(5,2) * t315 + t253 * t314;
t521 = t144 * t315;
t245 = Icges(6,5) * t328 + Icges(6,6) * t327;
t136 = Icges(6,3) * t315 + t245 * t314;
t523 = t136 * t315;
t382 = t138 * t327 + t150 * t328;
t609 = t314 * t382;
t623 = t314 * t378 - t521 + t523 + t609;
t599 = -t140 * t315 - t314 * t375 + t623;
t147 = Icges(4,6) * t314 + t315 * t387;
t120 = t153 * t511;
t420 = t141 * t315 - t120;
t51 = -t147 * t512 - t420;
t137 = -Icges(6,3) * t314 + t245 * t315;
t135 = t315 * t137;
t622 = -t145 * t315 + t674 * t511 + t512 * t671 + t135;
t598 = t51 + t622;
t668 = t677 * t314 + t663 * t509 + t671 * t510;
t133 = t314 * t144;
t667 = -t314 * t140 + t664 * t509 - t672 * t510 - t133;
t666 = -t146 + t672;
t665 = -t147 + t671;
t661 = t659 * qJD(3);
t660 = t643 * qJD(3);
t655 = t327 * t659 + t643 * t328;
t524 = t136 * t314;
t597 = -t146 * t510 - t524 - t667;
t596 = -t137 * t314 - t147 * t510 + t668;
t658 = t676 * qJD(3);
t657 = t673 * qJD(3);
t656 = t245 - t249 - t253;
t654 = -t327 * t643 + t328 * t659;
t592 = t662 * t315;
t593 = t662 * t314;
t546 = -rSges(6,3) - qJ(5);
t620 = pkin(4) * t509 + t314 * t546;
t486 = rSges(6,1) * t509 + rSges(6,2) * t510 + t620;
t653 = -t661 * t315 + (t314 * t387 - t533 - t672) * qJD(2);
t652 = -t661 * t314 + (t315 * t669 - t147 + t678) * qJD(2);
t651 = -t660 * t315 + (-t261 * t314 + t539 + t670) * qJD(2);
t650 = -qJD(2) * t663 + t314 * t660;
t613 = t314 * t655 + t592;
t649 = -t327 * t663 + t328 * t665;
t595 = t327 * t664 + t328 * t666;
t612 = t315 * t655 - t593;
t648 = t662 * qJD(3);
t519 = t147 * t327;
t647 = t327 * t671 + t328 * t663 - t519;
t646 = t375 - t378 - t382;
t644 = (t136 - t144) * qJD(2);
t642 = -qJD(2) * t662 + t654 * qJD(3) + t658 * t327 + t657 * t328;
t641 = t596 * t314 - t315 * t597;
t640 = t598 * t314 - t599 * t315;
t639 = (-t137 + t677) * qJD(2);
t638 = t655 * qJD(2) + qJD(3) * t656;
t637 = t315 ^ 2;
t636 = t612 * qJD(2);
t635 = qJD(3) * t649 + t327 * t653 + t328 * t651 + t639;
t634 = -qJD(2) * t140 - qJD(3) * t595 - t327 * t652 + t328 * t650 + t644;
t633 = t613 * qJD(2);
t632 = t643 - t676;
t631 = t659 + t673;
t630 = t315 * t659 + t663;
t629 = -t260 * t315 - t510 * t683 + t665 + t679;
t628 = qJD(2) * t646 + t314 * t648 + t639;
t627 = t644 + t648 * t315 + (-t249 * t314 + t528 - t647) * qJD(2);
t446 = qJD(2) * qJD(3);
t215 = -qJDD(3) * t315 + t314 * t446;
t451 = qJD(4) * t328;
t317 = t327 * qJ(4);
t605 = t328 * pkin(3) + t317;
t206 = qJD(3) * t605 - t451;
t322 = t327 * rSges(6,2);
t606 = t328 * rSges(6,1) + t322;
t480 = -qJD(3) * t606 - t206;
t561 = pkin(4) * t328;
t348 = (-t561 * qJD(3) + t480) * qJD(3);
t316 = qJD(4) * t327;
t430 = t314 * t316;
t449 = qJD(5) * t315;
t361 = -t430 - t449;
t199 = t605 * t314;
t309 = t315 * pkin(6);
t220 = pkin(2) * t314 - t309;
t482 = -t199 - t220;
t577 = -rSges(6,3) * t315 - t314 * t606;
t488 = pkin(4) * t511 + qJ(5) * t315 - t577;
t413 = t482 - t488;
t555 = rSges(6,1) * t327;
t268 = -rSges(6,2) * t328 + t555;
t562 = pkin(4) * t327;
t421 = t268 + t562;
t267 = pkin(3) * t327 - qJ(4) * t328;
t445 = qJD(3) * qJD(4);
t621 = qJDD(4) * t327 + t328 * t445;
t441 = t215 * t267 + t315 * t621;
t453 = qJD(3) * t327;
t434 = t314 * t453;
t241 = pkin(4) * t434;
t417 = -t241 + t449;
t455 = qJD(3) * t314;
t506 = -t268 * t455 + t417 + (t315 * t606 + t620) * qJD(2);
t221 = t315 * pkin(2) + t314 * pkin(6);
t211 = t221 * qJD(2);
t452 = qJD(3) * t328;
t456 = qJD(2) * t327;
t175 = t314 * t452 + t315 * t456;
t294 = pkin(3) * t509;
t242 = pkin(3) * t434;
t405 = -t242 + t430;
t80 = qJ(4) * t175 + qJD(2) * t294 + t405;
t610 = -t211 - t80;
t2 = -qJDD(5) * t314 + t421 * t215 + t348 * t315 + t413 * qJDD(2) + (t361 - t506 + t610) * qJD(2) + t441;
t626 = t2 - g(1);
t214 = qJDD(3) * t314 + t315 * t446;
t204 = t315 * t317 + t294;
t265 = t315 * t316;
t457 = qJD(2) * t315;
t297 = pkin(6) * t457;
t458 = qJD(2) * t314;
t481 = qJD(2) * (-pkin(2) * t458 + t297) + qJDD(2) * t221;
t432 = t315 * t453;
t349 = -t328 * t458 - t432;
t436 = t314 * t456;
t431 = t315 * t452;
t237 = qJ(4) * t431;
t478 = t237 + t265;
t79 = pkin(3) * t349 - qJ(4) * t436 + t478;
t351 = qJDD(2) * t204 + t481 + t621 * t314 + (t265 + t79) * qJD(2);
t411 = -t267 - t421;
t239 = rSges(6,2) * t431;
t450 = qJD(5) * t314;
t507 = -rSges(6,1) * t432 + pkin(4) * t349 - qJ(5) * t457 + qJD(2) * t577 + t239 - t450;
t3 = qJDD(5) * t315 + t486 * qJDD(2) + t507 * qJD(2) + t411 * t214 + (-qJD(2) * qJD(5) + t348) * t314 + t351;
t625 = t3 - g(2);
t443 = -rSges(6,1) - pkin(3) - pkin(4);
t624 = t443 * t328 - pkin(2);
t403 = t328 * rSges(5,1) + t327 * rSges(5,3);
t468 = -t605 - t403;
t619 = qJD(3) * t640 + t633;
t618 = qJD(3) * t641 + t636;
t617 = qJD(3) * t646 + t327 * t650 + t328 * t652;
t616 = qJD(3) * t647 + t327 * t651 - t328 * t653;
t615 = -t314 * t638 + t315 * t642;
t614 = t314 * t642 + t315 * t638;
t611 = t521 + t668;
t217 = qJD(2) * t220;
t608 = -qJD(2) * t199 - t217;
t416 = t265 - t450;
t454 = qJD(3) * t315;
t604 = t411 * t454 + t416;
t603 = -t317 - t322 + t624;
t305 = t314 * rSges(5,2);
t166 = rSges(5,1) * t509 + rSges(5,3) * t510 + t305;
t412 = t204 + t221;
t74 = t412 + t166;
t600 = qJD(2) * t488;
t591 = t605 + t606 + t561;
t590 = (-t327 * t632 + t328 * t631) * qJD(2);
t589 = -t327 * t630 + t328 * t629;
t588 = t628 * t637 + (t635 * t314 + (-t627 + t634) * t315) * t314;
t587 = t634 * t637 + (t627 * t314 + (-t628 + t635) * t315) * t314;
t586 = t656 * qJD(2);
t585 = -Icges(4,2) * t511 + t314 * t645 - t282 - t664;
t584 = t314 * t643 - t666;
t573 = -pkin(4) * t452 + qJD(3) * t591 + t480;
t572 = t327 * t585 + t328 * t584;
t571 = m(5) / 0.2e1;
t570 = m(6) / 0.2e1;
t569 = m(2) + m(3);
t568 = -m(5) - m(6);
t567 = t214 / 0.2e1;
t566 = t215 / 0.2e1;
t563 = -rSges(5,1) - pkin(3);
t560 = g(2) * t314;
t557 = rSges(4,1) * t328;
t556 = rSges(5,1) * t327;
t270 = rSges(4,1) * t327 + rSges(4,2) * t328;
t203 = t270 * t315;
t304 = t314 * rSges(4,3);
t167 = rSges(4,1) * t509 - rSges(4,2) * t510 + t304;
t112 = t167 + t221;
t62 = qJD(2) * t112 - t270 * t455;
t552 = t203 * t62;
t435 = t270 * t454;
t467 = rSges(4,2) * t512 + t315 * rSges(4,3);
t164 = rSges(4,1) * t511 - t467;
t487 = -t164 - t220;
t61 = qJD(2) * t487 - t435;
t551 = t314 * t61;
t550 = t315 * t61;
t548 = -rSges(6,2) - qJ(4);
t547 = -rSges(5,3) - qJ(4);
t485 = -t166 - t204;
t484 = t314 * t199 + t315 * t204;
t277 = qJ(4) * t509;
t200 = -pkin(3) * t510 + t277;
t483 = qJD(2) * t200 + t314 * t451;
t479 = -t403 * qJD(3) - t206;
t477 = rSges(5,2) * t457 + rSges(5,3) * t431;
t476 = rSges(4,2) * t436 + rSges(4,3) * t457;
t269 = -rSges(5,3) * t328 + t556;
t469 = -t267 - t269;
t465 = t314 ^ 2 + t637;
t442 = t199 * t457 + t314 * t80 + t315 * t79;
t275 = qJ(4) * t511;
t195 = -pkin(3) * t512 + t275;
t440 = t195 * t455 + t200 * t454 + t316;
t307 = t315 * rSges(5,2);
t163 = t314 * t403 - t307;
t439 = -t163 + t482;
t438 = -t204 - t486;
t437 = t315 * t563;
t429 = -pkin(2) - t557;
t425 = -t455 / 0.2e1;
t424 = t455 / 0.2e1;
t423 = -t454 / 0.2e1;
t422 = t454 / 0.2e1;
t419 = -t140 + t519;
t418 = qJD(3) * t479;
t415 = t315 * t443;
t219 = rSges(3,1) * t315 - rSges(3,2) * t314;
t218 = rSges(3,1) * t314 + rSges(3,2) * t315;
t274 = -rSges(4,2) * t327 + t557;
t391 = -t314 * t62 - t550;
t107 = rSges(4,1) * t349 - rSges(4,2) * t431 + t476;
t198 = t270 * t314;
t110 = -qJD(3) * t198 + (t274 * t315 + t304) * qJD(2);
t384 = t107 * t315 + t110 * t314;
t373 = t164 * t314 + t167 * t315;
t363 = t199 * t455 + t204 * t454 + qJD(1) - t451;
t362 = t454 * t469 + t265;
t347 = -qJDD(4) * t328 + t214 * t199 + t327 * t445 + t79 * t454 + t80 * t455 + qJDD(1);
t343 = t327 * t547 + t328 * t563 - pkin(2);
t289 = rSges(6,2) * t509;
t288 = rSges(5,3) * t509;
t285 = rSges(6,2) * t511;
t284 = rSges(5,3) * t511;
t266 = t315 * t451;
t233 = t274 * qJD(3);
t213 = t267 * t458;
t212 = t267 * t455;
t202 = -rSges(5,1) * t510 + t288;
t201 = -rSges(6,1) * t510 + t289;
t197 = -rSges(5,1) * t512 + t284;
t196 = -rSges(6,1) * t512 + t285;
t176 = t431 - t436;
t174 = t465 * t453;
t109 = -t269 * t455 + (t315 * t403 + t305) * qJD(2);
t106 = rSges(5,1) * t349 - rSges(5,3) * t436 + t477;
t60 = qJD(3) * t373 + qJD(1);
t45 = -t212 + (-qJD(3) * t269 + t316) * t314 + t74 * qJD(2);
t44 = qJD(2) * t439 + t362;
t43 = (t163 * t314 + t166 * t315) * qJD(3) + t363;
t42 = -t212 + (-qJD(3) * t268 + t316) * t314 + (t221 - t438) * qJD(2) + t417;
t41 = qJD(2) * t413 + t604;
t34 = qJD(2) * t107 + qJDD(2) * t167 - t214 * t270 - t233 * t455 + t481;
t33 = -t233 * t454 + t215 * t270 + t487 * qJDD(2) + (-t110 - t211) * qJD(2);
t32 = (t314 * t488 + t315 * t486) * qJD(3) + t363;
t25 = qJD(3) * t384 + t164 * t214 - t167 * t215 + qJDD(1);
t18 = qJD(2) * t106 + qJDD(2) * t166 + t214 * t469 + t314 * t418 + t351;
t17 = t215 * t269 + t315 * t418 + t439 * qJDD(2) + (-t109 - t430 + t610) * qJD(2) + t441;
t4 = t163 * t214 + t485 * t215 + (t106 * t315 + t109 * t314) * qJD(3) + t347;
t1 = t488 * t214 + t438 * t215 + (t314 * t506 + t315 * t507) * qJD(3) + t347;
t5 = [t569 * qJDD(1) + m(4) * t25 + m(5) * t4 + m(6) * t1 + (-m(4) + t568 - t569) * g(3); -m(3) * (-g(1) * t218 + g(2) * t219) + (((t51 - t120 + (t141 + t520) * t315 + t667) * t315 + (-t609 + (-t137 - t378) * t314 + t611 + t623) * t314) * qJD(3) + t636) * t422 + (t655 * qJD(3) + t657 * t327 - t658 * t328) * qJD(2) + (-(-t41 - t600 + t604 + t608) * t42 + t41 * (t241 + t242 + t361) + t42 * (t237 + t239 + t297 + t416) + (t42 * t327 * t415 + t41 * (t328 * t548 + t555) * t314) * qJD(3) + ((t41 * t603 + t42 * t546) * t315 + (t41 * (-pkin(6) - t546) + t603 * t42) * t314) * qJD(2) + t625 * (t412 + t486) + t626 * (t309 + t546 * t315 + (t327 * t548 + t624) * t314)) * m(6) + (-(-qJD(2) * t163 + t362 - t44 + t608) * t45 - t44 * t405 + t45 * (t297 + t477 + t478) + (t45 * t327 * t437 + t44 * (t328 * t547 + t556) * t314) * qJD(3) + (t44 * t343 * t315 + (t44 * (-rSges(5,2) - pkin(6)) + t45 * (-pkin(2) + t468)) * t314) * qJD(2) + (t18 - g(2)) * t74 + (t17 - g(1)) * (t314 * t343 + t307 + t309)) * m(5) + (-(-qJD(2) * t164 - t217 - t435 - t61) * t62 + t62 * (t297 + t476) + (t270 * t551 - t552) * qJD(3) + ((-pkin(2) - t274) * t550 + (t61 * (-rSges(4,3) - pkin(6)) + t62 * t429) * t314) * qJD(2) + (t34 - g(2)) * t112 + (t33 - g(1)) * (t314 * t429 + t309 + t467)) * m(4) + (m(3) * (t218 ^ 2 + t219 ^ 2) + Icges(3,3) - t654) * qJDD(2) + (t612 - t649) * t567 + (-t595 + t613) * t566 + (((t315 * t419 + t523 + t596 - t611) * t315 + (t314 * t419 - t133 + t135 + t420 + t524 + t597 - t622) * t314) * qJD(3) + t619 - t633) * t425 + (t615 + t616) * t424 + (t614 - t617 + t618) * t423; t641 * t567 + t640 * t566 + (qJD(2) * t615 + t587 * qJD(3) + qJDD(2) * t612 + t596 * t214 + t597 * t215) * t314 / 0.2e1 - (qJD(2) * t614 + t588 * qJD(3) + qJDD(2) * t613 + t598 * t214 + t599 * t215) * t315 / 0.2e1 - (((t314 * t630 - t585 * t315) * t328 + (t314 * t629 + t584 * t315) * t327) * qJD(3) + (t327 * t631 + t328 * t632) * qJD(2)) * qJD(2) / 0.2e1 + (t617 * t315 + t616 * t314 + (-t314 * t595 - t315 * t649) * qJD(2)) * qJD(2) / 0.2e1 + (-t314 * t649 + t315 * t595) * qJDD(2) / 0.2e1 + t619 * t458 / 0.2e1 + t618 * t457 / 0.2e1 + ((t455 * t592 - t586) * t314 + ((t572 * t315 + (t589 - t593) * t314) * qJD(3) + t590) * t315) * t425 + ((t314 * t597 + t315 * t596) * qJD(2) + t587) * t424 + ((t314 * t599 + t315 * t598) * qJD(2) + t588) * t423 + ((t454 * t593 + t586) * t315 + ((t589 * t314 + (t572 - t592) * t315) * qJD(3) + t590) * t314) * t422 + (t1 * t484 + (t1 * t488 + t3 * t411) * t314 + (t1 * t486 + t2 * t411) * t315 - g(1) * (t277 + t289) - g(2) * (t275 + t285) - g(3) * t591 - (g(1) * t415 + t443 * t560) * t327 + (-t483 + t573 * t314 + (pkin(4) * t510 + t315 * t411 - t201) * qJD(2)) * t42 + (t213 - t266 + t573 * t315 + (t268 * t314 + t195 + t196) * qJD(2)) * t41 + (t442 + (qJD(2) * t438 + t506) * t314 + (t507 + t600) * t315 - t440 - (t196 * t314 + t201 * t315 - t465 * t562) * qJD(3)) * t32) * m(6) + (t44 * t213 + t4 * t484 + t43 * t442 + (t17 * t469 + t44 * t479 + t4 * t166 + t43 * t106 + (t43 * t163 + t45 * t469) * qJD(2)) * t315 + (t18 * t469 + t45 * t479 + t4 * t163 + t43 * t109 + (t44 * t269 + t43 * t485) * qJD(2)) * t314 - t44 * (t266 + (-t195 - t197) * qJD(2)) - t45 * (qJD(2) * t202 + t483) - t43 * t440 - ((t43 * t202 + t44 * t468) * t315 + (t43 * t197 + t45 * t468) * t314) * qJD(3) - g(1) * (t277 + t288) - g(2) * (t275 + t284) + g(3) * t468 - (g(1) * t437 + t560 * t563) * t327) * m(5) + (t25 * t373 + t60 * ((t164 * t315 - t167 * t314) * qJD(2) + t384) + t391 * t233 + (-t34 * t314 - t33 * t315 + (-t315 * t62 + t551) * qJD(2)) * t270 - (t198 * t61 - t552) * qJD(2) - (t60 * (-t198 * t314 - t203 * t315) + t391 * t274) * qJD(3) + g(1) * t203 + g(2) * t198 - g(3) * t274) * m(4); t568 * (-g(3) * t328 + (g(1) * t315 + t560) * t327) - m(5) * (t174 * t43 + t175 * t45 + t176 * t44) - m(6) * (t174 * t32 + t175 * t42 + t176 * t41) + 0.2e1 * ((t44 * t454 + t45 * t455 - t4) * t571 + (t41 * t454 + t42 * t455 - t1) * t570) * t328 + 0.2e1 * ((qJD(3) * t43 + t17 * t315 + t18 * t314 - t44 * t458 + t45 * t457) * t571 + (qJD(3) * t32 + t2 * t315 + t3 * t314 - t41 * t458 + t42 * t457) * t570) * t327; (-t314 * t626 + t625 * t315) * m(6);];
tau = t5;
