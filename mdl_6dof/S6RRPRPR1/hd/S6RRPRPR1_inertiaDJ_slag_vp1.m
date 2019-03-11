% Calculate time derivative of joint inertia matrix for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:07:55
% EndTime: 2019-03-09 10:08:30
% DurationCPUTime: 20.07s
% Computational Cost: add. (41921->1034), mult. (35067->1385), div. (0->0), fcn. (32299->12), ass. (0->523)
t395 = qJ(2) + pkin(10);
t384 = qJ(4) + t395;
t374 = cos(t384);
t399 = cos(pkin(11));
t403 = sin(qJ(1));
t579 = t403 * t399;
t398 = sin(pkin(11));
t405 = cos(qJ(1));
t584 = t398 * t405;
t318 = t374 * t579 - t584;
t580 = t403 * t398;
t583 = t399 * t405;
t440 = t374 * t580 + t583;
t373 = sin(t384);
t595 = t373 * t403;
t204 = Icges(6,4) * t318 - Icges(6,2) * t440 + Icges(6,6) * t595;
t206 = Icges(6,1) * t318 - Icges(6,4) * t440 + Icges(6,5) * t595;
t474 = Icges(5,5) * t374 - Icges(5,6) * t373;
t273 = -Icges(5,3) * t405 + t403 * t474;
t319 = -t374 * t584 + t579;
t320 = t374 * t583 + t580;
t202 = Icges(6,5) * t318 - Icges(6,6) * t440 + Icges(6,3) * t595;
t594 = t373 * t405;
t543 = t202 * t594;
t626 = Icges(5,4) * t374;
t479 = -Icges(5,2) * t373 + t626;
t275 = -Icges(5,6) * t405 + t403 * t479;
t627 = Icges(5,4) * t373;
t486 = Icges(5,1) * t374 - t627;
t277 = -Icges(5,5) * t405 + t403 * t486;
t458 = t275 * t373 - t277 * t374;
t677 = t405 * t458;
t689 = -t319 * t204 - t320 * t206 - t403 * t273 - t543 + t677;
t205 = Icges(6,4) * t320 + Icges(6,2) * t319 + Icges(6,6) * t594;
t207 = Icges(6,1) * t320 + Icges(6,4) * t319 + Icges(6,5) * t594;
t274 = Icges(5,3) * t403 + t405 * t474;
t276 = Icges(5,6) * t403 + t405 * t479;
t278 = Icges(5,5) * t403 + t405 * t486;
t457 = t276 * t373 - t278 * t374;
t203 = Icges(6,5) * t320 + Icges(6,6) * t319 + Icges(6,3) * t594;
t541 = t203 * t594;
t688 = t319 * t205 + t320 * t207 + t403 * t274 - t405 * t457 + t541;
t393 = pkin(11) + qJ(6);
t382 = cos(t393);
t394 = qJD(2) + qJD(4);
t380 = sin(t393);
t624 = Icges(7,4) * t382;
t477 = -Icges(7,2) * t380 + t624;
t593 = t374 * t394;
t625 = Icges(7,4) * t380;
t154 = t477 * t593 + (Icges(7,6) * t394 + (-Icges(7,2) * t382 - t625) * qJD(6)) * t373;
t257 = -Icges(7,6) * t374 + t373 * t477;
t484 = Icges(7,1) * t382 - t625;
t258 = -Icges(7,5) * t374 + t373 * t484;
t687 = -t380 * t154 + (-t257 * t382 - t258 * t380) * qJD(6);
t332 = rSges(5,1) * t373 + rSges(5,2) * t374;
t438 = t332 * t394;
t375 = pkin(5) * t399 + pkin(4);
t643 = pkin(4) - t375;
t686 = t373 * t643;
t401 = -pkin(9) - qJ(5);
t578 = qJ(5) + t401;
t685 = t374 * t578;
t684 = t374 * t643;
t402 = sin(qJ(2));
t404 = cos(qJ(2));
t630 = Icges(3,4) * t404;
t483 = -Icges(3,2) * t402 + t630;
t314 = Icges(3,6) * t403 + t405 * t483;
t631 = Icges(3,4) * t402;
t490 = Icges(3,1) * t404 - t631;
t316 = Icges(3,5) * t403 + t405 * t490;
t453 = t314 * t402 - t316 * t404;
t683 = t403 * t453;
t381 = sin(t395);
t383 = cos(t395);
t628 = Icges(4,4) * t383;
t481 = -Icges(4,2) * t381 + t628;
t293 = Icges(4,6) * t403 + t405 * t481;
t629 = Icges(4,4) * t381;
t488 = Icges(4,1) * t383 - t629;
t295 = Icges(4,5) * t403 + t405 * t488;
t455 = t293 * t381 - t295 * t383;
t682 = t403 * t455;
t681 = t403 * t457;
t390 = t404 * pkin(2);
t376 = t390 + pkin(1);
t644 = pkin(1) - t376;
t680 = t403 * t644;
t313 = -Icges(3,6) * t405 + t403 * t483;
t315 = -Icges(3,5) * t405 + t403 * t490;
t454 = t313 * t402 - t315 * t404;
t679 = t405 * t454;
t292 = -Icges(4,6) * t405 + t403 * t481;
t294 = -Icges(4,5) * t405 + t403 * t488;
t456 = t292 * t381 - t294 * t383;
t678 = t405 * t456;
t472 = Icges(7,5) * t382 - Icges(7,6) * t380;
t153 = t472 * t593 + (Icges(7,3) * t394 + (-Icges(7,5) * t380 - Icges(7,6) * t382) * qJD(6)) * t373;
t611 = t257 * t380;
t676 = -t394 * t611 - t153;
t495 = rSges(7,1) * t382 - rSges(7,2) * t380;
t160 = t495 * t593 + (rSges(7,3) * t394 + (-rSges(7,1) * t380 - rSges(7,2) * t382) * qJD(6)) * t373;
t415 = -t373 * t578 - t684;
t675 = -t415 * t394 - t160;
t582 = t403 * t380;
t588 = t382 * t405;
t302 = -t374 * t582 - t588;
t581 = t403 * t382;
t589 = t380 * t405;
t303 = t374 * t581 - t589;
t496 = -t303 * rSges(7,1) - t302 * rSges(7,2);
t183 = rSges(7,3) * t595 - t496;
t304 = -t374 * t589 + t581;
t305 = t374 * t588 + t582;
t184 = t305 * rSges(7,1) + t304 * rSges(7,2) + rSges(7,3) * t594;
t674 = -t403 * t183 - t405 * t184;
t591 = t374 * t405;
t358 = pkin(4) * t591;
t307 = qJ(5) * t594 + t358;
t371 = pkin(5) * t580;
t446 = t375 * t591 - t401 * t594 + t371;
t198 = t446 - t307;
t673 = -t184 - t198;
t243 = t685 - t686;
t259 = -rSges(7,3) * t374 + t373 * t495;
t557 = qJD(1) * t403;
t246 = t259 * t557;
t672 = t243 * t557 + t246;
t671 = -t243 - t259;
t388 = t403 * rSges(4,3);
t637 = rSges(4,2) * t381;
t670 = -t405 * t637 + t388;
t649 = pkin(2) * t402;
t504 = -pkin(3) * t381 - t649;
t447 = -t332 + t504;
t253 = t447 * t403;
t254 = t447 * t405;
t669 = t253 * t403 + t254 * t405;
t668 = qJD(1) * t273;
t499 = -t318 * rSges(6,1) + rSges(6,2) * t440;
t210 = rSges(6,3) * t595 - t499;
t556 = qJD(1) * t405;
t586 = t394 * t403;
t423 = t373 * t556 + t374 * t586;
t539 = t373 * t586;
t251 = qJD(1) * t319 + t398 * t539;
t252 = qJD(1) * t320 - t399 * t539;
t500 = t252 * rSges(6,1) + t251 * rSges(6,2);
t522 = t373 * t557;
t585 = t394 * t405;
t538 = t373 * t585;
t249 = qJD(1) * t440 + t398 * t538;
t250 = -qJD(1) * t318 - t399 * t538;
t537 = t374 * t585;
t528 = t250 * rSges(6,1) + t249 * rSges(6,2) + rSges(6,3) * t537;
t667 = t403 * (rSges(6,3) * t423 + t500) + t405 * (-rSges(6,3) * t522 + t528) + t210 * t556;
t475 = Icges(4,5) * t383 - Icges(4,6) * t381;
t290 = -Icges(4,3) * t405 + t403 * t475;
t476 = Icges(3,5) * t404 - Icges(3,6) * t402;
t311 = -Icges(3,3) * t405 + t403 * t476;
t506 = qJD(1) * t374 - qJD(6);
t666 = t403 * t506 + t538;
t422 = -t522 + t537;
t122 = Icges(6,5) * t250 + Icges(6,6) * t249 + Icges(6,3) * t422;
t123 = Icges(6,5) * t252 + Icges(6,6) * t251 + Icges(6,3) * t423;
t665 = -(t123 * t373 + t202 * t593) * t405 + (t122 * t373 + t203 * t593) * t403;
t664 = t405 * t506 - t539;
t329 = Icges(5,2) * t374 + t627;
t330 = Icges(5,1) * t373 + t626;
t452 = t329 * t373 - t330 * t374;
t662 = qJD(1) * t452 + t474 * t394;
t400 = -qJ(3) - pkin(7);
t392 = -pkin(8) + t400;
t562 = t392 - t400;
t648 = pkin(3) * t383;
t348 = t376 + t648;
t569 = t348 - t376;
t247 = t403 * t569 + t405 * t562;
t197 = -pkin(5) * t584 + t403 * t415;
t339 = qJ(5) * t537;
t343 = pkin(4) * t539;
t646 = pkin(5) * t398;
t545 = qJD(1) * t646;
t570 = t401 * t522 + t405 * t545;
t592 = t374 * t401;
t597 = t373 * t375;
t614 = qJ(5) * t373;
t507 = -qJD(6) * t374 + qJD(1);
t450 = t507 * t405;
t165 = t380 * t666 + t382 * t450;
t166 = t380 * t450 - t382 * t666;
t532 = t166 * rSges(7,1) + t165 * rSges(7,2) + rSges(7,3) * t537;
t98 = -rSges(7,3) * t522 + t532;
t451 = t403 * t507;
t167 = -t380 * t664 + t382 * t451;
t168 = t380 * t451 + t382 * t664;
t497 = t168 * rSges(7,1) + t167 * rSges(7,2);
t99 = rSges(7,3) * t423 + t497;
t661 = (t183 + t197) * t556 + (t98 - t339 + (-t592 + t686) * t585 + (t614 + t684) * t557 + t570) * t405 + (t99 + t343 + (-t597 - t685) * t586 + (t405 * t415 + t371) * qJD(1)) * t403;
t660 = 2 * m(3);
t659 = 2 * m(4);
t658 = 2 * m(5);
t657 = 2 * m(6);
t656 = 2 * m(7);
t396 = t403 ^ 2;
t397 = t405 ^ 2;
t655 = m(6) / 0.2e1;
t654 = m(7) / 0.2e1;
t653 = t403 / 0.2e1;
t652 = -t405 / 0.2e1;
t362 = rSges(3,1) * t402 + rSges(3,2) * t404;
t651 = m(3) * t362;
t650 = m(5) * t332;
t647 = pkin(4) * t374;
t645 = t403 * pkin(7);
t391 = t405 * pkin(7);
t642 = -pkin(7) - t400;
t641 = rSges(3,1) * t404;
t640 = rSges(4,1) * t383;
t639 = rSges(5,1) * t374;
t638 = rSges(3,2) * t402;
t636 = rSges(3,3) * t405;
t174 = Icges(7,5) * t303 + Icges(7,6) * t302 + Icges(7,3) * t595;
t176 = Icges(7,4) * t303 + Icges(7,2) * t302 + Icges(7,6) * t595;
t178 = Icges(7,1) * t303 + Icges(7,4) * t302 + Icges(7,5) * t595;
t468 = -t176 * t380 + t178 * t382;
t93 = Icges(7,5) * t168 + Icges(7,6) * t167 + Icges(7,3) * t423;
t95 = Icges(7,4) * t168 + Icges(7,2) * t167 + Icges(7,6) * t423;
t97 = Icges(7,1) * t168 + Icges(7,4) * t167 + Icges(7,5) * t423;
t23 = (t394 * t468 - t93) * t374 + (t174 * t394 - t380 * t95 + t382 * t97 + (-t176 * t382 - t178 * t380) * qJD(6)) * t373;
t635 = t23 * t405;
t175 = Icges(7,5) * t305 + Icges(7,6) * t304 + Icges(7,3) * t594;
t177 = Icges(7,4) * t305 + Icges(7,2) * t304 + Icges(7,6) * t594;
t179 = Icges(7,1) * t305 + Icges(7,4) * t304 + Icges(7,5) * t594;
t467 = -t177 * t380 + t179 * t382;
t92 = Icges(7,5) * t166 + Icges(7,6) * t165 + Icges(7,3) * t422;
t94 = Icges(7,4) * t166 + Icges(7,2) * t165 + Icges(7,6) * t422;
t96 = Icges(7,1) * t166 + Icges(7,4) * t165 + Icges(7,5) * t422;
t24 = (t394 * t467 - t92) * t374 + (t175 * t394 - t380 * t94 + t382 * t96 + (-t177 * t382 - t179 * t380) * qJD(6)) * t373;
t634 = t24 * t403;
t389 = t403 * rSges(3,3);
t387 = t403 * rSges(5,3);
t633 = -rSges(6,3) - qJ(5);
t632 = -rSges(7,3) + t401;
t608 = t258 * t382;
t607 = t292 * t383;
t606 = t293 * t383;
t605 = t294 * t381;
t604 = t295 * t381;
t603 = t313 * t404;
t602 = t314 * t404;
t601 = t315 * t402;
t600 = t316 * t402;
t599 = t329 * t394;
t598 = t330 * t394;
t596 = t373 * t394;
t587 = t392 * t405;
t211 = t320 * rSges(6,1) + t319 * rSges(6,2) + rSges(6,3) * t594;
t577 = -t211 - t307;
t498 = rSges(6,1) * t399 - rSges(6,2) * t398;
t240 = (rSges(6,3) * t373 + t374 * t498) * t394;
t260 = qJ(5) * t596 + (pkin(4) * t394 - qJD(5)) * t374;
t576 = -t240 - t260;
t334 = t405 * t348;
t364 = t405 * t376;
t248 = -t403 * t562 + t334 - t364;
t289 = -pkin(1) * t405 + t403 * t642 + t364;
t575 = -t248 - t289;
t265 = -rSges(6,3) * t374 + t373 * t498;
t331 = pkin(4) * t373 - qJ(5) * t374;
t574 = -t265 - t331;
t501 = -rSges(5,2) * t373 + t639;
t280 = -rSges(5,3) * t405 + t403 * t501;
t281 = rSges(5,1) * t591 - rSges(5,2) * t594 + t387;
t187 = t403 * t280 + t405 * t281;
t288 = t400 * t405 + t391 - t680;
t573 = t403 * t288 + t405 * t289;
t306 = (t614 + t647) * t403;
t572 = t403 * t306 + t405 * t307;
t338 = t504 * qJD(2);
t323 = t405 * t338;
t385 = qJD(3) * t403;
t571 = t323 + t385;
t568 = rSges(5,2) * t522 + rSges(5,3) * t556;
t521 = t381 * t557;
t567 = rSges(4,2) * t521 + rSges(4,3) * t556;
t520 = t402 * t557;
t369 = pkin(2) * t520;
t566 = pkin(3) * t521 + t369;
t367 = t392 * t557;
t386 = qJD(3) * t405;
t565 = t367 + t386;
t553 = qJD(2) * t402;
t546 = pkin(2) * t553;
t564 = t400 * t557 + t403 * t546;
t563 = t405 * t641 + t389;
t561 = t396 + t397;
t560 = qJD(1) * t274;
t291 = Icges(4,3) * t403 + t405 * t475;
t559 = qJD(1) * t291;
t312 = Icges(3,3) * t403 + t405 * t476;
t558 = qJD(1) * t312;
t555 = qJD(2) * t381;
t554 = qJD(2) * t383;
t552 = qJD(2) * t403;
t551 = qJD(2) * t404;
t550 = qJD(5) * t373;
t548 = t405 * t638;
t544 = t202 * t595;
t542 = t203 * t595;
t256 = -Icges(7,3) * t374 + t373 * t472;
t104 = t256 * t594 + t304 * t257 + t305 * t258;
t77 = -t175 * t374 + t373 * t467;
t536 = t104 / 0.2e1 + t77 / 0.2e1;
t103 = t256 * t595 + t257 * t302 + t258 * t303;
t76 = -t174 * t374 + t373 * t468;
t535 = t76 / 0.2e1 + t103 / 0.2e1;
t155 = t484 * t593 + (Icges(7,5) * t394 + (-Icges(7,1) * t380 - t624) * qJD(6)) * t373;
t534 = t373 * t382 * t155 + t256 * t596 + t593 * t608;
t533 = -t260 + t675;
t351 = t405 * t550;
t424 = -t374 * t557 - t538;
t531 = t403 * (qJ(5) * t423 + qJD(1) * t358 + t403 * t550 - t343) + t405 * (pkin(4) * t424 - qJ(5) * t522 + t339 + t351) + t306 * t556;
t530 = t403 * (-t403 * t438 + (t405 * t501 + t387) * qJD(1)) + t405 * (rSges(5,1) * t424 - rSges(5,2) * t537 + t568) + t280 * t556;
t508 = t405 * t546;
t523 = -t386 - t564;
t529 = t403 * ((-t405 * t644 - t645) * qJD(1) + t523) + t405 * (-t508 + t385 + (t405 * t642 + t680) * qJD(1)) + t288 * t556;
t527 = -t331 + t671;
t526 = -t307 + t575;
t308 = t331 * t557;
t525 = t308 + t566;
t524 = t351 + t571;
t519 = t593 / 0.2e1;
t518 = t557 / 0.2e1;
t517 = t556 / 0.2e1;
t340 = rSges(4,1) * t381 + rSges(4,2) * t383;
t516 = -t340 - t649;
t209 = t574 * t405;
t190 = -qJD(1) * t275 - t405 * t599;
t515 = t278 * t394 + t190;
t191 = qJD(1) * t276 - t403 * t599;
t514 = t277 * t394 + t191;
t192 = -qJD(1) * t277 - t405 * t598;
t513 = -t276 * t394 + t192;
t193 = qJD(1) * t278 - t403 * t598;
t512 = t275 * t394 - t193;
t511 = -t403 * t392 + t334;
t510 = -t374 * t375 - t348;
t509 = -t338 - t550;
t505 = t403 * t247 + t405 * t248 + t573;
t102 = t403 * t210 + t405 * t211 + t572;
t135 = t527 * t405;
t503 = -t638 + t641;
t502 = -t637 + t640;
t64 = t174 * t595 + t176 * t302 + t178 * t303;
t65 = t175 * t595 + t177 * t302 + t179 * t303;
t45 = t65 * t403 - t405 * t64;
t494 = t403 * t64 + t405 * t65;
t66 = t174 * t594 + t304 * t176 + t305 * t178;
t67 = t175 * t594 + t304 * t177 + t305 * t179;
t46 = t67 * t403 - t405 * t66;
t493 = t403 * t66 + t405 * t67;
t492 = t77 * t403 - t405 * t76;
t491 = t403 * t76 + t405 * t77;
t489 = Icges(3,1) * t402 + t630;
t487 = Icges(4,1) * t381 + t628;
t485 = Icges(6,1) * t399 - Icges(6,4) * t398;
t482 = Icges(3,2) * t404 + t631;
t480 = Icges(4,2) * t383 + t629;
t478 = Icges(6,4) * t399 - Icges(6,2) * t398;
t328 = Icges(5,5) * t373 + Icges(5,6) * t374;
t473 = Icges(6,5) * t399 - Icges(6,6) * t398;
t421 = t373 * t632 + t510;
t114 = (-t392 + t646) * t405 + t421 * t403 + t496;
t115 = t446 + t511 + t184;
t471 = t114 * t405 + t115 * t403;
t118 = t183 * t374 + t259 * t595;
t119 = -t374 * t184 - t259 * t594;
t470 = t118 * t405 + t119 * t403;
t425 = t373 * t633 - t348 - t647;
t409 = t403 * t425 - t587;
t136 = t409 + t499;
t137 = t511 - t577;
t469 = t136 * t405 + t137 * t403;
t466 = t183 * t405 - t184 * t403;
t465 = -t204 * t398 + t206 * t399;
t464 = -t205 * t398 + t207 * t399;
t444 = -t348 - t501;
t222 = (rSges(5,3) - t392) * t405 + t444 * t403;
t223 = t281 + t511;
t461 = t222 * t405 + t223 * t403;
t300 = t405 * t640 + t670;
t449 = -pkin(1) - t503;
t448 = -t331 + t504;
t445 = -t376 - t502;
t442 = t403 * (t403 * t338 + t556 * t569 - t367 + t564) + t405 * (-qJD(1) * t247 + t323 + t508) + t247 * t556 + t529;
t61 = t403 * t197 + t405 * t198 + t572 - t674;
t439 = (-t390 - t648) * qJD(2);
t437 = qJD(2) * t362;
t436 = qJD(2) * t340;
t435 = -t265 + t448;
t432 = t394 * t328;
t431 = qJD(2) * t489;
t430 = qJD(2) * t487;
t429 = qJD(2) * t482;
t428 = qJD(2) * t480;
t427 = qJD(2) * (-Icges(3,5) * t402 - Icges(3,6) * t404);
t426 = qJD(2) * (-Icges(4,5) * t381 - Icges(4,6) * t383);
t124 = Icges(6,4) * t250 + Icges(6,2) * t249 + Icges(6,6) * t422;
t125 = Icges(6,4) * t252 + Icges(6,2) * t251 + Icges(6,6) * t423;
t126 = Icges(6,1) * t250 + Icges(6,4) * t249 + Icges(6,5) * t422;
t127 = Icges(6,1) * t252 + Icges(6,4) * t251 + Icges(6,5) * t423;
t130 = -t273 * t405 - t403 * t458;
t131 = -t274 * t405 - t681;
t188 = -t405 * t432 - t668;
t189 = -t403 * t432 + t560;
t78 = -t204 * t440 + t206 * t318 + t544;
t79 = -t205 * t440 + t207 * t318 + t542;
t17 = t165 * t176 + t166 * t178 + t174 * t422 + t304 * t95 + t305 * t97 + t594 * t93;
t18 = t165 * t177 + t166 * t179 + t175 * t422 + t304 * t94 + t305 * t96 + t594 * t92;
t8 = qJD(1) * t493 - t17 * t405 + t18 * t403;
t420 = t45 * t557 + t46 * t556 + ((-t130 - t78) * t557 + t689 * t556) * t405 + (t8 + (t319 * t124 + t320 * t126 + t403 * t188 + t249 * t205 + t250 * t207 + (-t542 + t681 - t689) * qJD(1)) * t403 + (t79 + t131) * t557 + t688 * t556 + (-t319 * t125 - t320 * t127 - t249 * t204 - t250 * t206 + t665 + (-t515 * t373 + t513 * t374 - t189) * t403 + (t191 * t373 - t193 * t374 + t275 * t593 + t277 * t596 - t668) * t405 + (t544 + (t274 - t458) * t403 + t688) * qJD(1)) * t405) * t403;
t419 = t448 + t671;
t301 = t501 * t394;
t417 = -t301 + t439;
t416 = t439 - t260;
t164 = t435 * t405;
t414 = -t240 + t416;
t413 = t442 + t531;
t117 = t419 * t405;
t28 = -t103 * t374 + t373 * t494;
t29 = -t104 * t374 + t373 * t493;
t33 = t153 * t594 + t304 * t154 + t305 * t155 + t165 * t257 + t166 * t258 + t256 * t422;
t3 = (t394 * t493 - t33) * t374 + (-qJD(1) * t46 + t104 * t394 + t17 * t403 + t18 * t405) * t373;
t19 = t167 * t176 + t168 * t178 + t174 * t423 + t302 * t95 + t303 * t97 + t595 * t93;
t20 = t167 * t177 + t168 * t179 + t175 * t423 + t302 * t94 + t303 * t96 + t595 * t92;
t34 = t153 * t595 + t302 * t154 + t303 * t155 + t167 * t257 + t168 * t258 + t256 * t423;
t4 = (t394 * t494 - t34) * t374 + (-qJD(1) * t45 + t394 * t103 + t19 * t403 + t20 * t405) * t373;
t9 = qJD(1) * t494 - t19 * t405 + t20 * t403;
t412 = t3 * t653 + t4 * t652 - t374 * (qJD(1) * t491 + t634 - t635) / 0.2e1 + t28 * t518 + t29 * t517 + t492 * t596 / 0.2e1 + t9 * t595 / 0.2e1 + t8 * t594 / 0.2e1 + (t405 * t519 - t522 / 0.2e1) * t46 + (t373 * t517 + t403 * t519) * t45;
t411 = t416 + t675;
t410 = rSges(3,2) * t520 + rSges(3,3) * t556 - t405 * t437;
t12 = (t440 * t125 - t318 * t127 - t251 * t204 - t252 * t206 + (t79 - t543) * qJD(1)) * t405 + (-t440 * t124 + t318 * t126 + t251 * t205 + t252 * t207 + (t78 + t541) * qJD(1) + t665) * t403;
t16 = (t405 * t189 + (t131 + t677) * qJD(1)) * t405 + (t130 * qJD(1) + (-t190 * t373 + t192 * t374 - t276 * t593 - t278 * t596 + t560) * t403 + (-t188 + t512 * t374 + t514 * t373 + (-t273 - t457) * qJD(1)) * t405) * t403;
t408 = (-t9 - t16 - t12) * t405 + t420;
t297 = t479 * t394;
t298 = t486 * t394;
t407 = qJD(1) * t328 + (t298 - t599) * t374 + (-t297 - t598) * t373;
t231 = (Icges(6,3) * t373 + t374 * t473) * t394;
t232 = (Icges(6,6) * t373 + t374 * t478) * t394;
t233 = (Icges(6,5) * t373 + t374 * t485) * t394;
t261 = -Icges(6,3) * t374 + t373 * t473;
t262 = -Icges(6,6) * t374 + t373 * t478;
t263 = -Icges(6,5) * t374 + t373 * t485;
t406 = -t635 / 0.2e1 + t634 / 0.2e1 + (t231 * t594 + t319 * t232 + t320 * t233 + t249 * t262 + t250 * t263 + t261 * t422 + t403 * t662 + t407 * t405 + t33 + (t394 * t464 - t122 + t515) * t374 + (-t124 * t398 + t126 * t399 + t203 * t394 + t513) * t373) * t653 + (t231 * t595 - t440 * t232 + t318 * t233 + t251 * t262 + t252 * t263 + t261 * t423 + t407 * t403 - t405 * t662 + t34 + (t394 * t465 - t123 + t514) * t374 + (-t125 * t398 + t127 * t399 + t202 * t394 - t512) * t373) * t652 + (t261 * t595 - t262 * t440 + t263 * t318 - t328 * t405 - t403 * t452 + t103 + t76 + (-t202 + t275) * t374 + (t277 + t465) * t373) * t518 + (t261 * t594 + t319 * t262 + t320 * t263 + t403 * t328 - t405 * t452 + t104 + t77 + (-t203 + t276) * t374 + (t278 + t464) * t373) * t517;
t349 = t503 * qJD(2);
t327 = t502 * qJD(2);
t322 = -t548 + t563;
t321 = t403 * t503 - t636;
t299 = -rSges(4,3) * t405 + t403 * t502;
t287 = t516 * t405;
t286 = t516 * t403;
t272 = t645 + (pkin(1) - t638) * t405 + t563;
t271 = t403 * t449 + t391 + t636;
t255 = t265 * t557;
t245 = -t403 * t400 + t300 + t364;
t244 = (rSges(4,3) - t400) * t405 + t445 * t403;
t235 = t403 * t427 + t558;
t234 = -qJD(1) * t311 + t405 * t427;
t216 = t403 * t426 + t559;
t215 = -qJD(1) * t290 + t405 * t426;
t214 = t362 * t552 + ((-rSges(3,3) - pkin(7)) * t403 + t449 * t405) * qJD(1);
t213 = (t391 + (-pkin(1) - t641) * t403) * qJD(1) + t410;
t208 = t574 * t403;
t200 = -t340 * t556 - t403 * t327 + (-t402 * t556 - t403 * t551) * pkin(2);
t199 = t340 * t557 + t369 + (-pkin(2) * t551 - t327) * t405;
t163 = t435 * t403;
t152 = t403 * t312 - t405 * t453;
t151 = t403 * t311 - t679;
t150 = -t312 * t405 - t683;
t149 = -t311 * t405 - t403 * t454;
t146 = t340 * t552 + (t405 * t445 - t388) * qJD(1) - t523;
t145 = t385 + (-t376 - t640) * t557 + (-qJD(1) * t400 + qJD(2) * t516) * t405 + t567;
t143 = qJD(1) * t254 + t403 * t417;
t142 = t332 * t557 + t405 * t417 + t566;
t141 = t403 * t291 - t405 * t455;
t140 = t403 * t290 - t678;
t139 = -t291 * t405 - t682;
t138 = -t290 * t405 - t403 * t456;
t134 = t527 * t403;
t129 = (-t338 + t438) * t403 + (t405 * t444 - t387) * qJD(1) + t565;
t128 = -t405 * t438 + (-t587 + (-t348 - t639) * t403) * qJD(1) + t568 + t571;
t116 = t419 * t403;
t111 = -t256 * t374 + (t608 - t611) * t373;
t110 = qJD(1) * t209 + t403 * t576;
t109 = t405 * t576 + t255 + t308;
t108 = t466 * t373;
t107 = t111 * t596;
t89 = qJD(1) * t164 + t403 * t414;
t88 = t405 * t414 + t255 + t525;
t87 = -t281 * t557 + t530;
t84 = t505 + t187;
t75 = t343 + (t593 * t633 + t509) * t403 + t425 * t556 - t500 + t565;
t74 = -pkin(4) * t538 + qJD(1) * t409 + t339 + t524 + t528;
t63 = qJD(1) * t135 + t403 * t533;
t62 = t405 * t533 + t308 + t672;
t60 = t421 * t556 + (-t545 + (t374 * t632 + t597) * t394 + t509) * t403 - t497 + t565;
t59 = (-t592 - t597) * t585 + (-t587 + (-rSges(7,3) * t373 + t510) * t403) * qJD(1) + t524 + t532 + t570;
t58 = t102 + t505;
t57 = qJD(1) * t117 + t403 * t411;
t56 = t405 * t411 + t525 + t672;
t53 = (t259 * t586 + t99) * t374 + (t403 * t160 - t183 * t394 + t259 * t556) * t373;
t52 = (-t259 * t585 - t98) * t374 + (-t160 * t405 + t184 * t394 + t246) * t373;
t49 = t61 + t505;
t47 = t373 * t687 + t676 * t374 + t534;
t38 = t557 * t577 + t531 + t667;
t35 = (-t281 + t575) * t557 + t442 + t530;
t30 = t466 * t593 + (qJD(1) * t674 - t403 * t98 + t405 * t99) * t373;
t25 = (-t211 + t526) * t557 + t413 + t667;
t14 = (-t307 + t673) * t557 + t531 + t661;
t13 = t413 + (t526 + t673) * t557 + t661;
t1 = [(t136 * t75 + t137 * t74) * t657 + (t128 * t223 + t129 * t222) * t658 + (t145 * t245 + t146 * t244) * t659 + (t213 * t272 + t214 * t271) * t660 + (t114 * t60 + t115 * t59) * t656 + t534 + (t261 - t329) * t596 + (-t480 + t488) * t555 + (t481 + t487) * t554 + (-t482 + t490) * t553 + (t483 + t489) * t551 + (-t262 * t398 + t263 * t399 + t330) * t593 + (t297 - t231 + t676) * t374 + (-t232 * t398 + t233 * t399 + t298 + t687) * t373; t406 + m(3) * ((-t213 * t403 - t214 * t405) * t362 + (-t271 * t405 - t272 * t403) * t349) + m(4) * (t145 * t286 + t146 * t287 + t199 * t244 + t200 * t245) + m(5) * (t128 * t253 + t129 * t254 + t142 * t222 + t143 * t223) + m(6) * (t136 * t88 + t137 * t89 + t163 * t74 + t164 * t75) + m(7) * (t114 * t56 + t115 * t57 + t116 * t59 + t117 * t60) + ((t606 / 0.2e1 + t604 / 0.2e1 + t602 / 0.2e1 + t600 / 0.2e1 - t272 * t651) * t405 + (t271 * t651 + t603 / 0.2e1 + t601 / 0.2e1 + t607 / 0.2e1 + t605 / 0.2e1) * t403) * qJD(1) + ((-qJD(1) * t292 - t405 * t428) * t383 + (-qJD(1) * t294 - t405 * t430) * t381 + (-qJD(1) * t313 - t405 * t429) * t404 + (-qJD(1) * t315 - t405 * t431) * t402 + (-t453 - t455) * qJD(2)) * t653 + ((qJD(1) * t293 - t403 * t428) * t383 + (qJD(1) * t295 - t403 * t430) * t381 + (qJD(1) * t314 - t403 * t429) * t404 + (qJD(1) * t316 - t403 * t431) * t402 + (-t454 - t456) * qJD(2)) * t652 + (t476 + t475) * qJD(2) * (t397 / 0.2e1 + t396 / 0.2e1); t420 + ((t403 * t321 + t322 * t405) * ((qJD(1) * t321 + t410) * t405 + (-t403 * t437 + (-t322 - t548 + t389) * qJD(1)) * t403) + t561 * t362 * t349) * t660 + t403 * ((t403 * t215 + (t140 + t682) * qJD(1)) * t403 + (t141 * qJD(1) + (t292 * t554 + t294 * t555) * t405 + (-t216 + (-t604 - t606) * qJD(2) + (t291 - t456) * qJD(1)) * t403) * t405) + t403 * ((t403 * t234 + (t151 + t683) * qJD(1)) * t403 + (t152 * qJD(1) + (t313 * t551 + t315 * t553) * t405 + (-t235 + (-t600 - t602) * qJD(2) + (t312 - t454) * qJD(1)) * t403) * t405) - t405 * ((t405 * t216 + (t139 + t678) * qJD(1)) * t405 + (t138 * qJD(1) + (-t293 * t554 - t295 * t555 + t559) * t403 + (-t215 + (t605 + t607) * qJD(2) - t455 * qJD(1)) * t405) * t403) - t405 * ((t405 * t235 + (t150 + t679) * qJD(1)) * t405 + (t149 * qJD(1) + (-t314 * t551 - t316 * t553 + t558) * t403 + (-t234 + (t601 + t603) * qJD(2) - t453 * qJD(1)) * t405) * t403) - t405 * t9 - t405 * t16 - t405 * t12 + (t142 * t254 + t143 * t253 + t35 * t84) * t658 + (t163 * t89 + t164 * t88 + t25 * t58) * t657 + (t116 * t57 + t117 * t56 + t13 * t49) * t656 + (t287 * t199 + t286 * t200 + (t403 * t299 + t300 * t405 + t573) * ((qJD(1) * t299 - t405 * t436 + t567) * t405 + (-t403 * t436 + (-t289 - t300 + t670) * qJD(1)) * t403 + t529)) * t659 + ((-t138 - t149) * t405 + (t139 + t150) * t403) * t557 + ((-t140 - t151) * t405 + (t141 + t152) * t403) * t556; m(7) * (qJD(1) * t471 + t403 * t60 - t405 * t59) + m(6) * (qJD(1) * t469 + t403 * t75 - t405 * t74) + m(5) * (qJD(1) * t461 - t128 * t405 + t403 * t129) + m(4) * (-t145 * t405 + t403 * t146 + (t244 * t405 + t245 * t403) * qJD(1)); m(7) * (t403 * t56 - t405 * t57 + (t116 * t403 + t117 * t405) * qJD(1)) + m(6) * (t403 * t88 - t405 * t89 + (t163 * t403 + t164 * t405) * qJD(1)) + m(5) * (qJD(1) * t669 + t403 * t142 - t143 * t405) + m(4) * (t403 * t199 - t200 * t405 + (t286 * t403 + t287 * t405) * qJD(1)); 0; t406 + m(6) * (t109 * t136 + t110 * t137 + t208 * t74 + t209 * t75) + m(7) * (t114 * t62 + t115 * t63 + t134 * t59 + t135 * t60) + (-t128 * t403 - t129 * t405 + (t222 * t403 - t223 * t405) * qJD(1)) * t650 - m(5) * t461 * t301; t408 + m(5) * (t187 * t35 - t301 * t669 + t87 * t84) + m(7) * (t116 * t63 + t117 * t62 + t13 * t61 + t134 * t57 + t135 * t56 + t14 * t49) + m(6) * (t102 * t25 + t109 * t164 + t110 * t163 + t208 * t89 + t209 * t88 + t38 * t58) + (-t142 * t405 - t143 * t403 + (-t253 * t405 + t254 * t403) * qJD(1)) * t650; m(6) * (t109 * t403 - t110 * t405 + (t208 * t403 + t209 * t405) * qJD(1)) + m(7) * (t62 * t403 - t405 * t63 + (t134 * t403 + t135 * t405) * qJD(1)); (t134 * t63 + t135 * t62 + t14 * t61) * t656 + (t102 * t38 + t109 * t209 + t110 * t208) * t657 + (t301 * t332 * t561 + t187 * t87) * t658 + t408; 0.2e1 * (t469 * t655 + t471 * t654) * t593 + 0.2e1 * ((-t114 * t557 + t115 * t556 + t403 * t59 + t405 * t60) * t654 + (-t136 * t557 + t137 * t556 + t403 * t74 + t405 * t75) * t655) * t373; 0.2e1 * ((t116 * t586 + t117 * t585 - t13) * t654 + (t163 * t586 + t164 * t585 - t25) * t655) * t374 + 0.2e1 * ((t116 * t556 - t117 * t557 + t394 * t49 + t403 * t57 + t405 * t56) * t654 + (t163 * t556 - t164 * t557 + t394 * t58 + t403 * t89 + t405 * t88) * t655) * t373; 0; 0.2e1 * ((t134 * t586 + t135 * t585 - t14) * t654 + (t208 * t586 + t209 * t585 - t38) * t655) * t374 + 0.2e1 * ((t134 * t556 - t135 * t557 + t394 * t61 + t403 * t63 + t405 * t62) * t654 + (t102 * t394 + t109 * t405 + t110 * t403 + t208 * t556 - t209 * t557) * t655) * t373; 0.4e1 * (t655 + t654) * (-0.1e1 + t561) * t373 * t593; t107 + m(7) * (t114 * t53 + t115 * t52 + t118 * t60 + t119 * t59) + (-t47 + (t403 * t535 + t405 * t536) * t394) * t374 + ((t33 / 0.2e1 + t24 / 0.2e1) * t405 + (t23 / 0.2e1 + t34 / 0.2e1) * t403 + (-t403 * t536 + t405 * t535) * qJD(1)) * t373; t412 + m(7) * (t108 * t13 + t116 * t52 + t117 * t53 + t118 * t56 + t119 * t57 + t30 * t49); m(7) * (qJD(1) * t470 + t53 * t403 - t405 * t52); t412 + m(7) * (t108 * t14 + t118 * t62 + t119 * t63 + t134 * t52 + t135 * t53 + t30 * t61); m(7) * ((t394 * t470 - t30) * t374 + (t108 * t394 + t403 * t52 + t405 * t53 + (-t118 * t403 + t119 * t405) * qJD(1)) * t373); (t108 * t30 + t118 * t53 + t119 * t52) * t656 + (t47 * t374 - t107 + (t403 * t28 + t405 * t29 - t374 * t491) * t394) * t374 + (t405 * t3 + t403 * t4 + t491 * t596 + (-t111 * t394 - t23 * t403 - t24 * t405) * t374 + (t405 * t28 - t403 * t29 + t374 * t492) * qJD(1)) * t373;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
