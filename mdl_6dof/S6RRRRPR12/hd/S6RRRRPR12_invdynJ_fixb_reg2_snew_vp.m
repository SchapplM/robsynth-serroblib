% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 00:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRRPR12_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_invdynJ_fixb_reg2_snew_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 00:31:40
% EndTime: 2019-05-08 00:33:39
% DurationCPUTime: 58.67s
% Computational Cost: add. (625921->965), mult. (1551951->1533), div. (0->0), fcn. (1327256->16), ass. (0->626)
t440 = sin(pkin(6));
t447 = sin(qJ(2));
t451 = cos(qJ(2));
t439 = sin(pkin(7));
t663 = pkin(10) * t439;
t670 = pkin(2) * t451;
t513 = -t447 * t663 - t670;
t494 = t513 * t440;
t442 = cos(pkin(7));
t443 = cos(pkin(6));
t603 = qJD(1) * t443;
t433 = qJD(2) + t603;
t596 = -qJD(2) + t433;
t559 = t442 * t596;
t552 = pkin(10) * t559;
t572 = pkin(10) * t442 + pkin(9);
t589 = qJDD(1) * t447;
t617 = t443 * t451;
t659 = t451 * g(3);
t701 = (-t572 * t589 - t659 + (t451 * t552 + qJD(1) * (pkin(9) * t617 - t447 * t494)) * qJD(1)) * t440;
t700 = t439 * t701;
t595 = qJD(2) + t433;
t438 = sin(pkin(13));
t446 = sin(qJ(3));
t450 = cos(qJ(3));
t604 = qJD(1) * t440;
t618 = t442 * t451;
t624 = t439 * t446;
t410 = t433 * t624 + (t446 * t618 + t447 * t450) * t604;
t601 = qJD(1) * t451;
t570 = t440 * t601;
t418 = -t442 * t433 + t439 * t570 - qJD(3);
t445 = sin(qJ(4));
t449 = cos(qJ(4));
t393 = -t445 * t410 - t418 * t449;
t394 = t410 * t449 - t418 * t445;
t441 = cos(pkin(13));
t361 = -t441 * t393 + t394 * t438;
t363 = t438 * t393 + t441 * t394;
t319 = t363 * t361;
t588 = qJDD(1) * t451;
t590 = qJD(1) * qJD(2);
t502 = -t447 * t590 + t588;
t482 = t502 * t442;
t566 = t451 * t590;
t503 = t566 + t589;
t484 = t503 * t440;
t555 = qJDD(1) * t443 + qJDD(2);
t512 = t439 * t555;
t471 = t446 * t484 + (-t440 * t482 - t512) * t450;
t382 = -t410 * qJD(3) - t471;
t381 = qJDD(4) - t382;
t684 = -t319 + t381;
t698 = t438 * t684;
t697 = t441 * t684;
t444 = sin(qJ(6));
t556 = t442 * t570;
t602 = qJD(1) * t447;
t571 = t440 * t602;
t623 = t439 * t450;
t408 = -t433 * t623 + t446 * t571 - t450 * t556;
t463 = (t446 * t482 + t450 * t503) * t440 + t446 * t512;
t383 = -t408 * qJD(3) + t463;
t483 = t502 * t440;
t587 = t442 * t555 + qJDD(3);
t470 = t439 * t483 - t587;
t332 = t393 * qJD(4) + t449 * t383 - t445 * t470;
t560 = t383 * t445 + t449 * t470;
t507 = -qJD(4) * t394 - t560;
t562 = t332 * t438 - t441 * t507;
t277 = qJDD(6) + t562;
t404 = qJD(4) + t408;
t448 = cos(qJ(6));
t337 = t363 * t444 - t448 * t404;
t339 = t363 * t448 + t404 * t444;
t291 = t339 * t337;
t685 = t277 - t291;
t696 = t444 * t685;
t365 = t393 * t394;
t683 = t365 + t381;
t695 = t445 * t683;
t694 = t448 * t685;
t693 = t449 * t683;
t671 = sin(qJ(1));
t672 = cos(qJ(1));
t500 = g(1) * t671 - g(2) * t672;
t506 = -pkin(1) + t513;
t558 = t451 * t595;
t661 = t443 * g(3);
t436 = t447 ^ 2;
t437 = t451 ^ 2;
t682 = -t436 - t437;
t459 = -t661 + (t506 * qJDD(1) + (-pkin(9) * t604 + (t442 * t604 * t682 - t439 * t558) * pkin(10) + t595 * t447 * pkin(2)) * qJD(1) - t500) * t440;
t432 = t433 ^ 2;
t478 = qJDD(1) * pkin(1) + t500;
t472 = t443 * t478;
t501 = g(1) * t672 + g(2) * t671;
t675 = qJD(1) ^ 2;
t479 = pkin(1) * t675 + t501;
t515 = t555 * pkin(2);
t462 = t432 * t663 + t447 * t479 + t451 * t472 + t515;
t692 = t442 * (t701 + t462) + t439 * t459;
t397 = t408 * t418;
t350 = t397 + t383;
t461 = t439 * t462 - t442 * t459;
t460 = -pkin(11) * t350 - t461;
t625 = t418 * t410;
t561 = -t382 - t625;
t691 = t561 * pkin(3) + t460;
t635 = t363 * t404;
t256 = t562 + t635;
t626 = t410 * t408;
t465 = -t470 - t626;
t687 = t446 * t465;
t686 = t450 * t465;
t279 = t441 * t332 + t438 * t507;
t342 = t404 * t361;
t259 = -t342 + t279;
t375 = t404 * t393;
t315 = -t375 + t332;
t314 = t375 + t332;
t357 = qJD(6) + t361;
t563 = t279 * t444 - t448 * t381;
t218 = (qJD(6) - t357) * t339 + t563;
t311 = (qJD(4) - t404) * t394 + t560;
t435 = t440 ^ 2;
t680 = (t603 - t595) * t435;
t335 = t337 ^ 2;
t336 = t339 ^ 2;
t356 = t357 ^ 2;
t358 = t361 ^ 2;
t359 = t363 ^ 2;
t678 = t393 ^ 2;
t392 = t394 ^ 2;
t677 = t404 ^ 2;
t405 = t408 ^ 2;
t406 = t410 ^ 2;
t676 = t418 ^ 2;
t674 = 2 * qJD(5);
t385 = pkin(3) * t408 - pkin(11) * t410;
t622 = t440 * t447;
t434 = g(3) * t622;
t664 = pkin(9) * t443;
t360 = -t451 * t479 + t447 * t472 - t434 + pkin(10) * t512 - pkin(2) * t432 + (t572 * t588 + ((t447 * t664 + t451 * t494) * qJD(1) + t447 * t552) * qJD(1)) * t440;
t607 = t450 * t360;
t520 = -t408 * t385 + t607;
t662 = t676 * pkin(3);
t453 = -pkin(11) * t470 + t446 * t692 + t520 - t662;
t455 = t691 - t700;
t196 = t445 * t455 + t449 * t453;
t369 = pkin(4) * t404 - qJ(5) * t394;
t163 = -pkin(4) * t678 + qJ(5) * t507 - t404 * t369 + t196;
t452 = pkin(4) * t683 - qJ(5) * t315 - t445 * t453 + t449 * t455;
t565 = t163 * t438 - t441 * t452;
t107 = t363 * t674 + t565;
t599 = qJD(5) * t361;
t354 = -0.2e1 * t599;
t605 = t441 * t163 + t438 * t452;
t108 = t354 + t605;
t66 = -t107 * t441 + t108 * t438;
t673 = pkin(4) * t66;
t669 = pkin(3) * t446;
t668 = pkin(3) * t450;
t257 = t562 - t635;
t260 = t342 + t279;
t191 = -t257 * t438 - t260 * t441;
t667 = pkin(4) * t191;
t666 = pkin(5) * t438;
t665 = pkin(9) * t440;
t317 = pkin(5) * t361 - pkin(12) * t363;
t103 = -t381 * pkin(5) - t677 * pkin(12) + (t674 + t317) * t363 + t565;
t104 = -pkin(5) * t677 + pkin(12) * t381 - t317 * t361 + t108;
t275 = t360 * t446 - t450 * t692;
t262 = t470 * pkin(3) - t676 * pkin(11) + t385 * t410 + t275;
t213 = -t507 * pkin(4) - t678 * qJ(5) + t369 * t394 + qJDD(5) + t262;
t143 = pkin(5) * t256 - pkin(12) * t259 + t213;
t71 = t104 * t444 - t448 * t143;
t72 = t104 * t448 + t143 * t444;
t51 = t444 * t71 + t448 * t72;
t35 = -t103 * t441 + t438 * t51;
t36 = t103 * t438 + t441 * t51;
t21 = t35 * t449 + t36 * t445;
t22 = -t35 * t445 + t36 * t449;
t50 = t444 * t72 - t448 * t71;
t548 = t22 * t446 - t450 * t50;
t9 = -t21 * t439 + t442 * t548;
t660 = t447 * t9;
t658 = t445 * t66;
t652 = t449 * t66;
t67 = t107 * t438 + t441 * t108;
t44 = t445 * t67 + t652;
t45 = t449 * t67 - t658;
t543 = -t213 * t450 + t446 * t45;
t25 = -t439 * t44 + t442 * t543;
t657 = t447 * t25;
t527 = -t279 * t448 - t381 * t444;
t239 = -qJD(6) * t337 - t527;
t304 = t357 * t337;
t222 = t239 + t304;
t151 = -t218 * t444 - t222 * t448;
t153 = -t218 * t448 + t222 * t444;
t268 = t335 + t336;
t126 = t153 * t438 + t268 * t441;
t127 = t153 * t441 - t268 * t438;
t87 = -t126 * t445 + t127 * t449;
t546 = -t151 * t450 + t446 * t87;
t86 = t126 * t449 + t127 * t445;
t49 = -t439 * t86 + t442 * t546;
t656 = t447 * t49;
t270 = -t356 - t335;
t177 = t270 * t444 + t694;
t178 = t270 * t448 - t696;
t591 = qJD(6) + t357;
t219 = -t339 * t591 - t563;
t137 = t178 * t438 + t219 * t441;
t138 = t178 * t441 - t219 * t438;
t96 = -t137 * t445 + t138 * t449;
t545 = -t177 * t450 + t446 * t96;
t95 = t137 * t449 + t138 * t445;
t54 = -t439 * t95 + t442 * t545;
t655 = t447 * t54;
t271 = -t336 - t356;
t231 = t277 + t291;
t647 = t231 * t444;
t182 = t271 * t448 - t647;
t646 = t231 * t448;
t183 = -t271 * t444 - t646;
t223 = t337 * t591 + t527;
t139 = t183 * t438 + t223 * t441;
t140 = t183 * t441 - t223 * t438;
t98 = -t139 * t445 + t140 * t449;
t544 = -t182 * t450 + t446 * t98;
t97 = t139 * t449 + t140 * t445;
t56 = -t439 * t97 + t442 * t544;
t654 = t447 * t56;
t193 = -t257 * t441 + t260 * t438;
t131 = t191 * t449 + t193 * t445;
t133 = -t191 * t445 + t193 * t449;
t274 = -t358 - t359;
t537 = t133 * t446 - t274 * t450;
t89 = -t131 * t439 + t442 * t537;
t653 = t447 * t89;
t651 = t103 * t444;
t650 = t103 * t448;
t649 = t213 * t438;
t648 = t213 * t441;
t645 = t262 * t445;
t644 = t262 * t449;
t286 = t319 + t381;
t643 = t286 * t438;
t642 = t286 * t441;
t320 = t461 + t700;
t641 = t320 * t446;
t640 = t320 * t450;
t324 = -t365 + t381;
t639 = t324 * t445;
t638 = t324 * t449;
t637 = t357 * t444;
t636 = t357 * t448;
t370 = t470 - t626;
t634 = t370 * t446;
t633 = t370 * t450;
t632 = t381 * t446;
t631 = t404 * t438;
t630 = t404 * t441;
t629 = t404 * t445;
t628 = t404 * t449;
t627 = t408 * t450;
t621 = t440 * t451;
t620 = t442 * t446;
t619 = t442 * t450;
t616 = t445 * t446;
t305 = -t677 - t358;
t240 = t305 * t438 + t697;
t241 = t305 * t441 - t698;
t173 = t240 * t449 + t241 * t445;
t174 = -t240 * t445 + t241 * t449;
t534 = t174 * t446 - t256 * t450;
t110 = -t173 * t439 + t442 * t534;
t615 = t447 * t110;
t326 = -t359 - t677;
t249 = t326 * t441 - t643;
t250 = -t326 * t438 - t642;
t184 = t249 * t449 + t250 * t445;
t185 = -t249 * t445 + t250 * t449;
t533 = t185 * t446 - t259 * t450;
t113 = -t184 * t439 + t442 * t533;
t614 = t447 * t113;
t246 = -t311 * t445 - t315 * t449;
t248 = -t311 * t449 + t315 * t445;
t333 = t392 + t678;
t529 = t248 * t446 + t333 * t450;
t168 = -t246 * t439 + t442 * t529;
t613 = t447 * t168;
t344 = -t677 - t678;
t280 = t344 * t445 + t693;
t281 = t344 * t449 - t695;
t312 = (-qJD(4) - t404) * t394 - t560;
t526 = t281 * t446 + t312 * t450;
t198 = -t280 * t439 + t442 * t526;
t612 = t447 * t198;
t353 = -t392 - t677;
t283 = t353 * t449 - t639;
t284 = -t353 * t445 - t638;
t525 = t284 * t446 - t314 * t450;
t210 = -t283 * t439 + t442 * t525;
t611 = t447 * t210;
t368 = -t405 - t406;
t351 = -t397 + t383;
t464 = (-qJD(3) - t418) * t410 - t471;
t521 = -t351 * t450 + t446 * t464;
t273 = -t368 * t439 + t442 * t521;
t610 = t447 * t273;
t597 = t440 * t675;
t468 = pkin(9) * t597 + t478;
t466 = t443 * t468;
t469 = qJDD(1) * t665 - t479;
t398 = t447 * t469 + (t440 * g(3) - t466) * t451;
t609 = t447 * t398;
t598 = t435 * t675;
t431 = t451 * t447 * t598;
t420 = t431 + t555;
t608 = t447 * t420;
t421 = -t431 + t555;
t606 = t451 * t421;
t594 = qJD(3) - t418;
t585 = pkin(10) * t589;
t584 = pkin(2) * t588;
t583 = t438 * t291;
t582 = t441 * t291;
t581 = t446 * t319;
t580 = t450 * t319;
t579 = t446 * t365;
t578 = t450 * t365;
t577 = t418 * t624;
t576 = t418 * t623;
t575 = t418 * t620;
t574 = t418 * t619;
t573 = -pkin(5) * t441 - pkin(4);
t569 = qJD(1) * t616;
t568 = t436 * t598;
t567 = t437 * t598;
t473 = -pkin(9) * t589 - t442 * t585 - t659;
t454 = (t445 * (t473 * t620 + (-t439 * t585 - t478 - t584) * t624 - t439 * pkin(11) * t588) + t449 * t439 * t473 + (((-t440 * t616 + t449 * t617) * qJD(1) * pkin(9) + (t445 * qJD(2) * pkin(11) + (t449 * t570 + t595 * t616) * pkin(2)) * t447 + (t436 * t449 * t604 - t558 * t616) * t663) * t439 + ((pkin(2) * t622 + t664) * t569 + (t559 * t616 + (t449 * t596 - t569 * t621) * t439) * pkin(10)) * t618) * qJD(1)) * t440;
t457 = pkin(11) * t587 + t462 * t620 - t624 * t661 + t520;
t195 = t445 * (t457 - t662) - t449 * t691 + t454;
t135 = t195 * t445 + t449 * t196;
t554 = t440 * t566;
t553 = pkin(4) * t249 - t605;
t551 = -pkin(11) * t446 - t668;
t380 = -t676 - t405;
t329 = t380 * t450 - t687;
t550 = pkin(10) * t329 + t640;
t384 = -t406 - t676;
t334 = -t384 * t446 + t633;
t549 = pkin(10) * t334 - t641;
t221 = t239 - t304;
t150 = t219 * t444 + t221 * t448;
t152 = t219 * t448 - t221 * t444;
t290 = t336 - t335;
t128 = t152 * t438 - t290 * t441;
t129 = t152 * t441 + t290 * t438;
t94 = -t128 * t445 + t129 * t449;
t547 = -t150 * t450 + t446 * t94;
t303 = -t336 + t356;
t207 = -t303 * t444 + t694;
t146 = t207 * t438 - t222 * t441;
t148 = t207 * t441 + t222 * t438;
t101 = -t146 * t445 + t148 * t449;
t205 = t303 * t448 + t696;
t542 = t101 * t446 - t205 * t450;
t302 = t335 - t356;
t208 = t302 * t448 - t647;
t147 = t208 * t438 + t218 * t441;
t149 = t208 * t441 - t218 * t438;
t102 = -t147 * t445 + t149 * t449;
t206 = t302 * t444 + t646;
t541 = t102 * t446 - t206 * t450;
t238 = -qJD(6) * t339 - t563;
t215 = -t238 * t444 + t337 * t636;
t169 = t215 * t438 + t582;
t171 = t215 * t441 - t583;
t121 = -t169 * t445 + t171 * t449;
t214 = -t238 * t448 - t337 * t637;
t540 = t121 * t446 + t214 * t450;
t217 = t239 * t448 - t339 * t637;
t170 = t217 * t438 - t582;
t172 = t217 * t441 + t583;
t122 = -t170 * t445 + t172 * t449;
t216 = t239 * t444 + t339 * t636;
t539 = t122 * t446 - t216 * t450;
t190 = -t256 * t438 + t259 * t441;
t192 = -t256 * t441 - t259 * t438;
t132 = -t190 * t445 + t192 * t449;
t318 = t359 - t358;
t538 = t132 * t446 - t318 * t450;
t536 = t135 * t446 - t262 * t450;
t244 = (-t337 * t448 + t339 * t444) * t357;
t211 = t244 * t438 - t277 * t441;
t212 = t244 * t441 + t277 * t438;
t145 = -t211 * t445 + t212 * t449;
t243 = (-t337 * t444 - t339 * t448) * t357;
t535 = t145 * t446 - t243 * t450;
t134 = -t195 * t449 + t196 * t445;
t341 = -t359 + t677;
t264 = t341 * t441 + t698;
t266 = -t341 * t438 + t697;
t201 = -t264 * t445 + t266 * t449;
t532 = t201 * t446 - t260 * t450;
t340 = t358 - t677;
t265 = t340 * t438 + t642;
t267 = t340 * t441 - t643;
t202 = -t265 * t445 + t267 * t449;
t531 = t202 * t446 + t257 * t450;
t247 = t312 * t449 - t314 * t445;
t364 = t392 - t678;
t530 = t247 * t446 - t364 * t450;
t415 = t440 * t468 + t661;
t491 = t433 * t439 + t556;
t276 = t607 + (t442 * (-g(3) * t621 + t451 * t466 + t515) + t439 * (-t440 * t584 - t415) + (t442 * (t433 * t491 - t442 * t554) + t439 * (-t439 * t554 - t491 * t570)) * pkin(10) + (t442 * t479 + ((-pkin(10) * t439 ^ 2 - t442 * t572) * qJDD(1) + (t439 * t595 + t556) * qJD(1) * pkin(2)) * t440) * t447) * t446;
t528 = t275 * t450 - t276 * t446;
t228 = t275 * t446 + t276 * t450;
t374 = -t392 + t677;
t298 = -t374 * t445 + t693;
t524 = t298 * t446 - t315 * t450;
t373 = -t677 + t678;
t299 = t373 * t449 - t639;
t523 = t299 * t446 + t311 * t450;
t346 = t410 * t594 + t471;
t522 = -t346 * t446 + t350 * t450;
t519 = t384 * t450 + t634;
t395 = t405 - t676;
t518 = t395 * t446 - t633;
t396 = -t406 + t676;
t517 = t396 * t450 + t687;
t516 = t380 * t446 + t686;
t251 = t361 * t631 - t441 * t562;
t252 = t361 * t630 + t438 * t562;
t188 = -t251 * t445 + t252 * t449;
t511 = t188 * t446 + t580;
t253 = t279 * t438 + t363 * t630;
t254 = t279 * t441 - t363 * t631;
t189 = -t253 * t445 + t254 * t449;
t510 = t189 * t446 - t580;
t308 = -t393 * t628 - t445 * t507;
t509 = t308 * t446 - t578;
t310 = t332 * t449 - t394 * t629;
t508 = t310 * t446 + t578;
t504 = pkin(4) * t35 - pkin(5) * t103 + pkin(12) * t51;
t10 = -pkin(3) * t21 - t504;
t15 = t22 * t450 + t446 * t50;
t11 = qJ(5) * t36 + (-pkin(12) * t438 + t573) * t50;
t16 = -qJ(5) * t35 + (-pkin(12) * t441 + t666) * t50;
t4 = -pkin(11) * t21 - t11 * t445 + t16 * t449;
t499 = pkin(10) * t15 + t10 * t450 + t4 * t446;
t43 = -pkin(12) * t151 - t50;
t29 = qJ(5) * t127 + t151 * t573 + t43 * t438;
t33 = -qJ(5) * t126 + t151 * t666 + t43 * t441;
t14 = -pkin(11) * t86 - t29 * t445 + t33 * t449;
t475 = pkin(4) * t126 + pkin(5) * t268 + pkin(12) * t153 + t51;
t27 = -pkin(3) * t86 - t475;
t65 = t151 * t446 + t450 * t87;
t498 = pkin(10) * t65 + t14 * t446 + t27 * t450;
t63 = -pkin(5) * t177 + t71;
t83 = -pkin(12) * t177 + t651;
t37 = -pkin(4) * t177 + qJ(5) * t138 + t438 * t83 + t441 * t63;
t40 = -qJ(5) * t137 - t438 * t63 + t441 * t83;
t19 = -pkin(11) * t95 - t37 * t445 + t40 * t449;
t477 = pkin(4) * t137 + pkin(5) * t219 + pkin(12) * t178 - t650;
t46 = -pkin(3) * t95 - t477;
t74 = t177 * t446 + t450 * t96;
t497 = pkin(10) * t74 + t19 * t446 + t450 * t46;
t64 = -pkin(5) * t182 + t72;
t84 = -pkin(12) * t182 + t650;
t38 = -pkin(4) * t182 + qJ(5) * t140 + t438 * t84 + t441 * t64;
t41 = -qJ(5) * t139 - t438 * t64 + t441 * t84;
t20 = -pkin(11) * t97 - t38 * t445 + t41 * t449;
t476 = pkin(4) * t139 + pkin(5) * t223 + pkin(12) * t183 + t651;
t47 = -pkin(3) * t97 - t476;
t75 = t182 * t446 + t450 * t98;
t496 = pkin(10) * t75 + t20 * t446 + t450 * t47;
t62 = -pkin(4) * t213 + qJ(5) * t67;
t26 = -pkin(11) * t44 - qJ(5) * t652 - t445 * t62;
t31 = -pkin(3) * t44 - t673;
t42 = t213 * t446 + t45 * t450;
t495 = pkin(10) * t42 + t26 * t446 + t31 * t450;
t154 = t174 * t450 + t256 * t446;
t481 = pkin(4) * t240 - t107;
t79 = -pkin(3) * t173 - t481;
t136 = -pkin(4) * t256 + qJ(5) * t241 - t648;
t157 = -qJ(5) * t240 + t649;
t80 = -pkin(11) * t173 - t136 * t445 + t157 * t449;
t493 = pkin(10) * t154 + t446 * t80 + t450 * t79;
t156 = t185 * t450 + t259 * t446;
t85 = -pkin(3) * t184 + t354 - t553;
t141 = -pkin(4) * t259 + qJ(5) * t250 + t649;
t160 = -qJ(5) * t249 + t648;
t90 = -pkin(11) * t184 - t141 * t445 + t160 * t449;
t492 = pkin(10) * t156 + t446 * t90 + t450 * t85;
t105 = -pkin(3) * t131 - t667;
t123 = t133 * t450 + t274 * t446;
t60 = -pkin(4) * t274 + qJ(5) * t193 + t67;
t61 = -qJ(5) * t191 - t66;
t34 = -pkin(11) * t131 - t445 * t60 + t449 * t61;
t490 = pkin(10) * t123 + t105 * t450 + t34 * t446;
t164 = t445 * t457 - t449 * t460 + t454 + (-t445 * t676 - t449 * t561 - t280) * pkin(3);
t226 = -pkin(11) * t280 + t645;
t237 = t281 * t450 - t312 * t446;
t489 = pkin(10) * t237 + t164 * t450 + t226 * t446;
t166 = -pkin(3) * t283 + t196;
t227 = -pkin(11) * t283 + t644;
t242 = t284 * t450 + t314 * t446;
t488 = pkin(10) * t242 + t166 * t450 + t227 * t446;
t306 = t351 * t446 + t450 * t464;
t487 = pkin(10) * t306 + t228;
t114 = -pkin(11) * t246 - t134;
t229 = t248 * t450 - t333 * t446;
t485 = pkin(10) * t229 + t114 * t446 - t246 * t668;
t425 = t433 * t570;
t424 = t433 * t571;
t423 = (t436 - t437) * t598;
t422 = -t432 - t567;
t419 = -t568 - t432;
t414 = -t424 + t483;
t413 = t424 + t483;
t412 = -t425 + t484;
t399 = t447 * t466 + t451 * t469 - t434;
t387 = t406 - t405;
t377 = t381 * t619;
t376 = t381 * t623;
t349 = -t408 * t594 + t463;
t345 = t408 * t577 + t410 * t576 - t442 * t470;
t328 = (t393 * t449 + t394 * t445) * t404;
t327 = (t393 * t445 - t394 * t449) * t404;
t322 = t383 * t624 + (t408 * t442 - t576) * t410;
t321 = t382 * t623 + (-t410 * t442 - t577) * t408;
t309 = t332 * t445 + t394 * t628;
t307 = -t393 * t629 + t449 * t507;
t301 = t439 * t518 + t442 * t464;
t300 = t351 * t442 + t439 * t517;
t297 = t373 * t445 + t638;
t296 = t374 * t449 + t695;
t295 = -t349 * t439 + t442 * t519;
t294 = t349 * t442 + t439 * t519;
t293 = (-t361 * t441 + t363 * t438) * t404;
t292 = (-t361 * t438 - t363 * t441) * t404;
t289 = -t346 * t439 + t442 * t516;
t288 = t346 * t442 + t439 * t516;
t282 = t387 * t442 + t439 * t522;
t269 = t327 * t442 + t328 * t624 - t376;
t245 = t312 * t445 + t314 * t449;
t236 = -t292 * t445 + t293 * t449;
t235 = t292 * t449 + t293 * t445;
t234 = t309 * t442 + t439 * t508;
t233 = t307 * t442 + t439 * t509;
t225 = t297 * t442 + t439 * t523;
t224 = t296 * t442 + t439 * t524;
t209 = t283 * t442 + t439 * t525;
t204 = t320 * t439 - t442 * t528;
t203 = -t320 * t442 - t439 * t528;
t200 = t265 * t449 + t267 * t445;
t199 = t264 * t449 + t266 * t445;
t197 = t280 * t442 + t439 * t526;
t187 = t253 * t449 + t254 * t445;
t186 = t251 * t449 + t252 * t445;
t181 = pkin(2) * t295 - t276 * t442 + t439 * t549;
t180 = -pkin(3) * t314 + pkin(11) * t284 + t645;
t179 = pkin(2) * t289 - t275 * t442 + t439 * t550;
t176 = pkin(3) * t312 + pkin(11) * t281 - t644;
t175 = t245 * t442 + t439 * t530;
t167 = t246 * t442 + t439 * t529;
t165 = t235 * t442 + t236 * t624 - t376;
t161 = pkin(2) * t273 + t439 * t487;
t155 = pkin(2) * t204 + t228 * t663;
t144 = t211 * t449 + t212 * t445;
t130 = t190 * t449 + t192 * t445;
t125 = t187 * t442 + t439 * t510;
t124 = t186 * t442 + t439 * t511;
t120 = t170 * t449 + t172 * t445;
t119 = t169 * t449 + t171 * t445;
t118 = t200 * t442 + t439 * t531;
t117 = t199 * t442 + t439 * t532;
t116 = -pkin(3) * t262 + pkin(11) * t135;
t115 = t135 * t450 + t262 * t446;
t112 = t184 * t442 + t439 * t533;
t111 = pkin(3) * t333 + pkin(11) * t248 + t135;
t109 = t173 * t442 + t439 * t534;
t100 = t147 * t449 + t149 * t445;
t99 = t146 * t449 + t148 * t445;
t93 = t128 * t449 + t129 * t445;
t92 = t144 * t442 + t439 * t535;
t91 = t130 * t442 + t439 * t538;
t88 = t131 * t442 + t439 * t537;
t82 = -t134 * t439 + t442 * t536;
t81 = t134 * t442 + t439 * t536;
t78 = pkin(2) * t210 + t180 * t442 + t439 * t488;
t77 = -pkin(3) * t259 + pkin(11) * t185 + t141 * t449 + t160 * t445;
t76 = pkin(2) * t198 + t176 * t442 + t439 * t489;
t73 = -pkin(3) * t256 + pkin(11) * t174 + t136 * t449 + t157 * t445;
t69 = t120 * t442 + t439 * t539;
t68 = t119 * t442 + t439 * t540;
t59 = pkin(2) * t168 + t111 * t442 + t439 * t485;
t58 = t100 * t442 + t439 * t541;
t57 = t439 * t542 + t442 * t99;
t55 = t439 * t544 + t442 * t97;
t53 = t439 * t545 + t442 * t95;
t52 = t439 * t547 + t442 * t93;
t48 = t439 * t546 + t442 * t86;
t39 = pkin(2) * t82 + t116 * t442 + (pkin(10) * t115 + t134 * t551) * t439;
t32 = -pkin(3) * t274 + pkin(11) * t133 + t445 * t61 + t449 * t60;
t30 = pkin(2) * t113 + t439 * t492 + t442 * t77;
t28 = pkin(2) * t110 + t439 * t493 + t442 * t73;
t24 = t439 * t543 + t44 * t442;
t23 = -pkin(3) * t213 + pkin(11) * t45 - qJ(5) * t658 + t449 * t62;
t18 = -pkin(3) * t182 + pkin(11) * t98 + t38 * t449 + t41 * t445;
t17 = -pkin(3) * t177 + pkin(11) * t96 + t37 * t449 + t40 * t445;
t13 = -pkin(3) * t151 + pkin(11) * t87 + t29 * t449 + t33 * t445;
t12 = pkin(2) * t89 + t32 * t442 + t439 * t490;
t8 = t21 * t442 + t439 * t548;
t7 = pkin(2) * t56 + t18 * t442 + t439 * t496;
t6 = pkin(2) * t54 + t17 * t442 + t439 * t497;
t5 = pkin(2) * t25 + t23 * t442 + t439 * t495;
t3 = pkin(2) * t49 + t13 * t442 + t439 * t498;
t2 = -pkin(3) * t50 + pkin(11) * t22 + t11 * t449 + t16 * t445;
t1 = pkin(2) * t9 + t2 * t442 + t439 * t499;
t70 = [0, 0, 0, 0, 0, qJDD(1), t500, t501, 0, 0, (t435 * t589 - t601 * t680) * t447, t443 * t423 + (t447 * t414 + (t425 + t484) * t451) * t440, t443 * t412 + (t608 + t451 * (t432 - t568)) * t440, (t435 * t588 + t602 * t680) * t451, t443 * t413 + (t447 * (-t432 + t567) + t606) * t440, t443 * t555, (-t398 + pkin(1) * (t420 * t451 + t422 * t447)) * t443 + (t451 * t415 + pkin(1) * t414 + pkin(9) * (t422 * t451 - t608)) * t440, -t415 * t622 - t443 * t399 + pkin(1) * ((t419 * t451 - t421 * t447) * t443 + (-qJD(1) * t558 - t589) * t435) + (-t447 * t419 - t606) * t665, t440 * t609 + t399 * t621 + pkin(1) * ((-t412 * t451 + t413 * t447) * t443 - t682 * t435 * t597) + (t447 * t412 + t413 * t451) * t665, pkin(1) * (t440 * t415 + (-t398 * t451 + t399 * t447) * t443) + (t399 * t451 + t609) * t665, t443 * t322 + (t447 * (t383 * t450 + t446 * t625) + t451 * (t383 * t620 + (-t408 * t439 - t574) * t410)) * t440, t443 * t282 + (t447 * (-t346 * t450 - t350 * t446) + t451 * (-t387 * t439 + t442 * t522)) * t440, t443 * t300 + (t447 * (-t396 * t446 + t686) + t451 * (-t351 * t439 + t442 * t517)) * t440, t443 * t321 + (t447 * (-t382 * t446 - t418 * t627) + t451 * (t382 * t619 + (t410 * t439 - t575) * t408)) * t440, t443 * t301 + (t447 * (t395 * t450 + t634) + t451 * (-t439 * t464 + t442 * t518)) * t440, (-t410 * t446 + t627) * t418 * t622 + (t408 * t575 + t410 * t574 + t439 * t470) * t621 + t443 * t345, (t179 + pkin(1) * (t289 * t451 + t329 * t447)) * t443 + (t447 * (-t641 + (-t288 * t439 - t289 * t442) * pkin(10)) + t451 * (-pkin(2) * t288 + t275 * t439 + t442 * t550) - pkin(1) * t288 + pkin(9) * (-t447 * t289 + t329 * t451)) * t440, (t181 + pkin(1) * (t295 * t451 + t334 * t447)) * t443 + (t447 * (-t640 + (-t294 * t439 - t295 * t442) * pkin(10)) + t451 * (-pkin(2) * t294 + t276 * t439 + t442 * t549) - pkin(1) * t294 + pkin(9) * (-t447 * t295 + t334 * t451)) * t440, (t161 + pkin(1) * (t273 * t451 + t306 * t447)) * t443 + (t447 * t528 + pkin(9) * (t306 * t451 - t610) + t506 * (t368 * t442 + t439 * t521) + (-pkin(10) * t610 + t451 * t487) * t442) * t440, (t155 + pkin(1) * (t204 * t451 + t228 * t447)) * t443 + (pkin(9) * (-t447 * t204 + t228 * t451) + (-pkin(1) - t670) * t203 + (t447 * (-t203 * t439 - t204 * t442) + t228 * t618) * pkin(10)) * t440, t443 * t234 + (t447 * (t310 * t450 - t579) + t451 * (-t309 * t439 + t442 * t508)) * t440, t443 * t175 + (t447 * (t247 * t450 + t364 * t446) + t451 * (-t245 * t439 + t442 * t530)) * t440, t443 * t224 + (t447 * (t298 * t450 + t315 * t446) + t451 * (-t296 * t439 + t442 * t524)) * t440, t443 * t233 + (t447 * (t308 * t450 + t579) + t451 * (-t307 * t439 + t442 * t509)) * t440, t443 * t225 + (t447 * (t299 * t450 - t311 * t446) + t451 * (-t297 * t439 + t442 * t523)) * t440, t443 * t269 + (t447 * (t328 * t450 + t632) + t451 * (-t327 * t439 + t328 * t620 - t377)) * t440, (t76 + pkin(1) * (t198 * t451 + t237 * t447)) * t443 + (t447 * (-t164 * t446 - t197 * t663 + t226 * t450) + t451 * (-pkin(2) * t197 - t176 * t439) - pkin(1) * t197 + pkin(9) * (t237 * t451 - t612) + (-pkin(10) * t612 + t451 * t489) * t442) * t440, (t78 + pkin(1) * (t210 * t451 + t242 * t447)) * t443 + (t447 * (-t166 * t446 - t209 * t663 + t227 * t450) + t451 * (-pkin(2) * t209 - t180 * t439) - pkin(1) * t209 + pkin(9) * (t242 * t451 - t611) + (-pkin(10) * t611 + t451 * t488) * t442) * t440, (t59 + pkin(1) * (t168 * t451 + t229 * t447)) * t443 + (t447 * (t114 * t450 - t167 * t663 + t246 * t669) + t451 * (-pkin(2) * t167 - t111 * t439) - pkin(1) * t167 + pkin(9) * (t229 * t451 - t613) + (-pkin(10) * t613 + t451 * t485) * t442) * t440, (t39 + pkin(1) * (t115 * t447 + t451 * t82)) * t443 + (t451 * (-pkin(2) * t81 - t116 * t439) - pkin(1) * t81 + pkin(9) * (t115 * t451 - t447 * t82) + (t447 * (-t439 * t81 - t442 * t82) + t115 * t618) * pkin(10) + (t447 * (-pkin(11) * t450 + t669) + t551 * t618) * t134) * t440, t443 * t125 + (t447 * (t189 * t450 + t581) + t451 * (-t187 * t439 + t442 * t510)) * t440, t443 * t91 + (t447 * (t132 * t450 + t318 * t446) + t451 * (-t130 * t439 + t442 * t538)) * t440, t443 * t117 + (t447 * (t201 * t450 + t260 * t446) + t451 * (-t199 * t439 + t442 * t532)) * t440, t443 * t124 + (t447 * (t188 * t450 - t581) + t451 * (-t186 * t439 + t442 * t511)) * t440, t443 * t118 + (t447 * (t202 * t450 - t257 * t446) + t451 * (-t200 * t439 + t442 * t531)) * t440, t443 * t165 + (t447 * (t236 * t450 + t632) + t451 * (-t235 * t439 + t236 * t620 - t377)) * t440, (t28 + pkin(1) * (t110 * t451 + t154 * t447)) * t443 + (t447 * (-t109 * t663 - t446 * t79 + t450 * t80) + t451 * (-pkin(2) * t109 - t439 * t73) - pkin(1) * t109 + pkin(9) * (t154 * t451 - t615) + (-pkin(10) * t615 + t451 * t493) * t442) * t440, (t30 + pkin(1) * (t113 * t451 + t156 * t447)) * t443 + (t447 * (-t112 * t663 - t446 * t85 + t450 * t90) + t451 * (-pkin(2) * t112 - t439 * t77) - pkin(1) * t112 + pkin(9) * (t156 * t451 - t614) + (-pkin(10) * t614 + t451 * t492) * t442) * t440, (t12 + pkin(1) * (t123 * t447 + t451 * t89)) * t443 + (t447 * (-t105 * t446 + t34 * t450 - t663 * t88) + t451 * (-pkin(2) * t88 - t32 * t439) - pkin(1) * t88 + pkin(9) * (t123 * t451 - t653) + (-pkin(10) * t653 + t451 * t490) * t442) * t440, (t5 + pkin(1) * (t25 * t451 + t42 * t447)) * t443 + (t447 * (-t24 * t663 + t26 * t450 - t31 * t446) + t451 * (-pkin(2) * t24 - t23 * t439) - pkin(1) * t24 + pkin(9) * (t42 * t451 - t657) + (-pkin(10) * t657 + t451 * t495) * t442) * t440, t443 * t69 + (t447 * (t122 * t450 + t216 * t446) + t451 * (-t120 * t439 + t442 * t539)) * t440, t443 * t52 + (t447 * (t150 * t446 + t450 * t94) + t451 * (-t439 * t93 + t442 * t547)) * t440, t443 * t57 + (t447 * (t101 * t450 + t205 * t446) + t451 * (-t439 * t99 + t442 * t542)) * t440, t443 * t68 + (t447 * (t121 * t450 - t214 * t446) + t451 * (-t119 * t439 + t442 * t540)) * t440, t443 * t58 + (t447 * (t102 * t450 + t206 * t446) + t451 * (-t100 * t439 + t442 * t541)) * t440, t443 * t92 + (t447 * (t145 * t450 + t243 * t446) + t451 * (-t144 * t439 + t442 * t535)) * t440, (t6 + pkin(1) * (t447 * t74 + t451 * t54)) * t443 + (t447 * (t19 * t450 - t446 * t46 - t53 * t663) + t451 * (-pkin(2) * t53 - t17 * t439) - pkin(1) * t53 + pkin(9) * (t451 * t74 - t655) + (-pkin(10) * t655 + t451 * t497) * t442) * t440, (t7 + pkin(1) * (t447 * t75 + t451 * t56)) * t443 + (t447 * (t20 * t450 - t446 * t47 - t55 * t663) + t451 * (-pkin(2) * t55 - t18 * t439) - pkin(1) * t55 + pkin(9) * (t451 * t75 - t654) + (-pkin(10) * t654 + t451 * t496) * t442) * t440, (t3 + pkin(1) * (t447 * t65 + t451 * t49)) * t443 + (t447 * (t14 * t450 - t27 * t446 - t48 * t663) + t451 * (-pkin(2) * t48 - t13 * t439) - pkin(1) * t48 + pkin(9) * (t451 * t65 - t656) + (-pkin(10) * t656 + t451 * t498) * t442) * t440, (t1 + pkin(1) * (t15 * t447 + t451 * t9)) * t443 + (t447 * (-t10 * t446 + t4 * t450 - t663 * t8) + t451 * (-pkin(2) * t8 - t2 * t439) - pkin(1) * t8 + pkin(9) * (t15 * t451 - t660) + (-pkin(10) * t660 + t451 * t499) * t442) * t440; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t431, t423, t412, t431, t413, t555, -t398, -t399, 0, 0, t322, t282, t300, t321, t301, t345, t179, t181, t161, t155, t234, t175, t224, t233, t225, t269, t76, t78, t59, t39, t125, t91, t117, t124, t118, t165, t28, t30, t12, t5, t69, t52, t57, t68, t58, t92, t6, t7, t3, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t626, t387, t351, -t626, t464, -t470, -t275, -t276, 0, 0, t309, t245, t296, t307, t297, t327, t176, t180, t111, t116, t187, t130, t199, t186, t200, t235, t73, t77, t32, t23, t120, t93, t99, t119, t100, t144, t17, t18, t13, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t365, t364, t315, t365, -t311, t381, -t195, -t196, 0, 0, t319, t318, t260, -t319, -t257, t381, t481, t553 + 0.2e1 * t599, t667, t673, t216, t150, t205, -t214, t206, t243, t477, t476, t475, t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, t259, t274, t213, 0, 0, 0, 0, 0, 0, t177, t182, t151, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, t290, t222, -t291, -t218, t277, -t71, -t72, 0, 0;];
tauJ_reg  = t70;