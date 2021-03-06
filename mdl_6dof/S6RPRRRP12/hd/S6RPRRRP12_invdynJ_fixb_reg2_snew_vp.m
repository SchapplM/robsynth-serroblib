% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRRP12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 02:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRRP12_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP12_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:22:33
% EndTime: 2019-05-06 02:23:43
% DurationCPUTime: 46.68s
% Computational Cost: add. (178436->825), mult. (555766->1252), div. (0->0), fcn. (477244->14), ass. (0->512)
t365 = sin(pkin(12));
t367 = sin(pkin(6));
t368 = cos(pkin(12));
t370 = cos(pkin(6));
t376 = cos(qJ(3));
t369 = cos(pkin(7));
t373 = sin(qJ(3));
t528 = t369 * t373;
t366 = sin(pkin(7));
t534 = t366 * t373;
t392 = t370 * t534 + (t365 * t376 + t368 * t528) * t367;
t335 = t392 * qJD(1);
t531 = t367 * t368;
t507 = t366 * t531;
t522 = qJD(1) * t370;
t342 = -qJD(1) * t507 + t369 * t522 + qJD(3);
t372 = sin(qJ(4));
t375 = cos(qJ(4));
t323 = t335 * t372 - t375 * t342;
t318 = qJD(5) + t323;
t602 = t318 ^ 2;
t325 = t335 * t375 + t342 * t372;
t371 = sin(qJ(5));
t374 = cos(qJ(5));
t523 = qJD(1) * t367;
t537 = t365 * t373;
t530 = t368 * t369;
t505 = t367 * t530;
t533 = t366 * t376;
t637 = -t370 * t533 - t376 * t505;
t333 = qJD(1) * t637 + t523 * t537;
t412 = qJD(4) + t333;
t299 = t325 * t371 - t374 * t412;
t603 = t299 ^ 2;
t266 = t603 - t602;
t389 = t392 * qJDD(1);
t313 = -t333 * qJD(3) + t389;
t514 = qJDD(1) * t367;
t497 = t368 * t514;
t513 = qJDD(1) * t370;
t341 = -t366 * t497 + t369 * t513 + qJDD(3);
t491 = t372 * t313 - t375 * t341;
t277 = -t325 * qJD(4) - t491;
t276 = qJDD(5) - t277;
t301 = t374 * t325 + t371 * t412;
t547 = t301 * t299;
t211 = -t547 - t276;
t559 = t211 * t371;
t169 = -t266 * t374 - t559;
t270 = t318 * t301;
t443 = -t375 * t313 - t372 * t341;
t278 = -t323 * qJD(4) - t443;
t498 = t365 * t514;
t427 = qJDD(1) * t637 + t373 * t498;
t404 = -qJD(3) * t335 - t427;
t397 = qJDD(4) - t404;
t492 = -t278 * t371 + t374 * t397;
t434 = qJD(5) * t301 - t492;
t189 = -t270 + t434;
t115 = t169 * t372 - t189 * t375;
t119 = t169 * t375 + t189 * t372;
t558 = t211 * t374;
t165 = -t266 * t371 + t558;
t463 = t119 * t373 - t165 * t376;
t730 = t370 * (t369 * t115 + t366 * t463) + (t365 * (t119 * t376 + t373 * t165) + t368 * (-t366 * t115 + t369 * t463)) * t367;
t386 = -t374 * t278 - t371 * t397;
t382 = -t299 * qJD(5) - t386;
t548 = t299 * t318;
t617 = -t548 + t382;
t561 = t617 * t371;
t625 = t270 + t434;
t131 = t625 * t374 + t561;
t297 = t301 ^ 2;
t622 = t297 - t603;
t101 = t131 * t375 - t372 * t622;
t129 = -t625 * t371 + t374 * t617;
t471 = t101 * t373 + t129 * t376;
t99 = t131 * t372 + t375 * t622;
t728 = t367 * (t368 * (-t366 * t99 + t369 * t471) + t365 * (t101 * t376 - t373 * t129)) + t370 * (t366 * t471 + t369 * t99);
t618 = -t547 + t276;
t557 = t618 * t371;
t615 = -t602 - t603;
t635 = t374 * t615 - t557;
t660 = t372 * t635 - t375 * t625;
t556 = t618 * t374;
t636 = t371 * t615 + t556;
t659 = t372 * t625 + t375 * t635;
t675 = t373 * t659 - t376 * t636;
t688 = -t366 * t660 + t369 * t675;
t700 = t365 * t688;
t715 = pkin(9) * t700;
t673 = t373 * t636 + t376 * t659;
t689 = t366 * t675 + t369 * t660;
t711 = qJ(2) * (t368 * t673 - t700) - pkin(1) * t689;
t709 = pkin(2) * t689;
t623 = -t297 - t602;
t158 = t374 * t623 + t559;
t708 = pkin(3) * t158;
t707 = pkin(4) * t158;
t706 = pkin(11) * t158;
t160 = -t371 * t623 + t558;
t705 = pkin(11) * t160;
t704 = t158 * t376;
t703 = t160 * t372;
t702 = t160 * t375;
t699 = t373 * t158;
t589 = pkin(9) * t366;
t697 = t589 * t689;
t696 = pkin(1) * (t365 * t673 + t368 * t688) + pkin(2) * t688;
t267 = -t297 + t602;
t663 = t267 * t374 + t557;
t616 = t548 + t382;
t662 = -t267 * t371 + t556;
t674 = t372 * t616 + t375 * t662;
t676 = t372 * t662 - t375 * t616;
t687 = t373 * t674 - t376 * t663;
t695 = t370 * (t366 * t687 + t369 * t676) + (t365 * (t373 * t663 + t376 * t674) + t368 * (-t366 * t676 + t369 * t687)) * t367;
t694 = pkin(9) * t673;
t685 = pkin(3) * t660;
t684 = pkin(10) * t660;
t677 = -pkin(3) * t636 + pkin(10) * t659;
t595 = sin(qJ(1));
t596 = cos(qJ(1));
t429 = t595 * g(1) - t596 * g(2);
t406 = qJDD(1) * pkin(1) + t429;
t493 = t370 * g(3) - qJDD(2);
t379 = -pkin(2) * t497 - t367 * t406 - t493;
t396 = t370 * t406;
t504 = pkin(2) * t513;
t380 = -g(3) * t531 + t368 * t396 + t504;
t430 = t596 * g(1) + t595 * g(2);
t398 = qJ(2) * t514 - t430;
t364 = t370 ^ 2;
t417 = (t365 * pkin(1) + t364 * t589) * qJD(1);
t539 = t365 * t366;
t593 = pkin(2) * t368;
t440 = -pkin(9) * t539 - t593;
t588 = pkin(9) * t369;
t594 = pkin(2) * t365;
t361 = t365 ^ 2;
t363 = t368 ^ 2;
t620 = -t361 - t363;
t672 = -2 * qJD(2);
t656 = t365 * t672;
t247 = qJD(1) * (t366 * t417 + t367 * (t522 * ((qJ(2) + 0.2e1 * t588) * t366 * t368 - t369 * t594) - (t440 * t539 + t369 * (t588 * t620 - qJ(2))) * t523 + t366 * t656)) + t366 * (-t365 * t398 + t380) - t369 * t379;
t670 = pkin(4) * t636;
t669 = pkin(11) * t635;
t668 = pkin(11) * t636;
t545 = t318 * t371;
t264 = t301 * t545;
t544 = t318 * t374;
t509 = t299 * t544;
t484 = t264 - t509;
t607 = t276 * t372 + t375 * t484;
t610 = -t375 * t276 + t372 * t484;
t626 = (t299 * t371 + t301 * t374) * t318;
t632 = t373 * t607 + t376 * t626;
t658 = t370 * (t366 * t632 + t369 * t610) + (t365 * (-t373 * t626 + t376 * t607) + t368 * (-t366 * t610 + t369 * t632)) * t367;
t183 = t299 * t545 - t374 * t434;
t438 = t371 * t434 + t509;
t511 = t372 * t547;
t608 = t375 * t438 - t511;
t510 = t375 * t547;
t609 = t372 * t438 + t510;
t633 = -t183 * t376 + t373 * t608;
t657 = t370 * (t366 * t633 + t369 * t609) + (t365 * (t373 * t183 + t376 * t608) + t368 * (-t366 * t609 + t369 * t633)) * t367;
t621 = t297 + t603;
t655 = pkin(4) * t621;
t654 = qJ(6) * t617;
t647 = t372 * t621;
t292 = t325 * t323;
t624 = -t292 + t397;
t645 = t372 * t624;
t641 = t375 * t621;
t639 = t375 * t624;
t408 = t412 ^ 2;
t535 = t366 * t370;
t435 = t505 + t535;
t629 = pkin(9) * t435;
t543 = t335 * t333;
t403 = t341 - t543;
t628 = t373 * t403;
t627 = t376 * t403;
t304 = t412 * t323;
t243 = t278 - t304;
t328 = t342 * t333;
t289 = -t328 + t313;
t529 = t369 * t370;
t619 = -t507 + t529;
t250 = pkin(5) * t299 - qJ(6) * t301;
t585 = t367 * g(3);
t388 = t396 - t585;
t418 = t440 * t367;
t495 = qJ(2) + t588;
t481 = t370 * t495;
t599 = 2 * qJD(2);
t282 = t368 * t398 + t388 * t365 + qJDD(1) * t629 + ((-pkin(1) * t368 - pkin(2) * t364) * qJD(1) + (t368 * t599 + (t365 * t481 + t368 * t418) * qJD(1)) * t367) * qJD(1);
t433 = -pkin(1) + t440;
t600 = qJD(1) ^ 2;
t519 = t367 * t600;
t527 = t369 * t376;
t532 = t367 * qJ(2);
t203 = t373 * t282 - (t504 + t388 * t368 + (-t495 * t514 + t430) * t365 + (t417 + (t656 + (-t365 * t418 + t368 * t481) * qJD(1)) * t367) * qJD(1)) * t527 - ((qJDD(1) * t433 - t429) * t367 + (t370 * t594 - t532 + (t367 * t369 * t620 - t368 * t535) * pkin(9)) * t519 - t493) * t533;
t312 = pkin(3) * t333 - pkin(10) * t335;
t601 = t342 ^ 2;
t172 = -t341 * pkin(3) - t601 * pkin(10) + t335 * t312 + t203;
t407 = t412 * t325;
t123 = -t243 * pkin(11) + (-t277 + t407) * pkin(4) + t172;
t204 = t376 * t282 + (t369 * (t600 * pkin(1) + t430) + ((-pkin(9) * t366 ^ 2 - t369 * t495) * qJDD(1) + (pkin(2) * qJD(1) * t435 + t369 * t672) * qJD(1)) * t367) * t537 + (t366 * t379 + t369 * t380 + ((-t366 * t367 + t368 * t529) * t532 + t619 * t629) * t600) * t373;
t173 = -t601 * pkin(3) + t341 * pkin(10) - t333 * t312 + t204;
t516 = qJD(3) + t342;
t285 = t516 * t335 + t427;
t377 = t285 * pkin(3) - t289 * pkin(10) - t247;
t126 = t375 * t173 + t372 * t377;
t283 = pkin(4) * t323 - pkin(11) * t325;
t96 = -pkin(4) * t408 + pkin(11) * t397 - t323 * t283 + t126;
t76 = t371 * t123 + t374 * t96;
t490 = -t276 * qJ(6) + t299 * t250 - t76;
t614 = -(t623 + t602) * pkin(5) - qJ(6) * t211 - t490;
t520 = qJD(6) * t318;
t315 = 0.2e1 * t520;
t480 = t315 - t490;
t57 = -pkin(5) * t602 + t480;
t75 = -t374 * t123 + t371 * t96;
t64 = -t276 * pkin(5) - qJ(6) * t602 + t250 * t301 + qJDD(6) + t75;
t35 = t371 * t64 + t374 * t57;
t496 = qJ(6) * t371 + pkin(4);
t590 = pkin(5) * t374;
t125 = t173 * t372 - t375 * t377;
t95 = -t397 * pkin(4) - t408 * pkin(11) + t283 * t325 + t125;
t393 = t434 * pkin(5) - t654 + t95;
t73 = (pkin(5) * t318 - 0.2e1 * qJD(6)) * t301 + t393;
t613 = -(t496 + t590) * t73 + pkin(11) * t35;
t391 = 0.2e1 * qJD(6) * t301 - t393;
t63 = (-t625 - t270) * pkin(5) + t391;
t612 = t374 * t63 - t496 * t625 + t669;
t62 = -pkin(5) * t270 + t391 + t654;
t611 = -t705 + (pkin(4) + t590) * t617 + t371 * t62;
t186 = t301 * t544 + t371 * t382;
t187 = t374 * t382 - t264;
t485 = t375 * t187 + t511;
t486 = t372 * t187 - t510;
t604 = t370 * (-t186 * t533 + t369 * t486 + t485 * t534) + (t365 * (t373 * t186 + t376 * t485) + t368 * (-t186 * t527 - t366 * t486 + t485 * t528)) * t367;
t321 = t323 ^ 2;
t322 = t325 ^ 2;
t331 = t333 ^ 2;
t332 = t335 ^ 2;
t80 = t125 * t372 + t375 * t126;
t70 = t373 * t172 + t376 * t80;
t598 = pkin(9) * t70;
t592 = pkin(3) * t376;
t591 = pkin(4) * t372;
t23 = t35 * t372 - t375 * t73;
t24 = t35 * t375 + t372 * t73;
t34 = t371 * t57 - t374 * t64;
t476 = t24 * t373 - t34 * t376;
t4 = -t366 * t23 + t369 * t476;
t587 = t365 * t4;
t40 = t371 * t75 + t374 * t76;
t31 = t372 * t40 - t375 * t95;
t32 = t372 * t95 + t375 * t40;
t39 = t371 * t76 - t374 * t75;
t475 = t32 * t373 - t376 * t39;
t8 = -t366 * t31 + t369 * t475;
t586 = t365 * t8;
t472 = -t172 * t376 + t373 * t80;
t79 = -t125 * t375 + t126 * t372;
t42 = -t366 * t79 + t369 * t472;
t584 = t365 * t42;
t562 = t616 * t374;
t128 = -t189 * t371 - t562;
t563 = t616 * t371;
t132 = -t189 * t374 + t563;
t92 = t132 * t375 - t647;
t474 = -t128 * t376 + t373 * t92;
t90 = t132 * t372 + t641;
t49 = -t366 * t90 + t369 * t474;
t583 = t365 * t49;
t191 = (-qJD(5) + t318) * t301 + t492;
t130 = t191 * t371 - t562;
t134 = t191 * t374 + t563;
t93 = t134 * t375 - t647;
t473 = -t130 * t376 + t373 * t93;
t91 = t134 * t372 + t641;
t50 = -t366 * t91 + t369 * t473;
t582 = t365 * t50;
t103 = t375 * t617 - t703;
t105 = -t372 * t617 - t702;
t469 = t105 * t373 + t704;
t60 = -t366 * t103 + t369 * t469;
t581 = t365 * t60;
t196 = (qJD(5) + t318) * t299 + t386;
t109 = t196 * t375 + t703;
t111 = -t196 * t372 + t702;
t467 = t111 * t373 - t704;
t67 = -t366 * t109 + t369 * t467;
t579 = t365 * t67;
t240 = -t325 * t333 + t491;
t244 = t278 + t304;
t177 = -t240 * t372 - t244 * t375;
t179 = -t240 * t375 + t244 * t372;
t262 = t321 + t322;
t458 = t179 * t373 + t262 * t376;
t98 = -t366 * t177 + t369 * t458;
t577 = t365 * t98;
t576 = t371 * t95;
t575 = t373 * t79;
t574 = t374 * t95;
t573 = t376 * t79;
t572 = qJ(6) * t374;
t279 = -t408 - t321;
t212 = t279 * t372 + t639;
t213 = t279 * t375 - t645;
t411 = 0.2e1 * qJD(4) + t333;
t241 = -t325 * t411 - t491;
t454 = t213 * t373 + t241 * t376;
t138 = -t366 * t212 + t369 * t454;
t571 = t138 * t365;
t281 = -t322 - t408;
t259 = t292 + t397;
t553 = t259 * t372;
t221 = t281 * t375 - t553;
t552 = t259 * t375;
t222 = -t281 * t372 - t552;
t245 = t323 * t411 + t443;
t453 = t222 * t373 + t245 * t376;
t140 = -t366 * t221 + t369 * t453;
t570 = t140 * t365;
t569 = t172 * t372;
t568 = t172 * t375;
t298 = -t331 - t332;
t290 = t328 + t313;
t394 = (-qJD(3) + t342) * t335 - t427;
t449 = -t290 * t376 + t373 * t394;
t206 = -t366 * t298 + t369 * t449;
t560 = t206 * t365;
t555 = t247 * t373;
t554 = t247 * t376;
t306 = -t341 - t543;
t546 = t306 * t376;
t542 = t335 * t373;
t362 = t367 ^ 2;
t541 = t361 * t362;
t540 = t362 * t363;
t538 = t365 * t367;
t525 = t373 * t306;
t518 = t368 * t600;
t517 = t370 * t600;
t515 = qJDD(1) * t362;
t508 = t376 * t292;
t503 = t373 * t292;
t502 = -pkin(4) * t375 - pkin(3);
t499 = t367 * t517;
t488 = -pkin(4) * t95 + pkin(11) * t40;
t487 = -pkin(5) * t64 + qJ(6) * t57;
t479 = -pkin(5) * t616 - qJ(6) * t189;
t305 = -t601 - t331;
t261 = t305 * t376 - t628;
t478 = pkin(9) * t261 + t554;
t309 = -t332 - t601;
t263 = -t373 * t309 + t546;
t477 = pkin(9) * t263 - t555;
t178 = t241 * t375 - t243 * t372;
t291 = t322 - t321;
t459 = t178 * t373 - t291 * t376;
t455 = t203 * t376 - t373 * t204;
t143 = t373 * t203 + t204 * t376;
t303 = -t322 + t408;
t232 = -t303 * t372 + t639;
t452 = t232 * t373 - t244 * t376;
t302 = t321 - t408;
t233 = t302 * t375 - t553;
t451 = t233 * t373 + t240 * t376;
t450 = -t285 * t373 + t289 * t376;
t447 = t305 * t373 + t627;
t446 = t309 * t376 + t525;
t326 = t331 - t601;
t445 = t326 * t373 - t546;
t327 = -t332 + t601;
t444 = t327 * t376 + t628;
t383 = (-pkin(1) * qJD(1) + t367 * t599) * qJD(1) + t398;
t390 = qJ(2) * t519 + t406;
t384 = t370 * t390 - t585;
t319 = t365 * t383 - t368 * t384;
t320 = t365 * t384 + t368 * t383;
t442 = t365 * t319 + t368 * t320;
t441 = -t333 * t373 - t335 * t376;
t400 = t375 * t304;
t237 = -t372 * t277 + t400;
t437 = t237 * t373 + t508;
t401 = t372 * t407;
t239 = t375 * t278 - t401;
t436 = t239 * t373 - t508;
t10 = t24 * t376 + t373 * t34;
t13 = -pkin(4) * t34 - t487;
t16 = -pkin(11) * t34 + (pkin(5) * t371 - t572) * t73;
t2 = -pkin(10) * t23 - t13 * t372 + t16 * t375;
t5 = -pkin(3) * t23 - t613;
t431 = pkin(9) * t10 + t2 * t373 + t376 * t5;
t11 = -pkin(3) * t31 - t488;
t15 = t32 * t376 + t373 * t39;
t9 = -pkin(10) * t31 + (-pkin(11) * t375 + t591) * t39;
t428 = pkin(9) * t15 + t11 * t376 + t373 * t9;
t426 = -pkin(4) * t625 - t574 + t669;
t425 = pkin(4) * t196 + t576 + t705;
t53 = (t621 - t602) * pkin(5) + t480;
t54 = qJ(6) * t621 + t64;
t26 = -pkin(11) * t128 - t371 * t53 + t374 * t54;
t81 = -pkin(4) * t128 - t479;
t14 = -pkin(10) * t90 + t26 * t375 - t372 * t81;
t405 = pkin(11) * t132 + t371 * t54 + t374 * t53 + t655;
t21 = -pkin(3) * t90 - t405;
t77 = t373 * t128 + t376 * t92;
t424 = pkin(9) * t77 + t14 * t373 + t21 * t376;
t43 = -0.2e1 * t520 - t614 + t707;
t44 = -pkin(5) * t561 + t374 * t62 + t706;
t19 = -pkin(10) * t103 - t372 * t43 + t375 * t44;
t36 = -pkin(3) * t103 - t611;
t84 = t105 * t376 - t699;
t423 = pkin(9) * t84 + t19 * t373 + t36 * t376;
t387 = pkin(5) * t618 + qJ(6) * t615 - t64;
t45 = -t387 - t670;
t46 = -t371 * t63 - t572 * t625 - t668;
t20 = -t372 * t45 + t375 * t46 - t684;
t38 = -t612 - t685;
t422 = t20 * t373 + t376 * t38 + t694;
t37 = -pkin(11) * t130 - t39;
t25 = -pkin(10) * t91 + t130 * t591 + t37 * t375;
t413 = pkin(11) * t134 + t40 + t655;
t27 = -pkin(3) * t91 - t413;
t78 = t373 * t130 + t376 * t93;
t421 = pkin(9) * t78 + t25 * t373 + t27 * t376;
t55 = t75 - t670;
t82 = t576 - t668;
t30 = -t372 * t55 + t375 * t82 - t684;
t51 = -t426 - t685;
t420 = t30 * t373 + t376 * t51 + t694;
t56 = t76 - t707;
t83 = t574 - t706;
t33 = -pkin(10) * t109 - t372 * t56 + t375 * t83;
t52 = -pkin(3) * t109 - t425;
t86 = t111 * t376 + t699;
t419 = pkin(9) * t86 + t33 * t373 + t376 * t52;
t141 = -pkin(10) * t212 + t569;
t157 = t213 * t376 - t373 * t241;
t88 = -pkin(3) * t212 + t125;
t416 = pkin(9) * t157 + t141 * t373 + t376 * t88;
t142 = -pkin(10) * t221 + t568;
t162 = t222 * t376 - t373 * t245;
t89 = -pkin(3) * t221 + t126;
t415 = pkin(9) * t162 + t142 * t373 + t376 * t89;
t246 = t373 * t290 + t376 * t394;
t414 = pkin(9) * t246 + t143;
t152 = t179 * t376 - t373 * t262;
t72 = -pkin(10) * t177 - t79;
t410 = pkin(9) * t152 - t177 * t592 + t373 * t72;
t402 = t372 * t304;
t399 = t375 * t407;
t395 = t376 * t404;
t257 = -t400 + t401;
t385 = t373 * t257 - t376 * t397;
t357 = t368 * t499;
t356 = t365 * t499;
t355 = t362 * t365 * t518;
t351 = (-t364 - t540) * t600;
t350 = (-t364 - t541) * t600;
t349 = t620 * t362 * t600;
t348 = -t355 + t513;
t347 = t355 + t513;
t346 = -t356 + t497;
t345 = t356 + t497;
t344 = -t357 + t498;
t343 = t357 + t498;
t336 = t367 * t390 + t493;
t314 = t332 - t331;
t288 = -t516 * t333 + t389;
t256 = -t402 - t399;
t238 = t372 * t278 + t399;
t236 = t375 * t277 + t402;
t231 = t302 * t372 + t552;
t230 = t303 * t375 + t645;
t229 = -t366 * t288 + t369 * t446;
t228 = t369 * t288 + t366 * t446;
t224 = -t366 * t285 + t369 * t447;
t223 = t369 * t285 + t366 * t447;
t205 = t369 * t298 + t366 * t449;
t176 = t241 * t372 + t243 * t375;
t139 = t369 * t221 + t366 * t453;
t137 = t369 * t212 + t366 * t454;
t136 = t366 * t247 - t369 * t455;
t135 = -t369 * t247 - t366 * t455;
t108 = pkin(3) * t245 + pkin(10) * t222 + t569;
t107 = pkin(3) * t241 + pkin(10) * t213 - t568;
t97 = t369 * t177 + t366 * t458;
t71 = -pkin(3) * t172 + pkin(10) * t80;
t69 = pkin(3) * t262 + pkin(10) * t179 + t80;
t65 = t369 * t109 + t366 * t467;
t58 = t369 * t103 + t366 * t469;
t48 = t366 * t473 + t369 * t91;
t47 = t366 * t474 + t369 * t90;
t41 = t366 * t472 + t369 * t79;
t29 = pkin(10) * t111 + t372 * t83 + t375 * t56 - t708;
t28 = t372 * t82 + t375 * t55 + t677;
t22 = pkin(10) * t93 + t502 * t130 + t37 * t372;
t18 = t372 * t46 + t375 * t45 + t677;
t17 = pkin(10) * t105 + t372 * t44 + t375 * t43 + t708;
t12 = -pkin(3) * t128 + pkin(10) * t92 + t26 * t372 + t375 * t81;
t7 = t369 * t31 + t366 * t475;
t6 = pkin(10) * t32 + (-pkin(11) * t372 + t502) * t39;
t3 = t369 * t23 + t366 * t476;
t1 = -pkin(3) * t34 + pkin(10) * t24 + t13 * t375 + t16 * t372;
t59 = [0, 0, 0, 0, 0, qJDD(1), t429, t430, 0, 0, t361 * t515, (t343 * t368 + t346 * t365) * t367 + (t361 - t363) * t362 * t517, t347 * t538 + t370 * t344 + t367 * (t364 - t541) * t518, t363 * t515, t348 * t531 + t370 * t345 + t365 * (-t364 + t540) * t519, t364 * qJDD(1), (-t319 + pkin(1) * (t347 * t368 + t351 * t365)) * t370 + (t368 * t336 + pkin(1) * t346 + qJ(2) * (-t347 * t365 + t351 * t368)) * t367, (-t320 + pkin(1) * (-t348 * t365 + t350 * t368)) * t370 + (-t365 * t336 - pkin(1) * t343 + qJ(2) * (-t348 * t368 - t350 * t365)) * t367, pkin(1) * (-t344 * t368 + t345 * t365) * t370 + (-pkin(1) * t349 + qJ(2) * (t344 * t365 + t345 * t368) + t442) * t367, pkin(1) * (t336 * t367 + (-t319 * t368 + t320 * t365) * t370) + t442 * t532, t370 * (t313 * t534 + (t333 * t369 + t342 * t533) * t335) + (t365 * (t313 * t376 - t342 * t542) + t368 * (t313 * t528 + (-t333 * t366 + t342 * t527) * t335)) * t367, t370 * (t369 * t314 + t366 * t450) + (t365 * (-t285 * t376 - t373 * t289) + t368 * (-t366 * t314 + t369 * t450)) * t367, t370 * (t369 * t290 + t366 * t444) + (t365 * (-t373 * t327 + t627) + t368 * (-t366 * t290 + t369 * t444)) * t367, (t376 * t328 - t373 * t404) * t538 + (t369 * t395 + (t335 * t366 + t342 * t528) * t333) * t531 + t370 * (t366 * t395 + (-t335 * t369 + t342 * t534) * t333), t370 * (t366 * t445 + t369 * t394) + (t365 * (t326 * t376 + t525) + t368 * (-t366 * t394 + t369 * t445)) * t367, t619 * t341 + ((t365 * (-t333 * t376 + t542) + t441 * t530) * t367 + t441 * t535) * t342, (pkin(2) * t224 - t369 * t203 + pkin(1) * (t224 * t368 + t261 * t365) + t478 * t366) * t370 + (t365 * (-t555 + (-t223 * t366 - t224 * t369) * pkin(9)) + t368 * (-pkin(2) * t223 + t366 * t203 + t369 * t478) - pkin(1) * t223 + qJ(2) * (-t224 * t365 + t261 * t368)) * t367, (pkin(2) * t229 - t204 * t369 + pkin(1) * (t229 * t368 + t263 * t365) + t477 * t366) * t370 + (t365 * (-t554 + (-t228 * t366 - t229 * t369) * pkin(9)) + t368 * (-pkin(2) * t228 + t204 * t366 + t369 * t477) - pkin(1) * t228 + qJ(2) * (-t229 * t365 + t263 * t368)) * t367, (pkin(2) * t206 + pkin(1) * (t206 * t368 + t246 * t365) + t414 * t366) * t370 + (t365 * t455 + qJ(2) * (t246 * t368 - t560) + t433 * t205 + (-pkin(9) * t560 + t368 * t414) * t369) * t367, (t143 * t589 + pkin(2) * t136 + pkin(1) * (t136 * t368 + t143 * t365)) * t370 + (qJ(2) * (-t136 * t365 + t143 * t368) + (-pkin(1) - t593) * t135 + (t365 * (-t135 * t366 - t136 * t369) + t143 * t530) * pkin(9)) * t367, t370 * (t369 * t238 + t366 * t436) + (t365 * (t239 * t376 + t503) + t368 * (-t366 * t238 + t369 * t436)) * t367, t370 * (t369 * t176 + t366 * t459) + (t365 * (t178 * t376 + t373 * t291) + t368 * (-t366 * t176 + t369 * t459)) * t367, t370 * (t369 * t230 + t366 * t452) + (t365 * (t232 * t376 + t373 * t244) + t368 * (-t366 * t230 + t369 * t452)) * t367, t370 * (t369 * t236 + t366 * t437) + (t365 * (t237 * t376 - t503) + t368 * (-t366 * t236 + t369 * t437)) * t367, t370 * (t369 * t231 + t366 * t451) + (t365 * (t233 * t376 - t373 * t240) + t368 * (-t366 * t231 + t369 * t451)) * t367, (t376 * t257 + t373 * t397) * t538 + (-t366 * t256 + t369 * t385) * t531 + t370 * (t369 * t256 + t366 * t385), (pkin(2) * t138 + t369 * t107 + pkin(1) * (t138 * t368 + t157 * t365) + t416 * t366) * t370 + (t365 * (-t137 * t589 + t141 * t376 - t373 * t88) + t368 * (-pkin(2) * t137 - t366 * t107) - pkin(1) * t137 + qJ(2) * (t157 * t368 - t571) + (-pkin(9) * t571 + t368 * t416) * t369) * t367, (pkin(2) * t140 + t369 * t108 + pkin(1) * (t140 * t368 + t162 * t365) + t415 * t366) * t370 + (t365 * (-t139 * t589 + t142 * t376 - t373 * t89) + t368 * (-pkin(2) * t139 - t366 * t108) - pkin(1) * t139 + qJ(2) * (t162 * t368 - t570) + (-pkin(9) * t570 + t368 * t415) * t369) * t367, (pkin(2) * t98 + t369 * t69 + pkin(1) * (t152 * t365 + t368 * t98) + t410 * t366) * t370 + (t365 * (t373 * pkin(3) * t177 + t376 * t72 - t97 * t589) + t368 * (-pkin(2) * t97 - t366 * t69) - pkin(1) * t97 + qJ(2) * (t152 * t368 - t577) + (-pkin(9) * t577 + t368 * t410) * t369) * t367, (pkin(2) * t42 + t369 * t71 + pkin(1) * (t365 * t70 + t368 * t42) + (t598 + (-pkin(10) * t373 - t592) * t79) * t366) * t370 + (t365 * (pkin(3) * t575 - pkin(10) * t573 - t41 * t589) + t368 * (-pkin(2) * t41 - t366 * t71) - pkin(1) * t41 + qJ(2) * (t368 * t70 - t584) + (-pkin(9) * t584 + t368 * (-pkin(3) * t573 - pkin(10) * t575 + t598)) * t369) * t367, t604, -t728, t695, t657, -t730, t658, (t369 * t28 + t420 * t366 + t696) * t370 + (t365 * (t30 * t376 - t373 * t51 - t697) + t368 * (-t366 * t28 - t709) + (t368 * t420 - t715) * t369 + t711) * t367, (pkin(2) * t67 + t369 * t29 + pkin(1) * (t365 * t86 + t368 * t67) + t419 * t366) * t370 + (t365 * (t33 * t376 - t373 * t52 - t65 * t589) + t368 * (-pkin(2) * t65 - t366 * t29) - pkin(1) * t65 + qJ(2) * (t368 * t86 - t579) + (-pkin(9) * t579 + t368 * t419) * t369) * t367, (pkin(2) * t50 + t369 * t22 + pkin(1) * (t365 * t78 + t368 * t50) + t421 * t366) * t370 + (t365 * (t25 * t376 - t373 * t27 - t48 * t589) + t368 * (-pkin(2) * t48 - t366 * t22) - pkin(1) * t48 + qJ(2) * (t368 * t78 - t582) + (-pkin(9) * t582 + t368 * t421) * t369) * t367, (pkin(2) * t8 + t369 * t6 + pkin(1) * (t15 * t365 + t368 * t8) + t428 * t366) * t370 + (t365 * (-t373 * t11 + t376 * t9 - t7 * t589) + t368 * (-pkin(2) * t7 - t366 * t6) - pkin(1) * t7 + qJ(2) * (t15 * t368 - t586) + (-pkin(9) * t586 + t368 * t428) * t369) * t367, t604, t695, t728, t658, t730, t657, (t369 * t18 + t422 * t366 + t696) * t370 + (t365 * (t20 * t376 - t373 * t38 - t697) + t368 * (-t366 * t18 - t709) + (t368 * t422 - t715) * t369 + t711) * t367, (pkin(2) * t49 + t369 * t12 + pkin(1) * (t365 * t77 + t368 * t49) + t424 * t366) * t370 + (t365 * (t14 * t376 - t373 * t21 - t47 * t589) + t368 * (-pkin(2) * t47 - t366 * t12) - pkin(1) * t47 + qJ(2) * (t368 * t77 - t583) + (-pkin(9) * t583 + t368 * t424) * t369) * t367, (pkin(2) * t60 + t369 * t17 + pkin(1) * (t365 * t84 + t368 * t60) + t423 * t366) * t370 + (t365 * (t19 * t376 - t373 * t36 - t58 * t589) + t368 * (-pkin(2) * t58 - t366 * t17) - pkin(1) * t58 + qJ(2) * (t368 * t84 - t581) + (-pkin(9) * t581 + t368 * t423) * t369) * t367, (pkin(2) * t4 + t369 * t1 + pkin(1) * (t10 * t365 + t368 * t4) + t431 * t366) * t370 + (t365 * (t2 * t376 - t3 * t589 - t373 * t5) + t368 * (-pkin(2) * t3 - t366 * t1) - pkin(1) * t3 + qJ(2) * (t10 * t368 - t587) + (-pkin(9) * t587 + t368 * t431) * t369) * t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t346, t343, t349, -t336, 0, 0, 0, 0, 0, 0, t223, t228, t205, t135, 0, 0, 0, 0, 0, 0, t137, t139, t97, t41, 0, 0, 0, 0, 0, 0, t689, t65, t48, t7, 0, 0, 0, 0, 0, 0, t689, t47, t58, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t543, t314, t290, -t543, t394, t341, -t203, -t204, 0, 0, t238, t176, t230, t236, t231, t256, t107, t108, t69, t71, t486, -t99, t676, t609, -t115, t610, t28, t29, t22, t6, t486, t676, t99, t610, t115, t609, t18, t12, t17, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292, t291, t244, -t292, -t240, t397, -t125, -t126, 0, 0, t186, t129, t663, t183, -t165, -t626, t426, t425, t413, t488, t186, t663, -t129, -t626, t165, t183, t612, t405, t611, t613; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t547, t622, t616, -t547, -t189, t276, -t75, -t76, 0, 0, t547, t616, -t622, t276, t189, -t547, t387, t479, t315 + t614, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t618, t616, t623, t64;];
tauJ_reg  = t59;
