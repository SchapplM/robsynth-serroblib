% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:35:56
% EndTime: 2019-03-10 03:36:30
% DurationCPUTime: 19.71s
% Computational Cost: add. (38751->848), mult. (90405->1066), div. (0->0), fcn. (69108->18), ass. (0->417)
t397 = sin(qJ(2));
t401 = cos(qJ(2));
t634 = sin(qJ(3));
t636 = cos(qJ(3));
t443 = t634 * t397 - t636 * t401;
t291 = t443 * qJD(1);
t310 = t397 * t636 + t401 * t634;
t292 = t310 * qJD(1);
t396 = sin(qJ(4));
t635 = cos(qJ(4));
t232 = t635 * t291 + t292 * t396;
t399 = cos(qJ(6));
t400 = cos(qJ(5));
t533 = qJD(5) + qJD(6);
t538 = qJD(6) * t399;
t540 = qJD(5) * t400;
t394 = sin(qJ(6));
t395 = sin(qJ(5));
t564 = t394 * t395;
t557 = t399 * t400;
t646 = t557 - t564;
t686 = -t646 * t232 - t399 * t540 - t400 * t538 + t533 * t564;
t563 = t394 * t400;
t309 = t395 * t399 + t563;
t253 = t533 * t309;
t685 = t309 * t232 + t253;
t667 = t232 * t400;
t684 = t540 + t667;
t509 = t636 * qJD(3);
t441 = t636 * qJD(2) + t509;
t388 = qJD(2) + qJD(3);
t379 = qJD(4) + t388;
t453 = -t396 * t291 + t292 * t635;
t213 = -t400 * t379 + t395 * t453;
t404 = -pkin(8) - pkin(7);
t344 = t404 * t401;
t320 = qJD(1) * t344;
t293 = t634 * t320;
t342 = t404 * t397;
t318 = qJD(1) * t342;
t612 = qJD(2) * pkin(2);
t300 = t318 + t612;
t241 = t636 * t300 + t293;
t285 = t292 * pkin(9);
t210 = -t285 + t241;
t198 = t388 * pkin(3) + t210;
t297 = t636 * t320;
t242 = t300 * t634 - t297;
t619 = t291 * pkin(9);
t211 = t242 - t619;
t205 = t635 * t211;
t137 = t396 * t198 + t205;
t134 = t379 * pkin(10) + t137;
t386 = t401 * pkin(2);
t374 = t386 + pkin(1);
t340 = t374 * qJD(1);
t257 = pkin(3) * t291 - t340;
t152 = pkin(4) * t232 - pkin(10) * t453 + t257;
t82 = t134 * t400 + t152 * t395;
t67 = -pkin(11) * t213 + t82;
t609 = t394 * t67;
t215 = t379 * t395 + t400 * t453;
t81 = -t134 * t395 + t400 * t152;
t66 = -pkin(11) * t215 + t81;
t662 = qJD(5) + t232;
t62 = pkin(5) * t662 + t66;
t26 = t399 * t62 - t609;
t607 = t399 * t67;
t27 = t394 * t62 + t607;
t387 = qJDD(2) + qJDD(3);
t378 = qJDD(4) + t387;
t536 = qJD(1) * qJD(2);
t505 = t401 * t536;
t535 = t397 * qJDD(1);
t255 = qJDD(2) * pkin(2) - t404 * (-t505 - t535);
t506 = t397 * t536;
t534 = t401 * qJDD(1);
t256 = t404 * (-t506 + t534);
t161 = -qJD(3) * t242 + t636 * t255 + t256 * t634;
t507 = t634 * qJD(3);
t440 = t634 * qJD(2) + t507;
t499 = t636 * qJDD(1);
t543 = qJD(1) * t397;
t498 = t634 * qJDD(1);
t661 = -qJD(1) * t441 - t498;
t212 = -t397 * t499 + t401 * t661 + t440 * t543;
t111 = t387 * pkin(3) + t212 * pkin(9) + t161;
t160 = t634 * t255 - t636 * t256 + t300 * t509 + t320 * t507;
t664 = t388 * t310;
t413 = t664 * qJD(1);
t436 = t443 * qJDD(1);
t410 = -t436 - t413;
t122 = pkin(9) * t410 + t160;
t511 = qJD(4) * t635;
t542 = qJD(4) * t396;
t489 = -t396 * t111 - t635 * t122 - t198 * t511 + t211 * t542;
t48 = pkin(10) * t378 - t489;
t541 = qJD(5) * t395;
t124 = qJD(4) * t453 - t396 * t212 - t635 * t410;
t362 = pkin(2) * t506;
t647 = -t212 * t635 - t291 * t511 - t292 * t542 - t396 * t413;
t408 = -t396 * t436 + t647;
t596 = qJDD(1) * pkin(1);
t57 = -pkin(2) * t534 - pkin(3) * t410 + t124 * pkin(4) - pkin(10) * t408 + t362 - t596;
t12 = -t134 * t541 + t152 * t540 + t395 * t57 + t400 * t48;
t525 = t379 * t541 + t395 * t408 + t453 * t540;
t458 = t378 * t400 - t525;
t10 = pkin(11) * t458 + t12;
t539 = qJD(6) * t394;
t123 = qJDD(5) + t124;
t56 = t400 * t57;
t13 = -qJD(5) * t82 - t395 * t48 + t56;
t91 = -t395 * t378 - t379 * t540 - t400 * t408 + t453 * t541;
t9 = pkin(5) * t123 + pkin(11) * t91 + t13;
t3 = (qJD(6) * t62 + t10) * t399 + t394 * t9 - t67 * t539;
t4 = -qJD(6) * t27 - t10 * t394 + t399 * t9;
t393 = qJ(2) + qJ(3);
t384 = qJ(4) + t393;
t368 = sin(t384);
t361 = g(3) * t368;
t369 = cos(t384);
t402 = cos(qJ(1));
t573 = t369 * t402;
t398 = sin(qJ(1));
t574 = t369 * t398;
t523 = g(1) * t573 + g(2) * t574 + t361;
t425 = t26 * t686 - t685 * t27 + t3 * t646 - t4 * t309 - t523;
t204 = t396 * t211;
t136 = t198 * t635 - t204;
t133 = -t379 * pkin(4) - t136;
t106 = t213 * pkin(5) + t133;
t488 = -t635 * t111 + t396 * t122 + t198 * t542 + t211 * t511;
t618 = t378 * pkin(4);
t49 = t488 - t618;
t23 = -pkin(5) * t458 + t49;
t392 = qJ(5) + qJ(6);
t380 = sin(t392);
t475 = g(1) * t402 + g(2) * t398;
t457 = t475 * t368;
t624 = g(3) * t380;
t418 = -t106 * t686 + t23 * t309 + t27 * t453 + t369 * t624 - t380 * t457;
t144 = t399 * t213 + t215 * t394;
t38 = t213 * t538 + t215 * t539 - t394 * t458 + t399 * t91;
t463 = t213 * t394 - t399 * t215;
t427 = qJD(6) * t463 + t394 * t91 + t399 * t458;
t7 = t144 * t686 + t309 * t427 - t38 * t646 + t685 * t463;
t678 = t662 * t395;
t85 = t395 * t458;
t18 = -t213 * t684 - t215 * t678 - t91 * t400 + t85;
t120 = qJDD(6) + t123;
t224 = qJD(6) + t662;
t42 = t309 * t120 - t224 * t686 + t453 * t463;
t15 = -t38 * t309 + t463 * t686;
t117 = t400 * t123;
t58 = t213 * t453 - t662 * t678 + t117;
t668 = t232 * t395;
t532 = pkin(11) * t668;
t490 = pkin(3) * t511;
t141 = t210 * t635 - t204;
t180 = pkin(4) * t453 + pkin(10) * t232;
t633 = pkin(3) * t292;
t162 = t180 + t633;
t92 = -t141 * t395 + t400 * t162;
t683 = -t395 * t490 - t92;
t373 = pkin(2) * t636 + pkin(3);
t516 = t634 * t396;
t237 = t373 * t511 + (-qJD(4) * t516 + (t635 * t636 - t516) * qJD(3)) * pkin(2);
t251 = t636 * t318 + t293;
t217 = -t285 + t251;
t250 = -t318 * t634 + t297;
t446 = t250 + t619;
t154 = t217 * t635 + t396 * t446;
t375 = pkin(2) * t543;
t155 = t162 + t375;
t94 = -t154 * t395 + t400 * t155;
t682 = -t395 * t237 - t94;
t93 = t400 * t141 + t395 * t162;
t681 = -t400 * t490 + t93;
t95 = t400 * t154 + t395 * t155;
t680 = -t400 * t237 + t95;
t648 = (t541 + t668) * pkin(5);
t482 = pkin(5) * t453 + pkin(11) * t667;
t382 = cos(t392);
t623 = g(3) * t382;
t576 = t368 * t402;
t577 = t368 * t398;
t663 = g(1) * t576 + g(2) * t577;
t423 = t106 * t685 - t23 * t646 - t26 * t453 - t369 * t623 + t382 * t663;
t16 = t144 * t685 + t427 * t646;
t88 = t91 * t395;
t60 = t215 * t684 - t88;
t43 = t120 * t646 + t144 * t453 - t224 * t685;
t116 = t395 * t123;
t59 = -t215 * t453 + t662 * t684 + t116;
t486 = t634 * t635;
t288 = pkin(2) * t486 + t396 * t373;
t283 = pkin(10) + t288;
t616 = -pkin(11) - t283;
t501 = qJD(5) * t616;
t674 = t395 * t501 - t532 - t680;
t673 = t400 * t501 - t482 + t682;
t631 = pkin(3) * t396;
t370 = pkin(10) + t631;
t615 = -pkin(11) - t370;
t500 = qJD(5) * t615;
t672 = t395 * t500 - t532 - t681;
t671 = t400 * t500 - t482 + t683;
t101 = t400 * t136 + t395 * t180;
t403 = -pkin(11) - pkin(10);
t519 = qJD(5) * t403;
t670 = t395 * t519 - t101 - t532;
t100 = -t136 * t395 + t400 * t180;
t669 = t400 * t519 - t100 - t482;
t595 = t133 * t232;
t594 = t144 * t463;
t666 = t232 * t453;
t115 = -t232 ^ 2 + t453 ^ 2;
t660 = -t144 ^ 2 + t463 ^ 2;
t659 = t144 * t224 - t38;
t568 = t382 * t398;
t571 = t380 * t402;
t268 = -t369 * t568 + t571;
t567 = t382 * t402;
t572 = t380 * t398;
t270 = t369 * t567 + t572;
t658 = g(1) * t270 - g(2) * t268 + t106 * t144 + t368 * t623 - t3;
t267 = t369 * t572 + t567;
t269 = -t369 * t571 + t568;
t657 = -g(1) * t269 + g(2) * t267 + t106 * t463 + t368 * t624 + t4;
t656 = -t224 * t463 + t427;
t432 = t232 * t257 + t489 + t523;
t102 = t232 * t379 + t408;
t590 = t224 * t453;
t589 = t662 * t453;
t651 = t241 * t388;
t140 = t210 * t396 + t205;
t478 = pkin(3) * t542 - t140;
t261 = t636 * t342 + t344 * t634;
t227 = -t310 * pkin(9) + t261;
t263 = t634 * t342 - t636 * t344;
t228 = -pkin(9) * t443 + t263;
t177 = t396 * t227 + t228 * t635;
t170 = t400 * t177;
t429 = t635 * t443;
t248 = t310 * t396 + t429;
t437 = t396 * t443;
t249 = t310 * t635 - t437;
t273 = pkin(3) * t443 - t374;
t175 = t248 * pkin(4) - t249 * pkin(10) + t273;
t105 = t395 * t175 + t170;
t550 = t237 - t154;
t549 = -t217 * t396 + t635 * t446 + t373 * t542 + (qJD(4) * t486 + (t396 * t636 + t486) * qJD(3)) * pkin(2);
t617 = t400 * pkin(5);
t371 = pkin(4) + t617;
t650 = -t368 * t403 + t369 * t371;
t649 = t369 * pkin(4) + t368 * pkin(10);
t645 = -t395 * t81 + t400 * t82;
t555 = t400 * t402;
t560 = t395 * t398;
t277 = t369 * t560 + t555;
t558 = t398 * t400;
t559 = t395 * t402;
t279 = -t369 * t559 + t558;
t644 = -g(1) * t279 + g(2) * t277;
t460 = t133 * t541 + t400 * t663 - t81 * t453;
t621 = g(3) * t395;
t641 = t133 * t540 + t369 * t621 + t49 * t395 + t453 * t82;
t625 = g(3) * t369;
t426 = -t257 * t453 - t488 - t625 + t663;
t103 = t379 * t453 - t124;
t381 = sin(t393);
t632 = pkin(3) * t381;
t383 = cos(t393);
t366 = pkin(3) * t383;
t546 = t366 + t386;
t317 = pkin(1) + t546;
t299 = t402 * t317;
t627 = g(2) * t299;
t622 = g(3) * t383;
t620 = g(3) * t401;
t258 = t616 * t395;
t385 = t400 * pkin(11);
t259 = t283 * t400 + t385;
t196 = t258 * t399 - t259 * t394;
t614 = qJD(6) * t196 + t394 * t673 + t399 * t674;
t197 = t258 * t394 + t259 * t399;
t613 = -qJD(6) * t197 - t394 * t674 + t399 * t673;
t11 = t12 * t400;
t611 = t662 * t81;
t610 = t662 * t82;
t305 = t615 * t395;
t306 = t370 * t400 + t385;
t239 = t305 * t399 - t306 * t394;
t603 = qJD(6) * t239 + t394 * t671 + t399 * t672;
t240 = t305 * t394 + t306 * t399;
t602 = -qJD(6) * t240 - t394 * t672 + t399 * t671;
t341 = t403 * t395;
t343 = pkin(10) * t400 + t385;
t260 = t341 * t399 - t343 * t394;
t601 = qJD(6) * t260 + t394 * t669 + t399 * t670;
t262 = t341 * t394 + t343 * t399;
t600 = -qJD(6) * t262 - t394 * t670 + t399 * t669;
t598 = pkin(7) * qJDD(1);
t597 = qJD(5) * t81;
t414 = t388 * t443;
t156 = qJD(4) * t429 + t310 * t542 + t396 * t664 + t414 * t635;
t593 = t156 * t400;
t592 = t213 * t395;
t591 = t215 * t213;
t580 = t249 * t395;
t579 = t249 * t400;
t578 = t292 * t291;
t570 = t381 * t398;
t569 = t381 * t402;
t566 = t383 * t398;
t565 = t383 * t402;
t562 = t395 * t156;
t551 = t648 + t549;
t548 = t648 + t478;
t390 = t397 ^ 2;
t391 = t401 ^ 2;
t545 = t390 - t391;
t544 = t390 + t391;
t530 = t635 * pkin(3);
t528 = pkin(10) * qJD(5) * t662;
t377 = t397 * t612;
t406 = qJD(1) ^ 2;
t526 = t397 * t406 * t401;
t522 = t366 + t649;
t520 = qJD(2) * t404;
t518 = t212 * t636;
t514 = t249 * t541;
t513 = -t49 - t625;
t389 = -pkin(9) + t404;
t504 = pkin(5) * t395 - t389;
t319 = t397 * t520;
t321 = t401 * t520;
t194 = t636 * t319 + t634 * t321 + t342 * t509 + t344 * t507;
t169 = -pkin(9) * t664 + t194;
t472 = -t319 * t634 + t636 * t321;
t409 = pkin(9) * t414 - t342 * t507 + t344 * t509 + t472;
t70 = t169 * t635 + t227 * t511 - t228 * t542 + t396 * t409;
t157 = -qJD(4) * t437 + t310 * t511 - t396 * t414 + t635 * t664;
t236 = pkin(3) * t664 + t377;
t79 = t157 * pkin(4) + t156 * pkin(10) + t236;
t502 = -t395 * t70 + t400 * t79;
t104 = t400 * t175 - t177 * t395;
t176 = -t635 * t227 + t228 * t396;
t491 = t11 - t523;
t372 = -t530 - pkin(4);
t485 = t397 * t505;
t484 = -pkin(4) * t577 + pkin(10) * t574;
t483 = -pkin(4) * t576 + pkin(10) * t573;
t481 = -g(1) * t577 + g(2) * t576;
t477 = -t137 + t648;
t474 = g(1) * t398 - g(2) * t402;
t473 = t366 + t650;
t287 = -pkin(2) * t516 + t373 * t635;
t471 = -pkin(10) * t123 + t595;
t76 = pkin(5) * t248 - pkin(11) * t579 + t104;
t83 = -pkin(11) * t580 + t105;
t44 = -t394 * t83 + t399 * t76;
t45 = t394 * t76 + t399 * t83;
t470 = t395 * t82 + t400 * t81;
t466 = -t123 * t283 + t595;
t465 = -t123 * t370 + t595;
t464 = -t136 * t232 + t137 * t453;
t461 = t368 * t371 + t369 * t403;
t282 = -pkin(4) - t287;
t459 = -t667 * t81 - t668 * t82 + t491;
t456 = t475 * t381;
t454 = -0.2e1 * pkin(1) * t536 - pkin(7) * qJDD(2);
t452 = t249 * t540 - t562;
t451 = -t514 - t593;
t19 = t175 * t540 - t177 * t541 + t395 * t79 + t400 * t70;
t450 = t461 * t398;
t449 = t461 * t402;
t90 = t458 * t400;
t439 = qJD(5) * t215 + t458;
t405 = qJD(2) ^ 2;
t435 = -pkin(7) * t405 + t474 + 0.2e1 * t596;
t434 = pkin(1) * t406 + t475 - t598;
t431 = g(1) * t565 + g(2) * t566 + g(3) * t381 - t340 * t291 - t160;
t424 = -qJD(5) * t470 - t13 * t395 + t11;
t71 = qJD(4) * t177 + t396 * t169 - t635 * t409;
t421 = -g(1) * t483 - g(2) * t484;
t420 = -t395 * t457 + t641;
t412 = g(1) * t569 + g(2) * t570 + t340 * t292 + t161 - t622;
t329 = t372 - t617;
t328 = -pkin(2) * t397 - t632;
t304 = t402 * t328;
t303 = t398 * t328;
t286 = -qJDD(1) * t374 + t362;
t280 = t369 * t555 + t560;
t278 = -t369 * t558 + t559;
t271 = t282 - t617;
t266 = t375 + t633;
t216 = -t291 ^ 2 + t292 ^ 2;
t195 = -qJD(3) * t263 + t472;
t187 = pkin(3) * t413 + qJDD(1) * t273 + t362;
t186 = t292 * t388 + t410;
t185 = t291 * t388 - t212;
t182 = t646 * t249;
t181 = t309 * t249;
t132 = pkin(5) * t580 + t176;
t61 = t592 * t662 + t90;
t54 = -t156 * t563 - t394 * t514 - t539 * t580 + (t533 * t579 - t562) * t399;
t53 = t156 * t557 + t249 * t253 - t394 * t562;
t52 = pkin(5) * t452 + t71;
t29 = t399 * t66 - t609;
t28 = -t394 * t66 - t607;
t20 = -qJD(5) * t105 + t502;
t17 = -pkin(11) * t452 + t19;
t14 = pkin(11) * t593 + pkin(5) * t157 + (-t170 + (pkin(11) * t249 - t175) * t395) * qJD(5) + t502;
t6 = -qJD(6) * t45 + t14 * t399 - t17 * t394;
t5 = qJD(6) * t44 + t14 * t394 + t17 * t399;
t1 = [0, 0, 0, 0, 0, qJDD(1), t474, t475, 0, 0, qJDD(1) * t390 + 0.2e1 * t485, 0.2e1 * t397 * t534 - 0.2e1 * t536 * t545, qJDD(2) * t397 + t401 * t405, qJDD(1) * t391 - 0.2e1 * t485, qJDD(2) * t401 - t397 * t405, 0, t397 * t454 + t401 * t435, -t397 * t435 + t401 * t454, 0.2e1 * t544 * t598 - t475, -g(1) * (-pkin(1) * t398 + pkin(7) * t402) - g(2) * (pkin(1) * t402 + pkin(7) * t398) + (pkin(7) ^ 2 * t544 + pkin(1) ^ 2) * qJDD(1), -t212 * t310 - t292 * t414 (-t441 * t291 + t310 * (-qJD(1) * t440 + t499) - t518 - t292 * t440) * t401 + (t212 * t634 + t291 * t440 - t292 * t441 + t310 * t661) * t397, t310 * t387 - t388 * t414, t291 * t664 - t410 * t443, -t387 * t443 - t388 * t664, 0, g(1) * t566 - g(2) * t565 + t195 * t388 + t261 * t387 + t286 * t443 + t291 * t377 - t340 * t664 + t374 * t410, -g(1) * t570 + g(2) * t569 - t194 * t388 + t374 * t212 - t263 * t387 + t286 * t310 + t292 * t377 + t340 * t414, -t161 * t310 - t194 * t291 - t195 * t292 + t261 * t212 - t242 * t664 + t263 * t410 - t475 + (-t160 + t651) * t443, t160 * t263 + t242 * t194 + t161 * t261 + t241 * t195 - t286 * t374 - t340 * t377 - g(1) * (-t374 * t398 - t402 * t404) - g(2) * (t374 * t402 - t398 * t404) -t453 * t156 + (t396 * (-t397 * t498 + t401 * t499) + t647) * t249, -t249 * t124 + t156 * t232 - t157 * t453 - t248 * t408, -t156 * t379 + t249 * t378, t124 * t248 + t157 * t232, -t157 * t379 - t248 * t378, 0, t124 * t273 + t157 * t257 - t176 * t378 + t187 * t248 + t232 * t236 + t369 * t474 - t379 * t71, -t257 * t156 - t177 * t378 + t187 * t249 + t236 * t453 + t273 * t408 - t70 * t379 + t481, -t177 * t124 + t136 * t156 - t137 * t157 + t176 * t408 - t70 * t232 + t248 * t489 + t249 * t488 + t453 * t71 - t475, -t489 * t177 + t137 * t70 + t488 * t176 - t136 * t71 + t187 * t273 + t257 * t236 - g(1) * (-t317 * t398 - t389 * t402) - g(2) * (-t389 * t398 + t299) t215 * t451 - t579 * t91 (t213 * t400 + t215 * t395) * t156 + (t90 + t88 + (-t215 * t400 + t592) * qJD(5)) * t249, t117 * t249 + t157 * t215 - t248 * t91 + t451 * t662, t213 * t452 - t249 * t85, -t116 * t249 - t213 * t157 + t248 * t458 - t452 * t662, t123 * t248 + t157 * t662, -g(1) * t278 - g(2) * t280 + t104 * t123 + t13 * t248 + t133 * t452 + t81 * t157 - t176 * t458 + t20 * t662 + t71 * t213 + t49 * t580, -g(1) * t277 - g(2) * t279 - t105 * t123 - t12 * t248 + t133 * t451 - t157 * t82 - t176 * t91 - t19 * t662 + t215 * t71 + t49 * t579, -t19 * t213 + t105 * t458 - t20 * t215 + t104 * t91 + t470 * t156 + (-qJD(5) * t645 - t12 * t395 - t13 * t400) * t249 - t481, -t627 + t13 * t104 + t12 * t105 + t133 * t71 + t49 * t176 + t82 * t19 + t81 * t20 + (g(1) * t389 - g(2) * t649) * t402 + (-g(1) * (-t317 - t649) + g(2) * t389) * t398, -t182 * t38 + t463 * t53, t144 * t53 + t181 * t38 + t182 * t427 + t463 * t54, t120 * t182 - t157 * t463 - t224 * t53 - t248 * t38, t144 * t54 - t181 * t427, -t120 * t181 - t144 * t157 - t224 * t54 + t248 * t427, t120 * t248 + t157 * t224, -g(1) * t268 - g(2) * t270 + t106 * t54 + t120 * t44 - t132 * t427 + t144 * t52 + t157 * t26 + t181 * t23 + t224 * t6 + t248 * t4, -g(1) * t267 - g(2) * t269 - t106 * t53 - t120 * t45 - t132 * t38 - t157 * t27 + t182 * t23 - t224 * t5 - t248 * t3 - t463 * t52, -t144 * t5 - t181 * t3 - t182 * t4 + t26 * t53 - t27 * t54 + t38 * t44 + t427 * t45 + t463 * t6 - t481, -t627 + t106 * t52 + t23 * t132 + t26 * t6 + t27 * t5 + t3 * t45 + t4 * t44 + (-g(1) * t504 - g(2) * t650) * t402 + (-g(1) * (-t317 - t650) - g(2) * t504) * t398; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t526, t545 * t406, t535, t526, t534, qJDD(2), t397 * t434 - t620, g(3) * t397 + t401 * t434, 0, 0, t578, t216, t185, -t578, t186, t387, -t250 * t388 + (-t291 * t543 + t387 * t636 - t388 * t507) * pkin(2) + t412, t251 * t388 + (-t292 * t543 - t387 * t634 - t388 * t509) * pkin(2) + t431 (t242 + t250) * t292 - (t241 - t251) * t291 + (-t291 * t509 + t292 * t507 + t410 * t634 + t518) * pkin(2), -t241 * t250 - t242 * t251 + (t636 * t161 + t160 * t634 - t620 + (-t241 * t634 + t242 * t636) * qJD(3) + (qJD(1) * t340 + t475) * t397) * pkin(2), t666, t115, t102, -t666, t103, t378, -t232 * t266 + t287 * t378 - t379 * t549 + t426, -t266 * t453 - t288 * t378 - t379 * t550 + t432, -t288 * t124 - t232 * t550 - t287 * t408 + t453 * t549 + t464, -g(3) * t546 - t136 * t549 + t137 * t550 - t257 * t266 - t287 * t488 - t288 * t489 - t328 * t475, t60, t18, t59, t61, t58, -t589, t282 * t525 + (-t282 * t378 + t513) * t400 + t466 * t395 + t549 * t213 + (-t283 * t540 + t682) * t662 + t460, -t282 * t91 + t466 * t400 + t549 * t215 + (t283 * t541 + t680) * t662 + t420, t95 * t213 + t94 * t215 + (-t237 * t213 + t283 * t439 - t597) * t400 + (t237 * t215 - t283 * t91 - t13 + (t213 * t283 - t82) * qJD(5)) * t395 + t459, t49 * t282 - t82 * t95 - t81 * t94 - g(1) * (t304 + t483) - g(2) * (t303 + t484) - g(3) * (t386 + t522) + t645 * t237 + t549 * t133 + t424 * t283, t15, t7, t42, t16, t43, -t590, t120 * t196 + t144 * t551 + t224 * t613 - t271 * t427 + t423, -t120 * t197 - t224 * t614 - t271 * t38 - t463 * t551 + t418, -t144 * t614 + t196 * t38 + t197 * t427 + t463 * t613 + t425, t3 * t197 + t4 * t196 + t23 * t271 - g(1) * (t304 - t449) - g(2) * (t303 - t450) - g(3) * (t386 + t473) + t614 * t27 + t613 * t26 + t551 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t578, t216, t185, -t578, t186, t387, t242 * t388 + t412, t431 + t651, 0, 0, t666, t115, t102, -t666, t103, t378, t140 * t379 + (-t232 * t292 + t378 * t635 - t379 * t542) * pkin(3) + t426, t141 * t379 + (-t292 * t453 - t378 * t396 - t379 * t511) * pkin(3) + t432, -t124 * t631 - t408 * t530 + t464 + t478 * t453 + (t141 - t490) * t232, t136 * t140 - t137 * t141 + (-t635 * t488 - t622 - t257 * t292 - t396 * t489 + t456 + (-t136 * t396 + t137 * t635) * qJD(4)) * pkin(3), t60, t18, t59, t61, t58, -t589, t372 * t525 + (-t372 * t378 + t513) * t400 + t465 * t395 + t478 * t213 + (-t370 * t540 + t683) * t662 + t460, -t372 * t91 + t465 * t400 + t478 * t215 + (t370 * t541 + t681) * t662 + t420, t93 * t213 + t92 * t215 + (-t213 * t490 + t370 * t439 - t597) * t400 + (t215 * t490 - t370 * t91 - t13 + (t213 * t370 - t82) * qJD(5)) * t395 + t459, t49 * t372 - t82 * t93 - t81 * t92 - t133 * t140 - g(3) * t522 + (t456 + (t133 * t396 + t635 * t645) * qJD(4)) * pkin(3) + t424 * t370 + t421, t15, t7, t42, t16, t43, -t590, t120 * t239 + t144 * t548 + t224 * t602 - t329 * t427 + t423, -t120 * t240 - t224 * t603 - t329 * t38 - t463 * t548 + t418, -t144 * t603 + t239 * t38 + t240 * t427 + t463 * t602 + t425, -g(3) * t473 + t106 * t548 + t23 * t329 + t4 * t239 + t3 * t240 + t26 * t602 + t27 * t603 + t475 * (t461 + t632); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t666, t115, t102, -t666, t103, t378, t137 * t379 + t426, t136 * t379 + t432, 0, 0, t60, t18, t59, t213 * t678 + t90, t58, -t589, -pkin(4) * t525 - t100 * t662 - t137 * t213 + t471 * t395 + (t513 - t528 + t618) * t400 + t460, pkin(4) * t91 + t101 * t662 - t137 * t215 + t471 * t400 + (-t457 + t528) * t395 + t641, t100 * t215 + t101 * t213 + (pkin(10) * t439 - t611) * t400 + (-t13 - t610 + (qJD(5) * t213 - t91) * pkin(10)) * t395 + t491, -t49 * pkin(4) + pkin(10) * t424 - g(3) * t649 - t81 * t100 - t82 * t101 - t133 * t137 + t421, t15, t7, t42, t16, t43, -t590, t120 * t260 + t144 * t477 + t224 * t600 + t371 * t427 + t423, -t120 * t262 - t224 * t601 + t371 * t38 - t463 * t477 + t418, -t144 * t601 + t260 * t38 + t262 * t427 + t463 * t600 + t425, g(1) * t449 + g(2) * t450 - g(3) * t650 + t106 * t477 - t23 * t371 + t26 * t600 + t4 * t260 + t3 * t262 + t27 * t601; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t591, -t213 ^ 2 + t215 ^ 2, t213 * t662 - t91, -t591, t215 * t662 + t458, t123, -t134 * t540 - t133 * t215 + t610 + t56 + (-qJD(5) * t152 + t361 - t48) * t395 + t644, g(1) * t280 - g(2) * t278 + t133 * t213 + t361 * t400 - t12 + t611, 0, 0, -t594, t660, t659, t594, t656, t120, -t224 * t28 + (t120 * t399 - t144 * t215 - t224 * t539) * pkin(5) + t657, t224 * t29 + (-t120 * t394 + t215 * t463 - t224 * t538) * pkin(5) + t658, -t463 * t27 + t144 * t29 - t144 * t26 - t463 * t28 + (t38 * t399 + t427 * t394 + (-t144 * t399 - t394 * t463) * qJD(6)) * pkin(5), -t26 * t28 - t27 * t29 + (t3 * t394 + t4 * t399 - t106 * t215 + t368 * t621 + (-t26 * t394 + t27 * t399) * qJD(6) + t644) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t594, t660, t659, t594, t656, t120, t224 * t27 + t657, t224 * t26 + t658, 0, 0;];
tau_reg  = t1;
