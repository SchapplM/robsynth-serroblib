% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x28]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPR13_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:36
% EndTime: 2019-12-31 21:47:03
% DurationCPUTime: 12.18s
% Computational Cost: add. (6218->658), mult. (15832->923), div. (0->0), fcn. (16263->8), ass. (0->501)
t352 = sin(qJ(3));
t355 = cos(qJ(3));
t350 = sin(pkin(5));
t353 = sin(qJ(2));
t591 = t350 * t353;
t625 = cos(pkin(5));
t275 = t352 * t625 + t355 * t591;
t646 = t275 / 0.2e1;
t340 = t625 * t355;
t273 = t352 * t591 - t340;
t655 = pkin(3) + pkin(9);
t670 = t273 * t655;
t600 = t275 * t355;
t467 = t600 / 0.2e1;
t579 = t352 * t273;
t165 = -t579 / 0.2e1 + t467;
t546 = qJD(1) * t275;
t184 = t273 * t546;
t665 = -t165 * qJD(2) + t184;
t508 = t355 * qJD(2);
t332 = t352 * t508;
t662 = t165 * qJD(1) + t332;
t580 = t352 * qJ(4);
t631 = t355 * pkin(3);
t423 = -t580 - t631;
t307 = -pkin(2) + t423;
t356 = cos(qJ(2));
t496 = pkin(1) * t625;
t281 = pkin(7) * t591 - t356 * t496;
t255 = -pkin(2) * t625 + t281;
t601 = t275 * qJ(4);
t361 = t255 - t601;
t632 = t273 * pkin(3);
t120 = t361 + t632;
t622 = t120 * t352;
t664 = t307 * t646 + t622 / 0.2e1;
t669 = t355 * t655 + t580;
t354 = cos(qJ(5));
t566 = t354 * t355;
t456 = -t566 / 0.2e1;
t351 = sin(qJ(5));
t657 = t275 ^ 2;
t668 = t351 * t657;
t667 = t354 * t657;
t157 = t165 * qJD(3);
t541 = qJD(2) * t356;
t488 = t350 * t541;
t312 = t355 * t488;
t666 = t352 * t312 + t157;
t463 = t591 / 0.2e1;
t269 = t352 * t463 - t340 / 0.2e1;
t663 = -t269 * qJD(5) - t665;
t635 = t355 / 0.2e1;
t661 = qJD(5) * t635 + t662;
t345 = t351 ^ 2;
t347 = t354 ^ 2;
t327 = t345 - t347;
t447 = t327 * qJD(3);
t536 = qJD(4) * t355;
t660 = qJD(3) * t669 - t536;
t590 = t350 * t356;
t283 = pkin(7) * t590 + t353 * t496;
t256 = pkin(8) * t625 + t283;
t633 = pkin(8) * t353;
t428 = -pkin(2) * t356 - t633;
t257 = (-pkin(1) + t428) * t350;
t151 = t352 * t256 - t355 * t257;
t152 = t355 * t256 + t352 * t257;
t123 = pkin(3) * t590 + t151;
t617 = t123 * t355;
t326 = qJ(4) * t590;
t122 = t152 - t326;
t618 = t122 * t352;
t636 = -t355 / 0.2e1;
t637 = t352 / 0.2e1;
t659 = t151 * t636 + t152 * t637 - t618 / 0.2e1 + t617 / 0.2e1;
t399 = t275 * pkin(4) + t151;
t461 = t590 / 0.2e1;
t103 = -t273 * pkin(4) + t152;
t89 = t103 - t326;
t656 = t89 / 0.2e1;
t498 = t656 - t103 / 0.2e1;
t658 = -t352 * t498 - t355 * (-t655 * t461 + t399 / 0.2e1);
t654 = pkin(4) + pkin(8);
t653 = -qJ(4) / 0.2e1;
t135 = t351 * t590;
t602 = t273 * t354;
t203 = t135 + t602;
t652 = -t203 / 0.2e1;
t499 = t354 * t590;
t301 = t352 * t499;
t240 = t351 * t591 - t301;
t651 = t240 / 0.2e1;
t575 = t353 * t354;
t577 = t352 * t356;
t241 = (t351 * t577 + t575) * t350;
t650 = t241 / 0.2e1;
t265 = t352 * t281;
t649 = -t265 / 0.2e1;
t648 = -t273 / 0.2e1;
t647 = -t275 / 0.2e1;
t343 = pkin(3) * t352;
t624 = qJ(4) * t355;
t422 = t352 * pkin(9) - t624;
t287 = t343 + t422;
t645 = t287 / 0.2e1;
t302 = t355 * t135;
t644 = -t302 / 0.2e1;
t643 = -t307 / 0.2e1;
t318 = t654 * t352;
t642 = -t318 / 0.2e1;
t319 = t654 * t355;
t641 = t319 / 0.2e1;
t336 = pkin(3) * t591;
t640 = -t336 / 0.2e1;
t639 = -t351 / 0.2e1;
t638 = -t352 / 0.2e1;
t634 = t655 / 0.2e1;
t362 = t590 * t655 + t399;
t87 = t361 + t670;
t39 = t351 * t87 - t354 * t362;
t604 = t273 * qJ(4);
t124 = t275 * t655 + t604;
t587 = t351 * t124;
t626 = t103 - t89;
t5 = t399 * t203 + t39 * t273 + (t354 * t626 - t587) * t275;
t630 = t5 * qJD(1);
t603 = t273 * t351;
t204 = t499 - t603;
t40 = t351 * t362 + t354 * t87;
t572 = t354 * t124;
t6 = t399 * t204 + t40 * t273 + (-t351 * t626 - t572) * t275;
t629 = t6 * qJD(1);
t282 = (pkin(2) * t353 - pkin(8) * t356) * t350;
t266 = t352 * t282;
t267 = t355 * t281;
t548 = t266 - t267;
t576 = t353 * qJ(4);
t121 = (-pkin(4) * t577 + t576) * t350 + t548;
t564 = t355 * t282;
t449 = t265 + t564;
t156 = -t336 - t449;
t563 = t355 * t356;
t105 = (pkin(4) * t563 - pkin(9) * t353) * t350 + t156;
t573 = t354 * t105;
t504 = t350 * t577;
t497 = pkin(3) * t504 + t283;
t155 = t422 * t590 + t497;
t586 = t351 * t155;
t50 = t573 - t586;
t502 = t350 * t563;
t7 = -t121 * t203 + t89 * t240 + t50 * t275 - t39 * t502;
t628 = t7 * qJD(1);
t571 = t354 * t155;
t588 = t351 * t105;
t51 = t571 + t588;
t8 = -t121 * t204 + t89 * t241 - t51 * t275 - t40 * t502;
t627 = t8 * qJD(1);
t623 = t120 * t275;
t621 = t120 * t355;
t620 = t121 * t351;
t619 = t121 * t354;
t179 = t275 * pkin(3) + t604;
t15 = t120 * t179 - t122 * t151 + t123 * t152;
t616 = t15 * qJD(1);
t503 = t350 * t576;
t154 = -t503 - t548;
t182 = -t326 * t355 + t497;
t16 = t120 * t182 - t122 * t154 + t123 * t156;
t615 = t16 * qJD(1);
t452 = t151 / 0.2e1 - t123 / 0.2e1;
t453 = t122 / 0.2e1 - t152 / 0.2e1;
t462 = -t590 / 0.2e1;
t17 = (pkin(3) * t462 + t452) * t355 + (qJ(4) * t462 + t453) * t352;
t614 = t17 * qJD(1);
t19 = t89 * t203 + t39 * t275;
t613 = t19 * qJD(1);
t20 = -t89 * t204 - t40 * t275;
t612 = t20 * qJD(1);
t611 = t204 * t351;
t610 = t204 * t352;
t21 = (-t122 + t152) * t275 + (-t123 + t151) * t273;
t609 = t21 * qJD(1);
t608 = t241 * t354;
t607 = t255 * t355;
t26 = t154 * t273 + t156 * t275 + (t617 - t618) * t590;
t606 = t26 * qJD(1);
t27 = -t182 * t275 + (t122 * t353 + (t154 - t621) * t356) * t350;
t605 = t27 * qJD(1);
t164 = t275 * t351;
t28 = -t182 * t273 + (t123 * t353 + (-t156 - t622) * t356) * t350;
t599 = t28 * qJD(1);
t317 = t343 - t624;
t597 = t317 * t273;
t596 = t319 * t351;
t299 = t319 * t354;
t505 = t152 * t590;
t33 = -t179 * t273 - t505 - t623;
t595 = t33 * qJD(1);
t506 = t151 * t590;
t34 = t120 * t273 - t179 * t275 + t506;
t594 = t34 * qJD(1);
t344 = t350 ^ 2;
t349 = t356 ^ 2;
t593 = t344 * t349;
t592 = t344 * t353;
t589 = t351 * t103;
t585 = t351 * t203;
t584 = t351 * t287;
t583 = t351 * t318;
t582 = t351 * t352;
t581 = t351 * t355;
t578 = t352 * t354;
t574 = t354 * t103;
t570 = t354 * t203;
t569 = t354 * t204;
t568 = t354 * t287;
t567 = t354 * t318;
t565 = t355 * t273;
t41 = t151 * t591 - t255 * t504 - t283 * t273 + t449 * t590;
t562 = t41 * qJD(1);
t42 = t283 * t275 + (-t152 * t353 + (t548 + t607) * t356) * t350;
t561 = t42 * qJD(1);
t221 = t566 * t164;
t470 = -t610 / 0.2e1;
t471 = t203 * t637;
t44 = -t221 + (t470 + t651) * t354 + (t471 + t650) * t351;
t560 = t44 * qJD(1);
t47 = -t122 * t590 - t623;
t559 = t47 * qJD(1);
t54 = t241 * t203 + t204 * t240;
t558 = t54 * qJD(1);
t56 = -t569 + t585;
t59 = t56 * t275;
t557 = t59 * qJD(1);
t68 = -t255 * t273 - t506;
t556 = t68 * qJD(1);
t69 = -t255 * t275 - t505;
t555 = t69 * qJD(1);
t468 = -t600 / 0.2e1;
t377 = t345 * t468 + t351 * t470;
t72 = -t608 / 0.2e1 + t377;
t554 = t72 * qJD(1);
t389 = t463 + t467;
t469 = t610 / 0.2e1;
t78 = t469 - t301 / 0.2e1 + t389 * t351;
t553 = t78 * qJD(1);
t459 = t582 / 0.2e1;
t376 = (t575 / 0.2e1 + t356 * t459) * t350;
t381 = t275 * t456 + t471;
t80 = t376 - t381;
t552 = t80 * qJD(1);
t501 = t203 * t590;
t82 = -t240 * t275 + t355 * t501;
t551 = t82 * qJD(1);
t500 = t204 * t590;
t83 = t241 * t275 - t355 * t500;
t550 = t83 * qJD(1);
t92 = -t501 - t668;
t549 = t92 * qJD(1);
t346 = t352 ^ 2;
t348 = t355 ^ 2;
t328 = t348 - t346;
t547 = qJD(1) * t204;
t545 = qJD(1) * t356;
t544 = qJD(2) * t350;
t543 = qJD(2) * t352;
t542 = qJD(2) * t354;
t540 = qJD(3) * t352;
t539 = qJD(3) * t355;
t538 = qJD(4) * t352;
t537 = qJD(4) * t354;
t535 = qJD(5) * t275;
t534 = qJD(5) * t351;
t533 = qJD(5) * t352;
t532 = qJD(5) * t354;
t531 = qJD(5) * t655;
t100 = t204 * t273 + t668;
t530 = t100 * qJD(1);
t101 = -t203 * t273 + t667;
t529 = t101 * qJD(1);
t236 = t275 * t582;
t153 = t236 / 0.2e1 + t275 * t459;
t118 = t355 * t499 - t153;
t528 = t118 * qJD(1);
t119 = t500 + t667;
t527 = t119 * qJD(1);
t107 = -t275 * t352 - t565;
t130 = t107 * t590;
t526 = t130 * qJD(1);
t525 = t151 * qJD(3);
t524 = t152 * qJD(3);
t237 = t275 * t578;
t183 = -t302 - t237;
t521 = t183 * qJD(1);
t195 = -t273 * t591 + t352 * t593;
t520 = t195 * qJD(1);
t196 = -t275 * t591 + t355 * t593;
t519 = t196 * qJD(1);
t518 = t203 * qJD(5);
t212 = pkin(1) * t592 + t283 * t625;
t517 = t212 * qJD(1);
t213 = t344 * pkin(1) * t356 - t281 * t625;
t516 = t213 * qJD(1);
t268 = t273 * qJD(3);
t515 = t273 * qJD(4);
t288 = (-t353 ^ 2 + t349) * t344;
t514 = t288 * qJD(1);
t513 = t346 * qJD(2);
t512 = t346 * qJD(4);
t511 = t350 * qJD(3);
t510 = t351 * qJD(3);
t509 = t354 * qJD(3);
t507 = pkin(8) * t540;
t495 = t351 * t508;
t494 = t354 * t508;
t493 = t351 * t533;
t492 = t352 * t532;
t491 = t351 * t546;
t490 = t344 * t545;
t489 = t353 * t544;
t487 = t350 * t545;
t486 = t356 * t511;
t485 = qJD(4) * t590;
t484 = t355 * t510;
t333 = t352 * t539;
t483 = t275 * t538;
t482 = t351 * t532;
t481 = t351 * t509;
t334 = t352 * t542;
t480 = t352 * t537;
t479 = qJ(4) * t203 / 0.2e1;
t478 = t204 * t653;
t477 = qJ(4) * t650;
t465 = t275 * t641;
t464 = -t591 / 0.2e1;
t460 = -t588 / 0.2e1;
t458 = -t581 / 0.2e1;
t457 = t573 / 0.2e1;
t455 = t565 / 0.2e1;
t454 = t275 * t634;
t451 = t266 / 0.2e1 - t267 / 0.2e1;
t450 = t625 * qJD(1);
t448 = t299 - t584;
t243 = t579 / 0.2e1;
t149 = t243 - t389;
t446 = t149 * qJD(1) - t332;
t331 = t351 * t543;
t445 = t164 * qJD(1) + t331;
t443 = pkin(8) * t461;
t442 = -pkin(8) * t577 / 0.2e1;
t441 = qJD(5) + t546;
t440 = qJD(5) + t543;
t439 = t541 * t592;
t438 = t353 * t490;
t437 = t273 * t487;
t311 = t355 * t487;
t436 = t351 * t494;
t435 = t355 * t481;
t434 = t355 * t462;
t433 = t354 * t461;
t431 = t350 * t450;
t430 = t625 * t544;
t429 = -t155 / 0.2e1 + t89 * t635;
t427 = 0.2e1 * t436;
t426 = t265 / 0.2e1 + t564 / 0.2e1;
t425 = t454 + t124 / 0.2e1;
t424 = t354 * t546 + t334;
t313 = -qJD(3) + t487;
t285 = -pkin(2) - t669;
t200 = t354 * t285 + t583;
t369 = t200 * t648 + t204 * t642 + t40 * t635;
t1 = t477 + (t121 / 0.2e1 + t275 * t645 + t124 * t637) * t354 + t658 * t351 + t369;
t414 = t568 + t596;
t58 = t414 * t352 - t319 * t582 + (t200 - t583) * t355;
t421 = -t1 * qJD(1) - t58 * qJD(2);
t199 = t351 * t285 - t567;
t358 = t199 * t648 + t203 * t642 + t39 * t635 + t448 * t647;
t391 = qJ(4) * t651 + t620 / 0.2e1;
t2 = t124 * t459 + (t465 - t658) * t354 + t358 + t391;
t57 = (-t199 - t567) * t355 + (t448 - t299) * t352;
t420 = -t2 * qJD(1) + t57 * qJD(2);
t419 = t355 * t433;
t418 = t601 - t670;
t132 = -t200 * t352 - t319 * t581;
t368 = t200 * t646 + t204 * t641 + t40 * t637;
t9 = t351 * t429 + t368 + t457;
t416 = -t9 * qJD(1) + t132 * qJD(2);
t415 = -t154 * t355 + t156 * t352;
t413 = t450 + qJD(2);
t370 = t199 * t647 + t319 * t652 + t39 * t638;
t10 = t354 * t429 + t370 + t460;
t131 = -t199 * t352 + t319 * t566;
t412 = -t10 * qJD(1) - t131 * qJD(2);
t210 = t307 * t355 + t317 * t352;
t359 = t621 / 0.2e1 + t179 * t637 + t273 * t643 + t317 * t646;
t29 = (t442 + t576) * t350 + t359 + t451;
t411 = t29 * qJD(1) + t210 * qJD(2);
t211 = -t307 * t352 + t317 * t355;
t407 = t443 - t282 / 0.2e1;
t31 = -t336 + t649 + t597 / 0.2e1 + (-t179 / 0.2e1 + t407) * t355 + t664;
t410 = t31 * qJD(1) - t211 * qJD(2);
t295 = t328 * t351;
t65 = -t236 + (t433 - t603 / 0.2e1 + t204 / 0.2e1) * t355;
t409 = -t65 * qJD(1) - t295 * qJD(2);
t297 = t328 * t354;
t64 = t644 - t237 + (-t602 / 0.2e1 + t652) * t355;
t408 = -t64 * qJD(1) - t297 * qJD(2);
t133 = t273 ^ 2 - t657;
t406 = t133 * qJD(1) + t107 * qJD(2);
t405 = t107 * qJD(1) + t328 * qJD(2);
t404 = t135 * qJD(1) - t510;
t403 = -t513 - t533;
t401 = -t268 + t437;
t400 = t312 - t437;
t398 = pkin(2) * t647 + t255 * t637;
t397 = t154 * t653 - t156 * pkin(3) / 0.2e1;
t396 = t269 * qJD(1) - t508 / 0.2e1;
t371 = t355 * t407 + t649;
t62 = t371 + t398;
t395 = pkin(2) * t543 - t62 * qJD(1);
t360 = pkin(2) * t273 / 0.2e1 + t607 / 0.2e1 + t350 * t442;
t60 = t360 + t451;
t394 = pkin(2) * t508 - t60 * qJD(1);
t393 = t350 * t428;
t214 = t275 * t311;
t392 = qJD(5) * t434 - t214;
t390 = -t120 * t317 / 0.2e1 + t179 * t643;
t388 = t585 / 0.2e1 - t569 / 0.2e1;
t387 = t352 * t634 - t624 / 0.2e1;
t13 = (t453 * t352 + t452 * t355) * pkin(8) + t390 + t397;
t386 = -t307 * t317 * qJD(2) + t13 * qJD(1);
t45 = t640 + t371 + t664;
t385 = t45 * qJD(1) + t307 * t543;
t384 = t350 * (-t307 * t356 + t633);
t383 = (t509 - t547) * t275;
t382 = -t354 * t441 - t334;
t380 = pkin(8) * t434 - t426 - t664;
t53 = (t570 + t611) * t355;
t70 = t203 ^ 2 - t204 ^ 2;
t379 = t70 * qJD(1) - t53 * qJD(2) - t56 * qJD(3);
t378 = t645 + t387;
t88 = t388 * t355;
t96 = t570 / 0.2e1 + t611 / 0.2e1;
t375 = t88 * qJD(2) - t96 * qJD(3) + t203 * t547;
t180 = t378 * t351;
t22 = t351 * t425 + t354 * t498 + t478;
t374 = -qJ(4) * t509 - t22 * qJD(1) - t180 * qJD(2);
t181 = t378 * t354;
t24 = -t351 * t498 + t354 * t425 + t479;
t373 = qJ(4) * t510 - t24 * qJD(1) - t181 * qJD(2);
t286 = (t347 / 0.2e1 - t345 / 0.2e1) * t355;
t372 = -t96 * qJD(1) + t286 * qJD(2) + t481;
t367 = qJD(3) * t423 + t536;
t366 = t348 * t351 * t542 - t88 * qJD(1) - t286 * qJD(3);
t296 = t327 * t348;
t365 = -t53 * qJD(1) - t296 * qJD(2) + 0.2e1 * t435;
t364 = -t56 * qJD(1) + t427 + t447;
t224 = t657 + t593;
t363 = t224 * qJD(1) + t275 * t543 - t486;
t342 = pkin(8) * t539;
t341 = t539 / 0.2e1;
t335 = t355 * t509;
t320 = qJD(2) * t463;
t306 = 0.2e1 * t355 * t482;
t289 = t313 * qJ(4);
t278 = -t311 + t539;
t277 = t313 * t354;
t272 = t286 * qJD(5);
t258 = (t490 - t511 / 0.2e1) * t353;
t215 = -t352 * t546 - t513;
t198 = t312 / 0.2e1 - t269 * qJD(3);
t148 = t243 + t468 + t463;
t147 = -t596 - t568 / 0.2e1 + t387 * t354;
t146 = t299 - t584 / 0.2e1 + t387 * t351;
t106 = (t312 - t268) * t275;
t104 = t107 * qJD(3);
t95 = t96 * qJD(5);
t86 = t88 * qJD(5);
t81 = t376 + t381;
t79 = t469 + t301 / 0.2e1 + (t467 + t464) * t351;
t71 = t608 / 0.2e1 + t377;
t67 = t203 * t635 + t354 * t455 + t237 + t644;
t66 = t204 * t636 + t351 * t455 + t236 + t419;
t63 = t355 * t443 + t398 + t426;
t61 = t360 - t451;
t55 = t56 * qJD(5);
t52 = t53 * qJD(5);
t46 = t640 + t380;
t43 = -t221 - t354 * t240 / 0.2e1 + t241 * t639 + t388 * t352;
t32 = -t597 / 0.2e1 + t179 * t635 - t336 + t380;
t30 = t352 * t443 - t359 + t451 + t503;
t25 = t354 * t454 + t479 + t89 * t639 - t572 / 0.2e1 - t589 / 0.2e1;
t23 = t351 * t454 + t478 + t354 * t656 - t587 / 0.2e1 + t574 / 0.2e1;
t18 = (-t580 / 0.2e1 - t631 / 0.2e1) * t590 + t659;
t14 = t659 * pkin(8) - t390 + t397;
t12 = t89 * t458 - t586 / 0.2e1 + t457 - t368;
t11 = t89 * t456 - t571 / 0.2e1 + t460 - t370;
t4 = t414 * t647 + (t572 + t589) * t638 - t399 * t458 + t89 * t459 + t477 + t619 / 0.2e1 - t369 + (-t434 * t655 + t465) * t351;
t3 = (t574 - t587) * t637 + t299 * t647 + t399 * t456 - t89 * t578 / 0.2e1 - t655 * t419 - t358 + t391;
t35 = [0, 0, 0, t439, t288 * qJD(2), t356 * t430, -t353 * t430, 0, -t212 * qJD(2), -t213 * qJD(2), t106, t130 * qJD(2) + t133 * qJD(3), -t196 * qJD(2) + t273 * t486, t195 * qJD(2) + t275 * t486, -t439, -t41 * qJD(2) - t69 * qJD(3), t42 * qJD(2) + t68 * qJD(3), t26 * qJD(2) + t21 * qJD(3) + t273 * t485, t28 * qJD(2) + t33 * qJD(3) + t275 * t515, t27 * qJD(2) + t34 * qJD(3) + t224 * qJD(4), t16 * qJD(2) + t15 * qJD(3) + t47 * qJD(4), (-qJD(2) * t241 - t275 * t510 - t518) * t204, t54 * qJD(2) + t59 * qJD(3) + t70 * qJD(5), t83 * qJD(2) + t100 * qJD(3) + t275 * t518, t82 * qJD(2) + t101 * qJD(3) + t204 * t535, t106, t7 * qJD(2) + t5 * qJD(3) - t92 * qJD(4) + t20 * qJD(5), t8 * qJD(2) + t6 * qJD(3) + t119 * qJD(4) + t19 * qJD(5); 0, 0, 0, t438, t514, t413 * t590, -t413 * t591, 0, -t283 * qJD(2) - t517, t281 * qJD(2) - t516, t214 + t666, t328 * t488 + t104 + t526, t352 * t489 - t519, t355 * t489 + t520, -t258, -t562 + (-t283 * t355 + t352 * t393) * qJD(2) + t63 * qJD(3), t561 + (t283 * t352 + t355 * t393) * qJD(2) + t61 * qJD(3), qJD(2) * t415 + t18 * qJD(3) + t606, t599 + (t182 * t355 + t352 * t384) * qJD(2) + t32 * qJD(3) + t148 * qJD(4), t605 + (-t182 * t352 + t355 * t384) * qJD(2) + t30 * qJD(3) + t483, t615 + (pkin(8) * t415 + t182 * t307) * qJD(2) + t14 * qJD(3) + t46 * qJD(4), t71 * qJD(3) - t86 + (-t495 - t547) * t241, t558 + t43 * qJD(3) - t52 + (t240 * t351 - t608) * t508, t550 + (-t135 * t348 + t241 * t352) * qJD(2) + t66 * qJD(3) + t81 * qJD(5), t551 + (-t240 * t352 - t348 * t499) * qJD(2) + t67 * qJD(3) + t79 * qJD(5), -t392 + t666, t628 + (t319 * t240 + t50 * t352 + (-t199 * t590 + t619) * t355) * qJD(2) + t3 * qJD(3) + t153 * qJD(4) + t12 * qJD(5), t627 + (t319 * t241 - t51 * t352 + (-t200 * t590 - t620) * t355) * qJD(2) + t4 * qJD(3) + t275 * t480 + t11 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t665, t406, t401, t313 * t275, t320, t63 * qJD(2) - t524 - t555, t61 * qJD(2) + t525 + t556, t609 + t18 * qJD(2) + (-t601 + t632) * qJD(3) - t515, t32 * qJD(2) + t524 + t595, t30 * qJD(2) - t485 - t525 + t594, t616 + t14 * qJD(2) + (-t152 * pkin(3) - t151 * qJ(4)) * qJD(3) + t122 * qJD(4), t71 * qJD(2) + t351 * t383 + t95, t43 * qJD(2) - t275 * t447 - t55 + t557, t66 * qJD(2) - t273 * t509 + t530, t67 * qJD(2) + t273 * t510 + t529, t663, t630 + t3 * qJD(2) + (-t351 * t399 - t354 * t418) * qJD(3) - t203 * qJD(4) + t23 * qJD(5), t629 + t4 * qJD(2) + (t351 * t418 - t354 * t399) * qJD(3) - qJD(4) * t204 + t25 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t401, t148 * qJD(2) + t184, t363, t46 * qJD(2) + t122 * qJD(3) + t559, 0, 0, 0, 0, 0, t153 * qJD(2) - qJD(3) * t203 - t549, -t204 * qJD(3) + t275 * t334 + t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t375, t379, t81 * qJD(2) + t203 * t441, t79 * qJD(2) + t204 * t441, t198, t12 * qJD(2) + t23 * qJD(3) - t40 * qJD(5) + t612, t11 * qJD(2) + t25 * qJD(3) + t39 * qJD(5) + t613; 0, 0, 0, -t438, -t514, -t356 * t431, t353 * t431, 0, t517, t516, -t214 + t157, t104 - t526, -t355 * t486 + t519, t352 * t486 - t520, t258, t62 * qJD(3) + t562, t60 * qJD(3) - t561, -t17 * qJD(3) - t355 * t485 - t606, -t31 * qJD(3) + t149 * qJD(4) - t599, -t29 * qJD(3) + t483 - t605, -t13 * qJD(3) - t45 * qJD(4) - t615, t72 * qJD(3) + t241 * t547 - t86, t44 * qJD(3) - t52 - t558, -t65 * qJD(3) - t80 * qJD(5) - t550, -t64 * qJD(3) + t78 * qJD(5) - t551, t157 + t392, -t2 * qJD(3) - t118 * qJD(4) - t9 * qJD(5) - t628, -t1 * qJD(3) - t183 * qJD(4) - t10 * qJD(5) - t627; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t333, t328 * qJD(3), 0, 0, 0, -pkin(2) * t540, -pkin(2) * t539, 0, t211 * qJD(3) - t352 * t536, -t210 * qJD(3) + t512, (qJD(3) * t317 - t538) * t307, -t333 * t345 + t348 * t482, -t296 * qJD(5) - 0.2e1 * t352 * t435, -t295 * qJD(3) - t355 * t492, -t297 * qJD(3) + t355 * t493, t333, t57 * qJD(3) + t132 * qJD(5) + t351 * t512, -t58 * qJD(3) - t131 * qJD(5) + t354 * t512; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t662, t405, t278, t313 * t352, qJD(1) * t464, -t342 - t395, -t394 + t507, t367 - t614, t342 - t410, -t411 - t507, pkin(8) * t367 - t386, t554 - t272 + (-t345 * t508 + t481) * t352, t560 + t306 + (-0.2e1 * t436 - t447) * t352, t335 + t409, t408 - t484, t661, t146 * qJD(5) - t318 * t510 - t660 * t354 + t420, t147 * qJD(5) - t318 * t509 + t660 * t351 + t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, t446, -t215, t342 - t385, 0, 0, 0, 0, 0, t351 * t513 + t335 - t528, t354 * t513 - t484 - t521; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t366, t365, -t440 * t566 - t552, t440 * t581 + t553, -t311 / 0.2e1 + t341, t146 * qJD(3) - t200 * qJD(5) + t416, t147 * qJD(3) + t199 * qJD(5) + t412; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t665, -t406, t400, (-t543 - t546) * t590, t320, -t62 * qJD(2) + t555, -t60 * qJD(2) - t556, t17 * qJD(2) - t609, t31 * qJD(2) - t595, t29 * qJD(2) - t485 - t594, -qJ(4) * t485 + t13 * qJD(2) - t616, -t72 * qJD(2) + t204 * t491 + t95, -t44 * qJD(2) - t55 - t557, t65 * qJD(2) - t275 * t534 - t530, t64 * qJD(2) - t275 * t532 - t529, -t663, t2 * qJD(2) - t135 * qJD(4) + t22 * qJD(5) - t630, t1 * qJD(2) + t24 * qJD(5) - t354 * t485 - t629; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t662, -t405, t311, -t352 * t487, qJD(1) * t463, t395, t394, t614, t410, t411, t386, t332 * t345 - t272 - t554, t352 * t427 + t306 - t560, -t409 - t493, -t408 - t492, -t661, t180 * qJD(5) - t420, t181 * qJD(5) - t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), -t482, t327 * qJD(5), 0, 0, 0, qJ(4) * t532 + qJD(4) * t351, -qJ(4) * t534 + t537; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t313, -t289, 0, 0, 0, 0, 0, -t404, -t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t372, t364, -t351 * t441 - t331, t382, t396, t351 * t531 - t374, t354 * t531 - t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t400, -t149 * qJD(2) - t184, -t363, qJ(4) * t486 + t45 * qJD(2) - t559, 0, 0, 0, 0, 0, t118 * qJD(2) + t135 * qJD(3) - t164 * qJD(5) + t549, -t527 + t183 * qJD(2) + (t486 - t535) * t354; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311, -t446, t215, t385, 0, 0, 0, 0, 0, t351 * t403 + t528, t354 * t403 + t521; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t313, t289, 0, 0, 0, 0, 0, t404, t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t445 - t534, t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t375, -t379, t80 * qJD(2) + (-qJD(1) * t203 + t510) * t275, -t78 * qJD(2) + t383, t198, t9 * qJD(2) - t22 * qJD(3) + t164 * qJD(4) - t612, t10 * qJD(2) - t24 * qJD(3) + t275 * t537 - t613; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t366, -t365, t552 + (t494 + t510) * t352, -t553 + (-t495 + t509) * t352, t311 / 0.2e1 + t341, -t180 * qJD(3) + t351 * t538 - t416, -t181 * qJD(3) - t412 + t480; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t372, -t364, t331 + t491, t424, -t396, t374, t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t445, t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t35;
