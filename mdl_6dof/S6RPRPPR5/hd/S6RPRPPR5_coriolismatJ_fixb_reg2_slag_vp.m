% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRPPR5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:45
% EndTime: 2019-03-09 02:52:01
% DurationCPUTime: 11.62s
% Computational Cost: add. (14679->521), mult. (27697->652), div. (0->0), fcn. (32632->8), ass. (0->406)
t380 = sin(pkin(10));
t610 = pkin(3) + qJ(5);
t608 = -pkin(8) - t610;
t357 = t608 * t380;
t385 = sin(qJ(6));
t382 = cos(pkin(10));
t450 = t608 * t382;
t617 = cos(qJ(6));
t254 = t357 * t617 + t385 * t450;
t627 = t254 / 0.2e1;
t471 = t617 * t380;
t544 = t385 * t382;
t349 = t471 + t544;
t383 = cos(pkin(9));
t618 = cos(qJ(3));
t370 = t618 * t383;
t381 = sin(pkin(9));
t386 = sin(qJ(3));
t542 = t386 * t381;
t351 = -t370 + t542;
t177 = t349 * t351;
t470 = t617 * t382;
t295 = t351 * t470;
t545 = t385 * t380;
t237 = t351 * t545 - t295;
t355 = t381 * t618 + t386 * t383;
t609 = pkin(7) + qJ(2);
t451 = t609 * t381;
t348 = t618 * t451;
t361 = t609 * t383;
t543 = t386 * t361;
t278 = t348 + t543;
t612 = t355 * pkin(4);
t420 = t278 + t612;
t403 = t382 * t420;
t371 = -t383 * pkin(2) - pkin(1);
t547 = t355 * qJ(4);
t421 = t371 - t547;
t387 = t355 * pkin(5) + t403 + (t351 * t608 - t421) * t380;
t661 = t351 * t610;
t394 = t421 + t661;
t660 = t380 * t420;
t118 = t382 * t394 + t660;
t250 = t382 * t351;
t96 = pkin(8) * t250 + t118;
t53 = t385 * t96 - t387 * t617;
t54 = t385 * t387 + t617 * t96;
t353 = t470 - t545;
t623 = t353 / 0.2e1;
t625 = -t349 / 0.2e1;
t253 = t357 * t385 - t450 * t617;
t629 = -t253 / 0.2e1;
t663 = t177 * t629 + t237 * t627 - t53 * t623 - t54 * t625;
t548 = t353 * t355;
t662 = t548 / 0.2e1;
t323 = t349 * qJD(6);
t402 = t471 / 0.2e1 + t544 / 0.2e1;
t624 = t349 / 0.2e1;
t392 = (t624 + t402) * t355;
t637 = t392 * qJD(1);
t659 = -t323 - t637;
t631 = t355 ^ 2;
t633 = t351 ^ 2;
t643 = t631 + t633;
t658 = qJD(2) * t643;
t456 = t177 * t624;
t552 = t353 * t237;
t109 = -t552 / 0.2e1 - t456;
t657 = t109 * qJD(6);
t656 = t177 * qJD(1);
t642 = t633 - t631;
t655 = t642 * qJD(1);
t654 = t642 * qJD(3);
t653 = t643 * qJD(1);
t484 = t353 * qJD(3);
t652 = -qJD(1) * t109 + t349 * t484;
t651 = qJD(3) * t109 - t237 * t656;
t399 = t349 * t355;
t630 = t399 / 0.2e1;
t650 = t399 * t53;
t622 = -t355 / 0.2e1;
t265 = t349 * t622;
t374 = t380 ^ 2;
t376 = t382 ^ 2;
t452 = t376 / 0.2e1 + t374 / 0.2e1;
t649 = t351 * t452;
t648 = t399 * t353;
t117 = -t380 * t394 + t403;
t432 = t117 * t380 - t118 * t382;
t64 = t432 * t355;
t647 = t484 + t656;
t398 = t625 + t402;
t159 = t398 * t351;
t646 = qJD(1) * t159 + t484;
t488 = t349 * qJD(3);
t526 = qJD(1) * t237;
t645 = t488 + t526;
t454 = t545 / 0.2e1;
t437 = t454 + t623;
t626 = -t295 / 0.2e1;
t154 = t351 * t437 + t626;
t644 = qJD(1) * t154 + t488;
t441 = t386 * t451;
t473 = t618 * t361;
t279 = t473 - t441;
t430 = t278 * t355 - t279 * t351;
t641 = qJD(2) * t430;
t640 = qJD(3) * (t547 - t661);
t391 = t355 * t402 + t265;
t639 = qJD(6) * t391;
t638 = qJD(6) * t392;
t636 = t430 * qJD(1);
t489 = t631 * qJD(1);
t342 = t542 / 0.2e1 - t370 / 0.2e1;
t635 = t177 ^ 2;
t634 = t349 ^ 2;
t632 = t353 ^ 2;
t628 = t253 / 0.2e1;
t368 = pkin(5) * t380 + qJ(4);
t621 = -t368 / 0.2e1;
t616 = pkin(3) * t351;
t615 = pkin(3) * t355;
t614 = pkin(8) * t355;
t613 = t351 * pkin(4);
t611 = t382 * pkin(5);
t550 = t353 * t177;
t554 = t349 * t237;
t110 = t550 + t554;
t607 = t110 * qJD(5);
t107 = -t554 / 0.2e1 - t550 / 0.2e1;
t606 = t107 * qJD(5);
t231 = t279 - t613;
t228 = t231 * t382;
t596 = qJ(4) * t351;
t229 = t355 * t610 + t596;
t91 = -t351 * pkin(5) + t228 + (-t229 - t614) * t380;
t604 = t385 * t91;
t296 = t355 * t470;
t239 = t355 * t545 - t296;
t475 = t617 * t91;
t227 = t231 * t380;
t124 = t382 * t229 + t227;
t100 = t382 * t614 + t124;
t546 = t385 * t100;
t59 = t475 - t546;
t472 = t617 * t100;
t60 = t472 + t604;
t4 = -t177 * t59 - t237 * t60 - t239 * t54 + t650;
t603 = t4 * qJD(1);
t600 = t59 * t353;
t599 = t60 * t349;
t474 = -pkin(4) - t611;
t171 = t355 * t474 - t278;
t172 = t351 * t474 + t279;
t7 = t171 * t172 - t53 * t59 + t54 * t60;
t598 = t7 * qJD(1);
t395 = t239 * t628 + t351 * t621 + t399 * t627;
t418 = t59 * t624 - t60 * t353 / 0.2e1;
t19 = t395 + t418;
t594 = qJD(1) * t19;
t30 = -t172 * t237 + t355 * t53;
t593 = qJD(1) * t30;
t31 = t172 * t177 - t355 * t54;
t592 = qJD(1) * t31;
t580 = t118 * t380;
t581 = t117 * t382;
t34 = -t231 * t351 + (t580 + t581) * t355;
t591 = qJD(1) * t34;
t63 = t432 * t351;
t590 = qJD(1) * t63;
t589 = qJD(1) * t64;
t551 = t353 * t239;
t202 = -t551 / 0.2e1;
t560 = t548 * t353;
t68 = t202 + t560 / 0.2e1 + 0.2e1 * t630 * t349;
t588 = qJD(1) * t68;
t564 = t399 * t349;
t69 = t551 - t564;
t587 = qJD(1) * t69;
t71 = (-t399 / 0.2e1 + t630) * t353 + (t662 + t239 / 0.2e1) * t349;
t586 = qJD(1) * t71;
t569 = t239 * t177;
t571 = t237 * t399;
t85 = t569 - t571;
t585 = qJD(1) * t85;
t87 = -t177 * t399 + t237 * t548;
t584 = qJD(1) * t87;
t412 = -t456 + t552 / 0.2e1;
t93 = t412 - t342;
t583 = qJD(1) * t93;
t123 = -t229 * t380 + t228;
t579 = t123 * t380;
t578 = t123 * t382;
t577 = t124 * t380;
t576 = t124 * t382;
t15 = t171 * t237 + t172 * t239 + t351 * t53 + t355 * t59;
t575 = t15 * qJD(1);
t16 = t171 * t177 + t172 * t399 + t351 * t54 - t355 * t60;
t574 = t16 * qJD(1);
t21 = -t172 * t351 + t239 * t53 + t399 * t54;
t573 = t21 * qJD(1);
t23 = -t54 * t548 - t650;
t572 = t23 * qJD(1);
t570 = t237 * t351;
t568 = t239 * t355;
t24 = t117 * t123 + t118 * t124 - t231 * t420;
t567 = t24 * qJD(1);
t565 = t177 * t351;
t562 = t399 * t355;
t561 = t548 * t349;
t25 = -t64 + (t576 - t579) * t351;
t559 = t25 * qJD(1);
t396 = t172 * t623 - t177 * t621 + t254 * t622;
t405 = -t546 / 0.2e1 + t475 / 0.2e1;
t26 = -t396 + t405;
t558 = t26 * qJD(1);
t397 = t172 * t625 + t237 * t621 + t355 * t628;
t404 = -t604 / 0.2e1 - t472 / 0.2e1;
t27 = -t397 + t404;
t557 = t27 * qJD(1);
t32 = (t123 - t228) * t355 + (-t117 + t403) * t351;
t556 = t32 * qJD(1);
t33 = (-t124 + t227) * t355 + (t118 - t660) * t351;
t555 = t33 * qJD(1);
t553 = t349 * t239;
t549 = t353 * t351;
t309 = t374 * t355;
t310 = t376 * t355;
t248 = t380 * t355;
t312 = -t596 / 0.2e1;
t439 = t452 * t355;
t406 = -t439 * t610 + t312;
t415 = t579 / 0.2e1 - t576 / 0.2e1;
t55 = t406 + t415;
t541 = t55 * qJD(1);
t246 = t421 + t616;
t255 = t596 + t615;
t79 = t246 * t255;
t540 = t79 * qJD(1);
t86 = -t569 - t571;
t539 = t86 * qJD(1);
t537 = t374 + t376;
t363 = t381 ^ 2 + t383 ^ 2;
t111 = -t568 + t570;
t536 = qJD(1) * t111;
t112 = -t568 - t570;
t535 = qJD(1) * t112;
t113 = t562 - t565;
t534 = qJD(1) * t113;
t114 = -t562 - t565;
t533 = qJD(1) * t114;
t125 = -t246 * t355 - t255 * t351;
t532 = qJD(1) * t125;
t126 = t246 * t351 - t255 * t355;
t531 = qJD(1) * t126;
t169 = t643 * t380;
t527 = qJD(1) * t169;
t524 = qJD(2) * t351;
t523 = qJD(3) * t250;
t522 = qJD(3) * t278;
t521 = qJD(3) * t368;
t520 = qJD(3) * t380;
t519 = qJD(4) * t631;
t518 = qJD(4) * t355;
t517 = qJD(5) * t237;
t516 = qJD(5) * t177;
t515 = qJD(5) * t355;
t514 = qJD(6) * t177;
t148 = t439 + t309 / 0.2e1 + t310 / 0.2e1;
t512 = t148 * qJD(1);
t155 = -t296 / 0.2e1 + t437 * t355;
t511 = t155 * qJD(1);
t455 = -t545 / 0.2e1;
t156 = t355 * t455 + t296 / 0.2e1 + t662;
t150 = t156 * qJD(1);
t161 = t398 * t355;
t508 = t161 * qJD(1);
t162 = 0.2e1 * t265;
t507 = t162 * qJD(1);
t401 = t470 / 0.2e1 + t455;
t163 = (t623 + t401) * t355;
t506 = t163 * qJD(1);
t313 = t596 / 0.2e1;
t167 = 0.2e1 * t313 + t615;
t505 = t167 * qJD(1);
t170 = t537 * t355 * t351;
t504 = t170 * qJD(1);
t268 = -t549 / 0.2e1;
t178 = 0.2e1 * t268;
t502 = t178 * qJD(1);
t208 = -t342 - t649;
t501 = t208 * qJD(1);
t213 = t642 * t380;
t500 = t213 * qJD(1);
t214 = t642 * t382;
t499 = t214 * qJD(1);
t215 = t643 * t382;
t498 = t215 * qJD(1);
t235 = t237 * qJD(6);
t495 = t250 * qJD(1);
t251 = -t309 - t310;
t494 = t251 * qJD(1);
t252 = t537 * t633;
t493 = t252 * qJD(1);
t256 = t279 * qJD(3);
t490 = t342 * qJD(1);
t487 = t351 * qJD(1);
t486 = t351 * qJD(3);
t485 = t351 * qJD(4);
t330 = t353 * qJD(6);
t483 = t355 * qJD(1);
t482 = t355 * qJD(3);
t359 = -0.1e1 / 0.2e1 - t452;
t481 = t359 * qJD(3);
t360 = t363 * qJ(2);
t480 = t360 * qJD(1);
t479 = t537 * qJD(3);
t478 = t363 * qJD(1);
t477 = t382 * qJD(3);
t468 = t237 * t483;
t467 = t399 * t483;
t466 = t548 * t483;
t465 = t246 * t483;
t464 = t380 * t489;
t462 = t380 * t486;
t461 = t380 * t477;
t460 = t351 * t483;
t276 = t351 * t482;
t459 = t349 * t330;
t458 = t380 * t487;
t457 = t382 * t489;
t378 = qJD(3) * qJ(4);
t393 = t473 / 0.2e1 - t441 / 0.2e1;
t389 = -t613 / 0.2e1 + t393;
t414 = -t578 / 0.2e1 - t577 / 0.2e1;
t62 = t389 + t414;
t449 = qJD(1) * t62 + t378;
t206 = t177 * t483;
t448 = qJD(3) * t156 + t206;
t447 = qJD(6) * t342 + t460;
t446 = qJD(1) * t371 + qJD(2);
t445 = qJD(5) + t521;
t444 = t382 * t460;
t443 = t355 * t458;
t442 = -pkin(4) / 0.2e1 - t611 / 0.2e1;
t440 = -t634 / 0.2e1 - t632 / 0.2e1;
t438 = -t348 / 0.2e1 - t543 / 0.2e1;
t436 = -t235 - t468;
t435 = t382 * t443;
t84 = -t550 + t554;
t236 = t237 ^ 2;
t95 = t236 - t635;
t433 = qJD(1) * t95 + qJD(3) * t84;
t431 = t577 + t578;
t127 = t253 * t353 - t254 * t349;
t390 = t355 * t442 + t438;
t13 = t390 + t663;
t429 = -qJD(1) * t13 + qJD(3) * t127;
t388 = t351 * t442 + t399 * t629 - t548 * t627 + t393;
t417 = -t600 / 0.2e1 - t599 / 0.2e1;
t18 = t388 + t417;
t428 = qJD(1) * t18 + t521;
t232 = -t632 + t634;
t427 = qJD(1) * t84 + qJD(3) * t232;
t347 = t537 * t610;
t400 = -t612 / 0.2e1 + t438;
t416 = -t581 / 0.2e1 - t580 / 0.2e1;
t57 = t400 - t416;
t426 = -qJD(1) * t57 + qJD(3) * t347;
t22 = t177 * t53 - t237 * t54;
t425 = -qJD(1) * t22 - qJD(4) * t107;
t175 = -0.1e1 / 0.2e1 + t440;
t424 = qJD(1) * t107 + qJD(3) * t175;
t135 = t236 + t635;
t423 = qJD(1) * t135 + qJD(3) * t110;
t257 = t632 + t634;
t422 = qJD(1) * t110 + qJD(3) * t257;
t407 = qJD(3) * t155 + t206 + t514;
t379 = qJ(4) * qJD(4);
t358 = 0.1e1 / 0.2e1 - t452;
t333 = t355 * qJD(2);
t311 = t342 * qJD(3);
t209 = -t342 + t649;
t180 = t549 / 0.2e1 + t268;
t174 = 0.1e1 / 0.2e1 + t440;
t168 = t313 + t312;
t166 = -t548 / 0.2e1 + t401 * t355;
t158 = t177 / 0.2e1 + t402 * t351;
t153 = t351 * t454 + t268 + t626;
t152 = t155 * qJD(6);
t151 = t156 * qJD(6);
t130 = -t330 - t150;
t94 = -t412 - t342;
t83 = t84 * qJD(6);
t74 = t202 + t551 / 0.2e1;
t72 = -t560 / 0.2e1 + t399 * t625 + t564 / 0.2e1 + t202;
t70 = t561 / 0.2e1 - t648 - t553 / 0.2e1;
t61 = t389 - t414;
t58 = t400 + t416;
t56 = t406 - t415;
t29 = t396 + t405;
t28 = t397 + t404;
t20 = t395 - t418;
t17 = t388 - t417;
t14 = t390 - t663;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t363 * qJD(2), t360 * qJD(2), -t276, t654, 0, t276, 0, 0, t371 * t482, -t371 * t486, t658, t641, 0, 0, 0, -t276, t654, t276, t658, qJD(3) * t125 + t355 * t485, qJD(3) * t126 + t519, qJD(3) * t79 - t246 * t518 + t641, t374 * t276, 0.2e1 * t380 * t382 * t276, -t213 * qJD(3), t376 * t276, -t214 * qJD(3), -t276, t215 * qJD(2) + t32 * qJD(3) + (-qJD(5) * t351 + t518) * t248, -t169 * qJD(2) + t33 * qJD(3) + (-t351 * t515 + t519) * t382, qJD(3) * t25 - qJD(4) * t170 + qJD(5) * t252, qJD(2) * t34 + qJD(3) * t24 + qJD(4) * t64 - qJD(5) * t63 (qJD(3) * t399 - t235) * t177, qJD(3) * t86 + qJD(6) * t95, qJD(3) * t113 - t235 * t355 (qJD(3) * t239 + t514) * t237, qJD(3) * t111 - t355 * t514, -t276, t112 * qJD(2) + t15 * qJD(3) + t31 * qJD(6) + (qJD(4) * t399 - t516) * t355, t114 * qJD(2) + t16 * qJD(3) + t30 * qJD(6) + (qJD(4) * t548 + t517) * t355, qJD(2) * t85 + qJD(3) * t4 + qJD(4) * t87 + qJD(5) * t135, qJD(2) * t21 + qJD(3) * t7 + qJD(4) * t23 + qJD(5) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t478, t480, 0, 0, 0, 0, 0, 0, 0, 0, t653, t636, 0, 0, 0, 0, 0, 0, t653, 0, 0, qJD(3) * t168 + t636, 0, 0, 0, 0, 0, 0, t498, -t527, 0, qJD(3) * t56 + qJD(5) * t209 + t591, 0, 0, 0, 0, 0, 0, -t152 + t535, qJD(3) * t180 - qJD(6) * t161 + t533, qJD(3) * t74 + t585, t573 + (t553 + t648) * qJD(2) + t20 * qJD(3) + t72 * qJD(4) + t94 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t460, t655, -t486, t460, -t482, 0, t371 * t483 - t256, -t371 * t487 + t522, 0, 0, 0, t486, t482, -t460, t655, t460 (-t547 + t616) * qJD(3) - t485, t256 + t532, -t522 + t531, t540 + t168 * qJD(2) + (-pkin(3) * t279 - qJ(4) * t278) * qJD(3) + t279 * qJD(4) (t374 * t487 + t461) * t355, 0.2e1 * t435 + (t310 - t309) * qJD(3), -t351 * t477 - t500 (t376 * t487 - t461) * t355, t462 - t499, -t460, -t420 * t520 + t556 - t250 * qJD(4) + (-t515 - t640) * t382, -t420 * t477 + t555 + t248 * qJD(5) + (t485 + t640) * t380, -qJD(3) * t431 + t559, t567 + t56 * qJD(2) + (-qJ(4) * t420 - t431 * t610) * qJD(3) + t61 * qJD(4) + t58 * qJD(5), t399 * t647 + t657, t539 + (-t551 - t564) * qJD(3) + t83, -t351 * t484 + t534 + t639, t239 * t645 - t657, t349 * t486 - t152 + t536, -t447, t575 + (t171 * t349 + t239 * t368 + t253 * t351) * qJD(3) + t153 * qJD(4) - t156 * qJD(5) + t29 * qJD(6), t574 + t180 * qJD(2) + (t171 * t353 + t254 * t351 + t368 * t399) * qJD(3) + t158 * qJD(4) - t162 * qJD(5) + t28 * qJD(6), t603 + t74 * qJD(2) + (-t239 * t254 + t253 * t399 - t599 - t600) * qJD(3) + t70 * qJD(4) + t607, t598 + t20 * qJD(2) + (t171 * t368 - t253 * t59 + t254 * t60) * qJD(3) + t17 * qJD(4) + t14 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t486, t460, t489, t256 - t465, 0, 0, 0, 0, 0, 0, t464 - t523, t457 + t462, -t504, qJD(3) * t61 + t589, 0, 0, 0, 0, 0, 0, qJD(3) * t153 + t467 + t639, qJD(3) * t158 + qJD(6) * t166 + t466, qJD(3) * t70 + t584, t572 + t72 * qJD(2) + t17 * qJD(3) + (-t561 + t648) * qJD(4) + t606; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t458 - t477) * t355, qJD(3) * t248 - t444, t493, qJD(2) * t209 + qJD(3) * t58 - t590, 0, 0, 0, 0, 0, 0, -t448, -qJD(3) * t162 + t468, t423, qJD(2) * t94 + qJD(3) * t14 - t425; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t651, t433, qJD(3) * t391 + t436, -t651, -t407, -t311, -qJD(2) * t155 + qJD(3) * t29 + qJD(4) * t391 - qJD(6) * t54 + t592, -qJD(2) * t161 + qJD(3) * t28 + qJD(4) * t166 + qJD(6) * t53 + t593, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t478, -t480, 0, 0, 0, 0, 0, 0, t482, -t486, -t653, -t636, 0, 0, 0, 0, 0, 0, -t653, -t482, t486, qJD(3) * t167 - t518 - t636, 0, 0, 0, 0, 0, 0, t462 - t498, t523 + t527, -t251 * qJD(3), -qJD(3) * t55 - qJD(4) * t148 - qJD(5) * t208 - t591, 0, 0, 0, 0, 0, 0, qJD(3) * t177 - t151 - t535, -qJD(3) * t178 - qJD(6) * t162 - t533, -qJD(3) * t69 - t585, -qJD(3) * t19 - qJD(4) * t68 - qJD(5) * t93 - t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t483, -t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t483, t487, t505, 0, 0, 0, 0, 0, 0, t458, t495, -t494, -t541, 0, 0, 0, 0, 0, 0, t656, -t502, -t587, -t594; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t483, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t512, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t588; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t501, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t583; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t323 - t507, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t460, -t655, 0, -t460, 0, 0, -t446 * t355, t446 * t351, 0, 0, 0, 0, 0, t460, -t655, -t460, 0, t333 - t532, -t524 - t531, -qJD(2) * t167 - t540, -t374 * t460, -0.2e1 * t435, t500, -t376 * t460, t499, t460, -t380 * t524 - t556, -qJD(2) * t250 - t555, qJD(2) * t251 - t559, qJD(2) * t55 + qJD(4) * t62 - qJD(5) * t57 - t567, -t399 * t656 + t657, t83 - t539, -t534 - t638, -t239 * t526 - t657, -t151 - t536, t447, -qJD(2) * t177 + qJD(4) * t154 - qJD(5) * t155 - qJD(6) * t26 - t575, qJD(2) * t178 + qJD(4) * t159 - qJD(5) * t161 - qJD(6) * t27 - t574, qJD(2) * t69 + qJD(4) * t71 - t603 + t607, qJD(2) * t19 + qJD(4) * t18 - qJD(5) * t13 - t598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t483, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t483, -t487, -t505, 0, 0, 0, 0, 0, 0, -t458, -t495, t494, t541, 0, 0, 0, 0, 0, 0, -t656, t502, t587, t594; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t379, 0, 0, 0, 0, 0, 0, qJD(4) * t380, qJD(4) * t382, t537 * qJD(5), qJD(5) * t347 + t379, -t459, t232 * qJD(6), 0, t459, 0, 0, qJD(4) * t349 + t330 * t368, qJD(4) * t353 - t323 * t368, qJD(5) * t257, qJD(4) * t368 + qJD(5) * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t378, 0, 0, 0, 0, 0, 0, t520, t477, 0, qJD(5) * t358 + t449, 0, 0, 0, 0, 0, 0, t644, t646, t586, qJD(5) * t174 + t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t479, qJD(4) * t358 + t426, 0, 0, 0, 0, 0, 0, -t511, -t508, t422, qJD(4) * t174 + t429; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t652, t427, t659, t652, t130, t490, -qJD(6) * t254 + t368 * t484 - t558, qJD(6) * t253 - t368 * t488 - t557, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t460, -t489, t333 + t465, 0, 0, 0, 0, 0, 0, -t464, -t457, t504, qJD(2) * t148 - qJD(3) * t62 - t589, 0, 0, 0, 0, 0, 0, -qJD(3) * t154 - t467 - t638, -qJD(3) * t159 - qJD(6) * t163 - t466, -qJD(3) * t71 - t584, qJD(2) * t68 - qJD(3) * t18 - t572 + t606; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t483, 0, 0, 0, 0, 0, 0, 0, 0, 0, t512, 0, 0, 0, 0, 0, 0, 0, 0, 0, t588; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t378, 0, 0, 0, 0, 0, 0, -t520, -t477, 0, qJD(5) * t359 - t449, 0, 0, 0, 0, 0, 0, -t644, -t646, -t586, qJD(5) * t175 - t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t481, 0, 0, 0, 0, 0, 0, 0, 0, 0, t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t659, -t330 - t506, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t443, t444, -t493, qJD(2) * t208 + qJD(3) * t57 + t590, 0, 0, 0, 0, 0, 0, t407, qJD(3) * t161 + t436, -t423, qJD(2) * t93 + qJD(3) * t13 + t425; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t501, 0, 0, 0, 0, 0, 0, 0, 0, 0, t583; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t479, -qJD(4) * t359 - t426, 0, 0, 0, 0, 0, 0, t330 + t511, -t323 + t508, -t422, -qJD(4) * t175 - t429; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t481, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t647, -t645, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t651, -t433, qJD(3) * t392 + t468, t651, t448, -t311, qJD(2) * t156 + qJD(3) * t26 + qJD(4) * t392 - t516 - t592, qJD(2) * t162 + qJD(3) * t27 + qJD(4) * t163 + t517 - t593, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t507, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t652, -t427, t637, -t652, t150, -t490, -t353 * t445 + t558, t349 * t445 + t557, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t637, t506, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t647, t645, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
