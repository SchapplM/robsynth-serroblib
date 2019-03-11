% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR15_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR15_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:38:44
% EndTime: 2019-03-09 20:39:20
% DurationCPUTime: 16.36s
% Computational Cost: add. (19074->827), mult. (54763->1143), div. (0->0), fcn. (46370->14), ass. (0->368)
t321 = sin(pkin(6));
t324 = sin(qJ(3));
t325 = sin(qJ(2));
t329 = cos(qJ(2));
t524 = cos(pkin(7));
t544 = cos(qJ(3));
t426 = t524 * t544;
t353 = t324 * t329 + t325 * t426;
t253 = t353 * t321;
t243 = qJD(1) * t253;
t320 = sin(pkin(7));
t494 = qJD(3) * t324;
t468 = t320 * t494;
t569 = t243 - t468;
t525 = cos(pkin(6));
t456 = t325 * t525;
t312 = pkin(1) * t456;
t395 = t321 * (-pkin(10) * t524 - pkin(9));
t341 = t329 * t395 - t312;
t218 = t341 * qJD(1);
t540 = pkin(10) * t320;
t399 = pkin(2) * t325 - t329 * t540;
t496 = qJD(1) * t321;
t255 = t399 * t496;
t141 = -t218 * t320 + t524 * t255;
t457 = t324 * t524;
t428 = t325 * t457;
t401 = t321 * t428;
t471 = t544 * t329;
t439 = t321 * t471;
t407 = qJD(1) * t439;
t244 = -qJD(1) * t401 + t407;
t568 = pkin(3) * t468 + qJ(4) * t244 - t141;
t326 = sin(qJ(1));
t330 = cos(qJ(1));
t452 = t330 * t525;
t266 = t325 * t326 - t329 * t452;
t267 = t325 * t452 + t326 * t329;
t474 = t320 * t544;
t440 = t321 * t474;
t148 = t266 * t426 + t267 * t324 + t330 * t440;
t458 = t321 * t524;
t222 = t266 * t320 - t330 * t458;
t323 = sin(qJ(5));
t328 = cos(qJ(5));
t110 = t148 * t323 + t222 * t328;
t512 = t321 * t330;
t517 = t320 * t324;
t149 = -t266 * t457 + t267 * t544 - t512 * t517;
t322 = sin(qJ(6));
t327 = cos(qJ(6));
t567 = t110 * t322 - t149 * t327;
t566 = t110 * t327 + t149 * t322;
t493 = qJD(4) * t324;
t545 = pkin(3) + pkin(11);
t564 = -t243 * t545 + (-t493 + (pkin(11) * t324 - qJ(4) * t544) * qJD(3)) * t320 + t568;
t311 = pkin(2) * t457;
t454 = t329 * t525;
t313 = pkin(1) * t454;
t305 = qJD(1) * t313;
t370 = t325 * t395;
t217 = qJD(1) * t370 + t305;
t368 = -t324 * t217 + t218 * t426;
t515 = t321 * t325;
t441 = t545 * t515;
t473 = t544 * t255;
t557 = pkin(4) + pkin(10);
t563 = -(t557 * t474 + t311) * qJD(3) + t244 * pkin(4) + (-qJD(1) * t441 - t473) * t320 - t368;
t397 = qJD(3) * t426;
t565 = pkin(2) * t397 - t544 * t217 - t218 * t457 - t255 * t517;
t513 = t321 * t329;
t357 = pkin(9) * t513 + t312;
t448 = t525 * qJD(1);
t405 = t448 + qJD(2);
t378 = t320 * t405;
t453 = t329 * t524;
t429 = t321 * t453;
t171 = t357 * qJD(1) + (qJD(1) * t429 + t378) * pkin(10);
t344 = pkin(2) * t525 + t370;
t177 = qJD(2) * pkin(2) + qJD(1) * t344 + t305;
t398 = pkin(2) * t329 + t325 * t540;
t383 = -pkin(1) - t398;
t247 = t383 * t321;
t235 = qJD(1) * t247;
t85 = t171 * t324 - t177 * t426 - t235 * t474;
t508 = -qJD(4) - t85;
t470 = t325 * t496;
t437 = t320 * t470;
t560 = qJ(4) * t437 - t524 * qJD(4) - t565;
t355 = -t323 * t524 - t328 * t474;
t504 = -qJD(5) * t355 + t569 * t323 + t328 * t437;
t265 = -t323 * t474 + t328 * t524;
t502 = t265 * qJD(5) - t323 * t437 + t569 * t328;
t382 = t426 * t513;
t559 = -qJD(1) * t382 - t544 * t378;
t467 = qJD(3) * t544;
t433 = t320 * t467;
t396 = t433 - t244;
t436 = t324 * t470;
t178 = t436 + t559;
t479 = t320 * t513;
t245 = qJD(1) * t479 - t405 * t524 - qJD(3);
t130 = -t328 * t178 - t245 * t323;
t129 = qJD(6) + t130;
t354 = t324 * t453 + t325 * t544;
t333 = qJD(2) * t353 + qJD(3) * t354;
t362 = t324 * t378;
t442 = t525 * qJDD(1);
t400 = t442 + qJDD(2);
t369 = t400 * t320;
t485 = qJDD(1) * t325;
t104 = qJD(3) * t362 - qJDD(1) * t382 + (qJD(1) * t333 + t324 * t485) * t321 - t544 * t369;
t558 = t148 * t328 - t222 * t323;
t180 = t354 * t496 + t362;
t487 = qJD(5) + t180;
t132 = t178 * t323 - t245 * t328;
t346 = t354 * t321;
t175 = qJD(1) * t346 + qJD(5) + t362;
t91 = t132 * t322 - t327 * t175;
t556 = t487 * t91;
t528 = pkin(4) * t243 - t557 * t468 - t560;
t445 = t328 * t487;
t377 = pkin(2) * t426 - pkin(10) * t517;
t259 = -pkin(3) * t524 - t377;
t208 = pkin(4) * t517 - pkin(11) * t524 + t259;
t381 = -pkin(3) * t544 - qJ(4) * t324 - pkin(2);
t232 = (-pkin(11) * t544 + t381) * t320;
t503 = t323 * t208 + t328 * t232;
t498 = pkin(10) * t474 + t311;
t555 = t498 * qJD(3) + t368;
t451 = t525 * t320;
t554 = -t324 * t515 + t451 * t544 + t382;
t491 = qJD(5) * t328;
t492 = qJD(5) * t323;
t553 = -t208 * t491 + t232 * t492 + t563 * t323 - t564 * t328;
t552 = (qJDD(2) + 0.2e1 * t442) * t321;
t507 = pkin(4) * t180 - t508;
t446 = qJDD(1) * t524;
t418 = t329 * t446;
t447 = t524 * qJD(2);
t464 = t321 * t485;
t103 = (t447 + qJD(3)) * t436 - qJD(2) * t407 - t544 * t464 + (-t321 * t418 - t369) * t324 + t559 * qJD(3);
t101 = -qJDD(5) + t103;
t240 = t245 * qJ(4);
t86 = t544 * t171 + t177 * t457 + t235 * t517;
t64 = -pkin(4) * t178 + t86;
t58 = -t240 + t64;
t551 = -t101 * t545 - t487 * t58;
t486 = qJD(1) * qJD(2);
t465 = t321 * t486;
t432 = t325 * t465;
t484 = qJDD(1) * t329;
t463 = t321 * t484;
t187 = t400 * t524 + qJDD(3) + (t432 - t463) * t320;
t434 = pkin(1) * qJD(2) * t525;
t402 = qJD(1) * t434;
t431 = pkin(1) * t442;
t475 = pkin(9) * t463 + t325 * t431 + t329 * t402;
t356 = -pkin(9) * t432 + t475;
t419 = qJD(1) * t447;
t120 = (t369 + (-t325 * t419 + t418) * t321) * pkin(10) + t356;
t358 = -t325 * t402 + t329 * t431;
t466 = t329 * t486;
t380 = -t466 - t485;
t359 = t380 * pkin(9);
t128 = t400 * pkin(2) + ((-t325 * t446 - t329 * t419) * pkin(10) + t359) * t321 + t358;
t367 = t399 * qJD(2);
t154 = (qJD(1) * t367 + qJDD(1) * t383) * t321;
t425 = qJD(3) * t457;
t393 = t324 * t120 - t128 * t426 - t154 * t474 + t171 * t467 + t177 * t425 + t235 * t468;
t366 = qJDD(4) + t393;
t16 = -pkin(4) * t103 - t187 * t545 + t366;
t83 = -t128 * t320 + t524 * t154;
t350 = qJ(4) * t103 - qJD(4) * t180 + t83;
t19 = t104 * t545 + t350;
t54 = t245 * t545 + t507;
t127 = -t177 * t320 + t524 * t235;
t387 = -qJ(4) * t180 + t127;
t57 = t178 * t545 + t387;
t25 = t323 * t54 + t328 * t57;
t6 = -qJD(5) * t25 + t328 * t16 - t19 * t323;
t50 = t323 * t104 + t178 * t491 + t328 * t187 + t245 * t492;
t93 = t132 * t327 + t175 * t322;
t21 = qJD(6) * t93 + t327 * t101 + t322 * t50;
t268 = -t330 * t325 - t326 * t454;
t269 = -t326 * t456 + t329 * t330;
t152 = -t268 * t426 + t269 * t324 - t326 * t440;
t223 = -t268 * t320 + t326 * t458;
t111 = t152 * t328 - t223 * t323;
t263 = -t524 * t525 + t479;
t408 = t263 * t323 - t328 * t554;
t376 = g(1) * t111 + g(2) * t558 + g(3) * t408;
t4 = pkin(5) * t101 - t6;
t550 = t129 * (pkin(5) * t132 + pkin(12) * t129) + t376 + t4;
t290 = pkin(5) * t323 - pkin(12) * t328 + qJ(4);
t423 = pkin(5) * t328 + pkin(12) * t323;
t450 = -t328 * t104 + t187 * t323;
t51 = qJD(5) * t132 + t450;
t49 = qJDD(6) + t51;
t549 = t129 * (qJD(5) * t423 - (-pkin(4) - t423) * t180 - t508) + t290 * t49;
t430 = t324 * t451;
t210 = t430 + t346;
t202 = (t429 + t451) * pkin(10) + t357;
t215 = t313 + t344;
t342 = t324 * t202 - t215 * t426 - t247 * t474;
t61 = t210 * pkin(4) + t263 * t545 + t342;
t138 = -t215 * t320 + t524 * t247;
t522 = qJ(4) * t210;
t386 = t138 - t522;
t66 = -t545 * t554 + t386;
t410 = t323 * t61 + t328 * t66;
t140 = t554 * qJD(3) + (-t401 + t439) * qJD(2);
t306 = t329 * t434;
t219 = qJD(2) * t370 + t306;
t220 = t341 * qJD(2);
t347 = -t202 * t467 - t215 * t425 - t324 * t219 + t220 * t426 - t247 * t468;
t256 = t321 * t367;
t472 = t544 * t256;
t42 = t140 * pkin(4) + (-qJD(2) * t441 - t472) * t320 - t347;
t139 = qJD(3) * t430 + t321 * t333;
t142 = -t220 * t320 + t524 * t256;
t349 = -qJ(4) * t140 - qJD(4) * t210 + t142;
t44 = t139 * t545 + t349;
t547 = -qJD(5) * t410 - t323 * t44 + t328 * t42;
t176 = t187 * qJ(4);
t239 = qJD(4) * t245;
t394 = -t544 * t120 - t128 * t457 - t154 * t517 + t171 * t494 - t177 * t397 - t235 * t433;
t26 = -t176 + t239 + t394;
t17 = -pkin(4) * t104 - t26;
t10 = pkin(5) * t51 - pkin(12) * t50 + t17;
t5 = t323 * t16 + t328 * t19 + t54 * t491 - t492 * t57;
t3 = -pkin(12) * t101 + t5;
t23 = pkin(12) * t175 + t25;
t33 = pkin(5) * t130 - pkin(12) * t132 + t58;
t415 = t23 * t322 - t327 * t33;
t1 = -t415 * qJD(6) + t322 * t10 + t327 * t3;
t546 = t180 ^ 2;
t332 = qJD(1) ^ 2;
t317 = t321 ^ 2;
t543 = pkin(1) * t317;
t542 = pkin(3) * t187;
t541 = pkin(3) * t554;
t539 = -t396 * pkin(5) + t503 * qJD(5) + t564 * t323 + t563 * t328;
t523 = qJ(4) * t178;
t90 = t180 * t545 + t523;
t537 = t323 * t64 + t328 * t90;
t536 = t129 * t91;
t535 = t129 * t93;
t488 = qJD(6) * t327;
t489 = qJD(6) * t322;
t20 = -t322 * t101 - t132 * t489 + t175 * t488 + t327 * t50;
t534 = t20 * t322;
t533 = t245 * t86;
t531 = t322 * t49;
t530 = t327 * t49;
t529 = t328 * t20;
t527 = -pkin(3) * t243 + (-qJ(4) * t467 - t493) * t320 + t568;
t521 = t180 * t178;
t520 = t180 * t323;
t519 = t317 * t332;
t518 = t320 * t323;
t516 = t320 * t328;
t514 = t321 * t326;
t511 = t322 * t545;
t510 = t327 * t545;
t388 = -t265 * t322 + t327 * t517;
t506 = qJD(6) * t388 + t396 * t322 - t504 * t327;
t227 = t265 * t327 + t322 * t517;
t505 = qJD(6) * t227 - t504 * t322 - t396 * t327;
t501 = pkin(10) * t468 + t560;
t500 = -(-pkin(3) * t470 - t473) * t320 + t555;
t318 = t325 ^ 2;
t497 = -t329 ^ 2 + t318;
t495 = qJD(2) * t321;
t490 = qJD(5) * t545;
t482 = t329 * t543;
t481 = t325 * t519;
t480 = t320 * t515;
t477 = t544 * t202 + t215 * t457 + t247 * t517;
t257 = -t524 * qJ(4) - t498;
t469 = t325 * t495;
t449 = t487 * t93;
t444 = t487 * t323;
t443 = t129 * t327;
t88 = t263 * qJ(4) - t477;
t231 = pkin(4) * t474 - t257;
t435 = t320 * t469;
t135 = -pkin(5) * t355 - pkin(12) * t265 + t231;
t427 = -pkin(12) * t396 - qJD(6) * t135 + t553;
t424 = t321 * t332 * t525;
t422 = -g(1) * t148 + g(2) * t152;
t153 = t269 * t544 + (t268 * t524 + t320 * t514) * t324;
t421 = g(1) * t149 - g(2) * t153;
t126 = pkin(12) * t517 + t503;
t420 = -t502 * pkin(5) - t504 * pkin(12) + qJD(6) * t126 - t528;
t114 = -t178 * t322 + t327 * t520;
t417 = -t327 * t492 - t114;
t12 = t23 * t327 + t322 * t33;
t30 = pkin(12) * t210 + t410;
t146 = -t263 * t328 - t323 * t554;
t67 = pkin(4) * t554 - t88;
t45 = -pkin(5) * t408 - pkin(12) * t146 + t67;
t414 = t30 * t327 + t322 * t45;
t413 = -t30 * t322 + t327 * t45;
t24 = -t323 * t57 + t328 * t54;
t411 = -t323 * t66 + t328 * t61;
t106 = t146 * t327 + t210 * t322;
t105 = t146 * t322 - t210 * t327;
t409 = t208 * t328 - t232 * t323;
t404 = 0.2e1 * t448 + qJD(2);
t391 = -t328 * t101 - t175 * t444;
t390 = -t129 * t488 - t531;
t389 = -t129 * t489 + t530;
t384 = t323 * t42 + t328 * t44 + t61 * t491 - t492 * t66;
t375 = g(1) * t152 + g(2) * t148 - g(3) * t554;
t374 = g(1) * t153 + g(2) * t149 + g(3) * t210;
t167 = -t266 * t324 + t267 * t426;
t169 = t268 * t324 + t269 * t426;
t373 = g(1) * t169 + g(2) * t167 + g(3) * t253;
t168 = -t266 * t544 - t267 * t457;
t170 = t268 * t544 - t269 * t457;
t254 = (-t428 + t471) * t321;
t372 = -g(1) * t170 - g(2) * t168 - g(3) * t254;
t363 = t17 - t374;
t361 = t323 * t101 - t175 * t445;
t360 = -g(1) * t269 - g(2) * t267 - g(3) * t515;
t22 = -pkin(5) * t175 - t24;
t352 = -pkin(12) * t49 + (t22 + t24) * t129;
t351 = -t202 * t494 + t215 * t397 + t544 * t219 + t220 * t457 + t247 * t433 + t256 * t517;
t2 = -qJD(6) * t12 + t327 * t10 - t322 * t3;
t345 = -qJD(6) * t129 * t545 + t374;
t47 = -qJ(4) * t435 + t263 * qJD(4) - t351;
t340 = (-pkin(12) * t178 - qJD(6) * t290 + t537) * t129 + t375;
t339 = -t374 - t394;
t338 = t375 - t393;
t40 = -pkin(4) * t139 - t47;
t336 = t405 * t357;
t68 = pkin(3) * t178 + t387;
t335 = t180 * t68 + qJDD(4) - t338;
t334 = -t245 * t178 - t103;
t258 = t381 * t320;
t183 = t253 * t323 + t328 * t480;
t134 = t169 * t323 + t269 * t516;
t133 = t167 * t323 + t267 * t516;
t125 = -pkin(5) * t517 - t409;
t115 = pkin(3) * t180 + t523;
t113 = t327 * t178 + t322 * t520;
t112 = t152 * t323 + t223 * t328;
t89 = t263 * pkin(3) + t342;
t84 = t386 - t541;
t82 = qJD(5) * t408 + t139 * t323 + t328 * t435;
t81 = qJD(5) * t146 - t139 * t328 + t323 * t435;
t73 = t240 - t86;
t71 = pkin(3) * t245 - t508;
t70 = t112 * t327 + t153 * t322;
t69 = -t112 * t322 + t153 * t327;
t56 = pkin(3) * t139 + t349;
t55 = (-pkin(3) * t469 - t472) * t320 - t347;
t39 = -qJD(6) * t105 + t140 * t322 + t327 * t82;
t38 = qJD(6) * t106 - t140 * t327 + t322 * t82;
t34 = pkin(5) * t178 + t323 * t90 - t328 * t64;
t29 = -pkin(5) * t210 - t411;
t28 = t366 - t542;
t27 = pkin(3) * t104 + t350;
t13 = pkin(5) * t81 - pkin(12) * t82 + t40;
t8 = -pkin(5) * t140 - t547;
t7 = pkin(12) * t140 + t384;
t9 = [qJDD(1), g(1) * t326 - g(2) * t330, g(1) * t330 + g(2) * t326 (qJDD(1) * t318 + 0.2e1 * t325 * t466) * t317, 0.2e1 * (t325 * t484 - t486 * t497) * t317, t329 * t404 * t495 + t552 * t325, t552 * t329 - t404 * t469, t400 * t525, 0.2e1 * qJDD(1) * t482 - 0.2e1 * t325 * t486 * t543 - qJD(2) * t336 + (-pkin(9) * t515 + t313) * t400 + (t321 * t359 + t358) * t525 + g(1) * t267 - g(2) * t269 -(-pkin(9) * t469 + t306) * t405 - t357 * t400 - t356 * t525 - g(1) * t266 - g(2) * t268 + 0.2e1 * t380 * t543, -t103 * t210 + t140 * t180, -t103 * t554 - t104 * t210 - t139 * t180 - t140 * t178, t103 * t263 - t140 * t245 + t180 * t435 + t187 * t210, t104 * t263 + t139 * t245 - t178 * t435 + t187 * t554, -t187 * t263 - t245 * t435 -(t320 * t472 + t347) * t245 - t342 * t187 + t393 * t263 - t85 * t435 + t142 * t178 + t138 * t104 - t83 * t554 + t127 * t139 + t421, -t138 * t103 + t127 * t140 + t142 * t180 - t187 * t477 + t83 * t210 + t245 * t351 - t263 * t394 - t435 * t86 + t422, g(1) * t222 - g(2) * t223 - t103 * t89 + t104 * t88 + t139 * t73 + t140 * t71 + t178 * t47 + t180 * t55 + t210 * t28 - t26 * t554, -t104 * t84 - t139 * t68 - t178 * t56 + t187 * t89 - t245 * t55 - t263 * t28 + t27 * t554 + t435 * t71 - t421, t103 * t84 - t140 * t68 - t180 * t56 - t187 * t88 - t210 * t27 + t245 * t47 + t26 * t263 - t435 * t73 - t422, t27 * t84 + t68 * t56 + t26 * t88 + t73 * t47 + t28 * t89 + t71 * t55 - g(1) * (-t326 * pkin(1) - t267 * pkin(2) - pkin(3) * t149 + pkin(9) * t512 - pkin(10) * t222 - qJ(4) * t148) - g(2) * (t330 * pkin(1) + t269 * pkin(2) + t153 * pkin(3) + pkin(9) * t514 + pkin(10) * t223 + t152 * qJ(4)) t132 * t82 + t146 * t50, -t130 * t82 - t132 * t81 - t146 * t51 + t408 * t50, -t101 * t146 + t132 * t140 + t175 * t82 + t210 * t50, -t101 * t408 - t130 * t140 - t175 * t81 - t210 * t51, -t101 * t210 + t140 * t175, g(1) * t110 - g(2) * t112 - t411 * t101 + t40 * t130 + t24 * t140 - t17 * t408 + t547 * t175 + t6 * t210 + t67 * t51 + t58 * t81, g(1) * t558 - g(2) * t111 + t410 * t101 + t40 * t132 - t25 * t140 + t17 * t146 - t384 * t175 - t5 * t210 + t67 * t50 + t58 * t82, t106 * t20 + t39 * t93, -t105 * t20 - t106 * t21 - t38 * t93 - t39 * t91, t106 * t49 + t129 * t39 - t20 * t408 + t81 * t93, -t105 * t49 - t129 * t38 + t21 * t408 - t81 * t91, t129 * t81 - t408 * t49 (-qJD(6) * t414 + t13 * t327 - t322 * t7) * t129 + t413 * t49 - t2 * t408 - t415 * t81 + t8 * t91 + t29 * t21 + t4 * t105 + t22 * t38 + g(1) * t566 - g(2) * t70 -(qJD(6) * t413 + t13 * t322 + t327 * t7) * t129 - t414 * t49 + t1 * t408 - t12 * t81 + t8 * t93 + t29 * t20 + t4 * t106 + t22 * t39 - g(1) * t567 - g(2) * t69; 0, 0, 0, -t329 * t481, t497 * t519, -t329 * t424 + t464, t325 * t424 + t463, t400, pkin(1) * t481 - g(1) * t268 + g(2) * t266 - g(3) * t513 + qJD(1) * t336 + t358 + (-t329 * t465 - t464) * pkin(9), t332 * t482 + (-pkin(9) * t470 + t305) * t448 + t305 * qJD(2) - t360 - t475, -t103 * t517 + t180 * t396, t244 * t178 + t180 * t243 + (-t544 * t103 - t104 * t324 + (-t178 * t544 - t180 * t324) * qJD(3)) * t320, -t103 * t524 + t244 * t245 + (-t180 * t470 + t187 * t324 - t245 * t467) * t320, -t104 * t524 - t243 * t245 + (t178 * t470 + t187 * t544 + t245 * t494) * t320, t187 * t524 + t245 * t437, t377 * t187 - t393 * t524 - t141 * t178 - t127 * t243 + t555 * t245 + (-pkin(2) * t104 + t127 * t494 + t245 * t473 + t470 * t85 - t544 * t83) * t320 + t372, -t498 * t187 + t394 * t524 - t141 * t180 - t127 * t244 + t565 * t245 + (t86 * t470 + pkin(2) * t103 + t83 * t324 + (-pkin(10) * t245 * t324 + t127 * t544) * qJD(3)) * t320 + t373, -t259 * t103 + t257 * t104 - t73 * t243 - t71 * t244 + t500 * t180 + t501 * t178 + (-t544 * t26 + t28 * t324 + (t324 * t73 + t544 * t71) * qJD(3) + t360) * t320, t28 * t524 - t258 * t104 + t259 * t187 + t68 * t243 - t500 * t245 - t527 * t178 + (t27 * t544 - t470 * t71 - t494 * t68) * t320 - t372, -t26 * t524 + t258 * t103 - t257 * t187 + t68 * t244 + t501 * t245 - t527 * t180 + (-t27 * t324 - t467 * t68 + t470 * t73) * t320 - t373, t27 * t258 + t26 * t257 + t28 * t259 - g(1) * (pkin(2) * t268 + pkin(3) * t170 + qJ(4) * t169 + t269 * t540) - g(2) * (-pkin(2) * t266 + pkin(3) * t168 + qJ(4) * t167 + t267 * t540) - g(3) * (pkin(3) * t254 + qJ(4) * t253 + t321 * t398) + t501 * t73 + t500 * t71 + t527 * t68, -t132 * t504 + t265 * t50, t130 * t504 - t132 * t502 - t265 * t51 + t355 * t50, -t265 * t101 + t132 * t396 - t175 * t504 + t50 * t517, -t101 * t355 - t130 * t396 - t175 * t502 - t51 * t517, -t101 * t517 + t175 * t396, -t409 * t101 + t6 * t517 + t231 * t51 - t17 * t355 - g(1) * t134 - g(2) * t133 - g(3) * t183 + t502 * t58 + t396 * t24 + ((-qJD(5) * t232 - t563) * t328 + (-qJD(5) * t208 - t564) * t323) * t175 + t528 * t130, t503 * t101 - t5 * t517 + t231 * t50 + t17 * t265 - g(1) * (t169 * t328 - t269 * t518) - g(2) * (t167 * t328 - t267 * t518) - g(3) * (t253 * t328 - t323 * t480) - t504 * t58 - t396 * t25 + t553 * t175 + t528 * t132, t20 * t227 + t506 * t93, t20 * t388 - t21 * t227 - t505 * t93 - t506 * t91, t129 * t506 - t20 * t355 + t227 * t49 + t502 * t93, -t129 * t505 + t21 * t355 + t388 * t49 - t502 * t91, t129 * t502 - t355 * t49 (-t126 * t322 + t135 * t327) * t49 - t2 * t355 + t125 * t21 - t4 * t388 - g(1) * (t134 * t327 + t170 * t322) - g(2) * (t133 * t327 + t168 * t322) - g(3) * (t183 * t327 + t254 * t322) + t539 * t91 + t505 * t22 + (t322 * t427 - t327 * t420) * t129 - t502 * t415 -(t126 * t327 + t135 * t322) * t49 + t1 * t355 + t125 * t20 + t4 * t227 - g(1) * (-t134 * t322 + t170 * t327) - g(2) * (-t133 * t322 + t168 * t327) - g(3) * (-t183 * t322 + t254 * t327) + t539 * t93 + t506 * t22 + (t322 * t420 + t327 * t427) * t129 - t502 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t521, -t178 ^ 2 + t546, t334, -t180 * t245 - t104, t187, -t127 * t180 + t338 - t533, t127 * t178 + t245 * t85 - t339, pkin(3) * t103 - qJ(4) * t104 + (-t73 - t86) * t180 + (t71 + t508) * t178, t115 * t178 + t335 + t533 - 0.2e1 * t542, t115 * t180 - t178 * t68 + t245 * t508 + 0.2e1 * t176 - t239 + t339, -t26 * qJ(4) - t28 * pkin(3) - t68 * t115 - t71 * t86 - g(1) * (-pkin(3) * t152 + qJ(4) * t153) - g(2) * (-pkin(3) * t148 + qJ(4) * t149) - g(3) * (t522 + t541) + t508 * t73, -t132 * t444 + t328 * t50 (-t132 * t487 - t51) * t328 + (t130 * t487 - t50) * t323, t132 * t178 + t391, -t130 * t178 + t361, t175 * t178, qJ(4) * t51 + t24 * t178 + t507 * t130 + (-t64 * t175 - t551) * t328 + ((t90 + t490) * t175 + t363) * t323, qJ(4) * t50 + t537 * t175 - t25 * t178 + t507 * t132 + t551 * t323 + (t175 * t490 + t363) * t328, t327 * t529 + (-t328 * t489 + t417) * t93, t113 * t93 + t114 * t91 + (t322 * t93 + t327 * t91) * t492 + (-t534 - t21 * t327 + (t322 * t91 - t327 * t93) * qJD(6)) * t328, t20 * t323 + t417 * t129 + (t449 + t389) * t328, -t21 * t323 + (t322 * t492 + t113) * t129 + (t390 - t556) * t328, t129 * t445 + t323 * t49, -t22 * t113 - t34 * t91 + t549 * t327 + t340 * t322 + (t49 * t511 + t2 + (-t22 * t322 - t545 * t91) * qJD(5) - t345 * t327) * t323 + (t22 * t488 - t415 * t180 + t545 * t21 + t4 * t322 + (t129 * t511 - t415) * qJD(5)) * t328, -t22 * t114 - t34 * t93 - t549 * t322 + t340 * t327 + (t49 * t510 - t1 + (-t22 * t327 - t545 * t93) * qJD(5) + t345 * t322) * t323 + (-t22 * t489 - t12 * t180 + t545 * t20 + t4 * t327 + (t129 * t510 - t12) * qJD(5)) * t328; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t334, t187 - t521, -t245 ^ 2 - t546, -t245 * t73 + t335 - t542, 0, 0, 0, 0, 0, t130 * t245 + t391, t132 * t245 + t361, 0, 0, 0, 0, 0, -t328 * t21 + (t327 * t245 - t322 * t445) * t129 + (t390 + t556) * t323, -t529 + (-t322 * t245 - t327 * t445) * t129 + (t449 - t389) * t323; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132 * t130, -t130 ^ 2 + t132 ^ 2, t130 * t175 + t50, -t450 + (-qJD(5) + t175) * t132, -t101, -t132 * t58 + t175 * t25 - t376 + t6, g(1) * t112 + g(2) * t110 + g(3) * t146 + t130 * t58 + t175 * t24 - t5, t443 * t93 + t534 (t20 - t536) * t327 + (-t21 - t535) * t322, t129 * t443 - t132 * t93 + t531, -t129 ^ 2 * t322 + t132 * t91 + t530, -t129 * t132, -pkin(5) * t21 + t132 * t415 - t25 * t91 + t352 * t322 - t327 * t550, -pkin(5) * t20 + t12 * t132 - t25 * t93 + t322 * t550 + t352 * t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93 * t91, -t91 ^ 2 + t93 ^ 2, t20 + t536, -t21 + t535, t49, -g(1) * t69 + g(2) * t567 + g(3) * t105 + t12 * t129 - t22 * t93 + t2, g(1) * t70 + g(2) * t566 + g(3) * t106 - t415 * t129 + t22 * t91 - t1;];
tau_reg  = t9;
