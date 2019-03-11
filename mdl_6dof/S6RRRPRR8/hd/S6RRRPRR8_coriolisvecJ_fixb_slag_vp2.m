% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:44:10
% EndTime: 2019-03-09 18:45:30
% DurationCPUTime: 42.32s
% Computational Cost: add. (30430->967), mult. (78501->1360), div. (0->0), fcn. (63308->12), ass. (0->408)
t367 = cos(qJ(2));
t363 = sin(qJ(2));
t357 = sin(pkin(6));
t421 = qJD(1) * t357;
t407 = t363 * t421;
t359 = cos(pkin(6));
t420 = qJD(1) * t359;
t412 = pkin(1) * t420;
t312 = -pkin(8) * t407 + t367 * t412;
t373 = t357 * (pkin(2) * t363 - pkin(9) * t367);
t313 = qJD(1) * t373;
t362 = sin(qJ(3));
t366 = cos(qJ(3));
t240 = -t362 * t312 + t366 * t313;
t472 = -qJ(4) - pkin(9);
t398 = qJD(3) * t472;
t566 = -(-qJ(4) * t366 * t367 + pkin(3) * t363) * t421 - t240 - qJD(4) * t362 + t366 * t398;
t241 = t366 * t312 + t362 * t313;
t406 = t367 * t421;
t392 = t362 * t406;
t565 = -qJ(4) * t392 - qJD(4) * t366 - t362 * t398 + t241;
t356 = sin(pkin(12));
t358 = cos(pkin(12));
t331 = t356 * t366 + t358 * t362;
t273 = t331 * t406;
t322 = t331 * qJD(3);
t423 = t273 - t322;
t330 = t356 * t362 - t358 * t366;
t274 = t330 * t406;
t323 = t330 * qJD(3);
t422 = -t274 + t323;
t537 = t356 * t566 - t565 * t358;
t315 = pkin(8) * t406 + t363 * t412;
t271 = pkin(3) * t392 + t315;
t417 = qJD(3) * t362;
t564 = pkin(3) * t417 - t271;
t563 = pkin(10) * t407 - t537;
t562 = -pkin(4) * t423 + t422 * pkin(10) + t564;
t361 = sin(qJ(5));
t365 = cos(qJ(5));
t245 = t274 * t361 + t365 * t407;
t431 = t361 * t323;
t394 = t245 - t431;
t414 = qJD(5) * t365;
t561 = t331 * t414 + t394;
t560 = t361 * t563 + t562 * t365;
t355 = -pkin(3) * t366 - pkin(2);
t262 = pkin(4) * t330 - pkin(10) * t331 + t355;
t341 = t472 * t362;
t342 = t472 * t366;
t279 = t341 * t356 - t342 * t358;
t415 = qJD(5) * t361;
t541 = t262 * t414 - t279 * t415 + t562 * t361 - t365 * t563;
t347 = qJD(2) + t420;
t418 = qJD(2) * t367;
t404 = t362 * t418;
t416 = qJD(3) * t366;
t259 = -t347 * t417 + (-t363 * t416 - t404) * t421;
t433 = t357 * t367;
t327 = t359 * t363 * pkin(1) + pkin(8) * t433;
t317 = t327 * qJD(2);
t304 = qJD(1) * t317;
t215 = -t259 * pkin(3) + t304;
t403 = t366 * t418;
t258 = t347 * t416 + (-t363 * t417 + t403) * t421;
t191 = t258 * t358 + t259 * t356;
t457 = t191 * Ifges(5,4);
t281 = pkin(9) * t347 + t315;
t308 = (-pkin(2) * t367 - pkin(9) * t363 - pkin(1)) * t357;
t289 = qJD(1) * t308;
t219 = t281 * t366 + t289 * t362;
t314 = qJD(2) * t373;
t302 = qJD(1) * t314;
t434 = t357 * t363;
t348 = pkin(8) * t434;
t477 = pkin(1) * t367;
t326 = t359 * t477 - t348;
t316 = t326 * qJD(2);
t303 = qJD(1) * t316;
t162 = -qJD(3) * t219 + t366 * t302 - t303 * t362;
t295 = t347 * t362 + t366 * t407;
t419 = qJD(2) * t357;
t400 = qJD(1) * t419;
t391 = t363 * t400;
t100 = pkin(3) * t391 - qJ(4) * t258 - qJD(4) * t295 + t162;
t161 = -t281 * t417 + t289 * t416 + t362 * t302 + t366 * t303;
t294 = t347 * t366 - t362 * t407;
t109 = qJ(4) * t259 + qJD(4) * t294 + t161;
t53 = t356 * t100 + t358 * t109;
t280 = -t347 * pkin(2) - t312;
t225 = -t294 * pkin(3) + qJD(4) + t280;
t375 = t294 * t356 + t358 * t295;
t525 = t358 * t294 - t356 * t295;
t119 = -pkin(4) * t525 - pkin(10) * t375 + t225;
t50 = pkin(10) * t391 + t53;
t190 = t258 * t356 - t358 * t259;
t81 = t190 * pkin(4) - t191 * pkin(10) + t215;
t340 = qJD(3) - t406;
t218 = -t281 * t362 + t366 * t289;
t188 = -qJ(4) * t295 + t218;
t168 = pkin(3) * t340 + t188;
t189 = qJ(4) * t294 + t219;
t432 = t358 * t189;
t99 = t356 * t168 + t432;
t93 = pkin(10) * t340 + t99;
t13 = t119 * t414 + t361 * t81 + t365 * t50 - t415 * t93;
t55 = t119 * t361 + t365 * t93;
t14 = -qJD(5) * t55 - t361 * t50 + t365 * t81;
t364 = cos(qJ(6));
t224 = qJD(5) - t525;
t202 = t340 * t361 + t365 * t375;
t54 = t365 * t119 - t361 * t93;
t39 = -pkin(11) * t202 + t54;
t37 = pkin(5) * t224 + t39;
t360 = sin(qJ(6));
t201 = t340 * t365 - t361 * t375;
t40 = pkin(11) * t201 + t55;
t439 = t360 * t40;
t16 = t364 * t37 - t439;
t112 = qJD(5) * t201 + t191 * t365 + t361 * t391;
t6 = pkin(5) * t190 - pkin(11) * t112 + t14;
t113 = -qJD(5) * t202 - t191 * t361 + t365 * t391;
t7 = pkin(11) * t113 + t13;
t2 = qJD(6) * t16 + t360 * t6 + t364 * t7;
t438 = t364 * t40;
t17 = t360 * t37 + t438;
t3 = -qJD(6) * t17 - t360 * t7 + t364 * t6;
t546 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t533 = t14 * mrSges(6,1) - t13 * mrSges(6,2) + t546;
t559 = t533 + t215 * mrSges(5,1) - t53 * mrSges(5,3) - t457 / 0.2e1;
t246 = -t274 * t365 + t361 * t407;
t266 = t365 * t279;
t430 = t365 * t323;
t558 = pkin(11) * t246 + pkin(11) * t430 + (-t266 + (pkin(11) * t331 - t262) * t361) * qJD(5) + t560 - t423 * pkin(5);
t557 = pkin(11) * t561 - t541;
t539 = t565 * t356 + t358 * t566;
t554 = Ifges(4,3) + Ifges(5,3);
t553 = t190 * Ifges(5,2);
t353 = pkin(3) * t356 + pkin(10);
t473 = pkin(11) + t353;
t397 = qJD(5) * t473;
t436 = t525 * t361;
t179 = t356 * t189;
t108 = t188 * t358 - t179;
t476 = pkin(3) * t295;
t143 = pkin(4) * t375 - pkin(10) * t525 + t476;
t63 = t365 * t108 + t361 * t143;
t552 = -pkin(11) * t436 + t361 * t397 + t63;
t475 = pkin(11) * t365;
t62 = -t108 * t361 + t365 * t143;
t551 = -pkin(5) * t375 - t365 * t397 + t475 * t525 - t62;
t538 = pkin(4) * t407 - t539;
t396 = t364 * t201 - t202 * t360;
t125 = Ifges(7,4) * t396;
t133 = t201 * t360 + t202 * t364;
t465 = Ifges(7,4) * t133;
t221 = qJD(6) + t224;
t498 = -t221 / 0.2e1;
t506 = -t133 / 0.2e1;
t508 = -t396 / 0.2e1;
t58 = Ifges(7,1) * t133 + Ifges(7,5) * t221 + t125;
t98 = t168 * t358 - t179;
t92 = -pkin(4) * t340 - t98;
t73 = -pkin(5) * t201 + t92;
t548 = (Ifges(7,5) * t396 - Ifges(7,6) * t133) * t498 + (t133 * t17 + t16 * t396) * mrSges(7,3) + (-Ifges(7,2) * t133 + t125 + t58) * t508 - t73 * (mrSges(7,1) * t133 + mrSges(7,2) * t396) + (Ifges(7,1) * t396 - t465) * t506;
t441 = t340 * Ifges(5,6);
t448 = t375 * Ifges(5,4);
t451 = t525 * Ifges(5,2);
t155 = t441 + t448 + t451;
t547 = t16 * mrSges(7,1) + t225 * mrSges(5,1) + t54 * mrSges(6,1) - t155 / 0.2e1 - t17 * mrSges(7,2) - t55 * mrSges(6,2) - t99 * mrSges(5,3);
t199 = t365 * t262 - t279 * t361;
t164 = pkin(5) * t330 - t331 * t475 + t199;
t200 = t361 * t262 + t266;
t435 = t331 * t361;
t175 = -pkin(11) * t435 + t200;
t87 = t164 * t364 - t175 * t360;
t545 = qJD(6) * t87 + t360 * t558 - t557 * t364;
t88 = t164 * t360 + t175 * t364;
t544 = -qJD(6) * t88 + t557 * t360 + t364 * t558;
t454 = t221 * Ifges(7,3);
t461 = t133 * Ifges(7,5);
t462 = t396 * Ifges(7,6);
t56 = t454 + t461 + t462;
t453 = t224 * Ifges(6,3);
t455 = t202 * Ifges(6,5);
t456 = t201 * Ifges(6,6);
t89 = t453 + t455 + t456;
t543 = t89 + t56;
t542 = -qJD(5) * t200 + t560;
t540 = pkin(5) * t561 + t538;
t35 = qJD(6) * t396 + t112 * t364 + t113 * t360;
t517 = t35 / 0.2e1;
t36 = -qJD(6) * t133 - t112 * t360 + t113 * t364;
t516 = t36 / 0.2e1;
t57 = Ifges(7,2) * t396 + Ifges(7,6) * t221 + t465;
t535 = t57 / 0.2e1;
t504 = t190 / 0.2e1;
t110 = Ifges(6,6) * t113;
t111 = Ifges(6,5) * t112;
t41 = Ifges(6,3) * t190 + t110 + t111;
t33 = Ifges(7,6) * t36;
t34 = Ifges(7,5) * t35;
t8 = Ifges(7,3) * t190 + t33 + t34;
t534 = t41 + t8;
t328 = t473 * t361;
t329 = t473 * t365;
t260 = -t328 * t364 - t329 * t360;
t528 = qJD(6) * t260 + t551 * t360 - t552 * t364;
t261 = -t328 * t360 + t329 * t364;
t527 = -qJD(6) * t261 + t552 * t360 + t551 * t364;
t69 = -mrSges(7,1) * t396 + mrSges(7,2) * t133;
t526 = m(7) * t73 + t69;
t374 = t360 * t361 - t364 * t365;
t252 = t374 * t331;
t307 = pkin(9) * t359 + t327;
t234 = -t362 * t307 + t366 * t308;
t325 = t359 * t362 + t366 * t434;
t196 = -pkin(3) * t433 - t325 * qJ(4) + t234;
t235 = t366 * t307 + t362 * t308;
t324 = t359 * t366 - t362 * t434;
t210 = qJ(4) * t324 + t235;
t129 = t356 * t196 + t358 * t210;
t123 = -pkin(10) * t433 + t129;
t248 = -t358 * t324 + t325 * t356;
t249 = t324 * t356 + t325 * t358;
t306 = t348 + (-pkin(2) - t477) * t359;
t257 = -t324 * pkin(3) + t306;
t153 = t248 * pkin(4) - t249 * pkin(10) + t257;
t71 = t365 * t123 + t361 * t153;
t158 = t374 * t525;
t522 = qJD(5) + qJD(6);
t267 = t522 * t374;
t425 = -t267 + t158;
t334 = t360 * t365 + t361 * t364;
t157 = t334 * t525;
t268 = t522 * t334;
t424 = -t268 + t157;
t524 = Ifges(4,5) * t258 + Ifges(5,5) * t191 + Ifges(4,6) * t259 - Ifges(5,6) * t190 + t554 * t391;
t523 = t13 * t365 - t14 * t361;
t445 = t295 * Ifges(4,4);
t213 = t294 * Ifges(4,2) + t340 * Ifges(4,6) + t445;
t290 = Ifges(4,4) * t294;
t214 = t295 * Ifges(4,1) + t340 * Ifges(4,5) + t290;
t376 = t218 * t366 + t219 * t362;
t469 = Ifges(4,4) * t366;
t470 = Ifges(4,4) * t362;
t478 = t366 / 0.2e1;
t481 = t340 / 0.2e1;
t485 = t295 / 0.2e1;
t486 = t294 / 0.2e1;
t521 = -t376 * mrSges(4,3) + t280 * (mrSges(4,1) * t362 + mrSges(4,2) * t366) + (-Ifges(4,2) * t362 + t469) * t486 + (Ifges(4,1) * t366 - t470) * t485 + (Ifges(4,5) * t366 - Ifges(4,6) * t362) * t481 - t362 * t213 / 0.2e1 + t214 * t478;
t380 = t361 * t55 + t365 * t54;
t383 = Ifges(6,5) * t365 - Ifges(6,6) * t361;
t466 = Ifges(6,4) * t365;
t385 = -Ifges(6,2) * t361 + t466;
t467 = Ifges(6,4) * t361;
t387 = Ifges(6,1) * t365 - t467;
t388 = mrSges(6,1) * t361 + mrSges(6,2) * t365;
t479 = t365 / 0.2e1;
t480 = -t361 / 0.2e1;
t495 = t224 / 0.2e1;
t500 = t202 / 0.2e1;
t502 = t201 / 0.2e1;
t468 = Ifges(6,4) * t202;
t90 = t201 * Ifges(6,2) + t224 * Ifges(6,6) + t468;
t198 = Ifges(6,4) * t201;
t91 = t202 * Ifges(6,1) + t224 * Ifges(6,5) + t198;
t520 = -t380 * mrSges(6,3) + t383 * t495 + t385 * t502 + t387 * t500 + t388 * t92 + t479 * t91 + t480 * t90;
t519 = Ifges(7,4) * t517 + Ifges(7,2) * t516 + Ifges(7,6) * t504;
t518 = Ifges(7,1) * t517 + Ifges(7,4) * t516 + Ifges(7,5) * t504;
t43 = t112 * Ifges(6,1) + t113 * Ifges(6,4) + t190 * Ifges(6,5);
t515 = t43 / 0.2e1;
t513 = pkin(1) * mrSges(3,1);
t512 = pkin(1) * mrSges(3,2);
t510 = t112 / 0.2e1;
t509 = t113 / 0.2e1;
t507 = t396 / 0.2e1;
t505 = t133 / 0.2e1;
t503 = -t201 / 0.2e1;
t501 = -t202 / 0.2e1;
t497 = t221 / 0.2e1;
t496 = -t224 / 0.2e1;
t494 = t525 / 0.2e1;
t493 = t375 / 0.2e1;
t492 = -t248 / 0.2e1;
t490 = t249 / 0.2e1;
t489 = t258 / 0.2e1;
t488 = t259 / 0.2e1;
t484 = t324 / 0.2e1;
t483 = t325 / 0.2e1;
t482 = -t340 / 0.2e1;
t471 = Ifges(3,4) * t363;
t464 = Ifges(3,5) * t367;
t459 = t190 * Ifges(5,4);
t458 = t191 * Ifges(5,1);
t452 = t525 * Ifges(5,4);
t450 = t525 * Ifges(5,6);
t449 = t375 * Ifges(5,1);
t447 = t375 * Ifges(5,5);
t446 = t294 * Ifges(4,6);
t444 = t295 * Ifges(4,5);
t443 = t303 * mrSges(3,2);
t442 = t340 * Ifges(5,5);
t440 = t347 * Ifges(3,5);
t170 = -qJD(3) * t235 + t366 * t314 - t316 * t362;
t270 = qJD(3) * t324 + t357 * t403;
t405 = t363 * t419;
t124 = pkin(3) * t405 - qJ(4) * t270 - qJD(4) * t325 + t170;
t169 = -t307 * t417 + t308 * t416 + t362 * t314 + t366 * t316;
t269 = -qJD(3) * t325 - t357 * t404;
t134 = qJ(4) * t269 + qJD(4) * t324 + t169;
t68 = t356 * t124 + t358 * t134;
t144 = -t268 * t331 + t323 * t374;
t172 = t245 * t360 + t246 * t364;
t429 = t144 - t172;
t145 = t252 * t522 + t334 * t323;
t171 = t245 * t364 - t246 * t360;
t428 = t145 - t171;
t135 = -mrSges(6,1) * t201 + mrSges(6,2) * t202;
t209 = mrSges(5,1) * t340 - mrSges(5,3) * t375;
t427 = t209 - t135;
t426 = -mrSges(3,1) * t347 - mrSges(4,1) * t294 + mrSges(4,2) * t295 + mrSges(3,3) * t407;
t410 = -Ifges(5,3) / 0.2e1 - Ifges(4,3) / 0.2e1;
t354 = -pkin(3) * t358 - pkin(4);
t114 = t190 * mrSges(5,1) + t191 * mrSges(5,2);
t52 = t100 * t358 - t356 * t109;
t70 = -t123 * t361 + t365 * t153;
t67 = t124 * t358 - t356 * t134;
t107 = t188 * t356 + t432;
t128 = t196 * t358 - t356 * t210;
t393 = -t246 - t430;
t278 = -t358 * t341 - t342 * t356;
t122 = pkin(4) * t433 - t128;
t389 = mrSges(6,1) * t365 - mrSges(6,2) * t361;
t386 = Ifges(6,1) * t361 + t466;
t384 = Ifges(6,2) * t365 + t467;
t382 = Ifges(6,5) * t361 + Ifges(6,6) * t365;
t381 = -t13 * t361 - t14 * t365;
t372 = -t365 * t249 + t361 * t433;
t47 = pkin(5) * t248 + pkin(11) * t372 + t70;
t222 = -t361 * t249 - t365 * t433;
t60 = pkin(11) * t222 + t71;
t24 = -t360 * t60 + t364 * t47;
t25 = t360 * t47 + t364 * t60;
t379 = t361 * t54 - t365 * t55;
t148 = -mrSges(6,2) * t224 + mrSges(6,3) * t201;
t149 = mrSges(6,1) * t224 - mrSges(6,3) * t202;
t378 = t148 * t365 - t149 * t361;
t377 = t161 * t366 - t162 * t362;
t159 = t222 * t364 + t360 * t372;
t160 = t222 * t360 - t364 * t372;
t65 = pkin(10) * t405 + t68;
t206 = -t358 * t269 + t270 * t356;
t207 = t269 * t356 + t270 * t358;
t232 = -t269 * pkin(3) + t317;
t96 = t206 * pkin(4) - t207 * pkin(10) + t232;
t22 = -t123 * t415 + t153 * t414 + t361 * t96 + t365 * t65;
t64 = -pkin(4) * t405 - t67;
t49 = -pkin(4) * t391 - t52;
t23 = -qJD(5) * t71 - t361 * t65 + t365 * t96;
t343 = Ifges(3,4) * t406;
t339 = t400 * t464;
t336 = -pkin(5) * t365 + t354;
t311 = -t347 * mrSges(3,2) + mrSges(3,3) * t406;
t276 = Ifges(3,1) * t407 + t343 + t440;
t275 = Ifges(3,6) * t347 + (Ifges(3,2) * t367 + t471) * t421;
t264 = mrSges(4,1) * t340 - mrSges(4,3) * t295;
t263 = -mrSges(4,2) * t340 + mrSges(4,3) * t294;
t251 = t334 * t331;
t239 = pkin(5) * t435 + t278;
t231 = -mrSges(4,2) * t391 + mrSges(4,3) * t259;
t230 = mrSges(4,1) * t391 - mrSges(4,3) * t258;
t212 = t340 * Ifges(4,3) + t444 + t446;
t208 = -mrSges(5,2) * t340 + mrSges(5,3) * t525;
t194 = -mrSges(4,1) * t259 + mrSges(4,2) * t258;
t178 = t258 * Ifges(4,1) + t259 * Ifges(4,4) + Ifges(4,5) * t391;
t177 = t258 * Ifges(4,4) + t259 * Ifges(4,2) + Ifges(4,6) * t391;
t167 = mrSges(5,1) * t391 - mrSges(5,3) * t191;
t166 = -mrSges(5,2) * t391 - mrSges(5,3) * t190;
t163 = -mrSges(5,1) * t525 + mrSges(5,2) * t375;
t156 = t442 + t449 + t452;
t154 = t340 * Ifges(5,3) + t447 + t450;
t138 = qJD(5) * t372 - t361 * t207 + t365 * t405;
t137 = qJD(5) * t222 + t365 * t207 + t361 * t405;
t103 = Ifges(5,5) * t391 + t458 - t459;
t102 = Ifges(5,6) * t391 + t457 - t553;
t86 = mrSges(7,1) * t221 - mrSges(7,3) * t133;
t85 = -mrSges(7,2) * t221 + mrSges(7,3) * t396;
t83 = -pkin(5) * t222 + t122;
t78 = pkin(5) * t436 + t107;
t75 = -mrSges(6,2) * t190 + mrSges(6,3) * t113;
t74 = mrSges(6,1) * t190 - mrSges(6,3) * t112;
t59 = -mrSges(6,1) * t113 + mrSges(6,2) * t112;
t46 = -qJD(6) * t160 - t137 * t360 + t138 * t364;
t45 = qJD(6) * t159 + t137 * t364 + t138 * t360;
t42 = t112 * Ifges(6,4) + t113 * Ifges(6,2) + t190 * Ifges(6,6);
t38 = -pkin(5) * t138 + t64;
t32 = -pkin(5) * t113 + t49;
t31 = -mrSges(7,2) * t190 + mrSges(7,3) * t36;
t30 = mrSges(7,1) * t190 - mrSges(7,3) * t35;
t19 = t364 * t39 - t439;
t18 = -t360 * t39 - t438;
t15 = pkin(11) * t138 + t22;
t12 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t11 = pkin(5) * t206 - pkin(11) * t137 + t23;
t5 = -qJD(6) * t25 + t11 * t364 - t15 * t360;
t4 = qJD(6) * t24 + t11 * t360 + t15 * t364;
t1 = [(t159 * t2 - t16 * t45 - t160 * t3 + t17 * t46) * mrSges(7,3) + (-Ifges(6,4) * t372 + Ifges(6,2) * t222) * t509 + (t13 * t222 - t137 * t54 + t138 * t55 + t14 * t372) * mrSges(6,3) + (Ifges(6,5) * t137 + Ifges(6,6) * t138) * t495 + (Ifges(7,4) * t45 + Ifges(7,2) * t46) * t507 + (Ifges(7,4) * t160 + Ifges(7,2) * t159) * t516 + (Ifges(6,4) * t137 + Ifges(6,2) * t138) * t502 - t190 * (Ifges(5,4) * t249 - Ifges(5,6) * t433) / 0.2e1 + (Ifges(5,4) * t207 + Ifges(5,6) * t405) * t494 + t191 * (Ifges(5,1) * t249 - Ifges(5,5) * t433) / 0.2e1 + (Ifges(5,1) * t207 + Ifges(5,5) * t405) * t493 + (t553 / 0.2e1 + (Ifges(6,3) + Ifges(7,3)) * t504 + Ifges(6,6) * t509 + Ifges(6,5) * t510 + t534 / 0.2e1 + Ifges(7,6) * t516 + Ifges(7,5) * t517 + t559) * t248 + (-Ifges(6,5) * t372 + Ifges(7,5) * t160 + Ifges(6,6) * t222 + Ifges(7,6) * t159) * t504 + (t303 * t367 + t304 * t363 + (-t312 * t367 - t315 * t363) * qJD(2)) * t357 * mrSges(3,3) + (Ifges(7,1) * t45 + Ifges(7,4) * t46) * t505 + (Ifges(7,1) * t160 + Ifges(7,4) * t159) * t517 + t426 * t317 + t98 * (mrSges(5,1) * t405 - mrSges(5,3) * t207) + t218 * (mrSges(4,1) * t405 - mrSges(4,3) * t270) + (-Ifges(6,1) * t372 + Ifges(6,4) * t222) * t510 + (Ifges(6,1) * t137 + Ifges(6,4) * t138) * t500 + ((t212 + t154) * t363 + t347 * (-Ifges(3,6) * t363 + t464) + t367 * t276) * t419 / 0.2e1 + t49 * (-mrSges(6,1) * t222 - mrSges(6,2) * t372) - t372 * t515 + t52 * (-mrSges(5,1) * t433 - t249 * mrSges(5,3)) + t162 * (-mrSges(4,1) * t433 - t325 * mrSges(4,3)) + t161 * (mrSges(4,2) * t433 + t324 * mrSges(4,3)) + (t207 * t225 + t215 * t249 - t405 * t99 + t433 * t53) * mrSges(5,2) + (-Ifges(5,6) * t481 - Ifges(5,4) * t493 - Ifges(5,2) * t494 + Ifges(6,3) * t495 + Ifges(7,3) * t497 + Ifges(6,5) * t500 + Ifges(6,6) * t502 + Ifges(7,5) * t505 + Ifges(7,6) * t507 + t547 + t543 / 0.2e1) * t206 - t275 * t405 / 0.2e1 + t219 * (-mrSges(4,2) * t405 + mrSges(4,3) * t269) + t304 * (-mrSges(4,1) * t324 + mrSges(4,2) * t325) + t316 * t311 + t306 * t194 + t280 * (-mrSges(4,1) * t269 + mrSges(4,2) * t270) + t269 * t213 / 0.2e1 + t270 * t214 / 0.2e1 + t169 * t263 + t170 * t264 + t257 * t114 + t234 * t230 + t235 * t231 + t232 * t163 + t222 * t42 / 0.2e1 + t207 * t156 / 0.2e1 + t68 * t208 + t67 * t209 + t129 * t166 + t128 * t167 + t32 * (-mrSges(7,1) * t159 + mrSges(7,2) * t160) + t22 * t148 + t23 * t149 + (Ifges(7,5) * t45 + Ifges(7,6) * t46) * t497 + m(3) * (t303 * t327 - t304 * t326 - t312 * t317 + t315 * t316) + m(7) * (t16 * t5 + t17 * t4 + t2 * t25 + t24 * t3 + t32 * t83 + t38 * t73) + m(6) * (t122 * t49 + t13 * t71 + t14 * t70 + t22 * t55 + t23 * t54 + t64 * t92) + m(5) * (t128 * t52 + t129 * t53 + t215 * t257 + t225 * t232 + t67 * t98 + t68 * t99) + m(4) * (t161 * t235 + t162 * t234 + t169 * t219 + t170 * t218 + t280 * t317 + t304 * t306) - t524 * t433 / 0.2e1 + (Ifges(4,5) * t270 + Ifges(5,5) * t207 + Ifges(4,6) * t269 + t554 * t405) * t481 + (t339 / 0.2e1 - t443 - t304 * mrSges(3,1)) * t359 + ((Ifges(3,5) * t359 / 0.2e1 - t326 * mrSges(3,3) + (-0.2e1 * t512 + 0.3e1 / 0.2e1 * Ifges(3,4) * t367) * t357) * t367 + (Ifges(4,5) * t483 + Ifges(4,6) * t484 + Ifges(5,5) * t490 + Ifges(5,6) * t492 - Ifges(3,6) * t359 - t327 * mrSges(3,3) + (-0.2e1 * t513 - 0.3e1 / 0.2e1 * t471) * t357 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + t410) * t433) * t363) * t400 + t24 * t30 + t25 * t31 + t45 * t58 / 0.2e1 + t38 * t69 + t73 * (-mrSges(7,1) * t46 + mrSges(7,2) * t45) + t70 * t74 + t71 * t75 + t83 * t12 + t4 * t85 + t5 * t86 + t178 * t483 + t177 * t484 + (Ifges(4,1) * t270 + Ifges(4,4) * t269 + Ifges(4,5) * t405) * t485 + (Ifges(4,4) * t270 + Ifges(4,2) * t269 + Ifges(4,6) * t405) * t486 + (Ifges(4,4) * t325 + Ifges(4,2) * t324 - Ifges(4,6) * t433) * t488 + (Ifges(4,1) * t325 + Ifges(4,4) * t324 - Ifges(4,5) * t433) * t489 + t103 * t490 + t102 * t492 + t160 * t518 + t159 * t519 + t46 * t535 + t122 * t59 + t64 * t135 + t137 * t91 / 0.2e1 + t138 * t90 / 0.2e1 + t92 * (-mrSges(6,1) * t138 + mrSges(6,2) * t137); t541 * t148 + (t13 * t200 + t14 * t199 + t278 * t49 + t538 * t92 + t54 * t542 + t541 * t55) * m(6) + t542 * t149 + t537 * t208 + t538 * t135 + t539 * t209 + t540 * t69 + (-t102 / 0.2e1 + t111 / 0.2e1 + t110 / 0.2e1 + t34 / 0.2e1 + t33 / 0.2e1 + t8 / 0.2e1 + t41 / 0.2e1 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t190 + t559) * t330 + (t362 * pkin(3) * t163 + (-t263 * t362 - t264 * t366) * pkin(9) + t521) * qJD(3) + ((t312 * mrSges(3,3) + t421 * t512 - t343 / 0.2e1 - t440 / 0.2e1 - t276 / 0.2e1 - t521) * t367 + (t315 * mrSges(3,3) + t99 * mrSges(5,2) - t450 / 0.2e1 - t447 / 0.2e1 - t98 * mrSges(5,1) - t444 / 0.2e1 - t446 / 0.2e1 - t218 * mrSges(4,1) + t219 * mrSges(4,2) + t275 / 0.2e1 - t154 / 0.2e1 - t212 / 0.2e1 + t410 * t340 + (t513 + t471 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t367) * t421 + (-qJD(2) + t347 / 0.2e1) * Ifges(3,6) + (Ifges(4,5) * t362 + Ifges(5,5) * t331 + Ifges(4,6) * t366 - Ifges(5,6) * t330) * qJD(2) / 0.2e1) * t363) * t421 + (t59 - t167) * t278 + (t215 * t355 + t225 * t564 - t278 * t52 + t279 * t53 + t537 * t99 + t539 * t98) * m(5) - t426 * t315 + (-t230 * t362 + t231 * t366) * pkin(9) + (-pkin(2) * t304 + (-qJD(3) * t376 + t377) * pkin(9) - t218 * t240 - t219 * t241 - t280 * t315) * m(4) - t525 * (-Ifges(5,4) * t274 - Ifges(5,2) * t273) / 0.2e1 + (mrSges(6,1) * t394 + mrSges(6,2) * t393) * t92 + t377 * mrSges(4,3) + (mrSges(6,2) * t423 - mrSges(6,3) * t394) * t55 + (t422 * t98 + t423 * t99) * mrSges(5,3) + (-mrSges(6,1) * t423 - mrSges(6,3) * t393) * t54 + (-mrSges(5,1) * t423 - mrSges(5,2) * t422) * t225 - t443 + (-mrSges(4,1) * t366 + mrSges(4,2) * t362 - mrSges(3,1)) * t304 + t339 + (-t245 / 0.2e1 + t431 / 0.2e1) * t90 + (-Ifges(6,5) * t430 + Ifges(6,6) * t431 + Ifges(6,3) * t322) * t495 + (-Ifges(6,1) * t430 + Ifges(6,4) * t431 + Ifges(6,5) * t322) * t500 + (-Ifges(6,4) * t430 + Ifges(6,2) * t431 + Ifges(6,6) * t322) * t502 + (-t430 / 0.2e1 - t246 / 0.2e1) * t91 + (-Ifges(7,1) * t252 - Ifges(7,4) * t251) * t517 - t375 * (-Ifges(5,1) * t274 - Ifges(5,4) * t273) / 0.2e1 + (-Ifges(5,5) * t274 - Ifges(5,6) * t273) * t482 + (-t323 / 0.2e1 + t274 / 0.2e1) * t156 + (-Ifges(5,5) * t323 - Ifges(5,6) * t322) * t481 + (-Ifges(5,1) * t323 - Ifges(5,4) * t322) * t493 + (-Ifges(5,4) * t323 - Ifges(5,2) * t322) * t494 + (-t16 * t423 + t251 * t32 - t428 * t73) * mrSges(7,1) + (t17 * t423 - t252 * t32 + t429 * t73) * mrSges(7,2) + (-t16 * t429 + t17 * t428 - t2 * t251 + t252 * t3) * mrSges(7,3) + (-Ifges(7,5) * t252 - Ifges(7,6) * t251) * t504 + (-Ifges(7,4) * t252 - Ifges(7,2) * t251) * t516 + t362 * t178 / 0.2e1 + t355 * t114 - t312 * t311 + t279 * t166 - t271 * t163 - t241 * t263 - t240 * t264 + t239 * t12 + t199 * t74 + t200 * t75 - pkin(2) * t194 + (-t155 + t543) * (t322 / 0.2e1 - t273 / 0.2e1) + t544 * t86 + t545 * t85 + (t16 * t544 + t17 * t545 + t2 * t88 + t239 * t32 + t3 * t87 + t540 * t73) * m(7) + (t49 * t388 + t383 * t504 + t387 * t510 + t385 * t509 + t43 * t479 + t458 / 0.2e1 + t103 / 0.2e1 + t215 * mrSges(5,2) - t459 / 0.2e1 - t52 * mrSges(5,3) + t42 * t480 + t381 * mrSges(6,3) + (t382 * t496 + t384 * t503 + t386 * t501 + t92 * t389 + t91 * t480 - t365 * t90 / 0.2e1 + t379 * mrSges(6,3)) * qJD(5)) * t331 + (-t171 / 0.2e1 + t145 / 0.2e1) * t57 + (-t172 / 0.2e1 + t144 / 0.2e1) * t58 + t87 * t30 + t88 * t31 + t177 * t478 + (Ifges(4,2) * t366 + t470) * t488 + (Ifges(4,1) * t362 + t469) * t489 + (Ifges(6,5) * t246 + Ifges(6,6) * t245 + Ifges(6,3) * t273) * t496 + (Ifges(7,5) * t144 + Ifges(7,6) * t145 + Ifges(7,3) * t322) * t497 + (Ifges(7,5) * t172 + Ifges(7,6) * t171 + Ifges(7,3) * t273) * t498 + (Ifges(6,1) * t246 + Ifges(6,4) * t245 + Ifges(6,5) * t273) * t501 + (Ifges(6,4) * t246 + Ifges(6,2) * t245 + Ifges(6,6) * t273) * t503 + (Ifges(7,1) * t144 + Ifges(7,4) * t145 + Ifges(7,5) * t322) * t505 + (Ifges(7,1) * t172 + Ifges(7,4) * t171 + Ifges(7,5) * t273) * t506 + (Ifges(7,4) * t144 + Ifges(7,2) * t145 + Ifges(7,6) * t322) * t507 + (Ifges(7,4) * t172 + Ifges(7,2) * t171 + Ifges(7,6) * t273) * t508 - t252 * t518 - t251 * t519; -t49 * t389 + (-mrSges(7,1) * t424 + mrSges(7,2) * t425) * t73 + t427 * t107 + (t218 * t294 + t219 * t295) * mrSges(4,3) + (-t442 / 0.2e1 - t452 / 0.2e1 - t225 * mrSges(5,2) - t156 / 0.2e1 - t449 / 0.2e1 + t98 * mrSges(5,3) - t520) * t525 + (-t16 * t425 + t17 * t424 - t2 * t374 - t3 * t334) * mrSges(7,3) + (Ifges(7,5) * t334 - Ifges(7,6) * t374 + t382) * t504 + t32 * (mrSges(7,1) * t374 + mrSges(7,2) * t334) + (Ifges(7,4) * t334 - Ifges(7,2) * t374) * t516 + (Ifges(7,1) * t334 - Ifges(7,4) * t374) * t517 - t374 * t519 + (-t267 / 0.2e1 + t158 / 0.2e1) * t58 + ((t356 * t53 + t358 * t52) * pkin(3) + t107 * t98 - t108 * t99 - t225 * t476) * m(5) + (-Ifges(7,5) * t267 - Ifges(7,6) * t268) * t497 + (-Ifges(7,1) * t267 - Ifges(7,4) * t268) * t505 + (-Ifges(7,4) * t267 - Ifges(7,2) * t268) * t507 + (-t268 / 0.2e1 + t157 / 0.2e1) * t57 + (-t455 / 0.2e1 - t456 / 0.2e1 - t56 / 0.2e1 - t462 / 0.2e1 - t453 / 0.2e1 - t454 / 0.2e1 - t461 / 0.2e1 + t448 / 0.2e1 + t441 / 0.2e1 - t547 + t451 / 0.2e1 - t89 / 0.2e1) * t375 + t524 + (-t163 * t295 + t166 * t356 + t167 * t358) * pkin(3) + (-Ifges(7,5) * t158 - Ifges(7,6) * t157) * t498 + (-Ifges(7,1) * t158 - Ifges(7,4) * t157) * t506 + (-Ifges(7,4) * t158 - Ifges(7,2) * t157) * t508 - (-Ifges(4,2) * t295 + t214 + t290) * t294 / 0.2e1 + t354 * t59 + t336 * t12 - t280 * (mrSges(4,1) * t295 + mrSges(4,2) * t294) + t260 * t30 + t261 * t31 - t218 * t263 + t219 * t264 - t108 * t208 - t161 * mrSges(4,2) + t162 * mrSges(4,1) - t63 * t148 - t62 * t149 + (-t361 * t74 + t365 * t75 + m(6) * t523 + (-m(6) * t380 - t361 * t148 - t365 * t149) * qJD(5)) * t353 + t523 * mrSges(6,3) - t295 * (Ifges(4,1) * t294 - t445) / 0.2e1 + (pkin(5) * t361 * t526 + t520) * qJD(5) + t527 * t86 + t528 * t85 + (t16 * t527 + t17 * t528 + t2 * t261 + t260 * t3 + t32 * t336 - t73 * t78) * m(7) + (-t107 * t92 + t354 * t49 - t54 * t62 - t55 * t63) * m(6) + t52 * mrSges(5,1) - t53 * mrSges(5,2) - t78 * t69 + t42 * t479 + (Ifges(4,5) * t294 - Ifges(4,6) * t295) * t482 + t213 * t485 + t384 * t509 + t386 * t510 + t361 * t515 + t334 * t518; -t374 * t30 + t334 * t31 + t361 * t75 + t365 * t74 + t424 * t86 + t425 * t85 + t378 * qJD(5) + (-t208 - t378) * t525 + (-t69 + t427) * t375 + t114 + (t16 * t424 + t17 * t425 + t2 * t334 - t3 * t374 - t375 * t73) * m(7) + (-t224 * t379 - t375 * t92 - t381) * m(6) + (t375 * t98 - t525 * t99 + t215) * m(5); -m(7) * (t16 * t18 + t17 * t19) + t133 * t535 + t533 + (t201 * t54 + t202 * t55) * mrSges(6,3) + t534 + (-Ifges(6,2) * t202 + t198 + t91) * t503 - t92 * (mrSges(6,1) * t202 + mrSges(6,2) * t201) - t54 * t148 + t55 * t149 + (t364 * t30 + t360 * t31 + m(7) * (t2 * t360 + t3 * t364) - t526 * t202 + (-t360 * t86 + t364 * t85 + m(7) * (-t16 * t360 + t17 * t364)) * qJD(6)) * pkin(5) - t19 * t85 - t18 * t86 + (Ifges(6,5) * t201 - Ifges(6,6) * t202) * t496 + t90 * t500 + (Ifges(6,1) * t201 - t468) * t501 + t548; -t16 * t85 + t17 * t86 + t57 * t505 + t546 + t548 + t8;];
tauc  = t1(:);
