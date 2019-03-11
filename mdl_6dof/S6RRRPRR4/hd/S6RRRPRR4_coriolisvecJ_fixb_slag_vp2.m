% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:15:35
% EndTime: 2019-03-09 18:16:05
% DurationCPUTime: 16.87s
% Computational Cost: add. (23655->730), mult. (56821->1000), div. (0->0), fcn. (42963->10), ass. (0->343)
t393 = sin(qJ(3));
t397 = cos(qJ(3));
t398 = cos(qJ(2));
t427 = qJD(1) * t398;
t394 = sin(qJ(2));
t428 = qJD(1) * t394;
t344 = t393 * t428 - t397 * t427;
t389 = sin(pkin(11));
t390 = cos(pkin(11));
t392 = sin(qJ(5));
t396 = cos(qJ(5));
t361 = t389 * t396 + t390 * t392;
t261 = t361 * t344;
t341 = t361 * qJD(5);
t546 = t261 + t341;
t431 = t390 * t396;
t402 = t389 * t392 - t431;
t262 = t402 * t344;
t340 = t402 * qJD(5);
t545 = t262 + t340;
t381 = pkin(2) * t393 + qJ(4);
t350 = (-pkin(9) - t381) * t389;
t387 = t390 * pkin(9);
t433 = t381 * t390;
t351 = t387 + t433;
t288 = t392 * t350 + t396 * t351;
t424 = qJD(3) * t397;
t377 = pkin(2) * t424 + qJD(4);
t222 = -qJD(5) * t288 - t361 * t377;
t363 = t393 * t398 + t397 * t394;
t345 = t363 * qJD(1);
t295 = pkin(3) * t345 + qJ(4) * t344;
t270 = pkin(2) * t428 + t295;
t510 = -pkin(8) - pkin(7);
t376 = t510 * t398;
t366 = qJD(1) * t376;
t346 = t393 * t366;
t375 = t510 * t394;
t365 = qJD(1) * t375;
t306 = t365 * t397 + t346;
t203 = t390 * t270 - t306 * t389;
t438 = t344 * t390;
t407 = t345 * pkin(4) + pkin(9) * t438;
t149 = t203 + t407;
t204 = t389 * t270 + t390 * t306;
t439 = t344 * t389;
t420 = pkin(9) * t439;
t169 = t420 + t204;
t86 = t396 * t149 - t169 * t392;
t567 = -t86 + t222;
t421 = qJD(5) * t396;
t435 = t377 * t389;
t221 = t350 * t421 + t377 * t431 + (-qJD(5) * t351 - t435) * t392;
t87 = t392 * t149 + t396 * t169;
t566 = -t87 + t221;
t565 = t546 * pkin(10);
t564 = -t345 * pkin(5) + pkin(10) * t545;
t563 = -t565 + t566;
t562 = t564 + t567;
t352 = qJD(2) * pkin(2) + t365;
t301 = t352 * t397 + t346;
t208 = t390 * t295 - t301 * t389;
t154 = t208 + t407;
t209 = t389 * t295 + t390 * t301;
t173 = t420 + t209;
t370 = (-pkin(9) - qJ(4)) * t389;
t371 = qJ(4) * t390 + t387;
t311 = t392 * t370 + t396 * t371;
t529 = -t361 * qJD(4) - qJD(5) * t311 - t396 * t154 + t173 * t392;
t423 = qJD(4) * t390;
t528 = t370 * t421 + (-t173 + t423) * t396 + (-qJD(4) * t389 - qJD(5) * t371 - t154) * t392;
t388 = qJD(2) + qJD(3);
t314 = t345 * t390 + t388 * t389;
t409 = -t345 * t389 + t390 * t388;
t234 = t314 * t396 + t392 * t409;
t391 = sin(qJ(6));
t395 = cos(qJ(6));
t541 = -t314 * t392 + t396 * t409;
t561 = -t234 * t391 + t395 * t541;
t136 = t234 * t395 + t391 * t541;
t489 = -t344 / 0.2e1;
t559 = t528 - t565;
t558 = t529 + t564;
t555 = t546 * pkin(5);
t269 = -pkin(3) * t388 + qJD(4) - t301;
t385 = -pkin(2) * t398 - pkin(1);
t374 = qJD(1) * t385;
t467 = Ifges(5,2) * t389;
t470 = Ifges(5,4) * t390;
t404 = -t467 + t470;
t471 = Ifges(5,4) * t389;
t405 = Ifges(5,1) * t390 - t471;
t406 = mrSges(5,1) * t389 + mrSges(5,2) * t390;
t455 = t388 * Ifges(4,5);
t484 = t390 / 0.2e1;
t485 = -t389 / 0.2e1;
t554 = (Ifges(5,1) * t314 + Ifges(5,4) * t409 + Ifges(5,5) * t344) * t484 + (Ifges(5,4) * t314 + Ifges(5,2) * t409 + Ifges(5,6) * t344) * t485 + t269 * t406 + t374 * mrSges(4,2) + t314 * t405 / 0.2e1 + t409 * t404 / 0.2e1 + t455 / 0.2e1;
t219 = -pkin(4) * t409 + t269;
t130 = -pkin(5) * t541 + t219;
t263 = t344 * pkin(3) - t345 * qJ(4) + t374;
t430 = t397 * t366;
t302 = t393 * t352 - t430;
t277 = qJ(4) * t388 + t302;
t184 = t390 * t263 - t277 * t389;
t129 = pkin(4) * t344 - pkin(9) * t314 + t184;
t185 = t389 * t263 + t390 * t277;
t144 = pkin(9) * t409 + t185;
t77 = t129 * t392 + t144 * t396;
t54 = pkin(10) * t541 + t77;
t453 = t391 * t54;
t339 = qJD(5) + t344;
t76 = t396 * t129 - t144 * t392;
t53 = -pkin(10) * t234 + t76;
t47 = pkin(5) * t339 + t53;
t20 = t395 * t47 - t453;
t452 = t395 * t54;
t21 = t391 * t47 + t452;
t308 = t388 * t363;
t294 = t308 * qJD(1);
t362 = t393 * t394 - t397 * t398;
t307 = t388 * t362;
t293 = t307 * qJD(1);
t121 = qJD(5) * t541 + t293 * t402;
t122 = -qJD(5) * t234 + t293 * t361;
t39 = qJD(6) * t561 + t121 * t395 + t122 * t391;
t40 = -qJD(6) * t136 - t121 * t391 + t122 * t395;
t419 = Ifges(7,5) * t39 + Ifges(7,6) * t40 + Ifges(7,3) * t294;
t468 = Ifges(7,4) * t136;
t329 = qJD(6) + t339;
t493 = -t329 / 0.2e1;
t500 = -t136 / 0.2e1;
t422 = qJD(5) * t392;
t426 = qJD(2) * t394;
t418 = pkin(2) * t426;
t543 = qJD(1) * t418;
t167 = pkin(3) * t294 + qJ(4) * t293 - qJD(4) * t345 + t543;
t416 = qJD(2) * t510;
t408 = qJD(1) * t416;
t355 = t394 * t408;
t401 = t398 * t408;
t425 = qJD(3) * t393;
t217 = t352 * t424 + t397 * t355 + t366 * t425 + t393 * t401;
t205 = qJD(4) * t388 + t217;
t103 = t390 * t167 - t205 * t389;
t442 = t293 * t390;
t73 = pkin(4) * t294 + pkin(9) * t442 + t103;
t104 = t389 * t167 + t390 * t205;
t443 = t293 * t389;
t85 = pkin(9) * t443 + t104;
t15 = t129 * t421 - t144 * t422 + t392 * t73 + t396 * t85;
t10 = pkin(10) * t122 + t15;
t16 = -qJD(5) * t77 - t392 * t85 + t396 * t73;
t7 = pkin(5) * t294 - pkin(10) * t121 + t16;
t3 = qJD(6) * t20 + t10 * t395 + t391 * t7;
t4 = -qJD(6) * t21 - t10 * t391 + t395 * t7;
t537 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t553 = t419 + t537 + (Ifges(7,5) * t561 - Ifges(7,6) * t136) * t493 + (t136 * t21 + t20 * t561) * mrSges(7,3) - t130 * (mrSges(7,1) * t136 + mrSges(7,2) * t561) + (Ifges(7,1) * t561 - t468) * t500;
t460 = t314 * Ifges(5,5);
t461 = t409 * Ifges(5,6);
t472 = Ifges(4,4) * t345;
t552 = t185 * mrSges(5,2) + t21 * mrSges(7,2) + t77 * mrSges(6,2) + Ifges(4,2) * t489 + t388 * Ifges(4,6) + t472 / 0.2e1 - t184 * mrSges(5,1) - t20 * mrSges(7,1) - t374 * mrSges(4,1) - t76 * mrSges(6,1) - t460 / 0.2e1 - t461 / 0.2e1;
t287 = t396 * t350 - t351 * t392;
t477 = pkin(10) * t361;
t245 = t287 - t477;
t354 = t402 * pkin(10);
t246 = -t354 + t288;
t157 = t245 * t391 + t246 * t395;
t550 = -qJD(6) * t157 - t563 * t391 + t562 * t395;
t156 = t245 * t395 - t246 * t391;
t549 = qJD(6) * t156 + t562 * t391 + t563 * t395;
t177 = t261 * t395 - t262 * t391;
t304 = t361 * t395 - t391 * t402;
t211 = -qJD(6) * t304 + t340 * t391 - t341 * t395;
t548 = t177 - t211;
t178 = t261 * t391 + t262 * t395;
t303 = -t361 * t391 - t395 * t402;
t210 = qJD(6) * t303 - t340 * t395 - t341 * t391;
t547 = t178 - t210;
t305 = t365 * t393 - t430;
t323 = pkin(4) * t439;
t248 = t305 - t323;
t544 = pkin(2) * t425 - t248 + t555;
t131 = Ifges(7,4) * t561;
t542 = -Ifges(7,2) * t136 + t131;
t518 = t39 / 0.2e1;
t517 = t40 / 0.2e1;
t506 = t121 / 0.2e1;
t505 = t122 / 0.2e1;
t494 = t294 / 0.2e1;
t310 = t396 * t370 - t371 * t392;
t266 = t310 - t477;
t267 = -t354 + t311;
t182 = t266 * t391 + t267 * t395;
t536 = -qJD(6) * t182 - t391 * t559 + t395 * t558;
t181 = t266 * t395 - t267 * t391;
t535 = qJD(6) * t181 + t391 * t558 + t395 * t559;
t534 = mrSges(4,2) * t345;
t242 = -t323 + t302;
t524 = -t242 + t555;
t300 = t362 * pkin(3) - t363 * qJ(4) + t385;
t317 = t375 * t393 - t376 * t397;
t223 = t390 * t300 - t317 * t389;
t172 = pkin(4) * t362 - t363 * t387 + t223;
t224 = t389 * t300 + t390 * t317;
t437 = t363 * t389;
t193 = -pkin(9) * t437 + t224;
t96 = t392 * t172 + t396 * t193;
t523 = t397 * t375 + t376 * t393;
t522 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t521 = m(4) / 0.2e1;
t520 = Ifges(7,4) * t518 + Ifges(7,2) * t517 + Ifges(7,6) * t494;
t519 = Ifges(7,1) * t518 + Ifges(7,4) * t517 + Ifges(7,5) * t494;
t516 = Ifges(6,4) * t506 + Ifges(6,2) * t505 + Ifges(6,6) * t494;
t515 = Ifges(6,1) * t506 + Ifges(6,4) * t505 + Ifges(6,5) * t494;
t66 = Ifges(7,2) * t561 + Ifges(7,6) * t329 + t468;
t514 = -t66 / 0.2e1;
t513 = t66 / 0.2e1;
t67 = Ifges(7,1) * t136 + Ifges(7,5) * t329 + t131;
t512 = -t67 / 0.2e1;
t511 = t67 / 0.2e1;
t509 = pkin(1) * mrSges(3,1);
t508 = pkin(1) * mrSges(3,2);
t469 = Ifges(6,4) * t234;
t124 = Ifges(6,2) * t541 + Ifges(6,6) * t339 + t469;
t504 = t124 / 0.2e1;
t230 = Ifges(6,4) * t541;
t125 = Ifges(6,1) * t234 + Ifges(6,5) * t339 + t230;
t503 = t125 / 0.2e1;
t502 = -t561 / 0.2e1;
t501 = t561 / 0.2e1;
t499 = t136 / 0.2e1;
t498 = -t541 / 0.2e1;
t497 = t541 / 0.2e1;
t496 = -t234 / 0.2e1;
t495 = t234 / 0.2e1;
t492 = t329 / 0.2e1;
t491 = -t339 / 0.2e1;
t490 = t339 / 0.2e1;
t488 = t344 / 0.2e1;
t487 = -t345 / 0.2e1;
t482 = m(4) * t374;
t481 = pkin(2) * t397;
t474 = mrSges(4,3) * t344;
t473 = Ifges(3,4) * t394;
t466 = t103 * mrSges(5,3);
t465 = t561 * Ifges(7,6);
t464 = t136 * Ifges(7,5);
t463 = t541 * Ifges(6,6);
t462 = t234 * Ifges(6,5);
t459 = t329 * Ifges(7,3);
t458 = t339 * Ifges(6,3);
t457 = t345 * mrSges(4,3);
t456 = t345 * Ifges(4,1);
t451 = Ifges(3,5) * qJD(2);
t450 = Ifges(3,6) * qJD(2);
t449 = qJD(2) * mrSges(3,1);
t448 = qJD(2) * mrSges(3,2);
t447 = t103 * t389;
t446 = t104 * t390;
t218 = qJD(3) * t302 + t355 * t393 - t397 * t401;
t445 = t218 * t523;
t440 = t307 * t389;
t434 = t377 * t390;
t206 = -mrSges(5,2) * t294 + mrSges(5,3) * t443;
t432 = t390 * t206;
t194 = pkin(3) * t308 + qJ(4) * t307 - qJD(4) * t363 + t418;
t367 = t394 * t416;
t368 = t398 * t416;
t236 = qJD(3) * t523 + t367 * t397 + t368 * t393;
t117 = t389 * t194 + t390 * t236;
t202 = -mrSges(5,1) * t443 - mrSges(5,2) * t442;
t429 = mrSges(4,1) * t388 + mrSges(5,1) * t409 - mrSges(5,2) * t314 - t457;
t417 = Ifges(6,5) * t121 + Ifges(6,6) * t122 + Ifges(6,3) * t294;
t382 = -pkin(4) * t390 - pkin(3);
t415 = t451 / 0.2e1;
t414 = -t450 / 0.2e1;
t13 = -t40 * mrSges(7,1) + t39 * mrSges(7,2);
t58 = -t122 * mrSges(6,1) + t121 * mrSges(6,2);
t95 = t396 * t172 - t193 * t392;
t116 = t390 * t194 - t236 * t389;
t265 = pkin(4) * t437 - t523;
t403 = Ifges(5,5) * t390 - Ifges(5,6) * t389;
t281 = t402 * t363;
t79 = pkin(5) * t362 + pkin(10) * t281 + t95;
t280 = t361 * t363;
t81 = -pkin(10) * t280 + t96;
t32 = -t391 * t81 + t395 * t79;
t33 = t391 * t79 + t395 * t81;
t195 = -t280 * t395 + t281 * t391;
t196 = -t280 * t391 - t281 * t395;
t322 = pkin(5) * t402 + t382;
t105 = pkin(9) * t440 + t117;
t92 = pkin(4) * t308 + t307 * t387 + t116;
t24 = t396 * t105 + t172 * t421 - t193 * t422 + t392 * t92;
t25 = -qJD(5) * t96 - t105 * t392 + t396 * t92;
t237 = qJD(3) * t317 + t367 * t393 - t397 * t368;
t155 = -pkin(4) * t443 + t218;
t179 = -pkin(4) * t440 + t237;
t123 = t458 + t462 + t463;
t145 = t294 * Ifges(5,6) - t293 * t404;
t146 = t294 * Ifges(5,5) - t293 * t405;
t213 = t344 * Ifges(5,3) + t460 + t461;
t336 = Ifges(4,4) * t344;
t275 = -t336 + t455 + t456;
t65 = t459 + t464 + t465;
t78 = -pkin(5) * t122 + t155;
t399 = (Ifges(5,5) * t389 + Ifges(6,5) * t361 + Ifges(7,5) * t304 + Ifges(5,6) * t390 - Ifges(6,6) * t402 + Ifges(7,6) * t303) * t494 - t402 * t516 + (Ifges(6,4) * t361 - Ifges(6,2) * t402) * t505 + (Ifges(6,1) * t361 - Ifges(6,4) * t402) * t506 + t155 * (mrSges(6,1) * t402 + mrSges(6,2) * t361) + (-t403 * t489 + t554) * t344 + (mrSges(6,1) * t546 - mrSges(6,2) * t545) * t219 + (-t15 * t402 - t16 * t361 + t545 * t76 - t546 * t77) * mrSges(6,3) + (Ifges(6,5) * t262 + Ifges(6,6) * t261) * t491 + (Ifges(7,5) * t178 + Ifges(7,6) * t177) * t493 + (Ifges(6,1) * t262 + Ifges(6,4) * t261) * t496 + (Ifges(7,4) * t178 + Ifges(7,2) * t177) * t502 - t301 * t474 + (-Ifges(4,1) * t344 + t123 + t213 - t472 + t65) * t487 + (-t184 * t438 - t185 * t439 + t446) * mrSges(5,3) + t145 * t484 + (-mrSges(5,1) * t390 + mrSges(5,2) * t389 - mrSges(4,1)) * t218 + (mrSges(7,1) * t548 - mrSges(7,2) * t547) * t130 + (t20 * t547 - t21 * t548 + t3 * t303 - t304 * t4) * mrSges(7,3) + (Ifges(6,5) * t496 + Ifges(7,5) * t500 - Ifges(4,2) * t488 + Ifges(6,6) * t498 + Ifges(7,6) * t502 + Ifges(5,3) * t489 + Ifges(6,3) * t491 + Ifges(7,3) * t493 + t552) * t345 + (Ifges(7,5) * t210 + Ifges(7,6) * t211) * t492 + (Ifges(7,1) * t210 + Ifges(7,4) * t211) * t499 + (Ifges(7,4) * t210 + Ifges(7,2) * t211) * t501 - t340 * t503 + (Ifges(6,4) * t262 + Ifges(6,2) * t261) * t498 + (-t336 + t275) * t488 + (Ifges(7,1) * t178 + Ifges(7,4) * t177) * t500 - t341 * t504 + t210 * t511 + t178 * t512 + t211 * t513 + t177 * t514 + t361 * t515 + (Ifges(7,4) * t304 + Ifges(7,2) * t303) * t517 + (Ifges(7,1) * t304 + Ifges(7,4) * t303) * t518 + t304 * t519 + t303 * t520 + t389 * t146 / 0.2e1 + t78 * (-mrSges(7,1) * t303 + mrSges(7,2) * t304) - Ifges(4,6) * t294 - Ifges(4,5) * t293 + (-Ifges(6,5) * t340 - Ifges(6,6) * t341) * t490 + (-Ifges(6,1) * t340 - Ifges(6,4) * t341) * t495 + (-Ifges(6,4) * t340 - Ifges(6,2) * t341) * t497 - (Ifges(5,1) * t389 + t470) * t442 / 0.2e1 + (Ifges(5,2) * t390 + t471) * t443 / 0.2e1 - t217 * mrSges(4,2) - t261 * t124 / 0.2e1 - t262 * t125 / 0.2e1;
t386 = Ifges(3,4) * t427;
t384 = -pkin(3) - t481;
t373 = mrSges(3,3) * t427 - t448;
t372 = -mrSges(3,3) * t428 + t449;
t369 = t382 - t481;
t343 = Ifges(3,1) * t428 + t386 + t451;
t342 = t450 + (Ifges(3,2) * t398 + t473) * qJD(1);
t320 = -mrSges(4,2) * t388 - t474;
t318 = t322 - t481;
t297 = mrSges(4,1) * t344 + t534;
t255 = mrSges(5,1) * t344 - mrSges(5,3) * t314;
t254 = -mrSges(5,2) * t344 + mrSges(5,3) * t409;
t207 = mrSges(5,1) * t294 + mrSges(5,3) * t442;
t200 = mrSges(6,1) * t339 - mrSges(6,3) * t234;
t199 = -mrSges(6,2) * t339 + mrSges(6,3) * t541;
t198 = pkin(5) * t280 + t265;
t153 = t307 * t361 + t340 * t363;
t152 = t307 * t402 - t341 * t363;
t139 = -mrSges(6,1) * t541 + mrSges(6,2) * t234;
t115 = mrSges(7,1) * t329 - mrSges(7,3) * t136;
t114 = -mrSges(7,2) * t329 + mrSges(7,3) * t561;
t102 = -mrSges(6,2) * t294 + mrSges(6,3) * t122;
t101 = mrSges(6,1) * t294 - mrSges(6,3) * t121;
t93 = -pkin(5) * t153 + t179;
t75 = -mrSges(7,1) * t561 + mrSges(7,2) * t136;
t56 = -qJD(6) * t196 - t152 * t391 + t153 * t395;
t55 = qJD(6) * t195 + t152 * t395 + t153 * t391;
t35 = -mrSges(7,2) * t294 + mrSges(7,3) * t40;
t34 = mrSges(7,1) * t294 - mrSges(7,3) * t39;
t23 = t395 * t53 - t453;
t22 = -t391 * t53 - t452;
t19 = pkin(10) * t153 + t24;
t18 = pkin(5) * t308 - pkin(10) * t152 + t25;
t6 = -qJD(6) * t33 + t18 * t395 - t19 * t391;
t5 = qJD(6) * t32 + t18 * t391 + t19 * t395;
t1 = [-(t385 * mrSges(4,2) + (Ifges(4,1) + Ifges(5,1) * t390 ^ 2 / 0.2e1 + (-t470 + t467 / 0.2e1) * t389) * t363) * t293 + (-t307 * t489 + t308 * t487) * Ifges(4,4) + m(5) * (t103 * t223 + t104 * t224 + t116 * t184 + t117 * t185 + t237 * t269 - t445) + (Ifges(7,1) * t196 + Ifges(7,4) * t195) * t518 - (t456 / 0.2e1 + t403 * t488 + t275 / 0.2e1 + (-t184 * t390 - t185 * t389) * mrSges(5,3) + t554) * t307 + (-t363 * Ifges(4,4) + (Ifges(5,3) + Ifges(4,2)) * t362 - t317 * mrSges(4,3) + t385 * mrSges(4,1)) * t294 + (t419 + t417) * t362 / 0.2e1 - t523 * t202 + (mrSges(4,1) * t543 + mrSges(5,1) * t103 - mrSges(5,2) * t104 - mrSges(4,3) * t217 + Ifges(6,5) * t506 + Ifges(7,5) * t518 + Ifges(6,6) * t505 + Ifges(7,6) * t517 + (Ifges(6,3) + Ifges(7,3)) * t494 - (-Ifges(4,4) + t403) * t293 + t522 + t537) * t362 + (t218 * t363 + t293 * t523 + t301 * t307 - t302 * t308) * mrSges(4,3) + m(6) * (t15 * t96 + t155 * t265 + t16 * t95 + t179 * t219 + t24 * t77 + t25 * t76) + m(7) * (t130 * t93 + t198 * t78 + t20 * t6 + t21 * t5 + t3 * t33 + t32 * t4) + (Ifges(7,4) * t196 + Ifges(7,2) * t195) * t517 + (t123 / 0.2e1 + t65 / 0.2e1 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t344 + t464 / 0.2e1 + t465 / 0.2e1 + t458 / 0.2e1 + t459 / 0.2e1 + t462 / 0.2e1 + t463 / 0.2e1 + t213 / 0.2e1 - t552) * t308 + (Ifges(6,5) * t152 + Ifges(6,6) * t153) * t490 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t492 + (Ifges(6,1) * t152 + Ifges(6,4) * t153) * t495 + (Ifges(6,4) * t152 + Ifges(6,2) * t153) * t497 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t499 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t501 + (t195 * t3 - t196 * t4 - t20 * t55 + t21 * t56) * mrSges(7,3) + t152 * t503 + t153 * t504 + t55 * t511 + t56 * t513 - t281 * t515 - t280 * t516 + t196 * t519 + t195 * t520 + t236 * t320 + m(4) * (t217 * t317 + t236 * t302 - t237 * t301 - t445) + (-Ifges(6,4) * t281 - Ifges(6,2) * t280) * t505 + (-Ifges(6,5) * t281 + Ifges(7,5) * t196 - Ifges(6,6) * t280 + Ifges(7,6) * t195 + t363 * t403) * t494 + (-t15 * t280 - t152 * t76 + t153 * t77 + t16 * t281) * mrSges(6,3) + (-Ifges(6,1) * t281 - Ifges(6,4) * t280) * t506 + t155 * (mrSges(6,1) * t280 - mrSges(6,2) * t281) + (t146 * t484 + t145 * t485 + t218 * t406 + (-t103 * t390 - t104 * t389) * mrSges(5,3)) * t363 + (-pkin(7) * t372 + t343 / 0.2e1 + t415 + (-0.2e1 * t508 + 0.3e1 / 0.2e1 * Ifges(3,4) * t398) * qJD(1)) * t398 * qJD(2) + t32 * t34 + t33 * t35 + (-pkin(7) * t373 - t342 / 0.2e1 + t414 + (-0.2e1 * t509 - 0.3e1 / 0.2e1 * t473 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t398) * qJD(1) + (t297 + 0.2e1 * t482 + t534) * pkin(2)) * t426 - t429 * t237 + t93 * t75 + t95 * t101 + t96 * t102 + t5 * t114 + t6 * t115 + t130 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t179 * t139 + t78 * (-mrSges(7,1) * t195 + mrSges(7,2) * t196) + t198 * t13 + t24 * t199 + t25 * t200 + t219 * (-mrSges(6,1) * t153 + mrSges(6,2) * t152) + t223 * t207 + t224 * t206 + t117 * t254 + t116 * t255 + t265 * t58; t381 * t432 + t549 * t114 + (t130 * t544 + t156 * t4 + t157 * t3 + t20 * t550 + t21 * t549 + t318 * t78) * m(7) + t550 * t115 + 0.2e1 * ((t217 * t393 - t218 * t397) * t521 + ((-t301 * t393 + t302 * t397) * t521 + (m(5) * t269 + m(6) * t219) * t393 / 0.2e1) * qJD(3)) * pkin(2) + t544 * t75 + ((t293 * t397 - t294 * t393) * mrSges(4,3) + (t320 * t397 + (t139 - t429) * t393) * qJD(3)) * pkin(2) + t399 - m(5) * (t184 * t203 + t185 * t204 + t269 * t305) - m(6) * (t219 * t248 + t76 * t86 + t77 * t87) + m(6) * (t15 * t288 + t155 * t369 + t16 * t287 + t221 * t77 + t222 * t76) + t384 * t202 + t369 * t58 - m(4) * (-t301 * t305 + t302 * t306) + t318 * t13 - t306 * t320 + t287 * t101 + t288 * t102 + t302 * t457 + t566 * t199 + t567 * t200 + t429 * t305 + (-t204 + t434) * t254 + m(5) * (t104 * t433 - t184 * t435 + t185 * t434 + t218 * t384 - t381 * t447) + (-t381 * t207 - t377 * t255 - t466) * t389 + ((t415 - t386 / 0.2e1 - t343 / 0.2e1 + qJD(1) * t508 + (t372 - t449) * pkin(7)) * t398 + (t414 + t342 / 0.2e1 + (t509 + t473 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t398) * qJD(1) + (t373 + t448) * pkin(7) + (-t297 - t482) * pkin(2)) * t394) * qJD(1) + t156 * t34 + t157 * t35 - t248 * t139 - t203 * t255; t535 * t114 + qJ(4) * t432 + t399 + t529 * t200 + t528 * t199 + t382 * t58 - t301 * t320 + t322 * t13 + t310 * t101 + t311 * t102 + (-t209 + t423) * t254 + t536 * t115 + (t429 + t457) * t302 + (-qJ(4) * t207 - qJD(4) * t255 - t466) * t389 + t524 * t75 + t181 * t34 + t182 * t35 - pkin(3) * t202 - t242 * t139 - t208 * t255 + (t130 * t524 + t181 * t4 + t182 * t3 + t20 * t536 + t21 * t535 + t322 * t78) * m(7) + (t15 * t311 + t155 * t382 + t16 * t310 - t219 * t242 + t528 * t77 + t529 * t76) * m(6) + (-t184 * t208 - t185 * t209 - t269 * t302 - pkin(3) * t218 + (-t184 * t389 + t185 * t390) * qJD(4) + (t446 - t447) * qJ(4)) * m(5); -t561 * t114 + t136 * t115 - t541 * t199 + t234 * t200 - t409 * t254 + t314 * t255 + t13 + t202 + t58 + (t136 * t20 - t21 * t561 + t78) * m(7) + (t234 * t76 - t541 * t77 + t155) * m(6) + (t184 * t314 - t185 * t409 + t218) * m(5); t522 - t136 * t514 + t561 * t512 + (-Ifges(6,2) * t234 + t125 + t230) * t498 - m(7) * (t20 * t22 + t21 * t23) + t417 + (Ifges(6,5) * t541 - Ifges(6,6) * t234) * t491 + t124 * t495 + (Ifges(6,1) * t541 - t469) * t496 + t542 * t502 + (-t234 * t75 + t395 * t34 + t391 * t35 + (t114 * t395 - t115 * t391) * qJD(6) + (-t130 * t234 + t3 * t391 + t395 * t4 + (-t20 * t391 + t21 * t395) * qJD(6)) * m(7)) * pkin(5) - t23 * t114 - t22 * t115 + (t234 * t77 + t541 * t76) * mrSges(6,3) - t76 * t199 + t77 * t200 - t219 * (mrSges(6,1) * t234 + mrSges(6,2) * t541) + t553; t66 * t499 - t20 * t114 + t21 * t115 + (t542 + t67) * t502 + t553;];
tauc  = t1(:);
