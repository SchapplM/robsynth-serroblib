% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR14_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:24
% EndTime: 2019-03-09 20:13:28
% DurationCPUTime: 34.12s
% Computational Cost: add. (18310->903), mult. (45925->1217), div. (0->0), fcn. (35078->10), ass. (0->400)
t334 = cos(qJ(2));
t326 = sin(pkin(6));
t416 = qJD(1) * t326;
t396 = t334 * t416;
t303 = -qJD(3) + t396;
t435 = cos(pkin(6));
t384 = t435 * qJD(1);
t376 = pkin(1) * t384;
t330 = sin(qJ(2));
t397 = t330 * t416;
t261 = -pkin(8) * t397 + t334 * t376;
t348 = (pkin(2) * t330 - pkin(9) * t334) * t326;
t262 = qJD(1) * t348;
t329 = sin(qJ(3));
t333 = cos(qJ(3));
t185 = -t329 * t261 + t262 * t333;
t325 = t333 * pkin(4);
t413 = qJD(3) * t333;
t494 = pkin(4) + pkin(9);
t495 = pkin(3) + pkin(10);
t580 = -(t325 * t334 - t330 * t495) * t416 + t185 + t494 * t413;
t414 = qJD(3) * t329;
t381 = pkin(3) * t414 - qJD(4) * t329;
t264 = pkin(8) * t396 + t330 * t376;
t379 = t329 * t396;
t401 = pkin(3) * t379 + t264;
t579 = -t381 + t401 + t303 * (pkin(10) * t329 - qJ(4) * t333);
t328 = sin(qJ(5));
t332 = cos(qJ(5));
t424 = t332 * t334;
t232 = (-t328 * t330 + t329 * t424) * t416;
t412 = qJD(5) * t328;
t578 = -t333 * t412 + t232;
t386 = -qJ(4) * t329 - pkin(2);
t278 = -t333 * t495 + t386;
t304 = t494 * t329;
t411 = qJD(5) * t332;
t569 = -t278 * t412 + t304 * t411 + t328 * t580 - t579 * t332;
t577 = t579 * t328 + t332 * t580;
t428 = t328 * t334;
t233 = (t329 * t428 + t330 * t332) * t416;
t284 = t328 * t304;
t378 = t333 * t396;
t385 = pkin(11) * t333 - t278;
t429 = t328 * t329;
t576 = -pkin(5) * t378 + pkin(11) * t233 + (pkin(5) * t333 - pkin(11) * t429) * qJD(3) + (t332 * t385 - t284) * qJD(5) + t577;
t575 = t569 + (t332 * t414 - t578) * pkin(11);
t351 = t384 + qJD(2);
t344 = t333 * t351;
t245 = t329 * t397 - t344;
t472 = -t245 / 0.2e1;
t574 = Ifges(5,6) * t472;
t246 = t329 * t351 + t333 * t397;
t469 = t246 / 0.2e1;
t467 = -t303 / 0.2e1;
t573 = t303 / 0.2e1;
t541 = Ifges(5,1) + Ifges(4,3);
t572 = Ifges(5,4) - Ifges(4,5);
t539 = Ifges(5,5) - Ifges(4,6);
t285 = t332 * t304;
t187 = pkin(5) * t329 + t328 * t385 + t285;
t208 = t332 * t278 + t284;
t425 = t332 * t333;
t193 = -pkin(11) * t425 + t208;
t327 = sin(qJ(6));
t331 = cos(qJ(6));
t109 = t187 * t331 - t193 * t327;
t571 = qJD(6) * t109 + t327 * t576 + t575 * t331;
t110 = t187 * t327 + t193 * t331;
t570 = -qJD(6) * t110 - t575 * t327 + t331 * t576;
t226 = pkin(9) * t351 + t264;
t258 = (-pkin(2) * t334 - pkin(9) * t330 - pkin(1)) * t326;
t237 = qJD(1) * t258;
t162 = t333 * t226 + t329 * t237;
t127 = -pkin(4) * t245 + t162;
t121 = t332 * t127;
t433 = qJ(4) * t245;
t144 = t246 * t495 + t433;
t458 = pkin(11) + t495;
t462 = pkin(11) * t246;
t568 = t458 * t412 + pkin(5) * t245 - t121 - (-t144 - t462) * t328;
t294 = t458 * t332;
t64 = t328 * t127 + t332 * t144;
t567 = qJD(5) * t294 + t332 * t462 + t64;
t566 = -qJD(5) * t208 + t577;
t186 = t333 * t261 + t329 * t262;
t173 = -qJ(4) * t397 - t186;
t152 = -pkin(4) * t379 - t173;
t399 = -pkin(5) * t332 - pkin(4);
t565 = -t152 + (-pkin(9) + t399) * t414 + t578 * pkin(5);
t564 = -t494 * t414 - t152;
t240 = qJD(5) + t246;
t192 = t245 * t328 - t303 * t332;
t161 = t226 * t329 - t333 * t237;
t349 = pkin(4) * t246 + t161;
t553 = qJD(4) + t349;
t94 = t303 * t495 + t553;
t225 = -pkin(2) * t351 - t261;
t339 = -t246 * qJ(4) + t225;
t99 = t245 * t495 + t339;
t46 = -t328 * t99 + t332 * t94;
t38 = -pkin(11) * t192 + t46;
t36 = pkin(5) * t240 + t38;
t191 = t245 * t332 + t303 * t328;
t47 = t328 * t94 + t332 * t99;
t39 = pkin(11) * t191 + t47;
t438 = t327 * t39;
t15 = t331 * t36 - t438;
t436 = t331 * t39;
t16 = t327 * t36 + t436;
t563 = t15 * mrSges(7,1) - t16 * mrSges(7,2);
t515 = -qJD(4) - t161;
t143 = pkin(3) * t303 - t515;
t562 = t143 * mrSges(5,1) + t161 * mrSges(4,3);
t409 = qJD(1) * qJD(2);
t389 = t334 * t409;
t374 = t326 * t389;
t201 = qJD(3) * t246 + t329 * t374;
t388 = t330 * t409;
t375 = t326 * t388;
t404 = Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1;
t431 = t326 * t330;
t403 = t329 * t431;
t377 = qJD(3) * t403;
t200 = qJD(1) * t377 - qJD(3) * t344 - t333 * t374;
t478 = -t200 / 0.2e1;
t513 = Ifges(5,2) * t478;
t415 = qJD(2) * t326;
t395 = t330 * t415;
t357 = t495 * t395;
t263 = qJD(2) * t348;
t253 = qJD(1) * t263;
t398 = pkin(1) * t435;
t276 = -pkin(8) * t431 + t334 * t398;
t265 = t276 * qJD(2);
t254 = qJD(1) * t265;
t87 = -t226 * t413 - t237 * t414 + t253 * t333 - t329 * t254;
t61 = -pkin(4) * t200 - qJD(1) * t357 - t87;
t430 = t326 * t334;
t277 = pkin(8) * t430 + t330 * t398;
t266 = t277 * qJD(2);
t255 = qJD(1) * t266;
t338 = t200 * qJ(4) - t246 * qJD(4) + t255;
t67 = t201 * t495 + t338;
t13 = t328 * t61 + t332 * t67 + t94 * t411 - t412 * t99;
t14 = -qJD(5) * t47 - t328 * t67 + t332 * t61;
t102 = qJD(5) * t191 + t201 * t328 + t332 * t375;
t6 = -pkin(5) * t200 - pkin(11) * t102 + t14;
t103 = -qJD(5) * t192 + t201 * t332 - t328 * t375;
t7 = pkin(11) * t103 + t13;
t2 = qJD(6) * t15 + t327 * t6 + t331 * t7;
t522 = qJD(6) * t16;
t3 = -t327 * t7 + t331 * t6 - t522;
t558 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t538 = t14 * mrSges(6,1) - t13 * mrSges(6,2) + t558;
t546 = -t201 / 0.2e1;
t548 = -Ifges(5,4) / 0.2e1;
t81 = t201 * pkin(3) + t338;
t82 = -pkin(3) * t375 - t87;
t561 = t82 * mrSges(5,1) - t87 * mrSges(4,3) - t81 * mrSges(5,3) + Ifges(5,6) * t546 - t404 * t201 + t375 * t548 + t513 + t538;
t426 = t331 * t332;
t352 = t327 * t328 - t426;
t353 = t327 * t332 + t331 * t328;
t514 = qJD(5) + qJD(6);
t210 = t514 * t353;
t346 = t353 * t246;
t419 = t346 + t210;
t410 = qJD(6) * t327;
t209 = -t327 * t412 - t328 * t410 + t426 * t514;
t521 = t352 * t246;
t557 = t209 - t521;
t560 = t15 * t419 - t16 * t557 - t2 * t353 + t3 * t352;
t383 = t331 * t191 - t192 * t327;
t115 = Ifges(7,4) * t383;
t119 = t191 * t327 + t192 * t331;
t452 = Ifges(7,4) * t119;
t228 = qJD(6) + t240;
t476 = -t228 / 0.2e1;
t487 = -t119 / 0.2e1;
t489 = -t383 / 0.2e1;
t52 = t119 * Ifges(7,1) + t228 * Ifges(7,5) + t115;
t291 = t303 * qJ(4);
t104 = t127 - t291;
t73 = -pkin(5) * t191 + t104;
t559 = (Ifges(7,5) * t383 - Ifges(7,6) * t119) * t476 + (t119 * t16 + t15 * t383) * mrSges(7,3) + (-Ifges(7,2) * t119 + t115 + t52) * t489 - t73 * (mrSges(7,1) * t119 + mrSges(7,2) * t383) + (Ifges(7,1) * t383 - t452) * t487;
t145 = t291 - t162;
t554 = t145 * mrSges(5,1) - t162 * mrSges(4,3);
t139 = t245 * pkin(3) + t339;
t552 = -t225 * mrSges(4,1) + t139 * mrSges(5,2) + Ifges(5,5) * t573 + Ifges(4,6) * t467 + (Ifges(4,2) + Ifges(5,3)) * t472 + (Ifges(4,4) + Ifges(5,6)) * t469;
t239 = Ifges(4,4) * t245;
t160 = t246 * Ifges(4,1) - t303 * Ifges(4,5) - t239;
t441 = t240 * Ifges(6,3);
t443 = t192 * Ifges(6,5);
t444 = t191 * Ifges(6,6);
t91 = t441 + t443 + t444;
t551 = t160 / 0.2e1 + t91 / 0.2e1 + t563;
t550 = t46 * mrSges(6,1) + t225 * mrSges(4,2) - t47 * mrSges(6,2) - t139 * mrSges(5,3) + Ifges(5,4) * t573 + Ifges(5,2) * t469 + t574;
t34 = qJD(6) * t383 + t102 * t331 + t103 * t327;
t500 = t34 / 0.2e1;
t35 = -qJD(6) * t119 - t102 * t327 + t103 * t331;
t499 = t35 / 0.2e1;
t51 = Ifges(7,2) * t383 + t228 * Ifges(7,6) + t452;
t547 = t51 / 0.2e1;
t100 = Ifges(6,6) * t103;
t101 = Ifges(6,5) * t102;
t40 = -Ifges(6,3) * t200 + t100 + t101;
t32 = Ifges(7,6) * t35;
t33 = Ifges(7,5) * t34;
t8 = -Ifges(7,3) * t200 + t32 + t33;
t545 = t40 + t8;
t542 = t87 * mrSges(4,1);
t526 = t254 * mrSges(3,2);
t293 = t458 * t328;
t216 = t293 * t327 - t294 * t331;
t525 = qJD(6) * t216 + t327 * t568 - t331 * t567;
t217 = -t293 * t331 - t294 * t327;
t524 = -qJD(6) * t217 + t327 * t567 + t331 * t568;
t523 = pkin(5) * t411 - t246 * t399 - t515;
t257 = pkin(9) * t435 + t277;
t183 = -t329 * t257 + t258 * t333;
t172 = pkin(3) * t430 - t183;
t273 = t329 * t435 + t333 * t431;
t128 = pkin(4) * t273 + pkin(10) * t430 + t172;
t272 = -t333 * t435 + t403;
t256 = -pkin(2) * t435 - t276;
t343 = -t273 * qJ(4) + t256;
t140 = t272 * t495 + t343;
t66 = t328 * t128 + t332 * t140;
t519 = mrSges(3,1) * t351 - mrSges(4,1) * t245 - mrSges(4,2) * t246 - mrSges(3,3) * t397;
t517 = t200 * t572 + t201 * t539 + t375 * t541;
t516 = t13 * t328 + t14 * t332;
t442 = t228 * Ifges(7,3);
t447 = t119 * Ifges(7,5);
t448 = t383 * Ifges(7,6);
t50 = t442 + t447 + t448;
t407 = Ifges(4,5) / 0.2e1 + t548;
t512 = t550 - t407 * t303 + t444 / 0.2e1 + t443 / 0.2e1 + t442 / 0.2e1 + t441 / 0.2e1 + t448 / 0.2e1 + t447 / 0.2e1 + t50 / 0.2e1 + t551 + t562;
t511 = t333 * t514;
t507 = -t161 * mrSges(4,1) - t162 * mrSges(4,2) + t143 * mrSges(5,2) - t264 * mrSges(3,3) - t145 * mrSges(5,3);
t359 = t328 * t46 - t332 * t47;
t454 = Ifges(6,4) * t328;
t364 = Ifges(6,2) * t332 + t454;
t453 = Ifges(6,4) * t332;
t366 = Ifges(6,1) * t328 + t453;
t369 = mrSges(6,1) * t332 - mrSges(6,2) * t328;
t449 = Ifges(6,6) * t332;
t451 = Ifges(6,5) * t328;
t465 = -t332 / 0.2e1;
t466 = -t328 / 0.2e1;
t474 = -t240 / 0.2e1;
t480 = -t192 / 0.2e1;
t482 = -t191 / 0.2e1;
t455 = Ifges(6,4) * t192;
t92 = t191 * Ifges(6,2) + t240 * Ifges(6,6) + t455;
t190 = Ifges(6,4) * t191;
t93 = t192 * Ifges(6,1) + t240 * Ifges(6,5) + t190;
t506 = mrSges(6,3) * t359 + t104 * t369 + (t449 + t451) * t474 + t364 * t482 + t366 * t480 + t465 * t92 + t466 * t93;
t470 = -t246 / 0.2e1;
t471 = t245 / 0.2e1;
t505 = -Ifges(4,4) * t469 - Ifges(4,2) * t472 + Ifges(5,6) * t470 + Ifges(5,3) * t471 + t467 * t539 - t552;
t473 = t240 / 0.2e1;
t479 = t192 / 0.2e1;
t481 = t191 / 0.2e1;
t504 = Ifges(4,1) * t469 + Ifges(4,4) * t472 + Ifges(6,5) * t479 - Ifges(5,2) * t470 - Ifges(5,6) * t471 + Ifges(6,6) * t481 + Ifges(6,3) * t473 - t467 * t572 + t550;
t502 = Ifges(7,4) * t500 + Ifges(7,2) * t499 + Ifges(7,6) * t478;
t501 = Ifges(7,1) * t500 + Ifges(7,4) * t499 + Ifges(7,5) * t478;
t42 = t102 * Ifges(6,1) + t103 * Ifges(6,4) - t200 * Ifges(6,5);
t498 = t42 / 0.2e1;
t497 = t92 / 0.2e1;
t496 = t93 / 0.2e1;
t492 = t102 / 0.2e1;
t491 = t103 / 0.2e1;
t488 = t383 / 0.2e1;
t486 = t119 / 0.2e1;
t475 = t228 / 0.2e1;
t463 = pkin(9) * t329;
t457 = Ifges(3,4) * t330;
t456 = Ifges(3,4) * t334;
t450 = Ifges(3,2) * t330;
t440 = t261 * mrSges(3,3);
t437 = t330 * Ifges(3,1);
t434 = Ifges(3,6) * qJD(2);
t427 = t329 * t332;
t125 = -mrSges(6,1) * t191 + mrSges(6,2) * t192;
t204 = mrSges(5,1) * t245 + mrSges(5,3) * t303;
t423 = t125 - t204;
t153 = t352 * t511 + t353 * t414;
t168 = t232 * t327 + t233 * t331;
t422 = t153 - t168;
t154 = -t352 * t414 + t353 * t511;
t167 = t232 * t331 - t233 * t327;
t421 = t154 - t167;
t202 = mrSges(4,2) * t303 - mrSges(4,3) * t245;
t418 = t202 - t204;
t203 = -mrSges(4,1) * t303 - mrSges(4,3) * t246;
t205 = mrSges(5,1) * t246 - mrSges(5,2) * t303;
t417 = t203 - t205;
t184 = t333 * t257 + t329 * t258;
t305 = t333 * pkin(9) + t325;
t408 = -Ifges(4,1) / 0.2e1 - Ifges(5,2) / 0.2e1;
t406 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t405 = -Ifges(4,2) / 0.2e1 - Ifges(5,3) / 0.2e1;
t394 = t334 * t415;
t391 = Ifges(3,5) * t435;
t390 = Ifges(3,6) * t435;
t387 = t415 / 0.2e1;
t177 = -t200 * mrSges(5,1) + mrSges(5,2) * t375;
t65 = t332 * t128 - t140 * t328;
t371 = qJD(1) * t387;
t368 = mrSges(6,1) * t328 + mrSges(6,2) * t332;
t367 = Ifges(6,1) * t332 - t454;
t365 = -Ifges(6,2) * t328 + t453;
t363 = Ifges(6,5) * t332 - Ifges(6,6) * t328;
t347 = -t272 * t328 + t326 * t424;
t49 = pkin(5) * t273 + pkin(11) * t347 + t65;
t214 = t272 * t332 + t326 * t428;
t54 = pkin(11) * t214 + t66;
t22 = -t327 * t54 + t331 * t49;
t23 = t327 * t49 + t331 * t54;
t74 = -mrSges(6,1) * t200 - mrSges(6,3) * t102;
t75 = mrSges(6,2) * t200 + mrSges(6,3) * t103;
t358 = t328 * t75 + t332 * t74;
t356 = t143 * t333 + t145 * t329;
t146 = -mrSges(6,2) * t240 + mrSges(6,3) * t191;
t147 = mrSges(6,1) * t240 - mrSges(6,3) * t192;
t355 = t332 * t146 - t328 * t147;
t354 = t161 * t333 - t162 * t329;
t150 = t214 * t331 + t327 * t347;
t151 = t214 * t327 - t331 * t347;
t350 = t330 * t371;
t171 = qJ(4) * t430 - t184;
t106 = -t257 * t413 - t258 * t414 + t263 * t333 - t329 * t265;
t213 = -t377 + (qJD(3) * t435 + t394) * t333;
t71 = pkin(4) * t213 - t106 - t357;
t212 = qJD(3) * t273 + t329 * t394;
t340 = -t213 * qJ(4) - t273 * qJD(4) + t266;
t80 = t212 * t495 + t340;
t20 = t128 * t411 - t140 * t412 + t328 * t71 + t332 * t80;
t86 = -t226 * t414 + t237 * t413 + t329 * t253 + t333 * t254;
t105 = -t257 * t414 + t258 * t413 + t329 * t263 + t333 * t265;
t141 = -pkin(4) * t272 - t171;
t77 = -qJ(4) * t375 + t303 * qJD(4) - t86;
t21 = -qJD(5) * t66 - t328 * t80 + t332 * t71;
t59 = -pkin(4) * t201 - t77;
t88 = -qJ(4) * t395 + qJD(4) * t430 - t105;
t342 = (t390 + (Ifges(3,2) * t334 + t457) * t326) * qJD(1);
t72 = -pkin(4) * t212 - t88;
t337 = t404 * t246 + t406 * t303 + t552 - t554;
t320 = pkin(5) * t328 + qJ(4);
t308 = Ifges(3,4) * t396;
t301 = Ifges(3,5) * t374;
t296 = -pkin(3) * t333 + t386;
t271 = pkin(5) * t425 + t305;
t270 = -qJ(4) * t413 + t381;
t268 = t353 * t333;
t267 = t352 * t333;
t260 = -mrSges(3,2) * t351 + mrSges(3,3) * t396;
t222 = Ifges(3,1) * t397 + Ifges(3,5) * t351 + t308;
t221 = t342 + t434;
t207 = -t278 * t328 + t285;
t189 = -qJ(4) * t378 + t401;
t182 = -mrSges(5,2) * t245 - mrSges(5,3) * t246;
t180 = pkin(3) * t246 + t433;
t179 = -mrSges(4,2) * t375 - mrSges(4,3) * t201;
t178 = mrSges(4,1) * t375 + mrSges(4,3) * t200;
t176 = mrSges(5,1) * t201 - mrSges(5,3) * t375;
t175 = -pkin(3) * t397 - t185;
t170 = t272 * pkin(3) + t343;
t159 = -t303 * Ifges(5,1) - t246 * Ifges(5,4) + t245 * Ifges(5,5);
t156 = t246 * Ifges(4,5) - t245 * Ifges(4,6) - t303 * Ifges(4,3);
t137 = qJD(5) * t214 + t212 * t328 + t332 * t395;
t136 = qJD(5) * t347 + t212 * t332 - t328 * t395;
t130 = mrSges(4,1) * t201 - mrSges(4,2) * t200;
t129 = -mrSges(5,2) * t201 + mrSges(5,3) * t200;
t114 = -t200 * Ifges(4,1) - t201 * Ifges(4,4) + Ifges(4,5) * t375;
t113 = -t200 * Ifges(4,4) - t201 * Ifges(4,2) + Ifges(4,6) * t375;
t111 = Ifges(5,5) * t375 + t200 * Ifges(5,6) + t201 * Ifges(5,3);
t97 = t212 * pkin(3) + t340;
t96 = -pkin(3) * t395 - t106;
t95 = -pkin(5) * t214 + t141;
t84 = mrSges(7,1) * t228 - mrSges(7,3) * t119;
t83 = -mrSges(7,2) * t228 + mrSges(7,3) * t383;
t63 = -t144 * t328 + t121;
t56 = -mrSges(7,1) * t383 + mrSges(7,2) * t119;
t53 = -mrSges(6,1) * t103 + mrSges(6,2) * t102;
t45 = -pkin(5) * t136 + t72;
t44 = -qJD(6) * t151 + t136 * t331 - t137 * t327;
t43 = qJD(6) * t150 + t136 * t327 + t137 * t331;
t41 = t102 * Ifges(6,4) + t103 * Ifges(6,2) - t200 * Ifges(6,6);
t37 = -pkin(5) * t103 + t59;
t29 = mrSges(7,2) * t200 + mrSges(7,3) * t35;
t28 = -mrSges(7,1) * t200 - mrSges(7,3) * t34;
t19 = t331 * t38 - t438;
t18 = -t327 * t38 - t436;
t17 = pkin(11) * t136 + t20;
t12 = pkin(5) * t213 - pkin(11) * t137 + t21;
t11 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t5 = -qJD(6) * t23 + t12 * t331 - t17 * t327;
t4 = qJD(6) * t22 + t12 * t327 + t17 * t331;
t1 = [(Ifges(6,5) * t137 + Ifges(6,6) * t136) * t473 + (-t15 * t43 + t150 * t2 - t151 * t3 + t16 * t44) * mrSges(7,3) + m(3) * (t254 * t277 + t264 * t265) + m(4) * (t105 * t162 - t106 * t161 + t183 * t87 + t184 * t86) + (-t272 * t81 - t430 * t82) * mrSges(5,2) - t394 * t440 + (Ifges(6,5) * t492 + Ifges(7,5) * t500 + Ifges(6,6) * t491 + Ifges(7,6) * t499 - t350 * t572 + t513 + t561) * t273 + t435 * (-Ifges(3,6) * t375 + t301) / 0.2e1 + (Ifges(7,5) * t43 + Ifges(7,6) * t44) * t475 + ((t159 + t156) * t330 + t351 * (Ifges(3,5) * t334 - Ifges(3,6) * t330) + t334 * t222) * t387 + m(5) * (t139 * t97 + t143 * t96 + t145 * t88 + t170 * t81 + t171 * t77 + t172 * t82) + m(6) * (t104 * t72 + t13 * t66 + t14 * t65 + t141 * t59 + t20 * t47 + t21 * t46) + m(7) * (t15 * t5 + t16 * t4 + t2 * t23 + t22 * t3 + t37 * t95 + t45 * t73) + (Ifges(6,4) * t137 + Ifges(6,2) * t136) * t481 + (Ifges(7,5) * t486 + Ifges(7,6) * t488 + Ifges(7,3) * t475 + t504 + t562 + t563) * t213 - (t342 + t221) * t395 / 0.2e1 + (t160 + t91 + t50) * t213 / 0.2e1 + (-0.2e1 * pkin(1) * (mrSges(3,1) * t330 + mrSges(3,2) * t334) * t409 + (-t450 + t456) * t389 + (Ifges(3,1) * t334 - t457) * t388) * t326 ^ 2 + t77 * (mrSges(5,1) * t272 + mrSges(5,3) * t430) + t86 * (mrSges(4,2) * t430 - mrSges(4,3) * t272) + t201 * (-Ifges(5,5) * t430 + Ifges(5,3) * t272) / 0.2e1 + (Ifges(7,4) * t43 + Ifges(7,2) * t44) * t488 + (Ifges(7,4) * t151 + Ifges(7,2) * t150) * t499 + (t505 + t554) * t212 + (t254 * t430 + t255 * t431 - t276 * t374 - t277 * t375) * mrSges(3,3) + t21 * t147 + t37 * (-mrSges(7,1) * t150 + mrSges(7,2) * t151) + t20 * t146 + t141 * t53 + t104 * (-mrSges(6,1) * t136 + mrSges(6,2) * t137) + t72 * t125 + (-Ifges(6,4) * t347 + Ifges(6,2) * t214) * t491 + (-Ifges(4,4) * t272 - Ifges(4,5) * t430 - Ifges(6,5) * t347 + Ifges(7,5) * t151 + Ifges(6,6) * t214 + Ifges(7,6) * t150 + (Ifges(4,1) + Ifges(6,3) + Ifges(7,3)) * t273) * t478 + (-Ifges(6,1) * t347 + Ifges(6,4) * t214) * t492 + t59 * (-mrSges(6,1) * t214 - mrSges(6,2) * t347) + (t13 * t214 + t136 * t47 - t137 * t46 + t14 * t347) * mrSges(6,3) - t347 * t498 + t95 * t11 + t4 * t83 + t5 * t84 + t73 * (-mrSges(7,1) * t44 + mrSges(7,2) * t43) + t65 * t74 + t66 * t75 + t45 * t56 + t43 * t52 / 0.2e1 + (-m(3) * t276 + m(4) * t256 - mrSges(3,1) * t435 + mrSges(4,1) * t272 + mrSges(4,2) * t273) * t255 + t22 * t28 + t23 * t29 + t200 * (-Ifges(5,4) * t430 + Ifges(5,6) * t272) / 0.2e1 + (Ifges(7,1) * t43 + Ifges(7,4) * t44) * t486 + (Ifges(7,1) * t151 + Ifges(7,4) * t150) * t500 + (Ifges(6,1) * t137 + Ifges(6,4) * t136) * t479 + t170 * t129 + t171 * t176 + t172 * t177 + t97 * t182 + t183 * t178 + t184 * t179 + t334 * (t391 + (t437 + t456) * t326) * t371 + t105 * t202 + t106 * t203 + t88 * t204 + t96 * t205 - t517 * t430 / 0.2e1 + (-m(3) * t261 + m(4) * t225 - t519) * t266 + t214 * t41 / 0.2e1 + t137 * t496 + t136 * t497 + t151 * t501 + t150 * t502 + (-Ifges(4,2) * t272 - Ifges(4,6) * t430) * t546 + t44 * t547 + t256 * t130 - t435 * t526 + t265 * t260 + t272 * t111 / 0.2e1 - t272 * t113 / 0.2e1 + (Ifges(5,4) * t470 + Ifges(4,5) * t469 + Ifges(5,5) * t471 + Ifges(4,6) * t472 + t467 * t541 + t507) * t395 + (t272 * t539 - t430 * t541) * t350 - t430 * t542 + (t114 + t545) * t273 / 0.2e1; (-t15 * t422 + t16 * t421 + t2 * t267 + t268 * t3) * mrSges(7,3) + (-Ifges(7,5) * t268 + Ifges(7,6) * t267) * t478 + (-Ifges(7,4) * t268 + Ifges(7,2) * t267) * t499 + (-Ifges(7,1) * t268 + Ifges(7,4) * t267) * t500 + t37 * (-mrSges(7,1) * t267 - mrSges(7,2) * t268) + (t354 * mrSges(4,3) + t356 * mrSges(5,1) + (m(4) * t354 + m(5) * t356) * pkin(9) + (Ifges(6,5) * t429 + Ifges(6,6) * t427) * t473 + t104 * (-mrSges(6,1) * t427 + mrSges(6,2) * t429) + (Ifges(6,1) * t429 + Ifges(6,4) * t427) * t479 + (Ifges(6,4) * t429 + Ifges(6,2) * t427) * t481 + t429 * t496 + t427 * t497 + (t427 * t47 - t429 * t46) * mrSges(6,3) + (-t418 * pkin(9) + t505) * t329 + (-t417 * pkin(9) + t50 + t504 + t551) * t333) * qJD(3) + t570 * t84 + t571 * t83 + (t109 * t3 + t110 * t2 + t15 * t570 + t16 * t571 + t271 * t37 + t565 * t73) * m(7) + (t255 * mrSges(4,2) + t100 / 0.2e1 + t101 / 0.2e1 + t114 / 0.2e1 + t40 / 0.2e1 + t8 / 0.2e1 + t33 / 0.2e1 + t32 / 0.2e1 + (t177 - t178) * pkin(9) + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t408) * t200 + t561) * t329 + ((-t434 / 0.2e1 - t156 / 0.2e1 - t159 / 0.2e1 + t221 / 0.2e1 + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t303 - t407 * t246 - t406 * t245 + (t390 / 0.2e1 + (pkin(1) * mrSges(3,1) + t457 / 0.2e1) * t326) * qJD(1) + (-t329 * t572 - t539 * t333) * qJD(2) / 0.2e1 - t507) * t330 + (t440 - qJD(2) * Ifges(3,5) / 0.2e1 + ((pkin(1) * mrSges(3,2) - t437 / 0.2e1 + t450 / 0.2e1) * t326 - t391 / 0.2e1) * qJD(1) - t222 / 0.2e1 - t308 / 0.2e1 + (t245 * t405 + t337) * t329 + (t245 * t404 + t246 * t408 - t512) * t333) * t334) * t416 + t566 * t147 + (t104 * t564 + t13 * t208 + t14 * t207 + t305 * t59 + t46 * t566 + t47 * t569) * m(6) + t569 * t146 + (-mrSges(7,1) * t421 + mrSges(7,2) * t422) * t73 + (Ifges(6,5) * t233 + Ifges(6,6) * t232) * t474 + (Ifges(7,5) * t153 + Ifges(7,6) * t154) * t475 + (Ifges(7,5) * t168 + Ifges(7,6) * t167) * t476 + (t154 / 0.2e1 - t167 / 0.2e1) * t51 + (-t189 + t270) * t182 + (-t232 * t47 + t233 * t46) * mrSges(6,3) + t564 * t125 + t565 * t56 - m(4) * (-t161 * t185 + t162 * t186 + t225 * t264) - m(5) * (t139 * t189 + t143 * t175 + t145 * t173) + (t153 / 0.2e1 - t168 / 0.2e1) * t52 - pkin(2) * t130 - t526 + t109 * t28 + t110 * t29 + (Ifges(6,1) * t233 + Ifges(6,4) * t232) * t480 + (Ifges(6,4) * t233 + Ifges(6,2) * t232) * t482 + m(4) * (-pkin(2) * t255 - t463 * t87) + m(5) * (t139 * t270 + t296 * t81 + t463 * t82) - t186 * t202 - t185 * t203 - t173 * t204 - t175 * t205 + t207 * t74 + t208 * t75 + t519 * t264 + (Ifges(7,1) * t153 + Ifges(7,4) * t154) * t486 + (Ifges(7,1) * t168 + Ifges(7,4) * t167) * t487 + t301 + (Ifges(7,4) * t153 + Ifges(7,2) * t154) * t488 + (Ifges(7,4) * t168 + Ifges(7,2) * t167) * t489 - t268 * t501 + t267 * t502 + (t59 * t369 - t103 * t364 / 0.2e1 - t102 * t366 / 0.2e1 + t42 * t466 + t41 * t465 - t255 * mrSges(4,1) + t81 * mrSges(5,2) - t77 * mrSges(5,1) + t86 * mrSges(4,3) + t113 / 0.2e1 - t111 / 0.2e1 + t405 * t201 + (t451 / 0.2e1 + t449 / 0.2e1 - t404) * t200 + (-t13 * t332 + t14 * t328) * mrSges(6,3) + (m(4) * t86 - m(5) * t77 - t176 + t179) * pkin(9) + (-t104 * t368 + t328 * t497 + t93 * t465 + t365 * t482 + t367 * t480 + t363 * t474 + (t328 * t47 + t332 * t46) * mrSges(6,3)) * qJD(5)) * t333 - t232 * t92 / 0.2e1 - t233 * t93 / 0.2e1 - t104 * (-mrSges(6,1) * t232 + mrSges(6,2) * t233) - t255 * mrSges(3,1) - t261 * t260 + t271 * t11 + t296 * t129 + t305 * t53; (-Ifges(7,5) * t210 - Ifges(7,6) * t209) * t475 + (-t346 / 0.2e1 - t210 / 0.2e1) * t52 + (-Ifges(7,1) * t210 - Ifges(7,4) * t209) * t486 + (-Ifges(7,4) * t210 - Ifges(7,2) * t209) * t488 + (-t239 / 0.2e1 + t574 + t512) * t245 + t418 * t161 + t423 * qJD(4) + t41 * t466 + t417 * t162 + (t53 - t176) * qJ(4) - t516 * mrSges(6,3) + (mrSges(7,1) * t557 - mrSges(7,2) * t419) * t73 + (t59 * qJ(4) + t553 * t104 - t46 * t63 - t47 * t64) * m(6) + t560 * mrSges(7,3) + t59 * t368 + t517 - t64 * t146 - t63 * t147 + t349 * t125 + (Ifges(7,5) * t346 - Ifges(7,6) * t521) * t476 + (Ifges(7,1) * t346 - Ifges(7,4) * t521) * t487 + (Ifges(7,4) * t346 - Ifges(7,2) * t521) * t489 + (t521 / 0.2e1 - t209 / 0.2e1) * t51 - t353 * t502 + (-Ifges(7,5) * t352 - Ifges(7,6) * t353 + t363) * t478 + (-Ifges(7,4) * t352 - Ifges(7,2) * t353) * t499 + (-Ifges(7,1) * t352 - Ifges(7,4) * t353) * t500 + t37 * (mrSges(7,1) * t353 - mrSges(7,2) * t352) - t352 * t501 - (m(6) * t516 + t358 + (-m(6) * t359 + t355) * qJD(5)) * t495 + (t337 + (t405 - t408) * t245 + t506) * t246 + t506 * qJD(5) - t86 * mrSges(4,2) + t82 * mrSges(5,2) - t77 * mrSges(5,3) + t542 - pkin(3) * t177 - t180 * t182 + (-pkin(3) * t82 - qJ(4) * t77 - t139 * t180 - t143 * t162 + t145 * t515) * m(5) + t216 * t28 + t217 * t29 + t365 * t491 + t367 * t492 + t332 * t498 + t523 * t56 + t524 * t84 + t525 * t83 + (t15 * t524 + t16 * t525 + t2 * t217 + t216 * t3 + t320 * t37 + t523 * t73) * m(7) + t320 * t11; -t352 * t28 + t353 * t29 - t419 * t84 + t557 * t83 + t355 * qJD(5) + (t56 + t423) * t303 + (t182 + t355) * t246 + t358 + t177 + (t303 * t73 - t560) * m(7) + (t104 * t303 - t240 * t359 + t516) * m(6) + (t139 * t246 - t145 * t303 + t82) * m(5); (Ifges(6,5) * t191 - Ifges(6,6) * t192) * t474 + t545 - t46 * t146 + t47 * t147 + (-Ifges(6,2) * t192 + t190 + t93) * t482 + t119 * t547 - m(7) * (t15 * t18 + t16 * t19) + (t191 * t46 + t192 * t47) * mrSges(6,3) - t19 * t83 - t18 * t84 + t538 + t92 * t479 + (Ifges(6,1) * t191 - t455) * t480 - t104 * (mrSges(6,1) * t192 + mrSges(6,2) * t191) + (-t192 * t56 + t331 * t28 + t327 * t29 + (-t327 * t84 + t331 * t83) * qJD(6) + (-t15 * t410 - t192 * t73 + t2 * t327 + (t3 + t522) * t331) * m(7)) * pkin(5) + t559; -t15 * t83 + t16 * t84 + t51 * t486 + t558 + t559 + t8;];
tauc  = t1(:);
