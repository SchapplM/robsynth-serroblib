% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:17
% EndTime: 2019-03-10 03:33:45
% DurationCPUTime: 14.55s
% Computational Cost: add. (30394->747), mult. (75991->1024), div. (0->0), fcn. (57618->10), ass. (0->364)
t362 = sin(qJ(6));
t363 = sin(qJ(5));
t367 = cos(qJ(6));
t368 = cos(qJ(5));
t385 = t362 * t363 - t367 * t368;
t523 = qJD(5) + qJD(6);
t278 = t523 * t385;
t365 = sin(qJ(3));
t366 = sin(qJ(2));
t370 = cos(qJ(3));
t371 = cos(qJ(2));
t424 = t370 * t371;
t325 = -t365 * t366 + t424;
t312 = t325 * qJD(1);
t327 = t365 * t371 + t370 * t366;
t313 = t327 * qJD(1);
t364 = sin(qJ(4));
t369 = cos(qJ(4));
t525 = t369 * t312 - t364 * t313;
t547 = t385 * t525;
t567 = t547 - t278;
t350 = pkin(3) * t364 + pkin(10);
t473 = -pkin(11) - t350;
t401 = qJD(5) * t473;
t415 = qJD(4) * t369;
t410 = pkin(3) * t415;
t546 = t525 * t363;
t412 = pkin(11) * t546;
t510 = -pkin(8) - pkin(7);
t345 = t510 * t371;
t332 = qJD(1) * t345;
t317 = t370 * t332;
t343 = t510 * t366;
t331 = qJD(1) * t343;
t320 = qJD(2) * pkin(2) + t331;
t272 = t320 * t365 - t317;
t477 = pkin(9) * t312;
t235 = t272 + t477;
t225 = t364 * t235;
t314 = t365 * t332;
t271 = t370 * t320 + t314;
t306 = t313 * pkin(9);
t234 = t271 - t306;
t158 = t234 * t369 - t225;
t387 = t312 * t364 + t369 * t313;
t206 = pkin(4) * t387 - pkin(10) * t525;
t480 = pkin(3) * t313;
t177 = t206 + t480;
t95 = t368 * t158 + t363 * t177;
t566 = t363 * t401 + t368 * t410 + t412 - t95;
t434 = t525 * t368;
t538 = pkin(5) * t387;
t398 = -pkin(11) * t434 + t538;
t94 = -t158 * t363 + t368 * t177;
t565 = -t363 * t410 + t368 * t401 - t398 - t94;
t352 = pkin(2) * t370 + pkin(3);
t416 = qJD(4) * t364;
t429 = t364 * t365;
t263 = t352 * t415 + (-t365 * t416 + (t369 * t370 - t429) * qJD(3)) * pkin(2);
t276 = -t331 * t365 + t317;
t238 = t276 - t477;
t277 = t370 * t331 + t314;
t239 = -t306 + t277;
t171 = t238 * t364 + t239 * t369;
t421 = qJD(1) * t366;
t356 = pkin(2) * t421;
t172 = t177 + t356;
t97 = t368 * t171 + t363 * t172;
t564 = t263 * t368 - t97;
t96 = -t171 * t363 + t368 * t172;
t563 = -t263 * t363 - t96;
t244 = pkin(5) * t546;
t428 = t365 * t369;
t308 = pkin(2) * t428 + t364 * t352;
t304 = pkin(10) + t308;
t474 = -pkin(11) - t304;
t402 = qJD(5) * t474;
t562 = t363 * t402 + t412 + t564;
t561 = t368 * t402 - t398 + t563;
t361 = qJD(2) + qJD(3);
t224 = pkin(3) * t361 + t234;
t150 = t224 * t369 - t225;
t101 = t368 * t150 + t363 * t206;
t509 = -pkin(11) - pkin(10);
t408 = qJD(5) * t509;
t560 = t363 * t408 - t101 + t412;
t100 = -t150 * t363 + t368 * t206;
t360 = t368 * pkin(11);
t559 = t360 * t525 + t368 * t408 - t100 - t538;
t414 = qJD(5) * t363;
t357 = pkin(5) * t414;
t558 = t357 - t244;
t326 = t362 * t368 + t363 * t367;
t548 = t326 * t525;
t557 = -Ifges(7,1) * t547 - Ifges(7,4) * t548;
t556 = -Ifges(7,4) * t547 - Ifges(7,2) * t548;
t555 = -Ifges(7,5) * t547 - Ifges(7,6) * t548;
t359 = qJD(4) + t361;
t139 = -pkin(4) * t359 - t150;
t236 = t359 * t368 - t363 * t387;
t118 = -pkin(5) * t236 + t139;
t280 = t361 * t325;
t265 = t280 * qJD(1);
t281 = t361 * t327;
t266 = t281 * qJD(1);
t148 = qJD(4) * t387 + t265 * t364 + t369 * t266;
t143 = Ifges(7,3) * t148;
t237 = t359 * t363 + t368 * t387;
t163 = t236 * t362 + t237 * t367;
t226 = t369 * t235;
t151 = t224 * t364 + t226;
t140 = pkin(10) * t359 + t151;
t354 = -pkin(2) * t371 - pkin(1);
t341 = qJD(1) * t354;
t282 = -t312 * pkin(3) + t341;
t168 = -pkin(4) * t525 - pkin(10) * t387 + t282;
t91 = t140 * t368 + t168 * t363;
t73 = pkin(11) * t236 + t91;
t449 = t362 * t73;
t254 = qJD(5) - t525;
t90 = -t140 * t363 + t368 * t168;
t72 = -pkin(11) * t237 + t90;
t61 = pkin(5) * t254 + t72;
t22 = t367 * t61 - t449;
t447 = t367 * t73;
t23 = t362 * t61 + t447;
t400 = t367 * t236 - t237 * t362;
t465 = Ifges(7,4) * t163;
t246 = qJD(6) + t254;
t493 = -t246 / 0.2e1;
t499 = -t163 / 0.2e1;
t554 = t143 + (Ifges(7,5) * t400 - Ifges(7,6) * t163) * t493 + (t163 * t23 + t22 * t400) * mrSges(7,3) - t118 * (mrSges(7,1) * t163 + mrSges(7,2) * t400) + (Ifges(7,1) * t400 - t465) * t499;
t147 = qJD(4) * t525 + t265 * t369 - t266 * t364;
t279 = t523 * t326;
t115 = -qJD(5) * t237 - t147 * t363;
t413 = qJD(5) * t368;
t409 = qJD(2) * t510;
t417 = qJD(3) * t370;
t418 = qJD(3) * t365;
t210 = t313 * t409 + t320 * t417 + t332 * t418;
t155 = -pkin(9) * t266 + t210;
t328 = t365 * t343;
t377 = (t424 * t510 - t328) * qJD(2) * qJD(1);
t476 = t265 * pkin(9);
t57 = t369 * t155 + t150 * qJD(4) + (-t320 * t418 + t332 * t417 + t377 - t476) * t364;
t242 = pkin(3) * t266 + qJD(2) * t356;
t71 = pkin(4) * t148 - pkin(10) * t147 + t242;
t16 = -t140 * t414 + t168 * t413 + t363 * t71 + t368 * t57;
t12 = pkin(11) * t115 + t16;
t114 = qJD(5) * t236 + t147 * t368;
t441 = qJD(5) * t91;
t17 = -t363 * t57 + t368 * t71 - t441;
t7 = pkin(5) * t148 - pkin(11) * t114 + t17;
t3 = qJD(6) * t22 + t12 * t367 + t362 * t7;
t211 = -qJD(3) * t272 + t377;
t58 = qJD(4) * t151 + t155 * t364 - t369 * (t211 - t476);
t35 = -pkin(5) * t115 + t58;
t390 = Ifges(6,5) * t363 + Ifges(6,6) * t368;
t467 = Ifges(6,4) * t363;
t392 = Ifges(6,2) * t368 + t467;
t466 = Ifges(6,4) * t368;
t394 = Ifges(6,1) * t363 + t466;
t397 = mrSges(6,1) * t368 - mrSges(6,2) * t363;
t42 = t114 * Ifges(6,4) + t115 * Ifges(6,2) + t148 * Ifges(6,6);
t43 = Ifges(6,1) * t114 + Ifges(6,4) * t115 + Ifges(6,5) * t148;
t463 = t16 * t368;
t483 = t368 / 0.2e1;
t492 = t246 / 0.2e1;
t498 = t163 / 0.2e1;
t500 = t400 / 0.2e1;
t502 = t148 / 0.2e1;
t503 = t115 / 0.2e1;
t504 = t114 / 0.2e1;
t156 = Ifges(7,4) * t400;
t81 = Ifges(7,1) * t163 + Ifges(7,5) * t246 + t156;
t511 = t81 / 0.2e1;
t512 = -t81 / 0.2e1;
t80 = Ifges(7,2) * t400 + Ifges(7,6) * t246 + t465;
t513 = t80 / 0.2e1;
t514 = -t80 / 0.2e1;
t47 = -qJD(6) * t163 - t114 * t362 + t115 * t367;
t515 = t47 / 0.2e1;
t46 = qJD(6) * t400 + t114 * t367 + t115 * t362;
t516 = t46 / 0.2e1;
t517 = Ifges(7,1) * t516 + Ifges(7,4) * t515 + Ifges(7,5) * t502;
t518 = Ifges(7,4) * t516 + Ifges(7,2) * t515 + Ifges(7,6) * t502;
t396 = mrSges(6,1) * t363 + mrSges(6,2) * t368;
t382 = t139 * t396;
t391 = Ifges(6,5) * t368 - Ifges(6,6) * t363;
t393 = -Ifges(6,2) * t363 + t466;
t395 = Ifges(6,1) * t368 - t467;
t231 = Ifges(6,4) * t236;
t130 = t237 * Ifges(6,1) + t254 * Ifges(6,5) + t231;
t426 = t368 * t130;
t468 = Ifges(6,4) * t237;
t129 = Ifges(6,2) * t236 + Ifges(6,6) * t254 + t468;
t430 = t363 * t129;
t494 = t237 / 0.2e1;
t541 = t382 + t426 / 0.2e1 - t430 / 0.2e1 + t254 * t391 / 0.2e1 + t395 * t494 + t236 * t393 / 0.2e1;
t553 = -t548 * t514 - t547 * t512 + (Ifges(7,4) * t326 - Ifges(7,2) * t385) * t515 + (Ifges(7,1) * t326 - Ifges(7,4) * t385) * t516 + t35 * (mrSges(7,1) * t385 + mrSges(7,2) * t326) - t385 * t518 + (Ifges(7,5) * t326 - Ifges(7,6) * t385 + t390) * t502 + (-Ifges(7,5) * t278 - Ifges(7,6) * t279) * t492 + (-Ifges(7,1) * t278 - Ifges(7,4) * t279) * t498 + (-Ifges(7,4) * t278 - Ifges(7,2) * t279) * t500 + (-t23 * t279 - t3 * t385) * mrSges(7,3) + (t567 * mrSges(7,2) + (t279 - t548) * mrSges(7,1)) * t118 + mrSges(6,3) * t463 + t541 * qJD(5) + (-t397 - mrSges(5,1)) * t58 - t57 * mrSges(5,2) + t42 * t483 + Ifges(5,5) * t147 - Ifges(5,6) * t148 + t392 * t503 + t394 * t504 - t278 * t511 - t279 * t513 + t326 * t517 + t363 * t43 / 0.2e1;
t250 = Ifges(5,4) * t525;
t321 = t473 * t363;
t322 = t350 * t368 + t360;
t268 = t321 * t362 + t322 * t367;
t551 = -qJD(6) * t268 - t362 * t566 + t565 * t367;
t267 = t321 * t367 - t322 * t362;
t550 = qJD(6) * t267 + t565 * t362 + t367 * t566;
t51 = -mrSges(6,1) * t115 + mrSges(6,2) * t114;
t549 = -m(6) * t58 - t51;
t157 = t234 * t364 + t226;
t545 = pkin(3) * t416 - t157 + t558;
t528 = -t369 * t238 + t239 * t364 - t352 * t416 - (t365 * t415 + (t364 * t370 + t428) * qJD(3)) * pkin(2);
t451 = t359 * Ifges(5,5);
t453 = t387 * Ifges(5,1);
t196 = t250 + t451 + t453;
t388 = t363 * t91 + t368 * t90;
t384 = t388 * mrSges(6,3);
t544 = t282 * mrSges(5,2) + t196 / 0.2e1 + t451 / 0.2e1 + t250 / 0.2e1 + t541 - t384;
t4 = -qJD(6) * t23 - t12 * t362 + t367 * t7;
t543 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + Ifges(7,5) * t46 + Ifges(7,6) * t47;
t542 = -Ifges(7,2) * t163 + t156;
t283 = t474 * t363;
t284 = t304 * t368 + t360;
t220 = t283 * t367 - t284 * t362;
t537 = qJD(6) * t220 + t362 * t561 + t367 * t562;
t221 = t283 * t362 + t284 * t367;
t536 = -qJD(6) * t221 - t362 * t562 + t367 * t561;
t469 = Ifges(5,4) * t387;
t342 = t509 * t363;
t344 = pkin(10) * t368 + t360;
t287 = t342 * t362 + t344 * t367;
t531 = -qJD(6) * t287 - t362 * t560 + t367 * t559;
t285 = t342 * t367 - t344 * t362;
t530 = qJD(6) * t285 + t362 * t559 + t367 * t560;
t275 = t325 * t364 + t327 * t369;
t208 = t385 * t275;
t529 = -t528 + t558;
t423 = -mrSges(5,1) * t359 - mrSges(6,1) * t236 + mrSges(6,2) * t237 + mrSges(5,3) * t387;
t527 = -t171 + t263;
t286 = t370 * t343 + t345 * t365;
t251 = -pkin(9) * t327 + t286;
t288 = -t370 * t345 + t328;
t252 = pkin(9) * t325 + t288;
t198 = t251 * t364 + t252 * t369;
t188 = t368 * t198;
t295 = -t325 * pkin(3) + t354;
t386 = t369 * t325 - t327 * t364;
t193 = -pkin(4) * t386 - t275 * pkin(10) + t295;
t108 = t363 * t193 + t188;
t526 = t369 * t251 - t252 * t364;
t524 = -t363 * t90 + t368 * t91;
t520 = t359 * Ifges(5,6) / 0.2e1 + t469 / 0.2e1;
t519 = t17 * mrSges(6,1) - t16 * mrSges(6,2) + Ifges(6,5) * t114 + Ifges(6,6) * t115 + t543;
t507 = pkin(1) * mrSges(3,1);
t506 = pkin(1) * mrSges(3,2);
t501 = -t400 / 0.2e1;
t454 = t525 * Ifges(5,2);
t497 = t454 / 0.2e1 + t520;
t496 = -t236 / 0.2e1;
t495 = -t237 / 0.2e1;
t491 = -t254 / 0.2e1;
t487 = t312 / 0.2e1;
t486 = -t313 / 0.2e1;
t485 = t313 / 0.2e1;
t484 = -t363 / 0.2e1;
t482 = m(4) * t341;
t481 = mrSges(6,3) * t90;
t479 = pkin(3) * t369;
t478 = pkin(5) * t368;
t475 = t4 * t326;
t472 = mrSges(4,3) * t312;
t471 = Ifges(3,4) * t366;
t470 = Ifges(4,4) * t313;
t464 = t150 * mrSges(5,3);
t462 = t400 * Ifges(7,6);
t461 = t163 * Ifges(7,5);
t460 = t17 * t363;
t459 = t526 * t58;
t458 = t236 * Ifges(6,6);
t457 = t237 * Ifges(6,5);
t456 = t246 * Ifges(7,3);
t455 = t254 * Ifges(6,3);
t452 = t313 * mrSges(4,3);
t445 = Ifges(3,5) * qJD(2);
t444 = Ifges(3,6) * qJD(2);
t443 = qJD(2) * mrSges(3,1);
t442 = qJD(2) * mrSges(3,2);
t440 = t151 * t387;
t431 = t275 * t363;
t420 = qJD(1) * t371;
t419 = qJD(2) * t366;
t411 = mrSges(6,3) * t460;
t353 = -pkin(4) - t478;
t407 = t445 / 0.2e1;
t406 = -t444 / 0.2e1;
t260 = pkin(2) * t419 + pkin(3) * t281;
t336 = t366 * t409;
t337 = t371 * t409;
t218 = t370 * t336 + t365 * t337 + t343 * t417 + t345 * t418;
t186 = -pkin(9) * t281 + t218;
t219 = -qJD(3) * t288 - t336 * t365 + t370 * t337;
t187 = -pkin(9) * t280 + t219;
t76 = qJD(4) * t526 + t186 * t369 + t187 * t364;
t173 = qJD(4) * t386 + t280 * t369 - t281 * t364;
t174 = qJD(4) * t275 + t280 * t364 + t369 * t281;
t87 = pkin(4) * t174 - pkin(10) * t173 + t260;
t403 = -t363 * t76 + t368 * t87;
t107 = t368 * t193 - t198 * t363;
t307 = -pkin(2) * t429 + t352 * t369;
t303 = -pkin(4) - t307;
t84 = -pkin(5) * t386 - t275 * t360 + t107;
t93 = -pkin(11) * t431 + t108;
t40 = -t362 * t93 + t367 * t84;
t41 = t362 * t84 + t367 * t93;
t63 = mrSges(6,1) * t148 - mrSges(6,3) * t114;
t64 = -mrSges(6,2) * t148 + mrSges(6,3) * t115;
t389 = -t363 * t63 + t368 * t64;
t383 = t173 * t363 + t275 * t413;
t19 = t193 * t413 - t198 * t414 + t363 * t87 + t368 * t76;
t77 = qJD(4) * t198 + t186 * t364 - t369 * t187;
t376 = m(6) * (-qJD(5) * t388 - t460 + t463);
t128 = t455 + t457 + t458;
t79 = t456 + t461 + t462;
t374 = t23 * mrSges(7,2) + t91 * mrSges(6,2) - t128 / 0.2e1 + t497 - t79 / 0.2e1 - t462 / 0.2e1 - t461 / 0.2e1 - t22 * mrSges(7,1) - t458 / 0.2e1 - t457 / 0.2e1 - t456 / 0.2e1 - t455 / 0.2e1 - t282 * mrSges(5,1) - t90 * mrSges(6,1) + t520;
t248 = Ifges(4,2) * t312 + Ifges(4,6) * t361 + t470;
t305 = Ifges(4,4) * t312;
t249 = Ifges(4,1) * t313 + Ifges(4,5) * t361 + t305;
t372 = (Ifges(7,3) * t387 + t555) * t493 + (Ifges(7,6) * t387 + t556) * t501 + (Ifges(7,5) * t387 + t557) * t499 + t553 - (-Ifges(5,2) * t387 + t196 + t250 + t426) * t525 / 0.2e1 - t23 * (-mrSges(7,2) * t387 - mrSges(7,3) * t548) - t22 * (mrSges(7,1) * t387 + mrSges(7,3) * t547) - t341 * (mrSges(4,1) * t313 + mrSges(4,2) * t312) - t90 * (mrSges(6,1) * t387 - mrSges(6,3) * t434) + t387 * t497 - (-Ifges(4,2) * t313 + t249 + t305) * t312 / 0.2e1 - t525 * t382 - t282 * (mrSges(5,1) * t387 + mrSges(5,2) * t525) + (Ifges(6,3) * t387 + t391 * t525) * t491 + (Ifges(6,5) * t387 + t395 * t525) * t495 + (Ifges(6,6) * t387 + t393 * t525) * t496 - t359 * (Ifges(5,5) * t525 - Ifges(5,6) * t387) / 0.2e1 - (Ifges(5,1) * t525 + t128 - t469 + t79) * t387 / 0.2e1 + t525 * t464 + t525 * t430 / 0.2e1 - t91 * (-mrSges(6,2) * t387 - mrSges(6,3) * t546) + (t22 * t278 - t475) * mrSges(7,3) + t271 * t472 + t248 * t485 + (Ifges(4,1) * t312 - t470) * t486 - t210 * mrSges(4,2) + t211 * mrSges(4,1) + Ifges(4,5) * t265 - Ifges(4,6) * t266 - t361 * (Ifges(4,5) * t312 - Ifges(4,6) * t313) / 0.2e1;
t355 = Ifges(3,4) * t420;
t340 = mrSges(3,3) * t420 - t442;
t339 = -mrSges(3,3) * t421 + t443;
t338 = t353 - t479;
t311 = Ifges(3,1) * t421 + t355 + t445;
t310 = t444 + (Ifges(3,2) * t371 + t471) * qJD(1);
t294 = t303 - t478;
t293 = mrSges(4,1) * t361 - t452;
t292 = -mrSges(4,2) * t361 + t472;
t291 = t356 + t480;
t270 = -mrSges(4,1) * t312 + mrSges(4,2) * t313;
t240 = -mrSges(5,2) * t359 + mrSges(5,3) * t525;
t207 = t326 * t275;
t205 = -mrSges(5,1) * t525 + mrSges(5,2) * t387;
t181 = mrSges(6,1) * t254 - mrSges(6,3) * t237;
t180 = -mrSges(6,2) * t254 + mrSges(6,3) * t236;
t144 = Ifges(6,3) * t148;
t138 = pkin(5) * t431 - t526;
t123 = mrSges(7,1) * t246 - mrSges(7,3) * t163;
t122 = -mrSges(7,2) * t246 + mrSges(7,3) * t400;
t120 = t151 + t244;
t92 = -mrSges(7,1) * t400 + mrSges(7,2) * t163;
t50 = -t326 * t173 + t208 * t523;
t49 = -t173 * t385 - t275 * t279;
t48 = pkin(5) * t383 + t77;
t27 = t367 * t72 - t449;
t26 = -t362 * t72 - t447;
t25 = -mrSges(7,2) * t148 + mrSges(7,3) * t47;
t24 = mrSges(7,1) * t148 - mrSges(7,3) * t46;
t20 = -qJD(5) * t108 + t403;
t18 = -pkin(11) * t383 + t19;
t14 = -mrSges(7,1) * t47 + mrSges(7,2) * t46;
t13 = -t173 * t360 + pkin(5) * t174 + (-t188 + (pkin(11) * t275 - t193) * t363) * qJD(5) + t403;
t6 = -qJD(6) * t41 + t13 * t367 - t18 * t362;
t5 = qJD(6) * t40 + t13 * t362 + t18 * t367;
t1 = [(t43 * t483 + t42 * t484 + Ifges(5,1) * t147 - Ifges(5,4) * t148 + t242 * mrSges(5,2) + t393 * t503 + t395 * t504 + t391 * t502 + (mrSges(5,3) + t396) * t58 + (-t16 * t363 - t17 * t368) * mrSges(6,3) + (-t368 * t129 / 0.2e1 + t130 * t484 + t139 * t397 + t392 * t496 + t394 * t495 + t390 * t491 - t524 * mrSges(6,3)) * qJD(5)) * t275 + (-t147 * t526 - t148 * t198 - t150 * t173 - t151 * t174) * mrSges(5,3) - t526 * t51 + (-t207 * t3 + t208 * t4 - t22 * t49 + t23 * t50) * mrSges(7,3) + (-Ifges(7,5) * t208 - Ifges(7,6) * t207) * t502 + (-Ifges(7,4) * t208 - Ifges(7,2) * t207) * t515 + (-Ifges(7,1) * t208 - Ifges(7,4) * t207) * t516 + t35 * (mrSges(7,1) * t207 - mrSges(7,2) * t208) + t341 * (mrSges(4,1) * t281 + mrSges(4,2) * t280) + t361 * (Ifges(4,5) * t280 - Ifges(4,6) * t281) / 0.2e1 + m(5) * (-t150 * t77 + t151 * t76 + t198 * t57 + t242 * t295 + t260 * t282 - t459) + m(6) * (t107 * t17 + t108 * t16 + t139 * t77 + t19 * t91 + t20 * t90 - t459) + (t210 * t325 - t211 * t327 - t265 * t286 - t266 * t288 - t271 * t280 - t272 * t281) * mrSges(4,3) + (-t325 * t266 - t281 * t487) * Ifges(4,2) + (t325 * t265 - t327 * t266 + t280 * t487 - t281 * t485) * Ifges(4,4) + t354 * (mrSges(4,1) * t266 + mrSges(4,2) * t265) - (-t57 * mrSges(5,3) + t143 / 0.2e1 + t144 / 0.2e1 - Ifges(5,4) * t147 + t242 * mrSges(5,1) + (Ifges(7,3) / 0.2e1 + Ifges(5,2) + Ifges(6,3) / 0.2e1) * t148 + t519) * t386 + m(7) * (t118 * t48 + t138 * t35 + t22 * t6 + t23 * t5 + t3 * t41 + t4 * t40) + m(4) * (t210 * t288 + t211 * t286 + t218 * t272 + t219 * t271) + t107 * t63 + t108 * t64 + t48 * t92 + t40 * t24 + t41 * t25 + (t453 / 0.2e1 + t544) * t173 + t118 * (-mrSges(7,1) * t50 + mrSges(7,2) * t49) + t5 * t122 + t6 * t123 + t138 * t14 + (-pkin(7) * t339 + t311 / 0.2e1 + t407 + (-0.2e1 * t506 + 0.3e1 / 0.2e1 * Ifges(3,4) * t371) * qJD(1)) * t371 * qJD(2) + (Ifges(7,5) * t49 + Ifges(7,6) * t50) * t492 + (Ifges(7,1) * t49 + Ifges(7,4) * t50) * t498 + (Ifges(7,4) * t49 + Ifges(7,2) * t50) * t500 + t49 * t511 + t50 * t513 - t208 * t517 - t207 * t518 + t19 * t180 + t20 * t181 + t76 * t240 + t260 * t205 + t423 * t77 + t280 * t249 / 0.2e1 - t281 * t248 / 0.2e1 + t218 * t292 + t219 * t293 + t295 * (mrSges(5,1) * t148 + mrSges(5,2) * t147) + (-t454 / 0.2e1 - t374) * t174 + (t327 * t265 + t280 * t485) * Ifges(4,1) + (-pkin(7) * t340 - t310 / 0.2e1 + t406 + (-0.2e1 * t507 - 0.3e1 / 0.2e1 * t471 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t371) * qJD(1) + (t270 + qJD(1) * (-mrSges(4,1) * t325 + mrSges(4,2) * t327) + 0.2e1 * t482) * pkin(2)) * t419; t527 * t240 + (t150 * t528 + t151 * t527 - t282 * t291 - t307 * t58 + t308 * t57) * m(5) - t423 * t528 + t529 * t92 + ((t292 * t370 - t293 * t365) * qJD(3) + (-t265 * t370 - t266 * t365) * mrSges(4,3)) * pkin(2) + t304 * t376 + (-t271 * t276 - t272 * t277 + (t210 * t365 + t211 * t370 + (-t271 * t365 + t272 * t370) * qJD(3)) * pkin(2)) * m(4) + t272 * t452 + t372 + (-t528 * t139 + t303 * t58 + t563 * t90 + t564 * t91) * m(6) + t536 * t123 + t537 * t122 + (t118 * t529 + t22 * t536 + t220 * t4 + t221 * t3 + t23 * t537 + t294 * t35) * m(7) - t97 * t180 - t96 * t181 + t220 * t24 + t221 * t25 - t291 * t205 - t277 * t292 - t276 * t293 + t294 * t14 + t303 * t51 + (-t147 * t307 - t148 * t308 + t440) * mrSges(5,3) + (-t263 * t181 + (-qJD(5) * t180 - t63) * t304 + (-t17 - t441) * mrSges(6,3)) * t363 + (t263 * t180 + t304 * t64 + (-t181 * t304 - t481) * qJD(5)) * t368 + ((t407 - t355 / 0.2e1 - t311 / 0.2e1 + qJD(1) * t506 + (t339 - t443) * pkin(7)) * t371 + (t406 + t310 / 0.2e1 + (t507 + t471 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t371) * qJD(1) + (t340 + t442) * pkin(7) + (-t270 - t482) * pkin(2)) * t366) * qJD(1); (-t313 * t205 + (-t147 * t369 - t148 * t364) * mrSges(5,3) + (t423 * t364 + (t368 * t180 - t363 * t181 + t240) * t369 + m(6) * (t139 * t364 + t369 * t524)) * qJD(4) + (t364 * t57 - t369 * t58 + 0.2e1 * t282 * t486 + (-t150 * t364 + t151 * t369) * qJD(4)) * m(5)) * pkin(3) - t384 * qJD(5) - m(6) * (t139 * t157 + t90 * t94 + t91 * t95) + mrSges(5,3) * t440 + t551 * t123 - m(5) * (-t150 * t157 + t151 * t158) + t372 + t545 * t92 + t550 * t122 - t411 - t95 * t180 - t94 * t181 - t158 * t240 + t267 * t24 + t268 * t25 - t423 * t157 - t271 * t292 + t338 * t14 + (t293 + t452) * t272 - t549 * (-pkin(4) - t479) + (t376 + t389 + (-t180 * t363 - t181 * t368) * qJD(5)) * t350 + (t118 * t545 + t22 * t551 + t23 * t550 + t267 * t4 + t268 * t3 + t338 * t35) * m(7); t530 * t122 + (t285 * t4 + t287 * t3 + t35 * t353 + t530 * t23 + t531 * t22 + (-t120 + t357) * t118) * m(7) + t531 * t123 + t553 + t555 * t493 + t557 * t499 + t556 * t501 + (-t22 * t567 + t23 * t548 - t475) * mrSges(7,3) + (t151 * mrSges(5,3) + t374) * t387 + pkin(10) * t376 + t389 * pkin(10) + (t464 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t387 - t544) * t525 - m(6) * (t100 * t90 + t101 * t91 + t139 * t151) - t411 - t120 * t92 + t549 * pkin(4) - t101 * t180 - t100 * t181 - t150 * t240 - t423 * t151 + t285 * t24 + t287 * t25 + t353 * t14 + ((-pkin(10) * t181 - t481) * t368 + (-mrSges(6,3) * t91 + pkin(5) * t92 - pkin(10) * t180) * t363) * qJD(5); t519 + t400 * t512 - t163 * t514 - m(7) * (t22 * t26 + t23 * t27) + t144 + (-t237 * t92 + t367 * t24 + t362 * t25 + (t122 * t367 - t123 * t362) * qJD(6) + (-t118 * t237 + t3 * t362 + t367 * t4 + (-t22 * t362 + t23 * t367) * qJD(6)) * m(7)) * pkin(5) - t27 * t122 + (t236 * t90 + t237 * t91) * mrSges(6,3) - t26 * t123 + (Ifges(6,5) * t236 - Ifges(6,6) * t237) * t491 + t129 * t494 + (Ifges(6,1) * t236 - t468) * t495 + t542 * t501 - t90 * t180 + t91 * t181 + (-Ifges(6,2) * t237 + t130 + t231) * t496 - t139 * (mrSges(6,1) * t237 + mrSges(6,2) * t236) + t554; t80 * t498 - t22 * t122 + t23 * t123 + (t542 + t81) * t501 + t543 + t554;];
tauc  = t1(:);
