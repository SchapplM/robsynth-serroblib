% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:33:52
% EndTime: 2018-11-23 18:34:23
% DurationCPUTime: 31.11s
% Computational Cost: add. (24484->935), mult. (62698->1251), div. (0->0), fcn. (48987->10), ass. (0->387)
t350 = sin(qJ(4));
t352 = sin(qJ(2));
t354 = cos(qJ(4));
t347 = sin(pkin(6));
t418 = qJD(1) * t347;
t355 = cos(qJ(3));
t356 = cos(qJ(2));
t428 = t355 * t356;
t264 = (-t350 * t428 + t352 * t354) * t418;
t351 = sin(qJ(3));
t414 = qJD(3) * t355;
t398 = t350 * t414;
t412 = qJD(4) * t354;
t583 = t351 * t412 + t264 + t398;
t402 = t352 * t418;
t348 = cos(pkin(6));
t480 = pkin(1) * t356;
t409 = t348 * t480;
t294 = -pkin(8) * t402 + qJD(1) * t409;
t368 = t347 * (pkin(2) * t352 - pkin(9) * t356);
t295 = qJD(1) * t368;
t218 = t355 * t294 + t351 * t295;
t203 = pkin(10) * t402 + t218;
t337 = t348 * t352 * pkin(1);
t386 = pkin(3) * t351 - pkin(10) * t355;
t435 = t347 * t356;
t222 = (t337 + (pkin(8) + t386) * t435) * qJD(1);
t132 = t354 * t203 + t350 * t222;
t316 = t386 * qJD(3);
t321 = -pkin(3) * t355 - pkin(10) * t351 - pkin(2);
t413 = qJD(4) * t350;
t415 = qJD(3) * t351;
t210 = t350 * t316 + t321 * t412 + (-t354 * t415 - t355 * t413) * pkin(9);
t582 = t210 - t132;
t131 = -t350 * t203 + t354 * t222;
t265 = (t350 * t352 + t354 * t428) * t418;
t429 = t354 * t355;
t339 = pkin(9) * t429;
t401 = t356 * t418;
t389 = t351 * t401;
t479 = pkin(9) * t350;
t420 = t354 * t316 + t415 * t479;
t581 = -pkin(4) * t389 + t265 * pkin(11) - t131 + (pkin(4) * t351 - pkin(11) * t429) * qJD(3) + (-t339 + (pkin(11) * t351 - t321) * t350) * qJD(4) + t420;
t580 = pkin(11) * t583 - t582;
t548 = Ifges(6,1) + Ifges(7,1);
t547 = Ifges(7,4) + Ifges(6,5);
t546 = Ifges(7,5) - Ifges(6,4);
t419 = pkin(8) * t435 + t337;
t297 = t419 * qJD(1);
t333 = qJD(1) * t348 + qJD(2);
t256 = t333 * pkin(9) + t297;
t291 = (-pkin(2) * t356 - pkin(9) * t352 - pkin(1)) * t347;
t269 = qJD(1) * t291;
t191 = -t351 * t256 + t355 * t269;
t276 = t333 * t355 - t351 * t402;
t277 = t333 * t351 + t355 * t402;
t212 = pkin(3) * t277 - pkin(10) * t276;
t120 = -t191 * t350 + t354 * t212;
t514 = -pkin(11) - pkin(10);
t403 = qJD(4) * t514;
t579 = -pkin(4) * t277 - t120 + (pkin(11) * t276 + t403) * t354;
t121 = t354 * t191 + t350 * t212;
t437 = t276 * t350;
t578 = -pkin(11) * t437 - t350 * t403 + t121;
t353 = cos(qJ(5));
t349 = sin(qJ(5));
t434 = t349 * t354;
t313 = t350 * t353 + t434;
t300 = t313 * t351;
t312 = t349 * t350 - t353 * t354;
t529 = qJD(4) + qJD(5);
t181 = -t300 * t529 - t312 * t414;
t197 = t264 * t349 + t265 * t353;
t425 = t181 - t197;
t431 = t351 * t354;
t432 = t350 * t351;
t182 = t414 * t434 + (t431 * t529 + t398) * t353 - t529 * t349 * t432;
t196 = -t353 * t264 + t265 * t349;
t424 = t182 - t196;
t577 = t389 - t415;
t311 = t354 * t321;
t237 = -pkin(11) * t431 + t311 + (-pkin(4) - t479) * t355;
t285 = t350 * t321 + t339;
t250 = -pkin(11) * t432 + t285;
t533 = t349 * t237 + t353 * t250;
t570 = -qJD(5) * t533 + t580 * t349 + t353 * t581;
t410 = qJD(5) * t353;
t411 = qJD(5) * t349;
t569 = t237 * t410 - t250 * t411 + t349 * t581 - t580 * t353;
t194 = t313 * t276;
t241 = t529 * t313;
t423 = t194 - t241;
t195 = t312 * t276;
t240 = t529 * t312;
t422 = -t195 + t240;
t217 = -t351 * t294 + t295 * t355;
t202 = -pkin(3) * t402 - t217;
t563 = pkin(4) * t583 + pkin(9) * t414 - t202;
t324 = qJD(3) - t401;
t223 = t277 * t354 + t324 * t350;
t370 = t277 * t350 - t324 * t354;
t362 = t353 * t223 - t349 * t370;
t416 = qJD(2) * t356;
t399 = t355 * t416;
t230 = t333 * t414 + (-t352 * t415 + t399) * t418;
t417 = qJD(2) * t347;
t391 = qJD(1) * t417;
t388 = t352 * t391;
t139 = -qJD(4) * t370 + t230 * t354 + t350 * t388;
t140 = -qJD(4) * t223 - t230 * t350 + t354 * t388;
t149 = t349 * t223 + t353 * t370;
t58 = -qJD(5) * t149 + t353 * t139 + t349 * t140;
t59 = qJD(5) * t362 + t349 * t139 - t353 * t140;
t296 = qJD(2) * t368;
t286 = qJD(1) * t296;
t436 = t347 * t352;
t334 = pkin(8) * t436;
t307 = -t334 + t409;
t298 = t307 * qJD(2);
t287 = qJD(1) * t298;
t124 = -t256 * t414 - t269 * t415 + t286 * t355 - t351 * t287;
t113 = -pkin(3) * t388 - t124;
t74 = -pkin(4) * t140 + t113;
t11 = pkin(5) * t59 - qJ(6) * t58 - qJD(6) * t362 + t74;
t400 = t351 * t416;
t231 = t333 * t415 + (t352 * t414 + t400) * t418;
t271 = qJD(4) - t276;
t260 = qJD(5) + t271;
t255 = -t333 * pkin(2) - t294;
t165 = -t276 * pkin(3) - t277 * pkin(10) + t255;
t192 = t355 * t256 + t351 * t269;
t171 = pkin(10) * t324 + t192;
t101 = t350 * t165 + t354 * t171;
t123 = -t256 * t415 + t269 * t414 + t351 * t286 + t355 * t287;
t112 = pkin(10) * t388 + t123;
t299 = t419 * qJD(2);
t288 = qJD(1) * t299;
t135 = t231 * pkin(3) - t230 * pkin(10) + t288;
t36 = -qJD(4) * t101 - t112 * t350 + t354 * t135;
t22 = pkin(4) * t231 - pkin(11) * t139 + t36;
t35 = t354 * t112 + t350 * t135 + t165 * t412 - t171 * t413;
t24 = pkin(11) * t140 + t35;
t100 = t354 * t165 - t171 * t350;
t86 = -pkin(11) * t223 + t100;
t72 = pkin(4) * t271 + t86;
t87 = -pkin(11) * t370 + t101;
t5 = t349 * t22 + t353 * t24 + t72 * t410 - t411 * t87;
t2 = qJ(6) * t231 + qJD(6) * t260 + t5;
t499 = t231 / 0.2e1;
t517 = t59 / 0.2e1;
t518 = -t59 / 0.2e1;
t519 = t58 / 0.2e1;
t556 = -t231 / 0.2e1;
t522 = mrSges(6,1) * t74 + mrSges(7,1) * t11 - mrSges(7,2) * t2 - mrSges(6,3) * t5 - t58 * Ifges(6,4) / 0.2e1 + Ifges(6,6) * t556 + Ifges(7,6) * t499 + 0.2e1 * Ifges(7,3) * t517 + (-t518 + t517) * Ifges(6,2) + (t546 + Ifges(7,5)) * t519;
t438 = t353 * t87;
t30 = t349 * t72 + t438;
t6 = -qJD(5) * t30 + t22 * t353 - t24 * t349;
t3 = -pkin(5) * t231 - t6;
t576 = -mrSges(6,2) * t74 - mrSges(7,2) * t3 + mrSges(6,3) * t6 + mrSges(7,3) * t11 - Ifges(6,4) * t518 - Ifges(7,5) * t517 - t547 * t499 - t548 * t519;
t444 = t324 * Ifges(4,6);
t447 = t277 * Ifges(4,4);
t450 = t276 * Ifges(4,2);
t187 = t444 + t447 + t450;
t575 = -t187 / 0.2e1;
t146 = Ifges(6,4) * t149;
t460 = Ifges(7,5) * t149;
t541 = t547 * t260 + t362 * t548 - t146 + t460;
t574 = t541 / 0.2e1;
t542 = t547 * t231 + t546 * t59 + t548 * t58;
t573 = t542 / 0.2e1;
t545 = -Ifges(6,6) + Ifges(7,6);
t572 = -qJ(6) * t577 - qJD(6) * t355 + t569;
t571 = pkin(5) * t577 - t570;
t145 = Ifges(7,5) * t362;
t77 = Ifges(7,6) * t260 + Ifges(7,3) * t149 + t145;
t462 = Ifges(6,4) * t362;
t80 = -Ifges(6,2) * t149 + Ifges(6,6) * t260 + t462;
t568 = t80 - t77;
t301 = t312 * t351;
t567 = pkin(5) * t424 - qJ(6) * t425 + qJD(6) * t301 + t563;
t566 = Ifges(4,5) * t230;
t564 = t287 * mrSges(3,2);
t325 = t514 * t350;
t326 = t514 * t354;
t369 = t353 * t325 + t326 * t349;
t537 = qJD(5) * t369 + t349 * t579 - t353 * t578;
t259 = t325 * t349 - t326 * t353;
t536 = -qJD(5) * t259 + t349 * t578 + t353 * t579;
t561 = -t124 * mrSges(4,1) + t123 * mrSges(4,2);
t445 = t324 * Ifges(4,5);
t560 = -t255 * mrSges(4,2) - t445 / 0.2e1;
t441 = t349 * t87;
t29 = t353 * t72 - t441;
t532 = -t29 + qJD(6);
t25 = -pkin(5) * t260 + t532;
t26 = qJ(6) * t260 + t30;
t365 = Ifges(5,6) * t370;
t559 = t255 * mrSges(4,1) + t100 * mrSges(5,1) + t29 * mrSges(6,1) - t25 * mrSges(7,1) - t101 * mrSges(5,2) - t30 * mrSges(6,2) + t26 * mrSges(7,3) - t444 / 0.2e1 - t365 / 0.2e1;
t90 = pkin(5) * t362 + qJ(6) * t149;
t554 = -t370 / 0.2e1;
t553 = t370 / 0.2e1;
t550 = -Ifges(3,6) * t333 / 0.2e1;
t544 = Ifges(6,3) + Ifges(7,2);
t14 = Ifges(6,5) * t58 - Ifges(6,6) * t59 + Ifges(6,3) * t231;
t15 = Ifges(7,4) * t58 + Ifges(7,2) * t231 + Ifges(7,6) * t59;
t543 = t14 + t15;
t153 = pkin(4) * t437 + t192;
t539 = pkin(4) * t413 - pkin(5) * t423 + qJ(6) * t422 - qJD(6) * t313 - t153;
t538 = -qJ(6) * t277 + t537;
t535 = pkin(5) * t277 - t536;
t32 = t353 * t86 - t441;
t534 = pkin(4) * t410 + qJD(6) - t32;
t289 = t334 + (-pkin(2) - t480) * t348;
t302 = -t348 * t355 + t351 * t436;
t303 = t348 * t351 + t355 * t436;
t199 = t302 * pkin(3) - t303 * pkin(10) + t289;
t290 = pkin(9) * t348 + t419;
t214 = t355 * t290 + t351 * t291;
t201 = -pkin(10) * t435 + t214;
t115 = t350 * t199 + t354 * t201;
t531 = t35 * t354 - t350 * t36;
t530 = -t36 * mrSges(5,1) + t35 * mrSges(5,2);
t528 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3);
t458 = Ifges(5,3) * t271;
t461 = Ifges(5,5) * t223;
t125 = -t365 + t458 + t461;
t78 = Ifges(6,5) * t362 - t149 * Ifges(6,6) + t260 * Ifges(6,3);
t79 = Ifges(7,4) * t362 + t260 * Ifges(7,2) + t149 * Ifges(7,6);
t527 = t125 + t79 + t78;
t244 = -qJD(3) * t302 + t347 * t399;
t246 = -t350 * t303 - t354 * t435;
t396 = t352 * t417;
t163 = qJD(4) * t246 + t354 * t244 + t350 * t396;
t245 = qJD(3) * t303 + t347 * t400;
t141 = -t290 * t415 + t291 * t414 + t351 * t296 + t355 * t298;
t129 = pkin(10) * t396 + t141;
t157 = t245 * pkin(3) - t244 * pkin(10) + t299;
t50 = -qJD(4) * t115 - t129 * t350 + t354 * t157;
t28 = pkin(4) * t245 - pkin(11) * t163 + t50;
t367 = -t354 * t303 + t350 * t435;
t164 = qJD(4) * t367 - t350 * t244 + t354 * t396;
t49 = t354 * t129 + t350 * t157 + t199 * t412 - t201 * t413;
t34 = pkin(11) * t164 + t49;
t102 = pkin(11) * t246 + t115;
t114 = t354 * t199 - t201 * t350;
t97 = pkin(4) * t302 + pkin(11) * t367 + t114;
t471 = t353 * t102 + t349 * t97;
t10 = -qJD(5) * t471 + t28 * t353 - t34 * t349;
t170 = -t324 * pkin(3) - t191;
t374 = t100 * t354 + t101 * t350;
t463 = Ifges(5,4) * t354;
t381 = -Ifges(5,2) * t350 + t463;
t384 = mrSges(5,1) * t350 + mrSges(5,2) * t354;
t526 = -t374 * mrSges(5,3) + t170 * t384 + t381 * t554;
t406 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t407 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t408 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t524 = t406 * t149 + t407 * t260 - t408 * t362 - t192 * mrSges(4,3) + t575 + t458 / 0.2e1 - t450 / 0.2e1 + t125 / 0.2e1 + t461 / 0.2e1 + t79 / 0.2e1 + t78 / 0.2e1 - t447 / 0.2e1 + t559;
t122 = pkin(4) * t370 + t170;
t51 = t149 * pkin(5) - qJ(6) * t362 + t122;
t523 = t122 * mrSges(6,1) + t51 * mrSges(7,1) + t77 / 0.2e1 - t80 / 0.2e1 - t26 * mrSges(7,2) - t30 * mrSges(6,3);
t65 = t139 * Ifges(5,1) + t140 * Ifges(5,4) + t231 * Ifges(5,5);
t516 = t65 / 0.2e1;
t513 = pkin(1) * mrSges(3,1);
t512 = pkin(1) * mrSges(3,2);
t465 = Ifges(5,4) * t223;
t126 = -Ifges(5,2) * t370 + Ifges(5,6) * t271 + t465;
t511 = -t126 / 0.2e1;
t221 = Ifges(5,4) * t370;
t127 = t223 * Ifges(5,1) + t271 * Ifges(5,5) - t221;
t510 = t127 / 0.2e1;
t509 = t139 / 0.2e1;
t508 = t140 / 0.2e1;
t507 = -t149 / 0.2e1;
t506 = t149 / 0.2e1;
t504 = -t362 / 0.2e1;
t503 = t362 / 0.2e1;
t501 = -t223 / 0.2e1;
t500 = t223 / 0.2e1;
t497 = -t260 / 0.2e1;
t496 = t260 / 0.2e1;
t495 = -t271 / 0.2e1;
t494 = t271 / 0.2e1;
t493 = -t276 / 0.2e1;
t492 = t276 / 0.2e1;
t491 = -t277 / 0.2e1;
t490 = t277 / 0.2e1;
t488 = -t302 / 0.2e1;
t486 = t303 / 0.2e1;
t484 = t348 / 0.2e1;
t483 = -t350 / 0.2e1;
t482 = t354 / 0.2e1;
t481 = m(6) * t122;
t478 = pkin(9) * t355;
t470 = mrSges(4,3) * t276;
t469 = mrSges(4,3) * t277;
t468 = mrSges(6,3) * t149;
t467 = mrSges(6,3) * t362;
t466 = Ifges(3,4) * t352;
t270 = Ifges(4,4) * t276;
t464 = Ifges(5,4) * t350;
t138 = Ifges(5,5) * t139;
t137 = Ifges(5,6) * t140;
t453 = t230 * Ifges(4,1);
t452 = t230 * Ifges(4,4);
t451 = t231 * Ifges(4,4);
t448 = t277 * Ifges(4,1);
t433 = t350 * t126;
t430 = t354 * t127;
t118 = mrSges(6,1) * t260 - t467;
t119 = -mrSges(7,1) * t260 + mrSges(7,2) * t362;
t427 = t118 - t119;
t156 = mrSges(5,1) * t370 + t223 * mrSges(5,2);
t233 = mrSges(4,1) * t324 - t469;
t426 = t156 - t233;
t421 = -mrSges(3,1) * t333 - mrSges(4,1) * t276 + mrSges(4,2) * t277 + mrSges(3,3) * t402;
t317 = pkin(4) * t432 + t351 * pkin(9);
t63 = Ifges(5,3) * t231 + t137 + t138;
t404 = -Ifges(4,6) * t231 + Ifges(4,3) * t388 + t566;
t344 = -pkin(4) * t354 - pkin(3);
t42 = -t231 * mrSges(7,1) + t58 * mrSges(7,2);
t213 = -t351 * t290 + t291 * t355;
t200 = pkin(3) * t435 - t213;
t385 = mrSges(5,1) * t354 - mrSges(5,2) * t350;
t383 = Ifges(5,1) * t354 - t464;
t382 = Ifges(5,1) * t350 + t463;
t380 = Ifges(5,2) * t354 + t464;
t379 = Ifges(5,5) * t354 - Ifges(5,6) * t350;
t378 = Ifges(5,5) * t350 + Ifges(5,6) * t354;
t45 = -t102 * t349 + t353 * t97;
t172 = t237 * t353 - t250 * t349;
t371 = t353 * t246 + t349 * t367;
t178 = t246 * t349 - t353 * t367;
t142 = -t290 * t414 - t291 * t415 + t296 * t355 - t351 * t298;
t9 = -t102 * t411 + t349 * t28 + t353 * t34 + t97 * t410;
t152 = -pkin(4) * t246 + t200;
t327 = Ifges(3,4) * t401;
t363 = -t294 * mrSges(3,3) + Ifges(3,1) * t402 / 0.2e1 + t327 / 0.2e1 + t333 * Ifges(3,5);
t130 = -pkin(3) * t396 - t142;
t361 = t528 + t543;
t89 = -pkin(4) * t164 + t130;
t188 = t270 + t445 + t448;
t360 = t191 * mrSges(4,3) - t188 / 0.2e1 - t448 / 0.2e1 - t270 / 0.2e1 + t560;
t359 = t191 * mrSges(4,1) + t324 * Ifges(4,3) + t277 * Ifges(4,5) + t276 * Ifges(4,6) + t550 - (Ifges(3,2) * t356 + t466) * t418 / 0.2e1 - t192 * mrSges(4,2) - t297 * mrSges(3,3);
t357 = t383 * t500 + t379 * t494 - t433 / 0.2e1 + t430 / 0.2e1 + t526;
t343 = -pkin(4) * t353 - pkin(5);
t341 = pkin(4) * t349 + qJ(6);
t320 = Ifges(3,5) * t356 * t391;
t293 = -t333 * mrSges(3,2) + mrSges(3,3) * t401;
t284 = -t350 * t478 + t311;
t235 = pkin(5) * t312 - qJ(6) * t313 + t344;
t232 = -mrSges(4,2) * t324 + t470;
t211 = -qJD(4) * t285 + t420;
t208 = -mrSges(4,2) * t388 - mrSges(4,3) * t231;
t207 = mrSges(4,1) * t388 - mrSges(4,3) * t230;
t204 = pkin(5) * t300 + qJ(6) * t301 + t317;
t176 = mrSges(5,1) * t271 - mrSges(5,3) * t223;
t175 = -t271 * mrSges(5,2) - mrSges(5,3) * t370;
t168 = pkin(5) * t355 - t172;
t167 = -qJ(6) * t355 + t533;
t161 = mrSges(4,1) * t231 + mrSges(4,2) * t230;
t144 = Ifges(4,5) * t388 - t451 + t453;
t143 = -t231 * Ifges(4,2) + Ifges(4,6) * t388 + t452;
t117 = -mrSges(6,2) * t260 - t468;
t116 = -mrSges(7,2) * t149 + mrSges(7,3) * t260;
t108 = -mrSges(5,2) * t231 + mrSges(5,3) * t140;
t107 = mrSges(5,1) * t231 - mrSges(5,3) * t139;
t92 = mrSges(6,1) * t149 + mrSges(6,2) * t362;
t91 = mrSges(7,1) * t149 - mrSges(7,3) * t362;
t88 = -mrSges(5,1) * t140 + mrSges(5,2) * t139;
t73 = pkin(4) * t223 + t90;
t70 = -pkin(5) * t371 - qJ(6) * t178 + t152;
t67 = qJD(5) * t178 + t163 * t349 - t353 * t164;
t66 = qJD(5) * t371 + t163 * t353 + t164 * t349;
t64 = t139 * Ifges(5,4) + t140 * Ifges(5,2) + t231 * Ifges(5,6);
t44 = -mrSges(7,2) * t59 + mrSges(7,3) * t231;
t43 = -mrSges(6,2) * t231 - mrSges(6,3) * t59;
t41 = mrSges(6,1) * t231 - mrSges(6,3) * t58;
t38 = -pkin(5) * t302 - t45;
t37 = qJ(6) * t302 + t471;
t31 = t349 * t86 + t438;
t20 = mrSges(6,1) * t59 + mrSges(6,2) * t58;
t19 = mrSges(7,1) * t59 - mrSges(7,3) * t58;
t12 = pkin(5) * t67 - qJ(6) * t66 - qJD(6) * t178 + t89;
t8 = -pkin(5) * t245 - t10;
t7 = qJ(6) * t245 + qJD(6) * t302 + t9;
t1 = [(-t404 / 0.2e1 - t566 / 0.2e1 + t287 * mrSges(3,3) - Ifges(4,6) * t556 + t561) * t435 + (t122 * t66 + t178 * t74 - t245 * t30 - t302 * t5) * mrSges(6,2) + (Ifges(6,4) * t66 + Ifges(6,6) * t245) * t507 + (Ifges(6,4) * t178 + Ifges(6,6) * t302) * t518 + t230 * (Ifges(4,1) * t303 - Ifges(4,4) * t302) / 0.2e1 - t522 * t371 + (-Ifges(5,5) * t367 + Ifges(5,6) * t246 + t547 * t178 - t545 * t371 + (Ifges(5,3) + t544) * t302) * t499 + t36 * (mrSges(5,1) * t302 + mrSges(5,3) * t367) + t113 * (-mrSges(5,1) * t246 - mrSges(5,2) * t367) + (-Ifges(5,4) * t367 + Ifges(5,2) * t246 + Ifges(5,6) * t302) * t508 + (-Ifges(5,1) * t367 + Ifges(5,4) * t246 + Ifges(5,5) * t302) * t509 - t367 * t516 + (-t11 * t178 + t2 * t302 + t245 * t26 - t51 * t66) * mrSges(7,3) - t348 * t564 + (Ifges(7,5) * t66 + Ifges(7,6) * t245) * t506 + (Ifges(7,5) * t178 + Ifges(7,6) * t302) * t517 + m(4) * (t123 * t214 + t124 * t213 + t141 * t192 + t142 * t191 + t255 * t299 + t288 * t289) + m(5) * (t100 * t50 + t101 * t49 + t113 * t200 + t114 * t36 + t115 * t35 + t130 * t170) + m(7) * (t11 * t70 + t12 * t51 + t2 * t37 + t25 * t8 + t26 * t7 + t3 * t38) + t144 * t486 + t143 * t488 + (Ifges(4,1) * t244 - Ifges(4,4) * t245) * t490 + (Ifges(4,4) * t244 - Ifges(4,2) * t245) * t492 + (Ifges(5,5) * t163 + Ifges(5,6) * t164 + Ifges(5,3) * t245) * t494 + (-t123 * t302 - t124 * t303 - t191 * t244 - t192 * t245) * mrSges(4,3) + t324 * (Ifges(4,5) * t244 - Ifges(4,6) * t245) / 0.2e1 + m(6) * (t10 * t29 + t122 * t89 + t152 * t74 + t30 * t9 + t45 * t6 + t471 * t5) + t471 * t43 + m(3) * (t287 * t419 - t288 * t307 - t294 * t299 + t297 * t298) + ((Ifges(3,5) * t484 - t307 * mrSges(3,3) + (-0.2e1 * t512 + 0.3e1 / 0.2e1 * Ifges(3,4) * t356) * t347) * t356 + (-Ifges(3,6) * t348 + Ifges(4,5) * t486 + Ifges(4,6) * t488 - t419 * mrSges(3,3) + (-0.2e1 * t513 - 0.3e1 / 0.2e1 * t466 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1) * t356) * t347) * t352) * t391 + (t363 * t356 + (t550 + t359) * t352) * t417 + (t245 * t544 + t547 * t66) * t496 + (t178 * t548 + t302 * t547) * t519 + (t245 * t547 + t548 * t66) * t503 + (-Ifges(6,2) * t507 + Ifges(7,3) * t506 + t496 * t545 + t503 * t546 + t523) * t67 + (t63 + t543) * t302 / 0.2e1 + (-mrSges(3,1) * t348 + mrSges(4,1) * t302 + mrSges(4,2) * t303 + mrSges(3,3) * t436) * t288 + t421 * t299 + t527 * t245 / 0.2e1 + t35 * (-mrSges(5,2) * t302 + mrSges(5,3) * t246) + t3 * (-mrSges(7,1) * t302 + mrSges(7,2) * t178) + t6 * (mrSges(6,1) * t302 - mrSges(6,3) * t178) + t298 * t293 + t289 * t161 + t255 * (mrSges(4,1) * t245 + mrSges(4,2) * t244) + t246 * t64 / 0.2e1 + t244 * t188 / 0.2e1 + t100 * (mrSges(5,1) * t245 - mrSges(5,3) * t163) + t101 * (-mrSges(5,2) * t245 + mrSges(5,3) * t164) + t25 * (-mrSges(7,1) * t245 + mrSges(7,2) * t66) + t29 * (mrSges(6,1) * t245 - mrSges(6,3) * t66) + t141 * t232 + t142 * t233 + t213 * t207 + t214 * t208 + t200 * t88 + t49 * t175 + t50 * t176 + t170 * (-mrSges(5,1) * t164 + mrSges(5,2) * t163) + t164 * t126 / 0.2e1 + t130 * t156 + t152 * t20 + t10 * t118 + t8 * t119 + t7 * t116 + t9 * t117 + t178 * t573 + t66 * t574 + t245 * t575 + t320 * t484 + t114 * t107 + t115 * t108 + (Ifges(5,1) * t163 + Ifges(5,4) * t164 + Ifges(5,5) * t245) * t500 + t163 * t510 + (Ifges(5,4) * t163 + Ifges(5,2) * t164 + Ifges(5,6) * t245) * t554 + t38 * t42 + t37 * t44 + t45 * t41 + (Ifges(4,4) * t303 - Ifges(4,2) * t302) * t556 + t70 * t19 + t12 * t91 + t89 * t92; t567 * t91 + t541 * (-t197 / 0.2e1 + t181 / 0.2e1) + t568 * (t196 / 0.2e1 - t182 / 0.2e1) + t569 * t117 + (t122 * t563 + t172 * t6 + t29 * t570 + t30 * t569 + t317 * t74 + t533 * t5) * m(6) + t570 * t118 + t571 * t119 + (t11 * t204 + t167 * t2 + t168 * t3 + t25 * t571 + t26 * t572 + t51 * t567) * m(7) + t572 * t116 + (Ifges(6,4) * t181 + Ifges(7,5) * t197 - Ifges(6,2) * t182 + Ifges(7,3) * t196) * t507 + (Ifges(6,4) * t197 + Ifges(7,5) * t181 - Ifges(6,2) * t196 + Ifges(7,3) * t182) * t506 + t533 * t43 + t563 * t92 + (t499 * t545 + t522) * t300 - m(4) * (t191 * t217 + t192 * t218 + t255 * t297) - m(5) * (t100 * t131 + t101 * t132 + t170 * t202) + (t196 * t545 + t197 * t547) * t497 + (t181 * t547 + t182 * t545) * t496 + (t181 * t548 + t182 * t546) * t503 + (t196 * t546 + t197 * t548) * t504 + (Ifges(5,5) * t265 + Ifges(5,6) * t264) * t495 + m(5) * (t100 * t211 + t101 * t210 + t284 * t36 + t285 * t35) + (t100 * t265 - t101 * t264) * mrSges(5,3) + t320 - t564 + (-t542 / 0.2e1 + t576) * t301 + (t211 - t131) * t176 + m(4) * (-pkin(2) * t288 + t123 * t478) + (t379 * t499 + t383 * t509 + t381 * t508 - t124 * mrSges(4,3) + t64 * t483 + t65 * t482 + t113 * t384 + t288 * mrSges(4,2) - t451 / 0.2e1 + t453 / 0.2e1 + t144 / 0.2e1 + (-t35 * t350 - t354 * t36) * mrSges(5,3) + (-m(4) * t124 + m(5) * t113 - t207 + t88) * pkin(9) + (t354 * t511 + t127 * t483 + t380 * t553 + t378 * t495 + t382 * t501 + t170 * t385 + (t100 * t350 - t101 * t354) * mrSges(5,3)) * qJD(4)) * t351 + (mrSges(6,1) * t424 + mrSges(6,2) * t425) * t122 + (mrSges(7,1) * t424 - mrSges(7,3) * t425) * t51 - t421 * t297 + t317 * t20 + (-t288 * mrSges(4,1) - t138 / 0.2e1 - t137 / 0.2e1 + pkin(9) * t208 + t123 * mrSges(4,3) + t143 / 0.2e1 - t14 / 0.2e1 - t15 / 0.2e1 - t63 / 0.2e1 + t452 / 0.2e1 - t406 * t59 + t408 * t58 + (-Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1 - t407) * t231 - t528 + t530) * t355 - t294 * t293 - t288 * mrSges(3,1) + t284 * t107 + t285 * t108 - t265 * t127 / 0.2e1 - t170 * (-mrSges(5,1) * t264 + mrSges(5,2) * t265) + (((-m(4) * t192 - t232) * pkin(9) + t524) * t351 + (t357 - t360 + (-m(4) * t191 + m(5) * t170 + t426) * pkin(9)) * t355) * qJD(3) + ((qJD(2) * (Ifges(4,5) * t351 + Ifges(4,6) * t355) / 0.2e1 + (t513 + t466 / 0.2e1) * t418 + (t333 / 0.2e1 - qJD(2)) * Ifges(3,6) - t359) * t352 + (-t327 / 0.2e1 + (t512 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t352) * t418 + t360 * t355 - t524 * t351 - t363) * t356) * t418 - t218 * t232 - t217 * t233 - t202 * t156 + t204 * t19 + t168 * t42 + t172 * t41 + t167 * t44 - pkin(2) * t161 + t582 * t175 + (Ifges(5,1) * t265 + Ifges(5,4) * t264) * t501 + t264 * t511 + (Ifges(5,4) * t265 + Ifges(5,2) * t264) * t553 + (t25 * t425 - t26 * t424) * mrSges(7,2) + (-t29 * t425 - t30 * t424) * mrSges(6,3); (-pkin(3) * t113 - t100 * t120 - t101 * t121 - t170 * t192) * m(5) + t541 * (-t240 / 0.2e1 + t195 / 0.2e1) + ((t92 + t481) * t350 * pkin(4) + t357) * qJD(4) - t568 * (t241 / 0.2e1 - t194 / 0.2e1) + (-t25 * t422 + t26 * t423) * mrSges(7,2) + (t29 * t422 + t30 * t423) * mrSges(6,3) - (t42 - t41) * t369 + (-t122 * t153 + t259 * t5 + t29 * t536 + t30 * t537 + t344 * t74 + t369 * t6) * m(6) + (t11 * t235 + t2 * t259 + t25 * t535 + t26 * t538 - t3 * t369 + t51 * t539) * m(7) + (-Ifges(6,4) * t240 - Ifges(7,5) * t195 - Ifges(6,2) * t241 + Ifges(7,3) * t194) * t507 + (-Ifges(6,4) * t195 - Ifges(7,5) * t240 - Ifges(6,2) * t194 + Ifges(7,3) * t241) * t506 + (-t240 * t548 + t241 * t546) * t503 + (-t240 * t547 + t241 * t545) * t496 + (Ifges(5,5) * t501 - Ifges(4,2) * t493 + Ifges(6,6) * t506 + Ifges(7,6) * t507 + Ifges(5,3) * t495 + t497 * t544 + t504 * t547 - t559) * t277 + (Ifges(4,1) * t491 + t379 * t495 + t383 * t501 - t526 + t560) * t276 + (t270 + t430 + t188) * t493 + (t312 * t545 + t378) * t499 + t187 * t490 + t433 * t492 - t561 + (t573 - t576) * t313 + (-t426 + t469) * t192 + (-t232 + t470) * t191 + (t194 * t545 - t195 * t547) * t497 + (t194 * t546 - t195 * t548) * t504 + t535 * t119 + t536 * t118 + t537 * t117 + t538 * t116 + t539 * t91 + (-mrSges(7,1) * t423 + mrSges(7,3) * t422) * t51 + (-mrSges(6,1) * t423 - mrSges(6,2) * t422) * t122 + t344 * t20 + t531 * mrSges(5,3) + (-t107 * t350 + t108 * t354 + (-m(5) * t374 - t350 * t175 - t354 * t176) * qJD(4) + m(5) * t531) * pkin(10) + (-t447 + t527) * t491 + t235 * t19 + t522 * t312 - t121 * t175 - t120 * t176 - t153 * t92 + t64 * t482 + t404 + t380 * t508 + t382 * t509 + t350 * t516 - t113 * t385 - pkin(3) * t88 + (t43 + t44) * t259; (t2 * t341 - t25 * t31 + t26 * t534 + t3 * t343 - t51 * t73) * m(7) + (m(7) * t25 * t411 + 0.2e1 * t481 * t501 + m(6) * (-t29 * t411 + t30 * t410 + t349 * t5 + t353 * t6) - t223 * t92 + t349 * t43 + t353 * t41 + (t117 * t353 - t349 * t427) * qJD(5)) * pkin(4) + t361 + t534 * t116 + (-Ifges(5,5) * t370 - Ifges(5,6) * t223) * t495 - t530 - m(6) * (-t29 * t31 + t30 * t32) + t63 + (-Ifges(5,2) * t223 - t221) * t553 + t427 * t31 + t341 * t44 + t343 * t42 + (-Ifges(6,2) * t506 + Ifges(7,3) * t507 + t497 * t545 + t504 * t546 - t523) * t362 - t100 * t175 + t101 * t176 - t32 * t117 + t126 * t500 + (-Ifges(5,1) * t370 - t465) * t501 + t370 * t510 + (-t100 * t370 + t101 * t223) * mrSges(5,3) - t170 * (t223 * mrSges(5,1) - mrSges(5,2) * t370) - t73 * t91 + (t122 * mrSges(6,2) + t25 * mrSges(7,2) - t29 * mrSges(6,3) - t51 * mrSges(7,3) - Ifges(6,4) * t506 - Ifges(7,5) * t507 - t547 * t497 - t548 * t504 + t574) * t149; t361 + (t427 + t467) * t30 + (-t116 - t117 - t468) * t29 + (t149 * t25 + t26 * t362) * mrSges(7,2) - t122 * (mrSges(6,1) * t362 - mrSges(6,2) * t149) + t80 * t503 + (Ifges(7,3) * t362 - t460) * t507 - t51 * (mrSges(7,1) * t362 + mrSges(7,3) * t149) + qJD(6) * t116 - pkin(5) * t42 + qJ(6) * t44 - t90 * t91 + (-t149 * t547 + t362 * t545) * t497 + (-pkin(5) * t3 + qJ(6) * t2 - t25 * t30 + t26 * t532 - t51 * t90) * m(7) + (-Ifges(6,2) * t362 - t146 + t541) * t506 + (-t149 * t548 + t145 - t462 + t77) * t504; -t260 * t116 + t362 * t91 + 0.2e1 * (t3 / 0.2e1 + t51 * t503 + t26 * t497) * m(7) + t42;];
tauc  = t1(:);
