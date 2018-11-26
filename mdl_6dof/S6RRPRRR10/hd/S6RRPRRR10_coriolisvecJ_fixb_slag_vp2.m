% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:27:30
% EndTime: 2018-11-23 17:28:03
% DurationCPUTime: 33.54s
% Computational Cost: add. (29553->926), mult. (76598->1309), div. (0->0), fcn. (63007->12), ass. (0->410)
t352 = sin(pkin(12));
t354 = cos(pkin(12));
t358 = sin(qJ(4));
t362 = cos(qJ(4));
t321 = t352 * t362 + t354 * t358;
t353 = sin(pkin(6));
t363 = cos(qJ(2));
t430 = t353 * t363;
t369 = t321 * t430;
t275 = qJD(1) * t369;
t316 = t321 * qJD(4);
t420 = t275 - t316;
t320 = t352 * t358 - t362 * t354;
t368 = t320 * t430;
t276 = qJD(1) * t368;
t315 = t320 * qJD(4);
t419 = -t276 + t315;
t359 = sin(qJ(2));
t384 = pkin(2) * t359 - qJ(3) * t363;
t418 = qJD(1) * t353;
t308 = t384 * t418;
t407 = t359 * t418;
t355 = cos(pkin(6));
t417 = qJD(1) * t355;
t410 = pkin(1) * t417;
t309 = -pkin(8) * t407 + t363 * t410;
t236 = t354 * t308 - t352 * t309;
t429 = t354 * t363;
t370 = (pkin(3) * t359 - pkin(9) * t429) * t353;
t206 = qJD(1) * t370 + t236;
t237 = t352 * t308 + t354 * t309;
t406 = t363 * t418;
t398 = t352 * t406;
t218 = -pkin(9) * t398 + t237;
t471 = pkin(9) + qJ(3);
t331 = t471 * t352;
t332 = t471 * t354;
t523 = -t362 * t331 - t332 * t358;
t540 = -t320 * qJD(3) + qJD(4) * t523 - t358 * t206 - t362 * t218;
t565 = pkin(10) * t407 - t540;
t310 = pkin(8) * t406 + t359 * t410;
t268 = pkin(3) * t398 + t310;
t564 = -t420 * pkin(4) + t419 * pkin(10) - t268;
t357 = sin(qJ(5));
t361 = cos(qJ(5));
t238 = t276 * t357 + t361 * t407;
t428 = t357 * t315;
t400 = t238 - t428;
t412 = qJD(5) * t361;
t563 = t321 * t412 + t400;
t562 = t565 * t357 + t564 * t361;
t350 = -pkin(3) * t354 - pkin(2);
t253 = pkin(4) * t320 - pkin(10) * t321 + t350;
t274 = -t331 * t358 + t332 * t362;
t413 = qJD(5) * t357;
t542 = t253 * t412 - t274 * t413 + t564 * t357 - t565 * t361;
t318 = t355 * t359 * pkin(1) + pkin(8) * t430;
t312 = t318 * qJD(2);
t301 = qJD(1) * t312;
t416 = qJD(2) * t353;
t403 = qJD(1) * t416;
t395 = t363 * t403;
t376 = t352 * t395;
t256 = pkin(3) * t376 + t301;
t365 = qJD(2) * t368;
t344 = qJD(2) + t417;
t295 = t344 * t354 - t352 * t407;
t296 = t344 * t352 + t354 * t407;
t524 = t362 * t295 - t358 * t296;
t176 = -qJD(1) * t365 + qJD(4) * t524;
t450 = t176 * Ifges(5,4);
t279 = qJ(3) * t344 + t310;
t303 = (-pkin(2) * t363 - qJ(3) * t359 - pkin(1)) * t353;
t288 = qJD(1) * t303;
t207 = -t352 * t279 + t354 * t288;
t158 = -pkin(3) * t406 - t296 * pkin(9) + t207;
t282 = (qJD(2) * t384 - qJD(3) * t359) * t353;
t264 = qJD(1) * t282;
t477 = pkin(1) * t363;
t411 = t355 * t477;
t343 = qJD(2) * t411;
t396 = t359 * t403;
t300 = -pkin(8) * t396 + qJD(1) * t343;
t265 = qJD(3) * t344 + t300;
t200 = t354 * t264 - t352 * t265;
t367 = qJD(2) * t370;
t159 = qJD(1) * t367 + t200;
t208 = t354 * t279 + t352 * t288;
t178 = pkin(9) * t295 + t208;
t201 = t352 * t264 + t354 * t265;
t179 = -pkin(9) * t376 + t201;
t414 = qJD(4) * t362;
t415 = qJD(4) * t358;
t51 = t158 * t414 + t358 * t159 - t178 * t415 + t362 * t179;
t270 = -t344 * pkin(2) + qJD(3) - t309;
t224 = -t295 * pkin(3) + t270;
t379 = t295 * t358 + t362 * t296;
t112 = -pkin(4) * t524 - pkin(10) * t379 + t224;
t49 = pkin(10) * t396 + t51;
t366 = qJD(2) * t369;
t177 = qJD(1) * t366 + qJD(4) * t379;
t83 = t177 * pkin(4) - t176 * pkin(10) + t256;
t336 = qJD(4) - t406;
t98 = t158 * t358 + t178 * t362;
t93 = pkin(10) * t336 + t98;
t13 = t112 * t412 + t357 * t83 + t361 * t49 - t413 * t93;
t54 = t112 * t357 + t361 * t93;
t14 = -qJD(5) * t54 - t357 * t49 + t361 * t83;
t360 = cos(qJ(6));
t223 = qJD(5) - t524;
t198 = t336 * t357 + t361 * t379;
t53 = t361 * t112 - t357 * t93;
t42 = -pkin(11) * t198 + t53;
t37 = pkin(5) * t223 + t42;
t356 = sin(qJ(6));
t197 = t336 * t361 - t357 * t379;
t43 = pkin(11) * t197 + t54;
t439 = t356 * t43;
t15 = t360 * t37 - t439;
t106 = qJD(5) * t197 + t176 * t361 + t357 * t396;
t6 = pkin(5) * t177 - pkin(11) * t106 + t14;
t107 = -qJD(5) * t198 - t176 * t357 + t361 * t396;
t7 = pkin(11) * t107 + t13;
t2 = qJD(6) * t15 + t356 * t6 + t360 * t7;
t438 = t360 * t43;
t16 = t356 * t37 + t438;
t3 = -qJD(6) * t16 - t356 * t7 + t360 * t6;
t547 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t534 = t14 * mrSges(6,1) - t13 * mrSges(6,2) + t547;
t561 = -t534 - t256 * mrSges(5,1) + t51 * mrSges(5,3) + t450 / 0.2e1;
t239 = -t276 * t361 + t357 * t407;
t260 = t361 * t274;
t427 = t361 * t315;
t560 = pkin(11) * t239 + pkin(11) * t427 + (-t260 + (pkin(11) * t321 - t253) * t357) * qJD(5) + t562 - t420 * pkin(5);
t559 = pkin(11) * t563 - t542;
t539 = -qJD(3) * t321 - qJD(4) * t274 - t206 * t362 + t358 * t218;
t52 = -t158 * t415 + t159 * t362 - t178 * t414 - t358 * t179;
t557 = t52 * mrSges(5,3);
t556 = t177 * Ifges(5,2);
t555 = t256 * mrSges(5,2);
t512 = -pkin(11) - pkin(10);
t408 = qJD(5) * t512;
t435 = t524 * t357;
t151 = pkin(4) * t379 - pkin(10) * t524;
t97 = t158 * t362 - t358 * t178;
t68 = t357 * t151 + t361 * t97;
t554 = pkin(11) * t435 + t357 * t408 - t68;
t476 = pkin(11) * t361;
t67 = t361 * t151 - t357 * t97;
t553 = -pkin(5) * t379 + t361 * t408 + t476 * t524 - t67;
t326 = t356 * t361 + t357 * t360;
t142 = t326 * t524;
t521 = qJD(5) + qJD(6);
t263 = t521 * t326;
t525 = t142 - t263;
t377 = t356 * t357 - t360 * t361;
t144 = t377 * t524;
t262 = t521 * t377;
t526 = t144 - t262;
t538 = pkin(4) * t407 - t539;
t402 = t360 * t197 - t198 * t356;
t116 = Ifges(7,4) * t402;
t122 = t197 * t356 + t198 * t360;
t462 = Ifges(7,4) * t122;
t219 = qJD(6) + t223;
t497 = -t219 / 0.2e1;
t505 = -t122 / 0.2e1;
t507 = -t402 / 0.2e1;
t58 = Ifges(7,1) * t122 + Ifges(7,5) * t219 + t116;
t92 = -pkin(4) * t336 - t97;
t73 = -pkin(5) * t197 + t92;
t550 = (Ifges(7,5) * t402 - Ifges(7,6) * t122) * t497 + (t122 * t16 + t15 * t402) * mrSges(7,3) + (-Ifges(7,2) * t122 + t116 + t58) * t507 - t73 * (mrSges(7,1) * t122 + mrSges(7,2) * t402) + (Ifges(7,1) * t402 - t462) * t505;
t222 = Ifges(5,4) * t524;
t459 = Ifges(5,5) * t336;
t141 = Ifges(5,1) * t379 + t222 + t459;
t549 = t224 * mrSges(5,2) + t141 / 0.2e1 - t97 * mrSges(5,3);
t442 = t336 * Ifges(5,6);
t466 = Ifges(5,4) * t379;
t140 = Ifges(5,2) * t524 + t442 + t466;
t548 = t16 * mrSges(7,2) + t54 * mrSges(6,2) + t98 * mrSges(5,3) + t140 / 0.2e1 - t15 * mrSges(7,1) - t224 * mrSges(5,1) - t53 * mrSges(6,1);
t189 = t361 * t253 - t274 * t357;
t152 = pkin(5) * t320 - t321 * t476 + t189;
t190 = t357 * t253 + t260;
t434 = t321 * t357;
t164 = -pkin(11) * t434 + t190;
t87 = t152 * t360 - t164 * t356;
t546 = qJD(6) * t87 + t356 * t560 - t559 * t360;
t88 = t152 * t356 + t164 * t360;
t545 = -qJD(6) * t88 + t559 * t356 + t360 * t560;
t446 = t219 * Ifges(7,3);
t454 = t122 * Ifges(7,5);
t455 = t402 * Ifges(7,6);
t56 = t446 + t454 + t455;
t445 = t223 * Ifges(6,3);
t447 = t198 * Ifges(6,5);
t448 = t197 * Ifges(6,6);
t89 = t445 + t447 + t448;
t544 = t89 + t56;
t543 = -qJD(5) * t190 + t562;
t541 = pkin(5) * t563 + t538;
t35 = qJD(6) * t402 + t106 * t360 + t107 * t356;
t515 = t35 / 0.2e1;
t36 = -qJD(6) * t122 - t106 * t356 + t107 * t360;
t514 = t36 / 0.2e1;
t57 = Ifges(7,2) * t402 + Ifges(7,6) * t219 + t462;
t536 = t57 / 0.2e1;
t503 = t177 / 0.2e1;
t103 = Ifges(6,6) * t107;
t104 = Ifges(6,5) * t106;
t39 = Ifges(6,3) * t177 + t103 + t104;
t33 = Ifges(7,6) * t36;
t34 = Ifges(7,5) * t35;
t8 = Ifges(7,3) * t177 + t33 + t34;
t535 = t39 + t8;
t337 = t512 * t357;
t338 = t512 * t361;
t283 = t337 * t360 + t338 * t356;
t529 = qJD(6) * t283 + t356 * t553 + t360 * t554;
t284 = t337 * t356 - t338 * t360;
t528 = -qJD(6) * t284 - t356 * t554 + t360 * t553;
t69 = -mrSges(7,1) * t402 + mrSges(7,2) * t122;
t527 = m(7) * t73 + t69;
t246 = t377 * t321;
t302 = qJ(3) * t355 + t318;
t231 = -t352 * t302 + t354 * t303;
t431 = t353 * t359;
t314 = t352 * t355 + t354 * t431;
t185 = -pkin(3) * t430 - t314 * pkin(9) + t231;
t232 = t354 * t302 + t352 * t303;
t313 = -t352 * t431 + t354 * t355;
t203 = pkin(9) * t313 + t232;
t118 = t358 * t185 + t362 * t203;
t115 = -pkin(10) * t430 + t118;
t241 = t313 * t358 + t314 * t362;
t345 = pkin(8) * t431;
t305 = t345 + (-pkin(2) - t477) * t355;
t249 = -t313 * pkin(3) + t305;
t378 = t362 * t313 - t314 * t358;
t138 = -pkin(4) * t378 - t241 * pkin(10) + t249;
t71 = t361 * t115 + t357 * t138;
t522 = t13 * t361 - t14 * t357;
t441 = t336 * Ifges(5,3);
t443 = t379 * Ifges(5,5);
t444 = t524 * Ifges(5,6);
t139 = t441 + t443 + t444;
t382 = t357 * t54 + t361 * t53;
t386 = Ifges(6,5) * t361 - Ifges(6,6) * t357;
t463 = Ifges(6,4) * t361;
t388 = -Ifges(6,2) * t357 + t463;
t464 = Ifges(6,4) * t357;
t391 = Ifges(6,1) * t361 - t464;
t393 = mrSges(6,1) * t357 + mrSges(6,2) * t361;
t479 = t361 / 0.2e1;
t482 = -t357 / 0.2e1;
t494 = t223 / 0.2e1;
t498 = t198 / 0.2e1;
t500 = t197 / 0.2e1;
t465 = Ifges(6,4) * t198;
t90 = Ifges(6,2) * t197 + Ifges(6,6) * t223 + t465;
t192 = Ifges(6,4) * t197;
t91 = Ifges(6,1) * t198 + Ifges(6,5) * t223 + t192;
t520 = -t382 * mrSges(6,3) + t386 * t494 + t388 * t500 + t391 * t498 + t393 * t92 + t479 * t91 + t482 * t90;
t519 = -0.2e1 * pkin(1);
t518 = Ifges(7,4) * t515 + Ifges(7,2) * t514 + Ifges(7,6) * t503;
t517 = Ifges(5,2) / 0.2e1;
t516 = Ifges(7,1) * t515 + Ifges(7,4) * t514 + Ifges(7,5) * t503;
t41 = t106 * Ifges(6,1) + t107 * Ifges(6,4) + t177 * Ifges(6,5);
t513 = t41 / 0.2e1;
t509 = t106 / 0.2e1;
t508 = t107 / 0.2e1;
t506 = t402 / 0.2e1;
t504 = t122 / 0.2e1;
t501 = -t197 / 0.2e1;
t499 = -t198 / 0.2e1;
t496 = t219 / 0.2e1;
t495 = -t223 / 0.2e1;
t493 = t524 / 0.2e1;
t492 = t379 / 0.2e1;
t491 = t378 / 0.2e1;
t489 = t241 / 0.2e1;
t488 = t313 / 0.2e1;
t487 = t314 / 0.2e1;
t486 = t336 / 0.2e1;
t485 = -t352 / 0.2e1;
t484 = t354 / 0.2e1;
t483 = t355 / 0.2e1;
t480 = t359 / 0.2e1;
t475 = t51 * mrSges(5,2);
t474 = t52 * mrSges(5,1);
t473 = t97 * mrSges(5,1);
t472 = t98 * mrSges(5,2);
t470 = mrSges(4,2) * t354;
t469 = Ifges(3,4) * t359;
t468 = Ifges(4,4) * t352;
t467 = Ifges(4,4) * t354;
t461 = Ifges(4,5) * t296;
t460 = Ifges(4,5) * t359;
t458 = Ifges(3,6) * t344;
t457 = Ifges(4,6) * t295;
t456 = Ifges(4,6) * t359;
t451 = t176 * Ifges(5,1);
t449 = t177 * Ifges(5,4);
t440 = t344 * Ifges(3,5);
t437 = t363 * Ifges(3,2);
t432 = t352 * t363;
t129 = -t263 * t321 + t315 * t377;
t161 = t238 * t356 + t239 * t360;
t426 = t129 - t161;
t130 = t246 * t521 + t326 * t315;
t160 = t238 * t360 - t239 * t356;
t425 = t130 - t160;
t123 = -mrSges(6,1) * t197 + mrSges(6,2) * t198;
t205 = mrSges(5,1) * t336 - mrSges(5,3) * t379;
t422 = t205 - t123;
t421 = -mrSges(3,1) * t344 - mrSges(4,1) * t295 + mrSges(4,2) * t296 + mrSges(3,3) * t407;
t405 = t359 * t416;
t311 = -pkin(8) * t405 + t343;
t294 = qJD(3) * t355 + t311;
t217 = t352 * t282 + t354 * t294;
t280 = mrSges(4,1) * t376 + t395 * t470;
t409 = Ifges(5,5) * t176 - Ifges(5,6) * t177 + Ifges(5,3) * t396;
t105 = t177 * mrSges(5,1) + t176 * mrSges(5,2);
t70 = -t115 * t357 + t361 * t138;
t117 = t185 * t362 - t358 * t203;
t399 = -t239 - t427;
t216 = t354 * t282 - t352 * t294;
t397 = t416 * t432;
t114 = pkin(4) * t430 - t117;
t394 = mrSges(6,1) * t361 - mrSges(6,2) * t357;
t392 = Ifges(4,1) * t354 - t468;
t390 = Ifges(6,1) * t357 + t463;
t389 = -Ifges(4,2) * t352 + t467;
t387 = Ifges(6,2) * t361 + t464;
t385 = Ifges(6,5) * t357 + Ifges(6,6) * t361;
t383 = -t13 * t357 - t14 * t361;
t373 = -t361 * t241 + t357 * t430;
t47 = -pkin(5) * t378 + pkin(11) * t373 + t70;
t220 = -t357 * t241 - t361 * t430;
t60 = pkin(11) * t220 + t71;
t24 = -t356 * t60 + t360 * t47;
t25 = t356 * t47 + t360 * t60;
t381 = t357 * t53 - t361 * t54;
t134 = -mrSges(6,2) * t223 + mrSges(6,3) * t197;
t135 = mrSges(6,1) * t223 - mrSges(6,3) * t198;
t380 = t134 * t361 - t135 * t357;
t146 = t220 * t360 + t356 * t373;
t147 = t220 * t356 - t360 * t373;
t183 = t216 + t367;
t199 = -pkin(9) * t397 + t217;
t66 = t183 * t362 - t185 * t415 - t358 * t199 - t203 * t414;
t375 = mrSges(4,1) * t359 - mrSges(4,3) * t429;
t374 = -mrSges(4,2) * t359 - mrSges(4,3) * t432;
t193 = qJD(4) * t378 - t365;
t194 = qJD(4) * t241 + t366;
t269 = pkin(3) * t397 + t312;
t102 = t194 * pkin(4) - t193 * pkin(10) + t269;
t65 = t358 * t183 + t185 * t414 + t362 * t199 - t203 * t415;
t62 = pkin(10) * t405 + t65;
t20 = t357 * t102 - t115 * t413 + t138 * t412 + t361 * t62;
t63 = -pkin(4) * t405 - t66;
t21 = -qJD(5) * t71 + t361 * t102 - t357 * t62;
t50 = -pkin(4) * t396 - t52;
t351 = -pkin(5) * t361 - pkin(4);
t339 = Ifges(3,4) * t406;
t334 = Ifges(3,5) * t395;
t317 = -t345 + t411;
t307 = -t344 * mrSges(3,2) + mrSges(3,3) * t406;
t290 = t375 * t403;
t289 = t374 * t403;
t278 = Ifges(3,1) * t407 + t339 + t440;
t277 = t458 + (t437 + t469) * t418;
t259 = -mrSges(4,1) * t406 - t296 * mrSges(4,3);
t258 = mrSges(4,2) * t406 + t295 * mrSges(4,3);
t251 = (t363 * t392 + t460) * t403;
t250 = (t363 * t389 + t456) * t403;
t245 = t326 * t321;
t235 = pkin(5) * t434 - t523;
t215 = Ifges(4,1) * t296 + Ifges(4,4) * t295 - Ifges(4,5) * t406;
t214 = Ifges(4,4) * t296 + Ifges(4,2) * t295 - Ifges(4,6) * t406;
t213 = -Ifges(4,3) * t406 + t457 + t461;
t204 = -mrSges(5,2) * t336 + mrSges(5,3) * t524;
t154 = -mrSges(5,2) * t396 - mrSges(5,3) * t177;
t153 = mrSges(5,1) * t396 - mrSges(5,3) * t176;
t150 = -mrSges(5,1) * t524 + mrSges(5,2) * t379;
t125 = qJD(5) * t373 - t357 * t193 + t361 * t405;
t124 = qJD(5) * t220 + t361 * t193 + t357 * t405;
t96 = Ifges(5,5) * t396 - t449 + t451;
t95 = Ifges(5,6) * t396 + t450 - t556;
t85 = mrSges(7,1) * t219 - mrSges(7,3) * t122;
t84 = -mrSges(7,2) * t219 + mrSges(7,3) * t402;
t80 = -pkin(5) * t220 + t114;
t76 = pkin(5) * t435 + t98;
t75 = -mrSges(6,2) * t177 + mrSges(6,3) * t107;
t74 = mrSges(6,1) * t177 - mrSges(6,3) * t106;
t59 = -mrSges(6,1) * t107 + mrSges(6,2) * t106;
t46 = -qJD(6) * t147 - t124 * t356 + t125 * t360;
t45 = qJD(6) * t146 + t124 * t360 + t125 * t356;
t40 = Ifges(6,4) * t106 + Ifges(6,2) * t107 + Ifges(6,6) * t177;
t38 = -pkin(5) * t125 + t63;
t32 = -pkin(5) * t107 + t50;
t29 = -mrSges(7,2) * t177 + mrSges(7,3) * t36;
t28 = mrSges(7,1) * t177 - mrSges(7,3) * t35;
t19 = t360 * t42 - t439;
t18 = -t356 * t42 - t438;
t17 = pkin(11) * t125 + t20;
t12 = pkin(5) * t194 - pkin(11) * t124 + t21;
t11 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t5 = -qJD(6) * t25 + t12 * t360 - t17 * t356;
t4 = qJD(6) * t24 + t12 * t356 + t17 * t360;
t1 = [-t241 * t557 + (-Ifges(6,5) * t373 + Ifges(7,5) * t147 + Ifges(6,6) * t220 + Ifges(7,6) * t146) * t503 + (-t556 / 0.2e1 + (-Ifges(6,3) - Ifges(7,3)) * t503 - Ifges(6,6) * t508 - Ifges(6,5) * t509 - Ifges(7,6) * t514 - Ifges(7,5) * t515 - t535 / 0.2e1 + t561) * t378 + (Ifges(5,1) * t492 + Ifges(5,4) * t493 + Ifges(5,5) * t486 + t549) * t193 + (Ifges(7,4) * t147 + Ifges(7,2) * t146) * t514 + (Ifges(7,4) * t45 + Ifges(7,2) * t46) * t506 - t430 * t474 - t177 * (Ifges(5,4) * t241 - Ifges(5,6) * t430) / 0.2e1 + t430 * t475 + ((-Ifges(3,6) * t355 + Ifges(4,5) * t487 + Ifges(4,6) * t488 + Ifges(5,5) * t489 + Ifges(5,6) * t491 - t318 * mrSges(3,3) + (mrSges(3,1) * t519 - 0.3e1 / 0.2e1 * t469) * t353) * t359 + (Ifges(3,5) * t483 - t317 * mrSges(3,3) + (Ifges(4,4) * t314 + Ifges(4,2) * t313) * t485 + (Ifges(4,1) * t314 + Ifges(4,4) * t313) * t484 + (mrSges(3,2) * t519 + (0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * t354 * Ifges(4,5) + 0.3e1 / 0.2e1 * t352 * Ifges(4,6)) * t363) * t353 + (-0.3e1 / 0.2e1 * Ifges(4,3) - Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t431) * t363) * t403 + (-Ifges(6,1) * t373 + Ifges(6,4) * t220) * t509 + (-t548 - Ifges(5,6) * t486 - Ifges(5,4) * t492 - Ifges(5,2) * t493 + Ifges(6,3) * t494 + Ifges(7,3) * t496 + Ifges(6,5) * t498 + Ifges(6,6) * t500 + Ifges(7,5) * t504 + Ifges(7,6) * t506 + t544 / 0.2e1) * t194 + (Ifges(6,4) * t124 + Ifges(6,2) * t125) * t500 + ((t295 * t389 / 0.2e1 + t296 * t392 / 0.2e1 + t440 / 0.2e1 + t270 * (mrSges(4,1) * t352 + t470) + t214 * t485 - t309 * mrSges(3,3) + t215 * t484 + t278 / 0.2e1 + (-t207 * t354 - t208 * t352) * mrSges(4,3)) * t363 + (t457 / 0.2e1 + t461 / 0.2e1 + t207 * mrSges(4,1) - t208 * mrSges(4,2) - t458 / 0.2e1 + t441 / 0.2e1 + t443 / 0.2e1 + t444 / 0.2e1 - t472 + t473 - t310 * mrSges(3,3) - t277 / 0.2e1 + t213 / 0.2e1 + t139 / 0.2e1) * t359) * t416 + t50 * (-mrSges(6,1) * t220 - mrSges(6,2) * t373) - t373 * t513 + t176 * (Ifges(5,1) * t241 - Ifges(5,5) * t430) / 0.2e1 + t334 * t483 + t251 * t487 + t250 * t488 + t96 * t489 + t95 * t491 + (t146 * t2 - t147 * t3 - t15 * t45 + t16 * t46) * mrSges(7,3) + (Ifges(7,1) * t147 + Ifges(7,4) * t146) * t515 + (Ifges(7,1) * t45 + Ifges(7,4) * t46) * t504 + (Ifges(6,1) * t124 + Ifges(6,4) * t125) * t498 + m(7) * (t15 * t5 + t16 * t4 + t2 * t25 + t24 * t3 + t32 * t80 + t38 * t73) + m(6) * (t114 * t50 + t13 * t71 + t14 * t70 + t20 * t54 + t21 * t53 + t63 * t92) + m(5) * (t117 * t52 + t118 * t51 + t224 * t269 + t249 * t256 + t65 * t98 + t66 * t97) + m(4) * (t200 * t231 + t201 * t232 + t207 * t216 + t208 * t217 + t270 * t312 + t301 * t305) + m(3) * (t300 * t318 - t301 * t317 - t309 * t312 + t310 * t311) + t201 * (mrSges(4,2) * t430 + t313 * mrSges(4,3)) + (-mrSges(3,1) * t355 - mrSges(4,1) * t313 + mrSges(4,2) * t314 + mrSges(3,3) * t431) * t301 - t409 * t430 / 0.2e1 + t200 * (-mrSges(4,1) * t430 - t314 * mrSges(4,3)) + t300 * (-t355 * mrSges(3,2) + mrSges(3,3) * t430) + t421 * t312 + t241 * t555 + t311 * t307 + t305 * t280 + t232 * t289 + t231 * t290 + t269 * t150 + t216 * t259 + t217 * t258 + t249 * t105 + t220 * t40 / 0.2e1 + t65 * t204 + t66 * t205 + (-t124 * t53 + t125 * t54 + t13 * t220 + t14 * t373) * mrSges(6,3) + (Ifges(7,5) * t45 + Ifges(7,6) * t46) * t496 + t147 * t516 + t146 * t518 + t24 * t28 + t25 * t29 + (Ifges(6,5) * t124 + Ifges(6,6) * t125) * t494 + t45 * t58 / 0.2e1 + t38 * t69 + t73 * (-mrSges(7,1) * t46 + mrSges(7,2) * t45) + t70 * t74 + t71 * t75 + t80 * t11 + t4 * t84 + t5 * t85 + t114 * t59 + t63 * t123 + t124 * t91 / 0.2e1 + t125 * t90 / 0.2e1 + t92 * (-mrSges(6,1) * t125 + mrSges(6,2) * t124) + t20 * t134 + t21 * t135 + t32 * (-mrSges(7,1) * t146 + mrSges(7,2) * t147) + t117 * t153 + t118 * t154 + (-Ifges(6,4) * t373 + Ifges(6,2) * t220) * t508 + t46 * t536; (t391 * t509 + t388 * t508 + t386 * t503 + t50 * t393 + t555 + t451 / 0.2e1 - t449 / 0.2e1 + t96 / 0.2e1 - t557 + t40 * t482 + t41 * t479 + t383 * mrSges(6,3) + (t92 * t394 + t385 * t495 + t387 * t501 + t390 * t499 + t91 * t482 - t361 * t90 / 0.2e1 + t381 * mrSges(6,3)) * qJD(5)) * t321 + (t13 * t190 + t14 * t189 - t50 * t523 + t53 * t543 + t538 * t92 + t54 * t542) * m(6) + (-t224 * t268 + t256 * t350 + t274 * t51 + t52 * t523 + t539 * t97 + t540 * t98) * m(5) - (t59 - t153) * t523 + (t104 / 0.2e1 + t103 / 0.2e1 + t8 / 0.2e1 + t39 / 0.2e1 + t34 / 0.2e1 + t33 / 0.2e1 - t95 / 0.2e1 + (t517 + Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t177 - t561) * t320 + (-t301 * mrSges(4,1) + t250 / 0.2e1 + t201 * mrSges(4,3) + qJD(3) * t258 + qJ(3) * t289) * t354 + t538 * t123 + (t301 * mrSges(4,2) + t251 / 0.2e1 - t200 * mrSges(4,3) - qJD(3) * t259 - qJ(3) * t290) * t352 - t524 * (-Ifges(5,4) * t276 - Ifges(5,2) * t275) / 0.2e1 - t336 * (-Ifges(5,5) * t276 - Ifges(5,6) * t275) / 0.2e1 - t379 * (-Ifges(5,1) * t276 - Ifges(5,4) * t275) / 0.2e1 + (-t315 / 0.2e1 + t276 / 0.2e1) * t141 + (t130 / 0.2e1 - t160 / 0.2e1) * t57 + (-t295 * (Ifges(4,4) * t429 - Ifges(4,2) * t432 + t456) / 0.2e1 - t296 * (Ifges(4,1) * t429 - Ifges(4,4) * t432 + t460) / 0.2e1 - t207 * t375 - t208 * t374 - t344 * (Ifges(3,5) * t363 - Ifges(3,6) * t359) / 0.2e1 - t270 * (mrSges(4,1) * t432 + mrSges(4,2) * t429) + t277 * t480 + (pkin(1) * (mrSges(3,1) * t359 + mrSges(3,2) * t363) + t437 * t480 + t363 * (Ifges(4,5) * t429 - Ifges(4,6) * t432 + Ifges(4,3) * t359) / 0.2e1) * t418 - t215 * t429 / 0.2e1 + t359 * t472 - t359 * t473 + t214 * t432 / 0.2e1 + (t309 * t363 + t310 * t359) * mrSges(3,3) + ((Ifges(5,5) * t321 / 0.2e1 - Ifges(5,6) * t320 / 0.2e1 + Ifges(4,5) * t352 / 0.2e1 + Ifges(4,6) * t484 - Ifges(3,6)) * t359 + ((Ifges(4,1) * t352 + t467) * t484 + (Ifges(4,2) * t354 + t468) * t485) * t363) * qJD(2) - (t339 + t278) * t363 / 0.2e1 - (t213 + 0.2e1 * t139 + (Ifges(3,1) * t363 - t469) * t418) * t359 / 0.2e1) * t418 + (Ifges(6,5) * t239 + Ifges(6,6) * t238 + Ifges(6,3) * t275) * t495 + (-t140 + t544) * (t316 / 0.2e1 - t275 / 0.2e1) + t545 * t85 + t546 * t84 + (t15 * t545 + t16 * t546 + t2 * t88 + t235 * t32 + t3 * t87 + t541 * t73) * m(7) + t541 * t69 + t542 * t134 + t543 * t135 + (t129 / 0.2e1 - t161 / 0.2e1) * t58 + t539 * t205 + t540 * t204 - t421 * t310 + (t419 * t97 + t420 * t98) * mrSges(5,3) + (mrSges(6,2) * t420 - mrSges(6,3) * t400) * t54 + (-mrSges(5,1) * t420 - mrSges(5,2) * t419) * t224 + (-mrSges(6,1) * t420 - mrSges(6,3) * t399) * t53 + (mrSges(6,1) * t400 + mrSges(6,2) * t399) * t92 + (-t207 * t236 - t208 * t237 - t270 * t310 - pkin(2) * t301 + (-t207 * t352 + t208 * t354) * qJD(3) + (-t200 * t352 + t201 * t354) * qJ(3)) * m(4) + t350 * t105 - t309 * t307 - t300 * mrSges(3,2) - t301 * mrSges(3,1) - pkin(2) * t280 + t274 * t154 - t268 * t150 - t236 * t259 - t237 * t258 + t235 * t11 + (Ifges(7,5) * t129 + Ifges(7,6) * t130 + Ifges(7,3) * t316) * t496 + (Ifges(7,5) * t161 + Ifges(7,6) * t160 + Ifges(7,3) * t275) * t497 + (Ifges(6,1) * t239 + Ifges(6,4) * t238 + Ifges(6,5) * t275) * t499 + (Ifges(6,4) * t239 + Ifges(6,2) * t238 + Ifges(6,6) * t275) * t501 + (-t15 * t420 + t245 * t32 - t425 * t73) * mrSges(7,1) + (t16 * t420 - t246 * t32 + t426 * t73) * mrSges(7,2) + (-t15 * t426 + t16 * t425 - t2 * t245 + t246 * t3) * mrSges(7,3) + (-Ifges(7,5) * t246 - Ifges(7,6) * t245) * t503 + (-Ifges(7,4) * t246 - Ifges(7,2) * t245) * t514 + (-Ifges(7,1) * t246 - Ifges(7,4) * t245) * t515 + (-Ifges(5,5) * t315 - Ifges(5,6) * t316) * t486 + (-Ifges(5,1) * t315 - Ifges(5,4) * t316) * t492 + (-Ifges(5,4) * t315 - Ifges(5,2) * t316) * t493 + (-t239 / 0.2e1 - t427 / 0.2e1) * t91 + (-t238 / 0.2e1 + t428 / 0.2e1) * t90 + (-Ifges(6,5) * t427 + Ifges(6,6) * t428 + Ifges(6,3) * t316) * t494 + (-Ifges(6,1) * t427 + Ifges(6,4) * t428 + Ifges(6,5) * t316) * t498 + (-Ifges(6,4) * t427 + Ifges(6,2) * t428 + Ifges(6,6) * t316) * t500 + (Ifges(7,1) * t129 + Ifges(7,4) * t130 + Ifges(7,5) * t316) * t504 + (Ifges(7,1) * t161 + Ifges(7,4) * t160 + Ifges(7,5) * t275) * t505 + (Ifges(7,4) * t129 + Ifges(7,2) * t130 + Ifges(7,6) * t316) * t506 + (Ifges(7,4) * t161 + Ifges(7,2) * t160 + Ifges(7,6) * t275) * t507 - t246 * t516 - t245 * t518 + t334 + t87 * t28 + t88 * t29 + t189 * t74 + t190 * t75; t280 + t361 * t74 + t357 * t75 - t377 * t28 + t326 * t29 + t296 * t259 - t295 * t258 + t105 + t526 * t84 + t525 * t85 + (-t69 + t422) * t379 + (-t204 - t380) * t524 + t380 * qJD(5) + (t15 * t525 + t16 * t526 + t2 * t326 - t3 * t377 - t379 * t73) * m(7) + (-t223 * t381 - t379 * t92 - t383) * m(6) + (t379 * t97 - t524 * t98 + t256) * m(5) + (t207 * t296 - t208 * t295 + t301) * m(4); (-t459 / 0.2e1 - t222 / 0.2e1 + (t517 - Ifges(5,1) / 0.2e1) * t379 - t520 - t549) * t524 + ((-m(6) * t382 - t357 * t134 - t361 * t135) * qJD(5) + m(6) * t522 - t357 * t74 + t361 * t75) * pkin(10) + t522 * mrSges(6,3) + t40 * t479 + (-Ifges(7,5) * t262 - Ifges(7,6) * t263) * t496 + (-Ifges(7,1) * t262 - Ifges(7,4) * t263) * t504 + (-Ifges(7,4) * t262 - Ifges(7,2) * t263) * t506 + (-t263 / 0.2e1 + t142 / 0.2e1) * t57 + (t548 - t447 / 0.2e1 - t448 / 0.2e1 - t445 / 0.2e1 - t446 / 0.2e1 + t442 / 0.2e1 - t56 / 0.2e1 + t466 / 0.2e1 - t455 / 0.2e1 - t89 / 0.2e1 - t454 / 0.2e1) * t379 + t409 + (Ifges(7,5) * t326 - Ifges(7,6) * t377 + t385) * t503 + t32 * (mrSges(7,1) * t377 + mrSges(7,2) * t326) + (Ifges(7,4) * t326 - Ifges(7,2) * t377) * t514 + (Ifges(7,1) * t326 - Ifges(7,4) * t377) * t515 - t377 * t518 + (-t262 / 0.2e1 + t144 / 0.2e1) * t58 + t474 - t475 + (-t15 * t526 + t16 * t525 - t2 * t377 - t3 * t326) * mrSges(7,3) + (-mrSges(7,1) * t525 + mrSges(7,2) * t526) * t73 + t422 * t98 + (-pkin(4) * t50 - t53 * t67 - t54 * t68 - t92 * t98) * m(6) + t351 * t11 + t284 * t29 + t283 * t28 - t97 * t204 + (-Ifges(7,5) * t144 - Ifges(7,6) * t142) * t497 + (-Ifges(7,1) * t144 - Ifges(7,4) * t142) * t505 + (-Ifges(7,4) * t144 - Ifges(7,2) * t142) * t507 + (pkin(5) * t357 * t527 + t520) * qJD(5) + t528 * t85 + t529 * t84 + (t15 * t528 + t16 * t529 + t2 * t284 + t283 * t3 + t32 * t351 - t73 * t76) * m(7) + t387 * t508 + t390 * t509 + t357 * t513 + t326 * t516 - pkin(4) * t59 - t76 * t69 - t68 * t134 - t67 * t135 - t50 * t394; (t197 * t53 + t198 * t54) * mrSges(6,3) + t122 * t536 + (Ifges(6,5) * t197 - Ifges(6,6) * t198) * t495 - m(7) * (t15 * t18 + t16 * t19) + t534 + (-Ifges(6,2) * t198 + t192 + t91) * t501 - t92 * (mrSges(6,1) * t198 + mrSges(6,2) * t197) + t535 + t90 * t498 + (Ifges(6,1) * t197 - t465) * t499 + (t360 * t28 + t356 * t29 + m(7) * (t2 * t356 + t3 * t360) - t527 * t198 + (-t356 * t85 + t360 * t84 + m(7) * (-t15 * t356 + t16 * t360)) * qJD(6)) * pkin(5) - t19 * t84 - t18 * t85 - t53 * t134 + t54 * t135 + t550; -t15 * t84 + t16 * t85 + t57 * t504 + t547 + t550 + t8;];
tauc  = t1(:);
