% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRR7
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
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:32:50
% EndTime: 2019-03-09 18:34:01
% DurationCPUTime: 35.83s
% Computational Cost: add. (32804->928), mult. (85323->1305), div. (0->0), fcn. (69326->12), ass. (0->411)
t355 = sin(pkin(12));
t357 = cos(pkin(12));
t361 = sin(qJ(3));
t365 = cos(qJ(3));
t327 = t355 * t365 + t357 * t361;
t366 = cos(qJ(2));
t356 = sin(pkin(6));
t435 = qJD(1) * t356;
t417 = t366 * t435;
t264 = t327 * t417;
t320 = t327 * qJD(3);
t437 = -t264 + t320;
t326 = -t355 * t361 + t357 * t365;
t321 = t326 * qJD(3);
t360 = sin(qJ(5));
t364 = cos(qJ(5));
t381 = t364 * t326 - t327 * t360;
t200 = qJD(5) * t381 - t320 * t360 + t321 * t364;
t265 = t326 * t417;
t217 = -t264 * t360 + t265 * t364;
t440 = t200 - t217;
t258 = t326 * t360 + t327 * t364;
t201 = qJD(5) * t258 + t364 * t320 + t321 * t360;
t216 = t364 * t264 + t265 * t360;
t439 = t201 - t216;
t362 = sin(qJ(2));
t358 = cos(pkin(6));
t434 = qJD(1) * t358;
t422 = pkin(1) * t434;
t312 = pkin(8) * t417 + t362 * t422;
t405 = t361 * t417;
t262 = pkin(3) * t405 + t312;
t431 = qJD(3) * t361;
t558 = pkin(3) * t431 - t262;
t418 = t362 * t435;
t309 = -pkin(8) * t418 + t366 * t422;
t377 = (pkin(2) * t362 - pkin(9) * t366) * t356;
t310 = qJD(1) * t377;
t242 = -t361 * t309 + t365 * t310;
t214 = (-qJ(4) * t365 * t366 + pkin(3) * t362) * t435 + t242;
t243 = t365 * t309 + t361 * t310;
t227 = -qJ(4) * t405 + t243;
t148 = t357 * t214 - t227 * t355;
t132 = pkin(4) * t418 - pkin(10) * t265 + t148;
t149 = t355 * t214 + t357 * t227;
t134 = -pkin(10) * t264 + t149;
t475 = -qJ(4) - pkin(9);
t409 = qJD(3) * t475;
t317 = qJD(4) * t365 + t361 * t409;
t318 = -qJD(4) * t361 + t365 * t409;
t245 = t357 * t317 + t355 * t318;
t218 = -pkin(10) * t320 + t245;
t244 = -t317 * t355 + t357 * t318;
t379 = -pkin(10) * t321 + t244;
t338 = t475 * t361;
t339 = t475 * t365;
t269 = t357 * t338 + t339 * t355;
t246 = -pkin(10) * t327 + t269;
t270 = t355 * t338 - t357 * t339;
t247 = pkin(10) * t326 + t270;
t383 = t364 * t246 - t247 * t360;
t552 = qJD(5) * t383 + (-t134 + t218) * t364 + (-t132 + t379) * t360;
t549 = pkin(4) * t437 + t558;
t557 = -pkin(11) * t418 + t552;
t556 = t439 * pkin(5) - pkin(11) * t440 + t549;
t555 = Ifges(4,3) + Ifges(5,3);
t353 = -pkin(3) * t365 - pkin(2);
t292 = -pkin(4) * t326 + t353;
t174 = -pkin(5) * t381 - pkin(11) * t258 + t292;
t178 = t246 * t360 + t247 * t364;
t359 = sin(qJ(6));
t363 = cos(qJ(6));
t114 = t174 * t359 + t178 * t363;
t554 = -qJD(6) * t114 - t359 * t557 + t556 * t363;
t113 = t174 * t363 - t178 * t359;
t553 = qJD(6) * t113 + t556 * t359 + t363 * t557;
t551 = t244 - t148;
t550 = t245 - t149;
t347 = qJD(2) + t434;
t284 = t347 * t365 - t361 * t418;
t285 = t347 * t361 + t365 * t418;
t233 = t284 * t355 + t285 * t357;
t406 = t357 * t284 - t285 * t355;
t545 = -t233 * t360 + t364 * t406;
t158 = qJD(6) - t545;
t544 = t233 * t364 + t360 * t406;
t107 = pkin(5) * t544 - pkin(11) * t545;
t548 = pkin(10) * t233;
t547 = t233 * Ifges(5,4);
t444 = t356 * t366;
t325 = t358 * t362 * pkin(1) + pkin(8) * t444;
t272 = pkin(9) * t347 + t312;
t303 = (-pkin(2) * t366 - pkin(9) * t362 - pkin(1)) * t356;
t279 = qJD(1) * t303;
t225 = -t272 * t361 + t365 * t279;
t191 = -qJ(4) * t285 + t225;
t337 = qJD(3) - t417;
t173 = pkin(3) * t337 + t191;
t226 = t272 * t365 + t279 * t361;
t192 = qJ(4) * t284 + t226;
t443 = t357 * t192;
t121 = t355 * t173 + t443;
t541 = pkin(10) * t406;
t101 = t121 + t541;
t428 = qJD(5) * t364;
t429 = qJD(5) * t360;
t432 = qJD(2) * t366;
t415 = t365 * t432;
t430 = qJD(3) * t365;
t253 = t347 * t430 + (-t362 * t431 + t415) * t435;
t416 = t361 * t432;
t254 = -t347 * t431 + (-t362 * t430 - t416) * t435;
t194 = t253 * t357 + t254 * t355;
t425 = qJD(1) * qJD(2);
t411 = t356 * t425;
t404 = t362 * t411;
t311 = qJD(2) * t377;
t293 = qJD(1) * t311;
t445 = t356 * t362;
t348 = pkin(8) * t445;
t485 = pkin(1) * t366;
t324 = t358 * t485 - t348;
t313 = t324 * qJD(2);
t294 = qJD(1) * t313;
t164 = -qJD(3) * t226 + t365 * t293 - t294 * t361;
t122 = pkin(3) * t404 - qJ(4) * t253 - qJD(4) * t285 + t164;
t163 = -t272 * t431 + t279 * t430 + t361 * t293 + t365 * t294;
t129 = qJ(4) * t254 + qJD(4) * t284 + t163;
t71 = t357 * t122 - t129 * t355;
t49 = pkin(4) * t404 - pkin(10) * t194 + t71;
t193 = -t253 * t355 + t254 * t357;
t72 = t355 * t122 + t357 * t129;
t51 = pkin(10) * t193 + t72;
t183 = t355 * t192;
t120 = t357 * t173 - t183;
t95 = pkin(4) * t337 + t120 - t548;
t10 = -t101 * t429 + t360 * t49 + t364 * t51 + t95 * t428;
t295 = t325 * t425;
t222 = -t254 * pkin(3) + t295;
t144 = -t193 * pkin(4) + t222;
t331 = qJD(5) + t337;
t147 = t331 * t359 + t363 * t544;
t90 = qJD(5) * t545 + t193 * t360 + t194 * t364;
t57 = -qJD(6) * t147 - t359 * t90 + t363 * t404;
t54 = Ifges(7,6) * t57;
t146 = t331 * t363 - t359 * t544;
t56 = qJD(6) * t146 + t359 * t404 + t363 * t90;
t55 = Ifges(7,5) * t56;
t91 = qJD(5) * t544 - t364 * t193 + t194 * t360;
t16 = Ifges(7,3) * t91 + t54 + t55;
t45 = t101 * t364 + t360 * t95;
t43 = pkin(11) * t331 + t45;
t271 = -t347 * pkin(2) - t309;
t229 = -t284 * pkin(3) + qJD(4) + t271;
t166 = -pkin(4) * t406 + t229;
t77 = -pkin(5) * t545 - pkin(11) * t544 + t166;
t22 = -t359 * t43 + t363 * t77;
t30 = t91 * pkin(5) - t90 * pkin(11) + t144;
t7 = pkin(11) * t404 + t10;
t2 = qJD(6) * t22 + t30 * t359 + t363 * t7;
t23 = t359 * t77 + t363 * t43;
t3 = -qJD(6) * t23 + t30 * t363 - t359 * t7;
t402 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t479 = t90 * Ifges(6,4);
t546 = t402 + t144 * mrSges(6,1) - t10 * mrSges(6,3) + t16 / 0.2e1 - t479 / 0.2e1;
t154 = Ifges(5,2) * t406 + t337 * Ifges(5,6) + t547;
t543 = t154 / 0.2e1;
t506 = -t406 / 0.2e1;
t540 = t91 * Ifges(6,2);
t539 = Ifges(5,4) * t406;
t392 = t22 * t363 + t23 * t359;
t538 = t392 * mrSges(7,3);
t152 = mrSges(6,1) * t331 - mrSges(6,3) * t544;
t92 = -mrSges(7,1) * t146 + mrSges(7,2) * t147;
t447 = t152 - t92;
t128 = t357 * t191 - t183;
t106 = t128 - t548;
t352 = pkin(3) * t357 + pkin(4);
t483 = pkin(3) * t355;
t316 = t360 * t352 + t364 * t483;
t127 = -t191 * t355 - t443;
t380 = t127 - t541;
t537 = t316 * qJD(5) - t106 * t360 + t364 * t380;
t302 = pkin(9) * t358 + t325;
t238 = -t361 * t302 + t365 * t303;
t323 = t358 * t361 + t365 * t445;
t204 = -pkin(3) * t444 - t323 * qJ(4) + t238;
t239 = t365 * t302 + t361 * t303;
t322 = t358 * t365 - t361 * t445;
t213 = qJ(4) * t322 + t239;
t139 = t357 * t204 - t355 * t213;
t249 = t322 * t355 + t323 * t357;
t112 = -pkin(4) * t444 - t249 * pkin(10) + t139;
t140 = t355 * t204 + t357 * t213;
t248 = t322 * t357 - t323 * t355;
t116 = pkin(10) * t248 + t140;
t534 = t360 * t112 + t364 * t116;
t157 = Ifges(6,4) * t545;
t451 = t331 * Ifges(6,5);
t461 = t544 * Ifges(6,1);
t100 = t157 + t451 + t461;
t394 = Ifges(7,5) * t363 - Ifges(7,6) * t359;
t372 = t158 * t394;
t469 = Ifges(7,4) * t359;
t398 = Ifges(7,1) * t363 - t469;
t373 = t147 * t398;
t468 = Ifges(7,4) * t363;
t396 = -Ifges(7,2) * t359 + t468;
t374 = t146 * t396;
t399 = mrSges(7,1) * t359 + mrSges(7,2) * t363;
t44 = -t101 * t360 + t364 * t95;
t42 = -pkin(5) * t331 - t44;
t375 = t42 * t399;
t488 = -t363 / 0.2e1;
t489 = t359 / 0.2e1;
t470 = Ifges(7,4) * t147;
t67 = Ifges(7,2) * t146 + Ifges(7,6) * t158 + t470;
t145 = Ifges(7,4) * t146;
t68 = Ifges(7,1) * t147 + Ifges(7,5) * t158 + t145;
t533 = t44 * mrSges(6,3) + t538 - t100 / 0.2e1 - t374 / 0.2e1 - t373 / 0.2e1 - t372 / 0.2e1 - t166 * mrSges(6,2) - t451 / 0.2e1 + t67 * t489 + t68 * t488 - t375 - t157 / 0.2e1;
t532 = Ifges(4,5) * t253 + Ifges(5,5) * t194 + Ifges(4,6) * t254 + Ifges(5,6) * t193 + t404 * t555;
t531 = -t22 * t359 + t23 * t363;
t454 = t285 * Ifges(4,4);
t220 = t284 * Ifges(4,2) + t337 * Ifges(4,6) + t454;
t280 = Ifges(4,4) * t284;
t221 = t285 * Ifges(4,1) + t337 * Ifges(4,5) + t280;
t384 = t225 * t365 + t226 * t361;
t472 = Ifges(4,4) * t365;
t473 = Ifges(4,4) * t361;
t486 = t365 / 0.2e1;
t491 = t337 / 0.2e1;
t496 = t285 / 0.2e1;
t497 = t284 / 0.2e1;
t530 = -t384 * mrSges(4,3) + t271 * (mrSges(4,1) * t361 + mrSges(4,2) * t365) + (-Ifges(4,2) * t361 + t472) * t497 + (Ifges(4,1) * t365 - t473) * t496 + (Ifges(4,5) * t365 - Ifges(4,6) * t361) * t491 - t361 * t220 / 0.2e1 + t221 * t486;
t11 = -qJD(5) * t45 - t360 * t51 + t364 * t49;
t260 = -qJD(3) * t323 - t356 * t416;
t261 = qJD(3) * t322 + t356 * t415;
t210 = t260 * t355 + t261 * t357;
t433 = qJD(2) * t356;
t414 = t362 * t433;
t176 = -qJD(3) * t239 + t365 * t311 - t313 * t361;
t137 = pkin(3) * t414 - qJ(4) * t261 - qJD(4) * t323 + t176;
t175 = -t302 * t431 + t303 * t430 + t361 * t311 + t365 * t313;
t143 = qJ(4) * t260 + qJD(4) * t322 + t175;
t81 = t357 * t137 - t143 * t355;
t65 = pkin(4) * t414 - pkin(10) * t210 + t81;
t209 = t260 * t357 - t261 * t355;
t82 = t355 * t137 + t357 * t143;
t70 = pkin(10) * t209 + t82;
t20 = -qJD(5) * t534 - t360 * t70 + t364 * t65;
t464 = t158 * Ifges(7,3);
t465 = t147 * Ifges(7,5);
t466 = t146 * Ifges(7,6);
t66 = t464 + t465 + t466;
t450 = t331 * Ifges(6,6);
t463 = t545 * Ifges(6,2);
t471 = Ifges(6,4) * t544;
t99 = t450 + t463 + t471;
t528 = t23 * mrSges(7,2) + t45 * mrSges(6,3) - t66 / 0.2e1 + t99 / 0.2e1 - t166 * mrSges(6,1) - t22 * mrSges(7,1);
t527 = Ifges(6,2) / 0.2e1;
t525 = t56 / 0.2e1;
t524 = t57 / 0.2e1;
t523 = t68 / 0.2e1;
t522 = t91 / 0.2e1;
t521 = pkin(1) * mrSges(3,1);
t520 = pkin(1) * mrSges(3,2);
t519 = -t146 / 0.2e1;
t518 = t146 / 0.2e1;
t517 = -t147 / 0.2e1;
t516 = t147 / 0.2e1;
t514 = -t158 / 0.2e1;
t513 = t158 / 0.2e1;
t512 = t545 / 0.2e1;
t511 = t544 / 0.2e1;
t382 = t364 * t248 - t249 * t360;
t510 = t382 / 0.2e1;
t180 = t248 * t360 + t249 * t364;
t509 = t180 / 0.2e1;
t508 = t193 / 0.2e1;
t507 = t194 / 0.2e1;
t505 = t406 / 0.2e1;
t504 = -t233 / 0.2e1;
t503 = t233 / 0.2e1;
t502 = t248 / 0.2e1;
t501 = t249 / 0.2e1;
t500 = t253 / 0.2e1;
t499 = t254 / 0.2e1;
t495 = t322 / 0.2e1;
t494 = t323 / 0.2e1;
t493 = t331 / 0.2e1;
t492 = -t337 / 0.2e1;
t490 = -t359 / 0.2e1;
t487 = t363 / 0.2e1;
t484 = pkin(3) * t285;
t482 = t2 * t363;
t481 = t3 * t359;
t480 = t90 * Ifges(6,1);
t478 = t91 * Ifges(6,4);
t474 = Ifges(3,4) * t362;
t467 = Ifges(3,5) * t366;
t462 = t545 * Ifges(6,6);
t460 = t544 * Ifges(6,5);
t457 = t406 * Ifges(5,6);
t456 = t233 * Ifges(5,5);
t455 = t284 * Ifges(4,6);
t453 = t285 * Ifges(4,5);
t452 = t294 * mrSges(3,2);
t449 = t331 * Ifges(6,3);
t448 = t347 * Ifges(3,5);
t442 = t359 * t200;
t441 = t363 * t200;
t438 = -mrSges(3,1) * t347 - mrSges(4,1) * t284 + mrSges(4,2) * t285 + mrSges(3,3) * t418;
t436 = t265 - t321;
t314 = t325 * qJD(2);
t427 = qJD(6) * t359;
t426 = qJD(6) * t363;
t423 = Ifges(6,5) * t90 - Ifges(6,6) * t91 + Ifges(6,3) * t404;
t421 = -Ifges(5,3) / 0.2e1 - Ifges(4,3) / 0.2e1;
t38 = t91 * mrSges(6,1) + t90 * mrSges(6,2);
t131 = -t193 * mrSges(5,1) + t194 * mrSges(5,2);
t197 = -t217 * t359 + t363 * t418;
t408 = t197 + t442;
t198 = t217 * t363 + t359 * t418;
t407 = -t198 + t441;
t236 = -pkin(3) * t260 + t314;
t196 = pkin(4) * t233 + t484;
t401 = -t2 * t359 - t3 * t363;
t400 = mrSges(7,1) * t363 - mrSges(7,2) * t359;
t397 = Ifges(7,1) * t359 + t468;
t395 = Ifges(7,2) * t363 + t469;
t393 = Ifges(7,5) * t359 + Ifges(7,6) * t363;
t28 = mrSges(7,1) * t91 - mrSges(7,3) * t56;
t29 = -mrSges(7,2) * t91 + mrSges(7,3) * t57;
t390 = -t359 * t28 + t363 * t29;
t59 = -pkin(11) * t444 + t534;
t301 = t348 + (-pkin(2) - t485) * t358;
t252 = -t322 * pkin(3) + t301;
t195 = -t248 * pkin(4) + t252;
t93 = -pkin(5) * t382 - t180 * pkin(11) + t195;
t32 = t359 * t93 + t363 * t59;
t31 = -t359 * t59 + t363 * t93;
t96 = -mrSges(7,2) * t158 + mrSges(7,3) * t146;
t97 = mrSges(7,1) * t158 - mrSges(7,3) * t147;
t389 = -t359 * t97 + t363 * t96;
t60 = t364 * t112 - t360 * t116;
t78 = t132 * t364 - t134 * t360;
t385 = t163 * t365 - t164 * t361;
t315 = t352 * t364 - t360 * t483;
t151 = -mrSges(6,2) * t331 + mrSges(6,3) * t545;
t378 = -t151 - t389;
t167 = -t359 * t180 - t363 * t444;
t376 = -t363 * t180 + t359 * t444;
t19 = t112 * t428 - t116 * t429 + t360 * t65 + t364 * t70;
t156 = -pkin(4) * t209 + t236;
t371 = -qJD(6) * t392 - t481;
t17 = t56 * Ifges(7,4) + t57 * Ifges(7,2) + t91 * Ifges(7,6);
t18 = t56 * Ifges(7,1) + t57 * Ifges(7,4) + t91 * Ifges(7,5);
t8 = -pkin(5) * t404 - t11;
t370 = t11 * mrSges(6,1) - t10 * mrSges(6,2) + mrSges(7,3) * t482 + qJD(6) * t375 + t17 * t487 + t18 * t489 + t395 * t524 + t397 * t525 - t8 * t400 + t423 - t67 * t427 / 0.2e1 + t426 * t523 + t393 * t522 + (t374 + t373 + t372) * qJD(6) / 0.2e1;
t369 = -t464 / 0.2e1 - t465 / 0.2e1 - t466 / 0.2e1 + t450 / 0.2e1 + t471 / 0.2e1 + t528;
t341 = Ifges(3,4) * t417;
t336 = t411 * t467;
t308 = pkin(11) + t316;
t307 = -pkin(5) - t315;
t306 = -t347 * mrSges(3,2) + mrSges(3,3) * t417;
t299 = t315 * qJD(5);
t267 = Ifges(3,1) * t418 + t341 + t448;
t266 = Ifges(3,6) * t347 + (Ifges(3,2) * t366 + t474) * t435;
t256 = mrSges(4,1) * t337 - mrSges(4,3) * t285;
t255 = -mrSges(4,2) * t337 + mrSges(4,3) * t284;
t235 = -mrSges(4,2) * t404 + mrSges(4,3) * t254;
t234 = mrSges(4,1) * t404 - mrSges(4,3) * t253;
t219 = t337 * Ifges(4,3) + t453 + t455;
t212 = mrSges(5,1) * t337 - mrSges(5,3) * t233;
t211 = -mrSges(5,2) * t337 + mrSges(5,3) * t406;
t199 = -mrSges(4,1) * t254 + mrSges(4,2) * t253;
t182 = t253 * Ifges(4,1) + t254 * Ifges(4,4) + Ifges(4,5) * t404;
t181 = t253 * Ifges(4,4) + t254 * Ifges(4,2) + Ifges(4,6) * t404;
t172 = mrSges(5,1) * t404 - mrSges(5,3) * t194;
t171 = -mrSges(5,2) * t404 + mrSges(5,3) * t193;
t165 = -mrSges(5,1) * t406 + mrSges(5,2) * t233;
t155 = t233 * Ifges(5,1) + t337 * Ifges(5,5) + t539;
t153 = t337 * Ifges(5,3) + t456 + t457;
t125 = t194 * Ifges(5,1) + t193 * Ifges(5,4) + Ifges(5,5) * t404;
t124 = t194 * Ifges(5,4) + t193 * Ifges(5,2) + Ifges(5,6) * t404;
t109 = qJD(5) * t178 + t218 * t360 - t364 * t379;
t105 = -mrSges(6,1) * t545 + mrSges(6,2) * t544;
t103 = qJD(5) * t180 - t364 * t209 + t210 * t360;
t102 = qJD(5) * t382 + t209 * t360 + t210 * t364;
t98 = t449 + t460 + t462;
t84 = -mrSges(6,2) * t404 - mrSges(6,3) * t91;
t83 = mrSges(6,1) * t404 - mrSges(6,3) * t90;
t80 = t107 + t196;
t75 = -pkin(5) * t418 - t78;
t74 = qJD(6) * t376 - t359 * t102 + t363 * t414;
t73 = qJD(6) * t167 + t363 * t102 + t359 * t414;
t58 = pkin(5) * t444 - t60;
t47 = t364 * t106 + t360 * t380;
t37 = Ifges(6,5) * t404 - t478 + t480;
t36 = Ifges(6,6) * t404 + t479 - t540;
t35 = pkin(5) * t103 - pkin(11) * t102 + t156;
t27 = t107 * t359 + t363 * t44;
t26 = t107 * t363 - t359 * t44;
t25 = t359 * t80 + t363 * t47;
t24 = -t359 * t47 + t363 * t80;
t21 = -mrSges(7,1) * t57 + mrSges(7,2) * t56;
t13 = -pkin(5) * t414 - t20;
t12 = pkin(11) * t414 + t19;
t5 = -qJD(6) * t32 - t12 * t359 + t35 * t363;
t4 = qJD(6) * t31 + t12 * t363 + t35 * t359;
t1 = [(Ifges(7,5) * t73 + Ifges(7,6) * t74) * t513 + m(7) * (t13 * t42 + t2 * t32 + t22 * t5 + t23 * t4 + t3 * t31 + t58 * t8) + m(5) * (t120 * t81 + t121 * t82 + t139 * t71 + t140 * t72 + t222 * t252 + t229 * t236) + m(4) * (t163 * t239 + t164 * t238 + t175 * t226 + t176 * t225 + t271 * t314 + t295 * t301) + m(3) * (t294 * t325 - t295 * t324 - t309 * t314 + t312 * t313) + (-Ifges(6,4) * t511 + Ifges(7,5) * t516 - Ifges(6,2) * t512 - Ifges(6,6) * t493 + Ifges(7,6) * t518 + Ifges(7,3) * t513 - t528) * t103 + (-t452 - t295 * mrSges(3,1) + t336 / 0.2e1) * t358 + t209 * t543 - t91 * (Ifges(6,4) * t180 - Ifges(6,6) * t444) / 0.2e1 + (Ifges(6,4) * t102 + Ifges(6,6) * t414) * t512 + t90 * (Ifges(6,1) * t180 - Ifges(6,5) * t444) / 0.2e1 + (Ifges(6,1) * t102 + Ifges(6,5) * t414) * t511 + t438 * t314 + (t10 * t444 + t102 * t166 + t144 * t180 - t414 * t45) * mrSges(6,2) + t11 * (-mrSges(6,1) * t444 - t180 * mrSges(6,3)) + (Ifges(7,4) * t73 + Ifges(7,2) * t74) * t518 + (Ifges(4,5) * t261 + Ifges(5,5) * t210 + Ifges(4,6) * t260 + Ifges(5,6) * t209 + t414 * t555) * t491 + (Ifges(6,5) * t102 + Ifges(6,3) * t414) * t493 - t266 * t414 / 0.2e1 + t121 * (-mrSges(5,2) * t414 + mrSges(5,3) * t209) + t226 * (-mrSges(4,2) * t414 + mrSges(4,3) * t260) + t44 * (mrSges(6,1) * t414 - mrSges(6,3) * t102) - (t423 + t532) * t444 / 0.2e1 + (Ifges(7,1) * t73 + Ifges(7,4) * t74) * t516 + t120 * (mrSges(5,1) * t414 - mrSges(5,3) * t210) + t225 * (mrSges(4,1) * t414 - mrSges(4,3) * t261) + ((t219 + t153 + t98) * t362 + t366 * t267 + t347 * (-Ifges(3,6) * t362 + t467)) * t433 / 0.2e1 + m(6) * (t10 * t534 + t11 * t60 + t144 * t195 + t156 * t166 + t19 * t45 + t20 * t44) + t534 * t84 + (-Ifges(7,5) * t376 + Ifges(7,6) * t167) * t522 + (-Ifges(7,4) * t376 + Ifges(7,2) * t167) * t524 + (t167 * t2 - t22 * t73 + t23 * t74 + t3 * t376) * mrSges(7,3) + (-Ifges(7,1) * t376 + Ifges(7,4) * t167) * t525 + t8 * (-mrSges(7,1) * t167 - mrSges(7,2) * t376) - t376 * t18 / 0.2e1 + t295 * (-mrSges(4,1) * t322 + mrSges(4,2) * t323) + t313 * t306 + t301 * t199 + t271 * (-mrSges(4,1) * t260 + mrSges(4,2) * t261) + t260 * t220 / 0.2e1 + t261 * t221 / 0.2e1 + t175 * t255 + t176 * t256 + t252 * t131 + t222 * (-mrSges(5,1) * t248 + mrSges(5,2) * t249) + t236 * t165 + t238 * t234 + t239 * t235 + t229 * (-mrSges(5,1) * t209 + mrSges(5,2) * t210) + t82 * t211 + t81 * t212 + t210 * t155 / 0.2e1 + t195 * t38 + t140 * t171 + t139 * t172 + t167 * t17 / 0.2e1 + t156 * t105 + t19 * t151 + t20 * t152 + t102 * t100 / 0.2e1 + t4 * t96 + t5 * t97 + t13 * t92 + ((-t324 * mrSges(3,3) + Ifges(3,5) * t358 / 0.2e1 + (-0.2e1 * t520 + 0.3e1 / 0.2e1 * Ifges(3,4) * t366) * t356) * t366 + (-t325 * mrSges(3,3) - Ifges(3,6) * t358 + Ifges(5,5) * t501 + Ifges(5,6) * t502 + Ifges(4,5) * t494 + Ifges(4,6) * t495 + Ifges(6,5) * t509 + Ifges(6,6) * t510 + (-0.2e1 * t521 - 0.3e1 / 0.2e1 * t474) * t356 + (-Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + t421) * t444) * t362) * t411 + t60 * t83 + t42 * (-mrSges(7,1) * t74 + mrSges(7,2) * t73) + t74 * t67 / 0.2e1 - (t540 / 0.2e1 + Ifges(7,6) * t524 + Ifges(7,5) * t525 + Ifges(7,3) * t522 + t546) * t382 + t58 * t21 + t31 * t28 + t32 * t29 + t73 * t523 + (Ifges(5,4) * t249 + Ifges(5,2) * t248 - Ifges(5,6) * t444) * t508 + t37 * t509 + t36 * t510 + t182 * t494 + t181 * t495 + (Ifges(4,1) * t261 + Ifges(4,4) * t260 + Ifges(4,5) * t414) * t496 + (Ifges(4,4) * t261 + Ifges(4,2) * t260 + Ifges(4,6) * t414) * t497 + (Ifges(4,4) * t323 + Ifges(4,2) * t322 - Ifges(4,6) * t444) * t499 + (Ifges(4,1) * t323 + Ifges(4,4) * t322 - Ifges(4,5) * t444) * t500 + t125 * t501 + t124 * t502 + (Ifges(5,1) * t210 + Ifges(5,4) * t209 + Ifges(5,5) * t414) * t503 + (Ifges(5,4) * t210 + Ifges(5,2) * t209 + Ifges(5,6) * t414) * t505 + (Ifges(5,1) * t249 + Ifges(5,4) * t248 - Ifges(5,5) * t444) * t507 + t71 * (-mrSges(5,1) * t444 - t249 * mrSges(5,3)) + t164 * (-mrSges(4,1) * t444 - t323 * mrSges(4,3)) + t72 * (mrSges(5,2) * t444 + t248 * mrSges(5,3)) + t163 * (mrSges(4,2) * t444 + t322 * mrSges(4,3)) + (t294 * t366 + t295 * t362 + (-t309 * t366 - t312 * t362) * qJD(2)) * t356 * mrSges(3,3); (t361 * pkin(3) * t165 + (-t255 * t361 - t256 * t365) * pkin(9) + t530) * qJD(3) - t331 * (Ifges(6,5) * t217 - Ifges(6,6) * t216) / 0.2e1 - t545 * (Ifges(6,4) * t217 - Ifges(6,2) * t216) / 0.2e1 - t544 * (Ifges(6,1) * t217 - Ifges(6,4) * t216) / 0.2e1 + (t8 * t399 + t394 * t522 + t396 * t524 + t398 * t525 + t37 / 0.2e1 + t144 * mrSges(6,2) - t478 / 0.2e1 + t480 / 0.2e1 + t18 * t487 + t17 * t490 - t11 * mrSges(6,3) + t401 * mrSges(7,3) + (-mrSges(7,3) * t531 + t393 * t514 + t395 * t519 + t397 * t517 + t400 * t42 + t488 * t67 + t490 * t68) * qJD(6)) * t258 + t336 + (mrSges(5,1) * t437 - mrSges(5,2) * t436) * t229 + (t120 * t436 - t121 * t437 + t326 * t72 - t327 * t71) * mrSges(5,3) - t438 * t312 + (mrSges(7,1) * t439 - mrSges(7,3) * t407) * t22 + (-mrSges(7,2) * t439 - mrSges(7,3) * t408) * t23 + (mrSges(6,1) * t439 + mrSges(6,2) * t440) * t166 + (-t439 * t45 - t44 * t440) * mrSges(6,3) + (t551 * t120 + t550 * t121 + t222 * t353 + t229 * t558 + t269 * t71 + t270 * t72) * m(5) + (t99 - t66) * (t216 / 0.2e1 - t201 / 0.2e1) + (-t198 / 0.2e1 + t441 / 0.2e1) * t68 + (-t197 / 0.2e1 - t442 / 0.2e1) * t67 + (mrSges(7,1) * t408 + mrSges(7,2) * t407) * t42 + (Ifges(5,1) * t265 - Ifges(5,4) * t264) * t504 + (Ifges(5,4) * t265 - Ifges(5,2) * t264) * t506 + (Ifges(5,5) * t265 - Ifges(5,6) * t264) * t492 + (-t320 / 0.2e1 + t264 / 0.2e1) * t154 + (Ifges(5,1) * t321 - Ifges(5,4) * t320) * t503 + (Ifges(5,4) * t321 - Ifges(5,2) * t320) * t505 + (Ifges(5,5) * t321 - Ifges(5,6) * t320) * t491 + (-t234 * t361 + t235 * t365) * pkin(9) + (-t217 / 0.2e1 + t200 / 0.2e1) * t100 - t452 + (-pkin(2) * t295 - t225 * t242 - t226 * t243 - t271 * t312 + (-qJD(3) * t384 + t385) * pkin(9)) * m(4) + t385 * mrSges(4,3) - t447 * t109 + (-mrSges(4,1) * t365 + mrSges(4,2) * t361 - mrSges(3,1)) * t295 + t553 * t96 + (t113 * t3 + t114 * t2 - t383 * t8 + (t109 - t75) * t42 + t553 * t23 + t554 * t22) * m(7) + t554 * t97 + (t10 * t178 + t11 * t383 + t144 * t292 + t552 * t45 + (-t109 - t78) * t44 + t549 * t166) * m(6) + t552 * t151 + t549 * t105 + t550 * t211 + t551 * t212 + ((t309 * mrSges(3,3) + t435 * t520 - t341 / 0.2e1 - t267 / 0.2e1 - t448 / 0.2e1 - t530) * t366 + (-t219 / 0.2e1 + t266 / 0.2e1 - t98 / 0.2e1 - t456 / 0.2e1 - t457 / 0.2e1 + t121 * mrSges(5,2) - t120 * mrSges(5,1) - t460 / 0.2e1 - t462 / 0.2e1 - t44 * mrSges(6,1) + t45 * mrSges(6,2) - t449 / 0.2e1 + t312 * mrSges(3,3) - t453 / 0.2e1 - t225 * mrSges(4,1) + t226 * mrSges(4,2) - t455 / 0.2e1 - t153 / 0.2e1 + t421 * t337 + (t521 + t474 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t366) * t435 + (-qJD(2) + t347 / 0.2e1) * Ifges(3,6) + (Ifges(4,5) * t361 + Ifges(5,5) * t327 + Ifges(6,5) * t258 + Ifges(4,6) * t365 + Ifges(5,6) * t326 + Ifges(6,6) * t381) * qJD(2) / 0.2e1) * t362) * t435 - (t21 - t83) * t383 + t361 * t182 / 0.2e1 + t353 * t131 + t326 * t124 / 0.2e1 + t222 * (-mrSges(5,1) * t326 + mrSges(5,2) * t327) + t327 * t125 / 0.2e1 - t309 * t306 + t292 * t38 + t269 * t172 + t270 * t171 - t262 * t165 - t243 * t255 - t242 * t256 + (t321 / 0.2e1 - t265 / 0.2e1) * t155 - pkin(2) * t199 + t178 * t84 - t78 * t152 + t113 * t28 + t114 * t29 - t75 * t92 - (t54 / 0.2e1 + t55 / 0.2e1 - t36 / 0.2e1 + (t527 + Ifges(7,3) / 0.2e1) * t91 + t546) * t381 + (Ifges(5,1) * t327 + Ifges(5,4) * t326) * t507 + (Ifges(5,4) * t327 + Ifges(5,2) * t326) * t508 + (Ifges(6,1) * t200 - Ifges(6,4) * t201) * t511 + (Ifges(6,4) * t200 - Ifges(6,2) * t201) * t512 + (Ifges(7,5) * t441 - Ifges(7,6) * t442 + Ifges(7,3) * t201) * t513 + (Ifges(7,5) * t198 + Ifges(7,6) * t197 + Ifges(7,3) * t216) * t514 + (Ifges(7,1) * t441 - Ifges(7,4) * t442 + Ifges(7,5) * t201) * t516 + (Ifges(7,1) * t198 + Ifges(7,4) * t197 + Ifges(7,5) * t216) * t517 + (Ifges(7,4) * t441 - Ifges(7,2) * t442 + Ifges(7,6) * t201) * t518 + (Ifges(7,4) * t198 + Ifges(7,2) * t197 + Ifges(7,6) * t216) * t519 + (Ifges(6,5) * t200 - Ifges(6,6) * t201) * t493 + (Ifges(4,2) * t365 + t473) * t499 + (Ifges(4,1) * t361 + t472) * t500 + t181 * t486; (t225 * t284 + t226 * t285) * mrSges(4,3) + t370 - t285 * (Ifges(4,1) * t284 - t454) / 0.2e1 + t532 + (t463 / 0.2e1 + t369) * t544 + (-t461 / 0.2e1 + t533) * t545 + t390 * t308 - (-Ifges(4,2) * t285 + t221 + t280) * t284 / 0.2e1 - t378 * t299 - mrSges(7,3) * t481 + (-Ifges(5,2) * t233 + t155 + t539) * t506 + t233 * t543 + (t120 * t406 + t121 * t233) * mrSges(5,3) + (Ifges(4,5) * t284 + Ifges(5,5) * t406 - Ifges(4,6) * t285 - Ifges(5,6) * t233) * t492 - t229 * (mrSges(5,1) * t233 + mrSges(5,2) * t406) + t316 * t84 + (-t165 * t285 + t171 * t355 + t172 * t357) * pkin(3) + t315 * t83 + t307 * t21 - t271 * (mrSges(4,1) * t285 + mrSges(4,2) * t284) - t225 * t255 + t226 * t256 - t128 * t211 - t127 * t212 - t196 * t105 - t163 * mrSges(4,2) + t164 * mrSges(4,1) - t47 * t151 - t25 * t96 - t24 * t97 + (Ifges(5,1) * t406 - t547) * t504 - t72 * mrSges(5,2) + t71 * mrSges(5,1) + t220 * t496 + (t307 * t8 - t22 * t24 - t23 * t25 + (t371 + t482) * t308 + t537 * t42 + t531 * t299) * m(7) + (t10 * t316 + t11 * t315 - t166 * t196 + (t299 - t47) * t45 - t537 * t44) * m(6) - t447 * t537 + ((-t359 * t96 - t363 * t97) * t308 - t538) * qJD(6) + ((t355 * t72 + t357 * t71) * pkin(3) - t120 * t127 - t121 * t128 - t229 * t484) * m(5); t389 * qJD(6) + t447 * t544 + t378 * t545 - t406 * t211 + t233 * t212 + t363 * t28 + t359 * t29 + t131 + t38 + (t158 * t531 - t544 * t42 - t401) * m(7) + (t44 * t544 - t45 * t545 + t144) * m(6) + (t120 * t233 - t121 * t406 + t222) * m(5); t370 + t371 * mrSges(7,3) + t447 * t45 + ((t527 - Ifges(6,1) / 0.2e1) * t544 + t533) * t545 + t369 * t544 + (-t97 * t426 - t96 * t427 + t390) * pkin(11) - t44 * t151 - t27 * t96 - t26 * t97 - pkin(5) * t21 + (-pkin(5) * t8 + (-t22 * t426 - t23 * t427 - t481 + t482) * pkin(11) - t22 * t26 - t23 * t27 - t42 * t45) * m(7); -t42 * (mrSges(7,1) * t147 + mrSges(7,2) * t146) + (Ifges(7,1) * t146 - t470) * t517 + t67 * t516 + (Ifges(7,5) * t146 - Ifges(7,6) * t147) * t514 - t22 * t96 + t23 * t97 + (t146 * t22 + t147 * t23) * mrSges(7,3) + t402 + t16 + (-Ifges(7,2) * t147 + t145 + t68) * t519;];
tauc  = t1(:);
