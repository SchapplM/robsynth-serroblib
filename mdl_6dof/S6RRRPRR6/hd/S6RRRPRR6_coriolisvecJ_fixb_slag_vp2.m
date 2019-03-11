% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRR6
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:25:09
% EndTime: 2019-03-09 18:25:59
% DurationCPUTime: 26.43s
% Computational Cost: add. (25452->809), mult. (62840->1135), div. (0->0), fcn. (46678->10), ass. (0->363)
t356 = sin(qJ(2));
t360 = cos(qJ(2));
t379 = pkin(2) * t356 - pkin(8) * t360;
t316 = t379 * qJD(1);
t359 = cos(qJ(3));
t355 = sin(qJ(3));
t409 = qJD(1) * t356;
t391 = t355 * t409;
t266 = pkin(7) * t391 + t359 * t316;
t416 = t359 * t360;
t367 = pkin(3) * t356 - qJ(4) * t416;
t449 = -qJ(4) - pkin(8);
t384 = qJD(3) * t449;
t569 = -qJD(1) * t367 - qJD(4) * t355 + t359 * t384 - t266;
t298 = t355 * t316;
t401 = qJD(4) * t359;
t417 = t356 * t359;
t418 = t355 * t360;
t568 = t298 + (-pkin(7) * t417 - qJ(4) * t418) * qJD(1) - t355 * t384 - t401;
t351 = sin(pkin(11));
t352 = cos(pkin(11));
t308 = t351 * t359 + t352 * t355;
t366 = t360 * t308;
t272 = qJD(1) * t366;
t294 = t308 * qJD(3);
t413 = -t272 + t294;
t369 = t351 * t355 - t352 * t359;
t365 = t369 * t360;
t273 = qJD(1) * t365;
t295 = t369 * qJD(3);
t412 = -t273 + t295;
t406 = qJD(2) * t359;
t314 = -t391 + t406;
t390 = t359 * t409;
t315 = qJD(2) * t355 + t390;
t250 = t314 * t351 + t315 * t352;
t354 = sin(qJ(5));
t358 = cos(qJ(5));
t382 = t352 * t314 - t315 * t351;
t177 = t250 * t358 + t354 * t382;
t353 = sin(qJ(6));
t357 = cos(qJ(6));
t522 = -t250 * t354 + t358 * t382;
t106 = t177 * t357 + t353 * t522;
t327 = -qJD(2) * pkin(2) + pkin(7) * t409;
t263 = -pkin(3) * t314 + qJD(4) + t327;
t190 = -pkin(4) * t382 + t263;
t116 = -pkin(5) * t522 + t190;
t408 = qJD(1) * t360;
t340 = qJD(3) - t408;
t329 = qJD(5) + t340;
t319 = -pkin(2) * t360 - t356 * pkin(8) - pkin(1);
t302 = t319 * qJD(1);
t347 = pkin(7) * t408;
t328 = qJD(2) * pkin(8) + t347;
t255 = t359 * t302 - t328 * t355;
t215 = -qJ(4) * t315 + t255;
t203 = pkin(3) * t340 + t215;
t256 = t302 * t355 + t328 * t359;
t216 = qJ(4) * t314 + t256;
t208 = t351 * t216;
t130 = t352 * t203 - t208;
t535 = pkin(9) * t250;
t109 = pkin(4) * t340 + t130 - t535;
t420 = t352 * t216;
t131 = t351 * t203 + t420;
t519 = pkin(9) * t382;
t110 = t131 + t519;
t51 = t358 * t109 - t110 * t354;
t553 = pkin(10) * t177;
t43 = t51 - t553;
t40 = pkin(5) * t329 + t43;
t52 = t109 * t354 + t110 * t358;
t546 = pkin(10) * t522;
t44 = t52 + t546;
t428 = t353 * t44;
t16 = t357 * t40 - t428;
t427 = t357 * t44;
t17 = t353 * t40 + t427;
t549 = -t177 * t353 + t357 * t522;
t398 = qJD(2) * qJD(3);
t403 = qJD(3) * t356;
t405 = qJD(2) * t360;
t264 = t359 * t398 + (-t355 * t403 + t359 * t405) * qJD(1);
t402 = qJD(3) * t359;
t525 = t355 * t405 + t356 * t402;
t265 = -qJD(1) * t525 - t355 * t398;
t200 = -t264 * t351 + t265 * t352;
t201 = t264 * t352 + t265 * t351;
t88 = qJD(5) * t522 + t200 * t354 + t201 * t358;
t89 = -qJD(5) * t177 + t200 * t358 - t201 * t354;
t31 = qJD(6) * t549 + t353 * t89 + t357 * t88;
t32 = -qJD(6) * t106 - t353 * t88 + t357 * t89;
t407 = qJD(2) * t356;
t385 = qJD(1) * t407;
t396 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t385;
t444 = Ifges(7,4) * t106;
t320 = qJD(6) + t329;
t466 = -t320 / 0.2e1;
t48 = Ifges(7,2) * t549 + Ifges(7,6) * t320 + t444;
t488 = t106 / 0.2e1;
t489 = -t106 / 0.2e1;
t97 = Ifges(7,4) * t549;
t49 = Ifges(7,1) * t106 + Ifges(7,5) * t320 + t97;
t491 = -t549 / 0.2e1;
t399 = qJD(5) * t358;
t400 = qJD(5) * t354;
t317 = t379 * qJD(2);
t303 = qJD(1) * t317;
t381 = pkin(7) * t385;
t189 = -qJD(3) * t256 + t359 * t303 + t355 * t381;
t129 = pkin(3) * t385 - qJ(4) * t264 - qJD(4) * t315 + t189;
t404 = qJD(3) * t355;
t188 = t302 * t402 + t355 * t303 - t328 * t404 - t359 * t381;
t137 = qJ(4) * t265 + qJD(4) * t314 + t188;
t71 = t352 * t129 - t137 * t351;
t58 = pkin(4) * t385 - pkin(9) * t201 + t71;
t72 = t351 * t129 + t352 * t137;
t60 = pkin(9) * t200 + t72;
t12 = t109 * t399 - t110 * t400 + t354 * t58 + t358 * t60;
t10 = pkin(10) * t89 + t12;
t13 = -qJD(5) * t52 - t354 * t60 + t358 * t58;
t9 = pkin(5) * t385 - pkin(10) * t88 + t13;
t2 = qJD(6) * t16 + t10 * t357 + t353 * t9;
t3 = -qJD(6) * t17 - t10 * t353 + t357 * t9;
t518 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t567 = (Ifges(7,1) * t549 - t444) * t489 + (Ifges(7,5) * t549 - Ifges(7,6) * t106) * t466 + (t106 * t17 + t16 * t549) * mrSges(7,3) - t116 * (mrSges(7,1) * t106 + mrSges(7,2) * t549) + t48 * t488 + t396 + t518 + (-Ifges(7,2) * t106 + t49 + t97) * t491;
t528 = t568 * t351 + t352 * t569;
t527 = t351 * t569 - t568 * t352;
t168 = Ifges(6,4) * t522;
t395 = Ifges(6,5) * t88 + Ifges(6,6) * t89 + Ifges(6,3) * t385;
t445 = Ifges(6,4) * t177;
t464 = -t329 / 0.2e1;
t483 = -t177 / 0.2e1;
t485 = -t522 / 0.2e1;
t512 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t96 = Ifges(6,1) * t177 + Ifges(6,5) * t329 + t168;
t566 = t395 + t512 + (Ifges(6,5) * t522 - Ifges(6,6) * t177) * t464 + (t177 * t52 + t51 * t522) * mrSges(6,3) + (-Ifges(6,2) * t177 + t168 + t96) * t485 - t190 * (mrSges(6,1) * t177 + mrSges(6,2) * t522) + (Ifges(6,1) * t522 - t445) * t483 + t567;
t564 = -pkin(4) * t409 + t412 * pkin(9) + t528;
t563 = t413 * pkin(9) - t527;
t242 = -t308 * t354 - t358 * t369;
t169 = qJD(5) * t242 - t294 * t354 - t295 * t358;
t207 = -t272 * t354 - t273 * t358;
t415 = t169 - t207;
t243 = t308 * t358 - t354 * t369;
t170 = -qJD(5) * t243 - t294 * t358 + t295 * t354;
t206 = -t272 * t358 + t273 * t354;
t414 = t170 - t206;
t324 = t449 * t355;
t325 = t449 * t359;
t259 = t352 * t324 + t325 * t351;
t232 = -pkin(9) * t308 + t259;
t260 = t351 * t324 - t352 * t325;
t233 = -pkin(9) * t369 + t260;
t532 = t232 * t399 - t233 * t400 + t354 * t564 - t358 * t563;
t153 = t354 * t232 + t358 * t233;
t531 = -qJD(5) * t153 + t354 * t563 + t358 * t564;
t552 = pkin(10) * t414 + t532;
t551 = -pkin(5) * t409 - pkin(10) * t415 + t531;
t140 = -t215 * t351 - t420;
t117 = t140 - t519;
t141 = t352 * t215 - t208;
t118 = t141 - t535;
t343 = pkin(3) * t352 + pkin(4);
t456 = pkin(3) * t351;
t290 = t343 * t354 + t358 * t456;
t509 = -t290 * qJD(5) - t358 * t117 + t118 * t354;
t289 = t358 * t343 - t354 * t456;
t508 = t289 * qJD(5) - t354 * t117 - t358 * t118;
t455 = pkin(3) * t355;
t305 = t408 * t455 + t347;
t550 = pkin(3) * t404 - t305;
t95 = Ifges(6,2) * t522 + Ifges(6,6) * t329 + t445;
t547 = t95 / 0.2e1;
t387 = Ifges(3,5) * qJD(2) / 0.2e1;
t545 = Ifges(4,3) + Ifges(5,3);
t540 = t508 + t553;
t539 = t546 + t509;
t526 = pkin(4) * t413 + t550;
t152 = t358 * t232 - t233 * t354;
t119 = -pkin(10) * t243 + t152;
t120 = pkin(10) * t242 + t153;
t66 = t119 * t353 + t120 * t357;
t534 = -qJD(6) * t66 - t353 * t552 + t357 * t551;
t65 = t119 * t357 - t120 * t353;
t533 = qJD(6) * t65 + t353 * t551 + t357 * t552;
t530 = t250 * Ifges(5,4);
t529 = -pkin(5) * t414 + t526;
t345 = Ifges(3,4) * t408;
t432 = t315 * Ifges(4,4);
t236 = t314 * Ifges(4,2) + t340 * Ifges(4,6) + t432;
t306 = Ifges(4,4) * t314;
t237 = t315 * Ifges(4,1) + t340 * Ifges(4,5) + t306;
t370 = t255 * t359 + t256 * t355;
t446 = Ifges(4,4) * t359;
t374 = -Ifges(4,2) * t355 + t446;
t447 = Ifges(4,4) * t355;
t376 = Ifges(4,1) * t359 - t447;
t377 = mrSges(4,1) * t355 + mrSges(4,2) * t359;
t442 = Ifges(4,6) * t355;
t443 = Ifges(4,5) * t359;
t459 = t359 / 0.2e1;
t460 = -t355 / 0.2e1;
t461 = t340 / 0.2e1;
t467 = t315 / 0.2e1;
t362 = -t370 * mrSges(4,3) + t327 * t377 + t314 * t374 / 0.2e1 + t376 * t467 + (-t442 + t443) * t461 + t236 * t460 + t237 * t459;
t524 = t362 + Ifges(3,1) * t409 / 0.2e1 + t345 / 0.2e1 + t387;
t161 = Ifges(5,2) * t382 + t340 * Ifges(5,6) + t530;
t520 = t161 / 0.2e1;
t477 = -t382 / 0.2e1;
t386 = -Ifges(3,6) * qJD(2) / 0.2e1;
t515 = Ifges(5,4) * t382;
t285 = pkin(5) + t289;
t225 = t285 * t353 + t290 * t357;
t511 = -qJD(6) * t225 - t353 * t540 + t357 * t539;
t224 = t285 * t357 - t290 * t353;
t510 = qJD(6) * t224 + t353 * t539 + t357 * t540;
t310 = t359 * t319;
t454 = pkin(7) * t355;
t252 = -qJ(4) * t417 + t310 + (-pkin(3) - t454) * t360;
t342 = pkin(7) * t416;
t275 = t355 * t319 + t342;
t419 = t355 * t356;
t258 = -qJ(4) * t419 + t275;
t180 = t352 * t252 - t351 * t258;
t284 = t369 * t356;
t149 = -pkin(4) * t360 + t284 * pkin(9) + t180;
t181 = t351 * t252 + t352 * t258;
t283 = t308 * t356;
t155 = -pkin(9) * t283 + t181;
t91 = t354 * t149 + t358 * t155;
t505 = Ifges(4,5) * t264 + Ifges(4,6) * t265;
t504 = Ifges(5,5) * t201 + Ifges(5,6) * t200 + t385 * t545 + t505;
t503 = pkin(1) * mrSges(3,2) * qJD(1);
t394 = Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1;
t448 = Ifges(3,4) * t356;
t502 = -t394 * t340 - t106 * Ifges(7,5) - t386 + (t360 * Ifges(3,2) + t448) * qJD(1) / 0.2e1 - t315 * Ifges(4,5) - t314 * Ifges(4,6) - t549 * Ifges(7,6) - t329 * Ifges(6,3) - t177 * Ifges(6,5) - t522 * Ifges(6,6) - t250 * Ifges(5,5) - t382 * Ifges(5,6) - t320 * Ifges(7,3) - t51 * mrSges(6,1) + t52 * mrSges(6,2) - t16 * mrSges(7,1) + t17 * mrSges(7,2) + t131 * mrSges(5,2) - t130 * mrSges(5,1) - t255 * mrSges(4,1) + t256 * mrSges(4,2) - t545 * t461;
t501 = -t189 * mrSges(4,1) - t71 * mrSges(5,1) + t188 * mrSges(4,2) + t72 * mrSges(5,2);
t500 = t31 / 0.2e1;
t499 = t32 / 0.2e1;
t496 = t88 / 0.2e1;
t495 = t89 / 0.2e1;
t494 = pkin(1) * mrSges(3,1);
t490 = t549 / 0.2e1;
t222 = -t283 * t358 + t284 * t354;
t223 = -t283 * t354 - t284 * t358;
t143 = t222 * t357 - t223 * t353;
t487 = t143 / 0.2e1;
t144 = t222 * t353 + t223 * t357;
t486 = t144 / 0.2e1;
t484 = t522 / 0.2e1;
t482 = t177 / 0.2e1;
t481 = t200 / 0.2e1;
t480 = t201 / 0.2e1;
t479 = t222 / 0.2e1;
t478 = t223 / 0.2e1;
t476 = t382 / 0.2e1;
t475 = -t250 / 0.2e1;
t474 = t250 / 0.2e1;
t473 = t264 / 0.2e1;
t472 = t265 / 0.2e1;
t471 = -t283 / 0.2e1;
t470 = -t284 / 0.2e1;
t469 = -t314 / 0.2e1;
t468 = -t315 / 0.2e1;
t465 = t320 / 0.2e1;
t463 = t329 / 0.2e1;
t462 = -t340 / 0.2e1;
t457 = pkin(3) * t315;
t133 = t206 * t357 - t207 * t353;
t172 = t242 * t353 + t243 * t357;
t69 = -qJD(6) * t172 - t169 * t353 + t170 * t357;
t426 = t133 - t69;
t134 = t206 * t353 + t207 * t357;
t171 = t242 * t357 - t243 * t353;
t68 = qJD(6) * t171 + t169 * t357 + t170 * t353;
t425 = t134 - t68;
t422 = qJD(2) * mrSges(3,2);
t410 = t359 * t317 + t407 * t454;
t167 = -t356 * t401 + t367 * qJD(2) + (-t342 + (qJ(4) * t356 - t319) * t355) * qJD(3) + t410;
t411 = t355 * t317 + t319 * t402;
t178 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t417 + (-qJD(4) * t356 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t360) * t355 + t411;
t100 = t351 * t167 + t352 * t178;
t318 = pkin(3) * t419 + t356 * pkin(7);
t271 = pkin(3) * t525 + pkin(7) * t405;
t344 = -pkin(3) * t359 - pkin(2);
t39 = -t89 * mrSges(6,1) + t88 * mrSges(6,2);
t8 = -t32 * mrSges(7,1) + t31 * mrSges(7,2);
t246 = -pkin(3) * t265 + qJD(2) * t347;
t124 = -t200 * mrSges(5,1) + t201 * mrSges(5,2);
t90 = t358 * t149 - t354 * t155;
t99 = t352 * t167 - t178 * t351;
t380 = m(4) * t327 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t314 + mrSges(4,2) * t315 + mrSges(3,3) * t409;
t253 = pkin(4) * t283 + t318;
t205 = pkin(4) * t250 + t457;
t378 = mrSges(4,1) * t359 - mrSges(4,2) * t355;
t375 = Ifges(4,1) * t355 + t446;
t373 = Ifges(4,2) * t359 + t447;
t372 = Ifges(4,5) * t355 + Ifges(4,6) * t359;
t67 = -pkin(5) * t360 - t223 * pkin(10) + t90;
t70 = pkin(10) * t222 + t91;
t35 = -t353 * t70 + t357 * t67;
t36 = t353 * t67 + t357 * t70;
t371 = t188 * t359 - t189 * t355;
t276 = pkin(4) * t369 + t344;
t226 = -qJD(2) * t366 + t369 * t403;
t182 = -pkin(4) * t226 + t271;
t150 = -pkin(4) * t200 + t246;
t227 = -qJD(2) * t365 - t294 * t356;
t81 = pkin(4) * t407 - pkin(9) * t227 + t99;
t84 = pkin(9) * t226 + t100;
t24 = t149 * t399 - t155 * t400 + t354 * t81 + t358 * t84;
t25 = -qJD(5) * t91 - t354 * t84 + t358 * t81;
t322 = mrSges(3,3) * t408 - t422;
t274 = -pkin(7) * t418 + t310;
t270 = mrSges(4,1) * t340 - mrSges(4,3) * t315;
t269 = -mrSges(4,2) * t340 + mrSges(4,3) * t314;
t267 = -pkin(7) * t390 + t298;
t245 = -mrSges(4,2) * t385 + mrSges(4,3) * t265;
t244 = mrSges(4,1) * t385 - mrSges(4,3) * t264;
t221 = mrSges(5,1) * t340 - mrSges(5,3) * t250;
t220 = -mrSges(5,2) * t340 + mrSges(5,3) * t382;
t213 = -qJD(3) * t275 + t410;
t212 = (-t356 * t406 - t360 * t404) * pkin(7) + t411;
t204 = -mrSges(4,1) * t265 + mrSges(4,2) * t264;
t202 = -pkin(5) * t242 + t276;
t192 = t264 * Ifges(4,1) + t265 * Ifges(4,4) + Ifges(4,5) * t385;
t191 = t264 * Ifges(4,4) + t265 * Ifges(4,2) + Ifges(4,6) * t385;
t184 = mrSges(5,1) * t385 - mrSges(5,3) * t201;
t183 = -mrSges(5,2) * t385 + mrSges(5,3) * t200;
t179 = -mrSges(5,1) * t382 + mrSges(5,2) * t250;
t165 = -pkin(5) * t222 + t253;
t162 = t250 * Ifges(5,1) + t340 * Ifges(5,5) + t515;
t157 = mrSges(6,1) * t329 - mrSges(6,3) * t177;
t156 = -mrSges(6,2) * t329 + mrSges(6,3) * t522;
t123 = t201 * Ifges(5,1) + t200 * Ifges(5,4) + Ifges(5,5) * t385;
t122 = t201 * Ifges(5,4) + t200 * Ifges(5,2) + Ifges(5,6) * t385;
t121 = pkin(5) * t177 + t205;
t115 = -qJD(5) * t223 + t226 * t358 - t227 * t354;
t114 = qJD(5) * t222 + t226 * t354 + t227 * t358;
t108 = -mrSges(6,1) * t522 + mrSges(6,2) * t177;
t93 = mrSges(7,1) * t320 - mrSges(7,3) * t106;
t92 = -mrSges(7,2) * t320 + mrSges(7,3) * t549;
t82 = -pkin(5) * t115 + t182;
t80 = -mrSges(6,2) * t385 + mrSges(6,3) * t89;
t79 = mrSges(6,1) * t385 - mrSges(6,3) * t88;
t55 = -pkin(5) * t89 + t150;
t50 = -mrSges(7,1) * t549 + mrSges(7,2) * t106;
t42 = -qJD(6) * t144 - t114 * t353 + t115 * t357;
t41 = qJD(6) * t143 + t114 * t357 + t115 * t353;
t38 = t88 * Ifges(6,1) + t89 * Ifges(6,4) + Ifges(6,5) * t385;
t37 = t88 * Ifges(6,4) + t89 * Ifges(6,2) + Ifges(6,6) * t385;
t27 = -mrSges(7,2) * t385 + mrSges(7,3) * t32;
t26 = mrSges(7,1) * t385 - mrSges(7,3) * t31;
t21 = pkin(10) * t115 + t24;
t20 = pkin(5) * t407 - pkin(10) * t114 + t25;
t19 = t357 * t43 - t428;
t18 = -t353 * t43 - t427;
t7 = t31 * Ifges(7,1) + t32 * Ifges(7,4) + Ifges(7,5) * t385;
t6 = t31 * Ifges(7,4) + t32 * Ifges(7,2) + Ifges(7,6) * t385;
t5 = -qJD(6) * t36 + t20 * t357 - t21 * t353;
t4 = qJD(6) * t35 + t20 * t353 + t21 * t357;
t1 = [m(7) * (t116 * t82 + t16 * t5 + t165 * t55 + t17 * t4 + t2 * t36 + t3 * t35) + m(6) * (t12 * t91 + t13 * t90 + t150 * t253 + t182 * t190 + t24 * t52 + t25 * t51) + m(5) * (t100 * t131 + t130 * t99 + t180 * t71 + t181 * t72 + t246 * t318 + t263 * t271) + t114 * t96 / 0.2e1 + t116 * (-mrSges(7,1) * t42 + mrSges(7,2) * t41) + t226 * t520 + t55 * (-mrSges(7,1) * t143 + mrSges(7,2) * t144) + (t380 * pkin(7) + t387 - 0.2e1 * t503 + t524) * t405 + (-Ifges(5,4) * t284 - Ifges(5,2) * t283) * t481 + (-Ifges(5,1) * t284 - Ifges(5,4) * t283) * t480 + t246 * (mrSges(5,1) * t283 - mrSges(5,2) * t284) + (-t130 * t227 + t131 * t226 - t283 * t72 + t284 * t71) * mrSges(5,3) + t115 * t547 + (pkin(7) * t204 + t191 * t460 + t376 * t473 + t374 * t472 + t192 * t459 + (-t188 * t355 - t189 * t359) * mrSges(4,3) + (-t359 * t236 / 0.2e1 + t373 * t469 + t375 * t468 + t327 * t378 + t372 * t462 + t237 * t460 + (t255 * t355 - t256 * t359) * mrSges(4,3)) * qJD(3) + (t386 + (-0.2e1 * t494 + Ifges(5,5) * t470 + Ifges(5,6) * t471 + Ifges(6,5) * t478 + Ifges(6,6) * t479 + Ifges(7,5) * t486 + Ifges(7,6) * t487 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t443 / 0.2e1 - t442 / 0.2e1) * t356) * qJD(1) - pkin(7) * t322 - t502) * qJD(2)) * t356 + m(4) * (t188 * t275 + t189 * t274 + t256 * t212 + t255 * t213) + (Ifges(7,4) * t144 + Ifges(7,2) * t143) * t499 + (Ifges(6,4) * t223 + Ifges(6,2) * t222) * t495 + t4 * t92 + t5 * t93 + t90 * t79 + t91 * t80 + t82 * t50 + t42 * t48 / 0.2e1 + t41 * t49 / 0.2e1 + t35 * t26 + t36 * t27 + (Ifges(7,4) * t41 + Ifges(7,2) * t42) * t490 + (Ifges(5,1) * t227 + Ifges(5,4) * t226) * t474 + (Ifges(5,4) * t227 + Ifges(5,2) * t226) * t476 + t38 * t478 + t37 * t479 + (Ifges(6,1) * t114 + Ifges(6,4) * t115) * t482 + (Ifges(6,4) * t114 + Ifges(6,2) * t115) * t484 + t7 * t486 + t6 * t487 + (Ifges(7,1) * t41 + Ifges(7,4) * t42) * t488 + (Ifges(6,5) * t114 + Ifges(6,6) * t115) * t463 + (Ifges(7,5) * t41 + Ifges(7,6) * t42) * t465 + t123 * t470 + t122 * t471 + t24 * t156 + t25 * t157 + t165 * t8 - (t396 + t395 + t504 + t505) * t360 / 0.2e1 + t182 * t108 + t181 * t183 + t180 * t184 + t190 * (-mrSges(6,1) * t115 + mrSges(6,2) * t114) + t100 * t220 + t99 * t221 + (Ifges(5,5) * t227 + Ifges(5,6) * t226) * t461 + t150 * (-mrSges(6,1) * t222 + mrSges(6,2) * t223) + t227 * t162 / 0.2e1 + (-Ifges(5,5) * t480 - Ifges(6,5) * t496 - Ifges(7,5) * t500 - Ifges(5,6) * t481 - Ifges(6,6) * t495 - Ifges(7,6) * t499 + (0.3e1 / 0.2e1 * Ifges(3,4) * t405 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + (m(4) * pkin(7) + t377) * pkin(7) - t394) * t407) * qJD(1) + t501 - t512 - t518) * t360 + t253 * t39 + t263 * (-mrSges(5,1) * t226 + mrSges(5,2) * t227) + (Ifges(6,1) * t223 + Ifges(6,4) * t222) * t496 + t212 * t269 + t213 * t270 + t271 * t179 + t274 * t244 + t275 * t245 + t318 * t124 + (Ifges(7,1) * t144 + Ifges(7,4) * t143) * t500 + (-t114 * t51 + t115 * t52 + t12 * t222 - t13 * t223) * mrSges(6,3) + (t143 * t2 - t144 * t3 - t16 * t41 + t17 * t42) * mrSges(7,3); t526 * t108 + t527 * t220 - m(4) * (t255 * t266 + t256 * t267) + t528 * t221 + t529 * t50 + (t170 / 0.2e1 - t206 / 0.2e1) * t95 + (Ifges(5,1) * t308 - Ifges(5,4) * t369) * t480 + (Ifges(5,4) * t308 - Ifges(5,2) * t369) * t481 + t246 * (mrSges(5,1) * t369 + mrSges(5,2) * t308) + (t130 * t412 - t131 * t413 - t308 * t71 - t369 * t72) * mrSges(5,3) - t369 * t122 / 0.2e1 + (t272 / 0.2e1 - t294 / 0.2e1) * t161 + (-Ifges(5,1) * t273 - Ifges(5,4) * t272) * t475 + (-Ifges(5,4) * t273 - Ifges(5,2) * t272) * t477 + (-Ifges(5,5) * t273 - Ifges(5,6) * t272) * t462 + (t273 / 0.2e1 - t295 / 0.2e1) * t162 + (-t134 / 0.2e1 + t68 / 0.2e1) * t49 + ((t387 + t503 - t345 / 0.2e1 + ((-m(4) * pkin(2) - mrSges(3,1) - t378) * qJD(2) - t380) * pkin(7) - t524) * t360 + ((t494 + t448 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t360) * qJD(1) + t386 + (t322 + t422) * pkin(7) + t502) * t356 + (Ifges(5,5) * t308 + Ifges(6,5) * t243 + Ifges(7,5) * t172 - Ifges(5,6) * t369 + Ifges(6,6) * t242 + Ifges(7,6) * t171 + t372) * t407 / 0.2e1) * qJD(1) + (t169 / 0.2e1 - t207 / 0.2e1) * t96 + (-Ifges(5,4) * t295 - Ifges(5,2) * t294) * t476 + (-Ifges(5,5) * t295 - Ifges(5,6) * t294) * t461 + (-Ifges(5,1) * t295 - Ifges(5,4) * t294) * t474 + (t528 * t130 + t527 * t131 + t246 * t344 + t259 * t71 + t260 * t72 + t263 * t550) * m(5) + (-t133 / 0.2e1 + t69 / 0.2e1) * t48 + t531 * t157 + t532 * t156 + (t12 * t153 + t13 * t152 + t150 * t276 + t190 * t526 + t51 * t531 + t52 * t532) * m(6) + t533 * t92 + (t116 * t529 + t16 * t534 + t17 * t533 + t2 * t66 + t202 * t55 + t3 * t65) * m(7) + t534 * t93 + t65 * t26 + t66 * t27 + (Ifges(7,1) * t134 + Ifges(7,4) * t133) * t489 + (Ifges(7,4) * t68 + Ifges(7,2) * t69) * t490 + (Ifges(7,4) * t134 + Ifges(7,2) * t133) * t491 + (Ifges(6,4) * t243 + Ifges(6,2) * t242) * t495 + (Ifges(6,1) * t243 + Ifges(6,4) * t242) * t496 + t375 * t473 + (Ifges(6,1) * t169 + Ifges(6,4) * t170) * t482 + (Ifges(6,1) * t207 + Ifges(6,4) * t206) * t483 + (Ifges(6,4) * t169 + Ifges(6,2) * t170) * t484 + (Ifges(6,4) * t207 + Ifges(6,2) * t206) * t485 + (Ifges(7,1) * t68 + Ifges(7,4) * t69) * t488 + (Ifges(6,5) * t207 + Ifges(6,6) * t206) * t464 + (Ifges(7,5) * t68 + Ifges(7,6) * t69) * t465 + (Ifges(7,5) * t134 + Ifges(7,6) * t133) * t466 + t373 * t472 + (mrSges(7,1) * t426 - mrSges(7,2) * t425) * t116 + (t16 * t425 - t17 * t426 + t171 * t2 - t172 * t3) * mrSges(7,3) + (Ifges(7,4) * t172 + Ifges(7,2) * t171) * t499 + (Ifges(7,1) * t172 + Ifges(7,4) * t171) * t500 + t371 * mrSges(4,3) + t152 * t79 + t153 * t80 + t171 * t6 / 0.2e1 + t55 * (-mrSges(7,1) * t171 + mrSges(7,2) * t172) + t172 * t7 / 0.2e1 + t202 * t8 - pkin(2) * t204 + t191 * t459 + (Ifges(6,5) * t169 + Ifges(6,6) * t170) * t463 + t242 * t37 / 0.2e1 + t243 * t38 / 0.2e1 + t150 * (-mrSges(6,1) * t242 + mrSges(6,2) * t243) + t259 * t184 + t260 * t183 - t267 * t269 - t266 * t270 + t276 * t39 + (t179 * t455 + t362) * qJD(3) + (m(4) * (-qJD(3) * t370 + t371) - t244 * t355 + t245 * t359 + (-t269 * t355 - t270 * t359) * qJD(3)) * pkin(8) - t305 * t179 + t308 * t123 / 0.2e1 + t344 * t124 + t355 * t192 / 0.2e1 + (mrSges(5,1) * t413 - mrSges(5,2) * t412) * t263 + (-mrSges(6,1) * t414 + mrSges(6,2) * t415) * t190 + (t12 * t242 - t13 * t243 + t414 * t52 - t415 * t51) * mrSges(6,3); (-t130 * t140 - t131 * t141 - t263 * t457 + (t351 * t72 + t352 * t71) * pkin(3)) * m(5) + t504 + t566 + (-t179 * t315 + t183 * t351 + t184 * t352) * pkin(3) - t501 - t121 * t50 + (Ifges(5,1) * t382 - t530) * t475 + (-Ifges(5,2) * t250 + t162 + t515) * t477 + (Ifges(4,5) * t314 + Ifges(5,5) * t382 - Ifges(4,6) * t315 - Ifges(5,6) * t250) * t462 - t263 * (mrSges(5,1) * t250 + mrSges(5,2) * t382) + (t130 * t382 + t131 * t250) * mrSges(5,3) + t250 * t520 + (t255 * t314 + t256 * t315) * mrSges(4,3) + t236 * t467 + (Ifges(4,1) * t314 - t432) * t468 + t508 * t156 + t509 * t157 + (t12 * t290 + t13 * t289 - t190 * t205 + t508 * t52 + t509 * t51) * m(6) + t510 * t92 + t511 * t93 + (-t116 * t121 + t16 * t511 + t17 * t510 + t2 * t225 + t224 * t3) * m(7) - t205 * t108 - t141 * t220 - t140 * t221 + t224 * t26 + t225 * t27 - t255 * t269 + t256 * t270 + (-Ifges(4,2) * t315 + t237 + t306) * t469 + t289 * t79 + t290 * t80 + t177 * t547 - t327 * (mrSges(4,1) * t315 + mrSges(4,2) * t314); t106 * t93 - t549 * t92 - t522 * t156 + t177 * t157 - t382 * t220 + t250 * t221 + t124 + t39 + t8 + (t106 * t16 - t17 * t549 + t55) * m(7) + (t177 * t51 - t52 * t522 + t150) * m(6) + (t130 * t250 - t131 * t382 + t246) * m(5); -t19 * t92 - t18 * t93 - t51 * t156 + t52 * t157 - m(7) * (t16 * t18 + t17 * t19) + (-t177 * t50 + t357 * t26 + t353 * t27 + (-t353 * t93 + t357 * t92) * qJD(6) + (-t116 * t177 + t2 * t353 + t3 * t357 + (-t16 * t353 + t17 * t357) * qJD(6)) * m(7)) * pkin(5) + t95 * t482 + t566; -t16 * t92 + t17 * t93 + t567;];
tauc  = t1(:);
