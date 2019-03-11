% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:11
% EndTime: 2019-03-09 09:34:58
% DurationCPUTime: 31.48s
% Computational Cost: add. (12450->837), mult. (25925->1097), div. (0->0), fcn. (17613->14), ass. (0->361)
t317 = sin(qJ(2));
t295 = t317 * qJD(1);
t320 = cos(qJ(5));
t312 = cos(pkin(10));
t386 = t312 * t295;
t311 = sin(pkin(10));
t316 = sin(qJ(5));
t425 = t311 * t316;
t187 = -t295 * t425 + t320 * t386;
t402 = qJD(5) * t320;
t403 = qJD(5) * t316;
t216 = -t311 * t403 + t312 * t402;
t529 = t216 + t187;
t347 = t320 * t311 + t316 * t312;
t336 = t347 * t317;
t188 = qJD(1) * t336;
t217 = t347 * qJD(5);
t528 = -t217 - t188;
t562 = qJD(2) / 0.2e1;
t561 = -mrSges(6,3) - mrSges(7,3);
t286 = pkin(2) * t295;
t321 = cos(qJ(2));
t436 = qJ(3) * t321;
t351 = qJ(4) * t317 - t436;
t206 = qJD(1) * t351 + t286;
t409 = qJD(1) * t321;
t287 = pkin(7) * t409;
t246 = pkin(3) * t409 + t287;
t142 = -t206 * t311 + t312 * t246;
t424 = t311 * t317;
t344 = pkin(4) * t321 - pkin(8) * t424;
t110 = qJD(1) * t344 + t142;
t143 = t312 * t206 + t311 * t246;
t124 = pkin(8) * t386 + t143;
t461 = pkin(2) + qJ(4);
t460 = -pkin(8) - t461;
t251 = t460 * t311;
t252 = t460 * t312;
t534 = -qJD(4) * t347 - t316 * t110 - t320 * t124 - t251 * t403 + t252 * t402;
t168 = t320 * t251 + t316 * t252;
t404 = qJD(4) * t312;
t405 = qJD(4) * t311;
t533 = -qJD(5) * t168 + (-t110 - t404) * t320 + (t124 + t405) * t316;
t399 = qJD(1) * qJD(2);
t249 = -t321 * qJDD(1) + t317 * t399;
t190 = qJDD(2) * t312 + t249 * t311;
t250 = qJDD(1) * t317 + t321 * t399;
t406 = qJD(3) * t317;
t345 = -qJD(4) * t321 - t406;
t435 = qJDD(1) * pkin(1);
t352 = -qJ(3) * t250 - t435;
t106 = qJD(1) * t345 + t249 * t461 + t352;
t234 = t250 * pkin(7);
t379 = qJDD(3) + t234;
t141 = pkin(3) * t250 - qJD(2) * qJD(4) - qJDD(2) * t461 + t379;
t62 = -t106 * t311 + t312 * t141;
t49 = pkin(4) * t250 - pkin(8) * t190 + t62;
t189 = -qJDD(2) * t311 + t249 * t312;
t63 = t312 * t106 + t311 * t141;
t55 = pkin(8) * t189 + t63;
t297 = t317 * qJ(3);
t378 = -pkin(1) - t297;
t191 = (-t321 * t461 + t378) * qJD(1);
t285 = pkin(7) * t295;
t245 = -pkin(3) * t295 - t285;
t200 = -qJD(2) * t461 + qJD(3) - t245;
t118 = -t191 * t311 + t312 * t200;
t401 = t312 * qJD(2);
t235 = t311 * t409 - t401;
t92 = pkin(4) * t295 + pkin(8) * t235 + t118;
t119 = t312 * t191 + t311 * t200;
t236 = -qJD(2) * t311 - t312 * t409;
t99 = pkin(8) * t236 + t119;
t11 = t316 * t49 + t320 * t55 + t92 * t402 - t403 * t99;
t154 = t235 * t320 - t236 * t316;
t74 = qJD(5) * t154 + t189 * t320 - t190 * t316;
t10 = pkin(9) * t74 + t11;
t319 = cos(qJ(6));
t274 = t295 + qJD(5);
t47 = -t316 * t99 + t320 * t92;
t38 = pkin(9) * t154 + t47;
t36 = pkin(5) * t274 + t38;
t315 = sin(qJ(6));
t373 = t235 * t316 + t320 * t236;
t48 = t316 * t92 + t320 * t99;
t39 = pkin(9) * t373 + t48;
t442 = t315 * t39;
t13 = t319 * t36 - t442;
t12 = -qJD(5) * t48 - t316 * t55 + t320 * t49;
t237 = qJDD(5) + t250;
t73 = qJD(5) * t373 + t189 * t316 + t190 * t320;
t9 = pkin(5) * t237 - pkin(9) * t73 + t12;
t2 = qJD(6) * t13 + t10 * t319 + t315 * t9;
t560 = t2 * mrSges(7,2);
t441 = t319 * t39;
t14 = t315 * t36 + t441;
t3 = -qJD(6) * t14 - t10 * t315 + t319 * t9;
t559 = t3 * mrSges(7,1);
t558 = t11 * mrSges(6,2);
t557 = t12 * mrSges(6,1);
t519 = m(5) + m(6) + m(7);
t390 = m(4) + t519;
t541 = Ifges(3,1) + Ifges(5,3);
t540 = -Ifges(4,4) + Ifges(3,5);
t539 = Ifges(4,5) - Ifges(3,6);
t556 = -pkin(5) * t409 - pkin(9) * t528 + t533;
t555 = pkin(9) * t529 - t534;
t362 = t321 * mrSges(4,2) - t317 * mrSges(4,3);
t366 = mrSges(3,1) * t321 - mrSges(3,2) * t317;
t554 = t362 - t366;
t306 = pkin(10) + qJ(5);
t291 = sin(t306);
t296 = t311 * pkin(4);
t253 = pkin(5) * t291 + t296;
t294 = qJ(6) + t306;
t279 = sin(t294);
t280 = cos(t294);
t292 = cos(t306);
t363 = mrSges(5,1) * t311 + mrSges(5,2) * t312;
t553 = -m(6) * t296 - m(7) * t253 - t291 * mrSges(6,1) - t279 * mrSges(7,1) - t292 * mrSges(6,2) - t280 * mrSges(7,2) - t363;
t313 = -pkin(8) - qJ(4);
t305 = -pkin(9) + t313;
t552 = -m(7) * (-pkin(2) + t305) + m(5) * t461 + mrSges(5,3) - m(6) * (-pkin(2) + t313) - t561;
t450 = Ifges(4,6) * t321;
t355 = -t317 * Ifges(4,2) - t450;
t551 = t14 * mrSges(7,2) + t48 * mrSges(6,2) + Ifges(4,4) * t562 + qJD(1) * t355 / 0.2e1 - t13 * mrSges(7,1) - t47 * mrSges(6,1);
t318 = sin(qJ(1));
t550 = g(2) * t318;
t111 = t187 * t319 - t188 * t315;
t346 = -t312 * t320 + t425;
t530 = -t315 * t347 - t319 * t346;
t85 = qJD(6) * t530 + t216 * t319 - t217 * t315;
t549 = t111 + t85;
t281 = t312 * pkin(4) + pkin(3);
t364 = mrSges(5,1) * t312 - mrSges(5,2) * t311;
t466 = pkin(5) * t292;
t548 = -t364 - m(5) * pkin(3) - m(6) * t281 - m(7) * (t281 + t466) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t547 = -t11 * t347 + t12 * t346 - t47 * t528 - t48 * t529;
t546 = t321 * t561 + t554;
t545 = t154 * t315 + t319 * t373;
t80 = t154 * t319 - t315 * t373;
t159 = t315 * t346 - t319 * t347;
t374 = t187 * t315 + t319 * t188;
t544 = -t13 * t374 - t159 * t2 + t3 * t530;
t26 = qJD(6) * t545 + t315 * t74 + t319 * t73;
t502 = t26 / 0.2e1;
t27 = qJD(6) * t80 - t315 * t73 + t319 * t74;
t501 = t27 / 0.2e1;
t494 = t73 / 0.2e1;
t493 = t74 / 0.2e1;
t231 = qJDD(6) + t237;
t478 = t231 / 0.2e1;
t477 = t237 / 0.2e1;
t543 = -t250 / 0.2e1;
t167 = -t251 * t316 + t320 * t252;
t127 = pkin(9) * t346 + t167;
t128 = -pkin(9) * t347 + t168;
t64 = t127 * t319 - t128 * t315;
t538 = qJD(6) * t64 + t315 * t556 - t319 * t555;
t65 = t127 * t315 + t128 * t319;
t537 = -qJD(6) * t65 + t315 * t555 + t319 * t556;
t503 = m(7) * pkin(5);
t535 = -mrSges(6,1) - t503;
t276 = qJ(3) + t296;
t302 = t321 * pkin(2);
t411 = t302 + t297;
t367 = qJ(4) * t321 + t411;
t232 = -pkin(1) - t367;
t486 = pkin(3) + pkin(7);
t264 = t486 * t317;
t241 = t312 * t264;
t134 = pkin(4) * t317 + t241 + (pkin(8) * t321 - t232) * t311;
t163 = t312 * t232 + t311 * t264;
t421 = t312 * t321;
t140 = -pkin(8) * t421 + t163;
t67 = t316 * t134 + t320 * t140;
t193 = -t281 * t295 - t285;
t531 = pkin(5) * t529 + qJD(3) - t193;
t388 = mrSges(4,1) * t409;
t261 = -qJD(2) * mrSges(4,3) - t388;
t527 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t409 - t261;
t389 = mrSges(4,1) * t295;
t526 = -mrSges(3,3) * t295 - t389 + (mrSges(3,1) - mrSges(4,2)) * qJD(2);
t451 = Ifges(4,6) * t317;
t525 = t317 * (-Ifges(4,2) * t321 + t451) + t321 * (Ifges(4,3) * t317 - t450);
t524 = t317 * t539 + t321 * t540;
t145 = -mrSges(5,2) * t250 + mrSges(5,3) * t189;
t146 = mrSges(5,1) * t250 - mrSges(5,3) * t190;
t523 = t311 * t145 + t312 * t146;
t233 = t249 * pkin(7);
t522 = -t233 * t321 + t234 * t317;
t192 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t233;
t205 = -qJDD(2) * pkin(2) + t379;
t521 = -t192 * t321 + t205 * t317;
t353 = t311 * t63 + t312 * t62;
t322 = cos(qJ(1));
t520 = g(1) * t322 + t550;
t166 = -mrSges(5,1) * t236 - mrSges(5,2) * t235;
t90 = -mrSges(6,1) * t373 - mrSges(6,2) * t154;
t518 = t261 - t166 - t90;
t310 = qJD(2) * qJ(3);
t211 = qJD(4) + t310 + t246;
t164 = -pkin(4) * t236 + t211;
t93 = -pkin(5) * t373 + t164;
t515 = -mrSges(7,1) * t93 + mrSges(7,3) * t14;
t514 = mrSges(7,2) * t93 - mrSges(7,3) * t13;
t266 = qJD(6) + t274;
t284 = Ifges(3,4) * t409;
t513 = Ifges(3,5) * qJD(2) - t235 * Ifges(5,5) - Ifges(6,5) * t154 - Ifges(7,5) * t80 + t236 * Ifges(5,6) + Ifges(6,6) * t373 + Ifges(7,6) * t545 + t274 * Ifges(6,3) + t266 * Ifges(7,3) + t295 * t541 + t284;
t440 = t321 * mrSges(4,3);
t512 = -t440 + t553 * t321 + (m(4) * pkin(2) - mrSges(4,2) + t552) * t317;
t354 = -Ifges(4,3) * t321 - t451;
t509 = t312 * (-t235 * Ifges(5,4) + t236 * Ifges(5,2) + Ifges(5,6) * t295) + t311 * (-t235 * Ifges(5,1) + t236 * Ifges(5,4) + Ifges(5,5) * t295) + Ifges(4,5) * qJD(2) + qJD(1) * t354;
t343 = t378 - t302;
t226 = t343 * qJD(1);
t254 = -qJD(2) * pkin(2) + qJD(3) + t285;
t259 = -t287 - t310;
t454 = Ifges(5,4) * t311;
t358 = Ifges(5,2) * t312 + t454;
t453 = Ifges(5,4) * t312;
t361 = Ifges(5,1) * t311 + t453;
t507 = t236 * (Ifges(5,6) * t321 + t317 * t358) / 0.2e1 - t235 * (Ifges(5,5) * t321 + t317 * t361) / 0.2e1 - t211 * t317 * t364 + m(4) * (t254 * t321 + t259 * t317) * pkin(7) + t119 * (mrSges(5,3) * t312 * t317 - mrSges(5,2) * t321) + t118 * (mrSges(5,1) * t321 - mrSges(5,3) * t424) + t226 * (-mrSges(4,2) * t317 - t440);
t505 = Ifges(7,4) * t502 + Ifges(7,2) * t501 + Ifges(7,6) * t478;
t504 = Ifges(7,1) * t502 + Ifges(7,4) * t501 + Ifges(7,5) * t478;
t500 = Ifges(6,4) * t494 + Ifges(6,2) * t493 + Ifges(6,6) * t477;
t499 = Ifges(6,1) * t494 + Ifges(6,4) * t493 + Ifges(6,5) * t477;
t468 = Ifges(7,4) * t80;
t41 = Ifges(7,2) * t545 + Ifges(7,6) * t266 - t468;
t498 = -t41 / 0.2e1;
t497 = t41 / 0.2e1;
t78 = Ifges(7,4) * t545;
t42 = -Ifges(7,1) * t80 + Ifges(7,5) * t266 + t78;
t496 = -t42 / 0.2e1;
t495 = t42 / 0.2e1;
t152 = Ifges(6,4) * t373;
t77 = -Ifges(6,1) * t154 + t274 * Ifges(6,5) + t152;
t491 = t77 / 0.2e1;
t490 = -t545 / 0.2e1;
t489 = t545 / 0.2e1;
t488 = t80 / 0.2e1;
t487 = -t80 / 0.2e1;
t485 = -t190 * Ifges(5,4) / 0.2e1 - t189 * Ifges(5,2) / 0.2e1 + Ifges(5,6) * t543;
t484 = -t373 / 0.2e1;
t483 = t373 / 0.2e1;
t482 = t154 / 0.2e1;
t481 = -t154 / 0.2e1;
t480 = t189 / 0.2e1;
t479 = t190 / 0.2e1;
t476 = t250 / 0.2e1;
t475 = -t266 / 0.2e1;
t474 = t266 / 0.2e1;
t473 = -t274 / 0.2e1;
t472 = t274 / 0.2e1;
t467 = pkin(5) * t154;
t465 = pkin(7) * t317;
t462 = g(3) * t321;
t300 = t321 * pkin(7);
t459 = mrSges(7,1) * t280;
t458 = mrSges(6,3) * t373;
t457 = mrSges(6,3) * t154;
t456 = Ifges(3,4) * t317;
t455 = Ifges(3,4) * t321;
t452 = Ifges(6,4) * t154;
t427 = t305 * t321;
t423 = t311 * t321;
t420 = t313 * t321;
t419 = t317 * t322;
t418 = t318 * t279;
t417 = t318 * t280;
t416 = t318 * t291;
t415 = t318 * t292;
t414 = t321 * t322;
t408 = qJD(2) * t317;
t289 = pkin(2) * t408;
t169 = qJD(2) * t351 + t289 + t345;
t407 = qJD(2) * t321;
t248 = t486 * t407;
t122 = t312 * t169 + t311 * t248;
t179 = t280 * t419 - t418;
t180 = t279 * t419 + t417;
t413 = t179 * mrSges(7,1) - t180 * mrSges(7,2);
t181 = t279 * t322 + t317 * t417;
t182 = t280 * t322 - t317 * t418;
t412 = t181 * mrSges(7,1) + t182 * mrSges(7,2);
t265 = t321 * pkin(3) + t300;
t410 = t322 * pkin(1) + t318 * pkin(7);
t394 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t231;
t393 = Ifges(6,5) * t73 + Ifges(6,6) * t74 + Ifges(6,3) * t237;
t215 = pkin(4) * t421 + t265;
t37 = -t74 * mrSges(6,1) + t73 * mrSges(6,2);
t8 = -t27 * mrSges(7,1) + t26 * mrSges(7,2);
t377 = -t399 / 0.2e1;
t208 = t250 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t123 = -t189 * mrSges(5,1) + t190 * mrSges(5,2);
t66 = t320 * t134 - t140 * t316;
t121 = -t169 * t311 + t312 * t248;
t365 = mrSges(3,1) * t317 + mrSges(3,2) * t321;
t360 = Ifges(3,2) * t321 + t456;
t356 = Ifges(5,5) * t311 + Ifges(5,6) * t312;
t203 = t347 * t321;
t56 = pkin(5) * t317 + pkin(9) * t203 + t66;
t202 = t346 * t321;
t57 = pkin(9) * t202 + t67;
t30 = -t315 * t57 + t319 * t56;
t31 = t315 * t56 + t319 * t57;
t131 = t202 * t319 + t203 * t315;
t132 = t202 * t315 - t203 * t319;
t342 = t394 + t559 - t560;
t341 = pkin(1) * t365;
t195 = t292 * t419 - t416;
t197 = t291 * t322 + t317 * t415;
t339 = t317 * (Ifges(3,1) * t321 - t456);
t194 = (-pkin(7) - t281) * t408;
t100 = qJD(2) * t344 + t121;
t107 = pkin(8) * t317 * t401 + t122;
t32 = t316 * t100 + t320 * t107 + t134 * t402 - t140 * t403;
t33 = -qJD(5) * t67 + t320 * t100 - t107 * t316;
t327 = t317 * (Ifges(5,3) * t321 + t317 * t356);
t326 = qJD(6) * t159 - t216 * t315 - t319 * t217;
t153 = -pkin(3) * t249 + qJDD(4) - t192;
t105 = -pkin(4) * t189 + t153;
t256 = -pkin(1) - t411;
t255 = t321 * t279 * mrSges(7,2);
t247 = t486 * t408;
t244 = -qJ(3) * t409 + t286;
t243 = t362 * qJD(1);
t222 = Ifges(3,6) * qJD(2) + qJD(1) * t360;
t209 = -qJ(3) * t407 + t289 - t406;
t207 = mrSges(4,1) * t249 - qJDD(2) * mrSges(4,3);
t198 = t292 * t322 - t317 * t416;
t196 = t291 * t419 + t415;
t186 = pkin(5) * t347 + t276;
t185 = mrSges(5,1) * t295 + mrSges(5,3) * t235;
t184 = -mrSges(5,2) * t295 + mrSges(5,3) * t236;
t162 = -t232 * t311 + t241;
t151 = -pkin(5) * t202 + t215;
t147 = pkin(2) * t249 - qJD(3) * t295 + t352;
t139 = t217 * t321 - t346 * t408;
t138 = qJD(2) * t336 + qJD(5) * t202;
t130 = mrSges(6,1) * t274 + t457;
t129 = -mrSges(6,2) * t274 + t458;
t102 = t190 * Ifges(5,1) + t189 * Ifges(5,4) + t250 * Ifges(5,5);
t98 = -pkin(5) * t139 + t194;
t76 = Ifges(6,2) * t373 + t274 * Ifges(6,6) - t452;
t69 = mrSges(7,1) * t266 + mrSges(7,3) * t80;
t68 = -mrSges(7,2) * t266 + mrSges(7,3) * t545;
t61 = -mrSges(6,2) * t237 + mrSges(6,3) * t74;
t60 = mrSges(6,1) * t237 - mrSges(6,3) * t73;
t51 = -qJD(6) * t132 - t138 * t315 + t139 * t319;
t50 = qJD(6) * t131 + t138 * t319 + t139 * t315;
t44 = -pkin(5) * t74 + t105;
t43 = -mrSges(7,1) * t545 - mrSges(7,2) * t80;
t22 = pkin(9) * t139 + t32;
t21 = pkin(5) * t407 - pkin(9) * t138 + t33;
t18 = -mrSges(7,2) * t231 + mrSges(7,3) * t27;
t17 = mrSges(7,1) * t231 - mrSges(7,3) * t26;
t16 = t319 * t38 - t442;
t15 = -t315 * t38 - t441;
t5 = -qJD(6) * t31 + t21 * t319 - t22 * t315;
t4 = qJD(6) * t30 + t21 * t315 + t22 * t319;
t1 = [-t249 * t360 / 0.2e1 + t147 * t362 + t63 * (-mrSges(5,2) * t317 - mrSges(5,3) * t421) - t102 * t423 / 0.2e1 + t62 * (mrSges(5,1) * t317 + mrSges(5,3) * t423) - t341 * t399 + (Ifges(6,5) * t138 + Ifges(6,6) * t139) * t472 + m(5) * (t118 * t121 + t119 * t122 + t153 * t265 + t162 * t62 + t163 * t63 - t211 * t247) + (-m(3) * t410 - t196 * mrSges(6,1) - t180 * mrSges(7,1) - t195 * mrSges(6,2) - t179 * mrSges(7,2) - t363 * t419 + (-m(5) * qJ(4) - mrSges(5,3)) * t414 - t390 * (pkin(2) * t414 + qJ(3) * t419 + t410) + t548 * t318 + (-m(7) * (t253 * t317 - t427) - m(6) * (pkin(4) * t424 - t420) - mrSges(2,1) + t546) * t322) * g(2) + (-t198 * mrSges(6,1) - t182 * mrSges(7,1) + t197 * mrSges(6,2) + t181 * mrSges(7,2) + (-m(4) * t343 + mrSges(2,1) + (-m(7) * (-qJ(3) - t253) + m(5) * qJ(3) + t363 + m(6) * t276) * t317 + (m(3) + t519) * pkin(1) + t552 * t321 - t554) * t318 + ((-m(3) - t390) * pkin(7) + t548) * t322) * g(1) + t317 * t557 + t317 * t559 + t525 * t377 + (-qJDD(2) * mrSges(3,1) + t208) * t465 + (t317 * t540 - t321 * t539) * qJDD(2) / 0.2e1 + (-Ifges(3,4) * t249 + Ifges(3,5) * qJDD(2) + Ifges(5,5) * t190 + Ifges(5,6) * t189 + t250 * t541 + t393 + t394) * t317 / 0.2e1 + (t317 * t541 - t321 * t356 + t455) * t476 + (Ifges(7,5) * t50 + Ifges(7,6) * t51) * t474 + (-t13 * t50 + t131 * t2 - t132 * t3 + t14 * t51) * mrSges(7,3) + t521 * mrSges(4,1) + t355 * t543 + t249 * t354 / 0.2e1 + (t524 * t562 + t507) * qJD(2) + m(7) * (t13 * t5 + t14 * t4 + t151 * t44 + t2 * t31 + t3 * t30 + t93 * t98) + m(6) * (t105 * t215 + t11 * t67 + t12 * t66 + t164 * t194 + t32 * t48 + t33 * t47) + (-t222 / 0.2e1 + t509 / 0.2e1 - t527 * pkin(7) + t259 * mrSges(4,1)) * t408 + (-qJDD(2) * mrSges(3,2) - t207) * t300 + t215 * t37 + m(4) * (pkin(7) * t521 + t147 * t256 + t209 * t226) + (-t249 * t300 + t250 * t465 + t522) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t522) + (Ifges(7,4) * t50 + Ifges(7,2) * t51) * t489 + t132 * t504 + t131 * t505 + t138 * t491 + t50 * t495 + t51 * t497 - t203 * t499 + t202 * t500 + (Ifges(7,4) * t132 + Ifges(7,2) * t131 + Ifges(7,6) * t317) * t501 + (Ifges(7,1) * t132 + Ifges(7,4) * t131 + Ifges(7,5) * t317) * t502 + (Ifges(7,5) * t132 + Ifges(7,6) * t131 + Ifges(7,3) * t317) * t478 + (Ifges(5,5) * t317 - t321 * t361) * t479 + (Ifges(5,6) * t317 - t321 * t358) * t480 + t421 * t485 + (Ifges(6,1) * t138 + Ifges(6,4) * t139) * t481 + (t321 * (-Ifges(3,2) * t317 + t455) + t339 + t327) * t399 / 0.2e1 + t105 * (-mrSges(6,1) * t202 - mrSges(6,2) * t203) + (-Ifges(6,4) * t203 + Ifges(6,2) * t202 + Ifges(6,6) * t317) * t493 + (-Ifges(6,1) * t203 + Ifges(6,4) * t202 + Ifges(6,5) * t317) * t494 + (-Ifges(6,5) * t203 + Ifges(6,6) * t202 + Ifges(6,3) * t317) * t477 + t194 * t90 + (Ifges(6,4) * t138 + Ifges(6,2) * t139) * t483 + t122 * t184 + t121 * t185 + t162 * t146 + t163 * t145 + t164 * (-mrSges(6,1) * t139 + mrSges(6,2) * t138) + t151 * t8 + (t11 * t202 + t12 * t203 - t138 * t47 + t139 * t48) * mrSges(6,3) + t139 * t76 / 0.2e1 + t44 * (-mrSges(7,1) * t131 + mrSges(7,2) * t132) + t32 * t129 + t33 * t130 + Ifges(2,3) * qJDD(1) + t98 * t43 + t93 * (-mrSges(7,1) * t51 + mrSges(7,2) * t50) + t67 * t61 + t4 * t68 + t5 * t69 + t66 * t60 + t31 * t18 + t30 * t17 + t153 * t364 * t321 + t209 * t243 - t247 * t166 - pkin(1) * (mrSges(3,1) * t249 + mrSges(3,2) * t250) + t256 * (-mrSges(4,2) * t249 - mrSges(4,3) * t250) + t265 * t123 + (t513 / 0.2e1 - t526 * pkin(7) + t254 * mrSges(4,1) + Ifges(7,6) * t489 + Ifges(7,3) * t474 + Ifges(6,5) * t481 + Ifges(6,6) * t483 + Ifges(7,5) * t487 + Ifges(6,3) * t472 - t551) * t407 - t317 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t250 + Ifges(4,6) * t249) / 0.2e1 - t317 * t558 + t321 * (Ifges(3,4) * t250 - Ifges(3,2) * t249 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t321 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t250 + Ifges(4,3) * t249) / 0.2e1 - t317 * t560 + t366 * t435 + (Ifges(7,1) * t50 + Ifges(7,4) * t51) * t487; t153 * t363 + (Ifges(7,1) * t487 + Ifges(7,4) * t489 + Ifges(7,5) * t474 + t495 + t514) * t326 - (Ifges(7,4) * t487 + Ifges(7,2) * t489 + Ifges(7,6) * t474 + t497 + t515) * t85 - (-Ifges(3,2) * t295 + t284 + t513) * t409 / 0.2e1 - t509 * t295 / 0.2e1 + t222 * t295 / 0.2e1 + (-Ifges(6,4) * t217 - Ifges(6,2) * t216) * t483 + (-Ifges(6,5) * t217 - Ifges(6,6) * t216) * t472 + (-qJ(3) * t390 * t414 + t322 * t512) * g(1) + (-t390 * t436 + t512) * t550 + (Ifges(7,4) * t374 + Ifges(7,2) * t111) * t490 + (Ifges(7,5) * t374 + Ifges(7,6) * t111) * t475 + t374 * t496 + (Ifges(7,1) * t374 + Ifges(7,4) * t111) * t488 - t93 * (-mrSges(7,1) * t111 + mrSges(7,2) * t374) + t534 * t129 + (t105 * t276 + t11 * t168 + t12 * t167 - t164 * t193 + t533 * t47 + t534 * t48) * m(6) + t533 * t130 + t524 * t377 + t526 * t287 + t527 * t285 - t529 * t76 / 0.2e1 + (mrSges(6,1) * t529 + mrSges(6,2) * t528) * t164 + t531 * t43 + (Ifges(6,1) * t188 + Ifges(6,4) * t187) * t482 + t537 * t69 + t538 * t68 + (t13 * t537 + t14 * t538 + t186 * t44 + t2 * t65 + t3 * t64 + t531 * t93) * m(7) + t539 * t249 + t540 * t250 + t547 * mrSges(6,3) - t259 * t389 - t254 * t388 + (Ifges(6,4) * t188 + Ifges(6,2) * t187) * t484 + (Ifges(6,5) * t188 + Ifges(6,6) * t187) * t473 + (-t207 + t123) * qJ(3) + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + t233 * mrSges(3,2) - t234 * mrSges(3,1) + t205 * mrSges(4,2) - pkin(2) * t208 + t520 * t365 - t353 * mrSges(5,3) + (-t404 - t142) * t185 + t159 * t505 - t217 * t491 + t111 * t498 + (Ifges(5,5) * t312 - Ifges(5,6) * t311) * t476 + (Ifges(5,1) * t312 - t454) * t479 + (-Ifges(5,2) * t311 + t453) * t480 + t311 * t485 + (-t405 - t143) * t184 + (-Ifges(6,1) * t217 - Ifges(6,4) * t216) * t481 + (-t111 * t14 - t544) * mrSges(7,3) - t192 * mrSges(4,3) - t193 * t90 + t186 * t8 - t188 * t77 / 0.2e1 + t167 * t60 + t168 * t61 + t64 * t17 + t65 * t18 + (-m(4) * t259 + m(5) * t211 + m(6) * t164 - t518) * qJD(3) + (-t507 + (-t327 / 0.2e1 - t339 / 0.2e1 + t341 + t525 / 0.2e1) * qJD(1)) * qJD(1) - t244 * t243 - t245 * t166 + (-m(4) * t411 - m(5) * t367 - t321 * mrSges(5,3) - m(7) * (t411 - t427) - m(6) * (t411 - t420) + t553 * t317 + t546) * g(3) + t276 * t37 + t312 * t102 / 0.2e1 + (Ifges(6,5) * t482 + Ifges(7,5) * t488 + Ifges(6,6) * t484 + Ifges(7,6) * t490 + Ifges(6,3) * t473 + Ifges(7,3) * t475 + t551) * t409 + (-pkin(2) * t205 - qJ(3) * t192 - t226 * t244) * m(4) - t347 * t500 - t346 * t499 + (-Ifges(6,4) * t346 - Ifges(6,2) * t347) * t493 + (-Ifges(6,1) * t346 - Ifges(6,4) * t347) * t494 + (-Ifges(6,5) * t346 - Ifges(6,6) * t347) * t477 + t105 * (mrSges(6,1) * t347 - mrSges(6,2) * t346) + (-t118 * t142 - t119 * t143 - t211 * t245 + qJ(3) * t153 - t353 * t461 + (-t118 * t312 - t119 * t311) * qJD(4)) * m(5) - t523 * t461 + (Ifges(7,1) * t530 + Ifges(7,4) * t159) * t502 + t530 * t504 + (Ifges(7,4) * t530 + Ifges(7,2) * t159) * t501 + (Ifges(7,5) * t530 + Ifges(7,6) * t159) * t478 + t44 * (-mrSges(7,1) * t159 + mrSges(7,2) * t530); (-t374 + t326) * t69 + t530 * t17 + t549 * t68 + ((t184 * t312 - t185 * t311 + t243) * qJD(1) - t520 * t390) * t317 + t529 * t129 + t390 * t462 + t523 + t528 * t130 - t159 * t18 + (-t43 + t518) * qJD(2) + t347 * t61 - t346 * t60 + t208 + (-qJD(2) * t93 + t13 * t326 + t14 * t549 + t544) * m(7) + (-qJD(2) * t164 - t547) * m(6) + (-qJD(2) * t211 - (t118 * t311 - t119 * t312) * t295 + t353) * m(5) + (qJD(2) * t259 + t226 * t295 + t205) * m(4); -t373 * t129 - t154 * t130 - t236 * t184 - t235 * t185 - t545 * t68 - t80 * t69 + t123 + t37 + t8 + (-t13 * t80 - t14 * t545 + t44) * m(7) + (-t154 * t47 - t373 * t48 + t105) * m(6) + (-t118 * t235 - t119 * t236 + t153) * m(5) - (t317 * g(3) + t321 * t520) * t519; -(-mrSges(6,1) * t292 + mrSges(6,2) * t291) * t462 + t342 - g(3) * (t255 + (-m(7) * t466 - t459) * t321) + t43 * t467 - m(7) * (t13 * t15 + t14 * t16 - t467 * t93) + (t2 * t315 + t3 * t319 + (-t13 * t315 + t14 * t319) * qJD(6)) * t503 + t76 * t481 + (Ifges(6,1) * t373 + t452) * t482 + (Ifges(6,5) * t373 + Ifges(6,6) * t154) * t473 - t164 * (-mrSges(6,1) * t154 + mrSges(6,2) * t373) - t16 * t68 - t15 * t69 - t558 + t557 + t393 + (Ifges(7,1) * t488 + Ifges(7,4) * t490 + Ifges(7,5) * t475 + t496 - t514) * t545 + (Ifges(7,4) * t488 + Ifges(7,2) * t490 + Ifges(7,6) * t475 + t498 - t515) * t80 + (Ifges(6,2) * t154 + t152 + t77) * t484 + (-t457 + t130) * t48 + (t458 - t129) * t47 + (-mrSges(6,2) * t198 + t197 * t535 - t412) * g(2) + (mrSges(6,2) * t196 + t195 * t535 - t413) * g(1) + ((-t315 * t69 + t319 * t68) * qJD(6) + t17 * t319 + t18 * t315) * pkin(5); -t93 * (-mrSges(7,1) * t80 + mrSges(7,2) * t545) + (Ifges(7,1) * t545 + t468) * t488 + t41 * t487 + (Ifges(7,5) * t545 + Ifges(7,6) * t80) * t475 - t13 * t68 + t14 * t69 - g(1) * t413 - g(2) * t412 - g(3) * (-t321 * t459 + t255) + (t13 * t545 - t14 * t80) * mrSges(7,3) + t342 + (Ifges(7,2) * t80 + t42 + t78) * t490;];
tau  = t1;
