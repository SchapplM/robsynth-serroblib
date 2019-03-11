% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:24
% EndTime: 2019-03-08 22:17:12
% DurationCPUTime: 28.65s
% Computational Cost: add. (12520->806), mult. (28739->1107), div. (0->0), fcn. (22413->18), ass. (0->356)
t313 = sin(qJ(3));
t317 = cos(qJ(3));
t346 = pkin(3) * t313 - qJ(4) * t317;
t395 = qJD(4) * t313;
t226 = qJD(3) * t346 - t395;
t306 = sin(pkin(12));
t308 = cos(pkin(12));
t314 = sin(qJ(2));
t399 = qJD(3) * t313;
t386 = pkin(8) * t399;
t307 = sin(pkin(6));
t404 = qJD(1) * t307;
t318 = cos(qJ(2));
t409 = t317 * t318;
t500 = t308 * t226 + t306 * t386 - (-t306 * t409 + t308 * t314) * t404;
t536 = t306 * t226 - (t306 * t314 + t308 * t409) * t404;
t410 = t308 * t317;
t342 = pkin(4) * t313 - pkin(9) * t410;
t535 = qJD(3) * t342 + t500;
t411 = t308 * t313;
t414 = t306 * t317;
t534 = (-pkin(8) * t411 - pkin(9) * t414) * qJD(3) + t536;
t312 = sin(qJ(5));
t316 = cos(qJ(5));
t257 = t306 * t316 + t308 * t312;
t337 = t257 * t317;
t209 = qJD(2) * t337;
t234 = t257 * qJD(5);
t526 = t209 - t234;
t343 = t306 * t312 - t308 * t316;
t336 = t343 * t317;
t210 = qJD(2) * t336;
t233 = t343 * qJD(5);
t525 = t210 - t233;
t356 = mrSges(4,1) * t317 - mrSges(4,2) * t313;
t435 = pkin(9) + qJ(4);
t489 = m(7) * (-pkin(10) - t435) - mrSges(7,3) - m(6) * t435 - mrSges(6,3) - m(5) * qJ(4) - mrSges(5,3);
t305 = pkin(12) + qJ(5);
t302 = qJ(6) + t305;
t294 = sin(t302);
t295 = cos(t302);
t296 = pkin(4) * t308 + pkin(3);
t299 = sin(t305);
t300 = cos(t305);
t355 = -t308 * mrSges(5,1) + t306 * mrSges(5,2);
t524 = m(7) * (pkin(5) * t300 + t296) + t295 * mrSges(7,1) - t294 * mrSges(7,2) + m(6) * t296 + t300 * mrSges(6,1) - t299 * mrSges(6,2) + m(5) * pkin(3) - t355;
t533 = t313 * t489 - t317 * t524 - mrSges(3,1) - t356;
t267 = -pkin(3) * t317 - qJ(4) * t313 - pkin(2);
t249 = t308 * t267;
t172 = -pkin(9) * t411 + t249 + (-pkin(8) * t306 - pkin(4)) * t317;
t206 = pkin(8) * t410 + t306 * t267;
t415 = t306 * t313;
t188 = -pkin(9) * t415 + t206;
t101 = t312 * t172 + t316 * t188;
t508 = -qJD(5) * t101 - t534 * t312 + t535 * t316;
t393 = qJD(5) * t316;
t394 = qJD(5) * t312;
t507 = t172 * t393 - t188 * t394 + t535 * t312 + t534 * t316;
t382 = t314 * t404;
t265 = qJD(2) * pkin(8) + t382;
t246 = t313 * t265;
t309 = cos(pkin(6));
t403 = qJD(1) * t309;
t201 = t317 * t403 - t246;
t260 = t346 * qJD(2);
t145 = -t201 * t306 + t308 * t260;
t120 = qJD(2) * t342 + t145;
t146 = t308 * t201 + t306 * t260;
t392 = t317 * qJD(2);
t380 = t306 * t392;
t130 = -pkin(9) * t380 + t146;
t268 = t435 * t306;
t269 = t435 * t308;
t396 = qJD(4) * t308;
t397 = qJD(4) * t306;
t506 = -t268 * t393 + (-t130 + t396) * t316 + (-qJD(5) * t269 - t120 - t397) * t312;
t190 = -t312 * t268 + t316 * t269;
t505 = -t257 * qJD(4) - qJD(5) * t190 - t316 * t120 + t130 * t312;
t531 = pkin(10) * t526 + t506;
t401 = qJD(2) * t313;
t530 = -pkin(5) * t401 - pkin(10) * t525 + t505;
t150 = -qJD(3) * t336 - t234 * t313;
t529 = pkin(5) * t399 - pkin(10) * t150 + t508;
t151 = -qJD(3) * t337 + t233 * t313;
t528 = -pkin(10) * t151 - t507;
t292 = qJD(5) - t392;
t250 = qJD(3) * t308 - t306 * t401;
t251 = qJD(3) * t306 + t308 * t401;
t358 = t316 * t250 - t251 * t312;
t169 = t250 * t312 + t251 * t316;
t428 = Ifges(6,4) * t169;
t84 = Ifges(6,2) * t358 + t292 * Ifges(6,6) + t428;
t464 = t84 / 0.2e1;
t165 = Ifges(6,4) * t358;
t85 = t169 * Ifges(6,1) + t292 * Ifges(6,5) + t165;
t463 = t85 / 0.2e1;
t311 = sin(qJ(6));
t315 = cos(qJ(6));
t520 = -t169 * t311 + t315 * t358;
t95 = t169 * t315 + t311 * t358;
t53 = -mrSges(7,1) * t520 + mrSges(7,2) * t95;
t98 = -mrSges(6,1) * t358 + mrSges(6,2) * t169;
t523 = t53 + t98;
t402 = qJD(2) * t307;
t373 = qJD(1) * t402;
t277 = t318 * t373;
t390 = qJDD(1) * t307;
t223 = t314 * t390 + t277;
t522 = qJDD(2) * pkin(8) + qJD(3) * t403 + t223;
t491 = -m(6) - m(4) - m(7) - m(5);
t521 = -pkin(2) * t491 - t533;
t354 = t306 * mrSges(5,1) + t308 * mrSges(5,2);
t519 = -t299 * mrSges(6,1) - t294 * mrSges(7,1) - t300 * mrSges(6,2) - t295 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - t354;
t518 = m(4) * pkin(8);
t391 = qJD(2) * qJD(3);
t263 = qJDD(2) * t313 + t317 * t391;
t211 = qJDD(3) * t308 - t263 * t306;
t212 = qJDD(3) * t306 + t263 * t308;
t79 = qJD(5) * t358 + t211 * t312 + t212 * t316;
t80 = -qJD(5) * t169 + t211 * t316 - t212 * t312;
t26 = qJD(6) * t520 + t311 * t80 + t315 * t79;
t474 = t26 / 0.2e1;
t27 = -qJD(6) * t95 - t311 * t79 + t315 * t80;
t473 = t27 / 0.2e1;
t466 = t79 / 0.2e1;
t465 = t80 / 0.2e1;
t453 = t211 / 0.2e1;
t452 = t212 / 0.2e1;
t262 = -t317 * qJDD(2) + t313 * t391;
t255 = qJDD(5) + t262;
t243 = qJDD(6) + t255;
t451 = t243 / 0.2e1;
t450 = t255 / 0.2e1;
t517 = -t262 / 0.2e1;
t449 = t262 / 0.2e1;
t516 = t263 / 0.2e1;
t100 = t316 * t172 - t188 * t312;
t219 = t343 * t313;
t70 = -pkin(5) * t317 + pkin(10) * t219 + t100;
t218 = t257 * t313;
t73 = -pkin(10) * t218 + t101;
t38 = t311 * t70 + t315 * t73;
t515 = -qJD(6) * t38 + t311 * t528 + t315 * t529;
t37 = -t311 * t73 + t315 * t70;
t514 = qJD(6) * t37 + t311 * t529 - t315 * t528;
t513 = qJD(3) / 0.2e1;
t189 = -t316 * t268 - t269 * t312;
t148 = -pkin(10) * t257 + t189;
t149 = -pkin(10) * t343 + t190;
t72 = t148 * t311 + t149 * t315;
t511 = -qJD(6) * t72 - t311 * t531 + t530 * t315;
t71 = t148 * t315 - t149 * t311;
t510 = qJD(6) * t71 + t530 * t311 + t315 * t531;
t475 = m(7) * pkin(5);
t504 = -t475 - mrSges(6,1);
t202 = t317 * t265 + t313 * t403;
t175 = pkin(4) * t380 + t202;
t502 = -pkin(5) * t526 - t175;
t385 = mrSges(4,3) * t401;
t501 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t250 + mrSges(5,2) * t251 + t385;
t499 = -t308 * t386 + t536;
t429 = Ifges(5,4) * t308;
t349 = -Ifges(5,2) * t306 + t429;
t430 = Ifges(5,4) * t306;
t351 = Ifges(5,1) * t308 - t430;
t496 = t250 * (Ifges(5,6) * t313 + t317 * t349) + t251 * (Ifges(5,5) * t313 + t317 * t351);
t389 = qJDD(1) * t309;
t383 = t313 * t389 + t317 * t522;
t118 = -t265 * t399 + t383;
t398 = qJD(3) * t317;
t119 = -t265 * t398 - t313 * t522 + t317 * t389;
t495 = t118 * t317 - t119 * t313;
t105 = qJDD(3) * qJ(4) + (qJD(4) - t246) * qJD(3) + t383;
t276 = t314 * t373;
t222 = t318 * t390 - t276;
t214 = -qJDD(2) * pkin(2) - t222;
t129 = pkin(3) * t262 - qJ(4) * t263 - qJD(2) * t395 + t214;
t62 = -t105 * t306 + t308 * t129;
t63 = t308 * t105 + t306 * t129;
t494 = -t306 * t62 + t308 * t63;
t10 = -t27 * mrSges(7,1) + t26 * mrSges(7,2);
t135 = -t211 * mrSges(5,1) + t212 * mrSges(5,2);
t39 = -t80 * mrSges(6,1) + t79 * mrSges(6,2);
t493 = t10 + t135 + t39;
t283 = qJD(6) + t292;
t492 = t251 * Ifges(5,5) + t169 * Ifges(6,5) + t95 * Ifges(7,5) + t250 * Ifges(5,6) + Ifges(6,6) * t358 + Ifges(7,6) * t520 - Ifges(5,3) * t392 + t292 * Ifges(6,3) + t283 * Ifges(7,3);
t297 = Ifges(4,4) * t392;
t488 = t308 * (t251 * Ifges(5,1) + t250 * Ifges(5,4) - Ifges(5,5) * t392) + Ifges(4,1) * t401 + Ifges(4,5) * qJD(3) + t297;
t487 = mrSges(4,1) + t524;
t486 = mrSges(4,2) + t489;
t196 = qJD(3) * qJ(4) + t202;
t381 = t318 * t404;
t204 = qJD(2) * t267 - t381;
t123 = -t196 * t306 + t308 * t204;
t124 = t308 * t196 + t306 * t204;
t192 = -qJD(3) * pkin(3) + qJD(4) - t201;
t266 = -qJD(2) * pkin(2) - t381;
t485 = -t192 * t317 * t354 - t124 * (-mrSges(5,2) * t313 - mrSges(5,3) * t414) - t123 * (mrSges(5,1) * t313 - mrSges(5,3) * t410) - t266 * (mrSges(4,1) * t313 + mrSges(4,2) * t317);
t88 = -pkin(4) * t392 - pkin(9) * t251 + t123;
t99 = pkin(9) * t250 + t124;
t51 = -t312 * t99 + t316 * t88;
t33 = -pkin(10) * t169 + t51;
t32 = pkin(5) * t292 + t33;
t52 = t312 * t88 + t316 * t99;
t34 = pkin(10) * t358 + t52;
t424 = t311 * t34;
t15 = t315 * t32 - t424;
t48 = pkin(4) * t262 - pkin(9) * t212 + t62;
t55 = pkin(9) * t211 + t63;
t12 = -qJD(5) * t52 - t312 * t55 + t316 * t48;
t6 = pkin(5) * t255 - pkin(10) * t79 + t12;
t11 = t312 * t48 + t316 * t55 + t88 * t393 - t394 * t99;
t7 = pkin(10) * t80 + t11;
t2 = qJD(6) * t15 + t311 * t6 + t315 * t7;
t423 = t315 * t34;
t16 = t311 * t32 + t423;
t3 = -qJD(6) * t16 - t311 * t7 + t315 * t6;
t483 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t482 = -m(5) * t192 - t501;
t481 = -t12 * mrSges(6,1) + t11 * mrSges(6,2);
t441 = pkin(4) * t306;
t264 = pkin(5) * t299 + t441;
t479 = -m(6) * (pkin(8) + t441) - m(7) * (pkin(8) + t264) - m(5) * pkin(8) - t518 + t519;
t432 = Ifges(4,4) * t313;
t350 = t317 * Ifges(4,2) + t432;
t478 = t15 * mrSges(7,1) + t51 * mrSges(6,1) - Ifges(4,6) * qJD(3) / 0.2e1 - qJD(2) * t350 / 0.2e1 - t16 * mrSges(7,2) - t52 * mrSges(6,2);
t319 = qJD(2) ^ 2;
t477 = Ifges(7,4) * t474 + Ifges(7,2) * t473 + Ifges(7,6) * t451;
t476 = Ifges(7,1) * t474 + Ifges(7,4) * t473 + Ifges(7,5) * t451;
t472 = Ifges(6,4) * t466 + Ifges(6,2) * t465 + Ifges(6,6) * t450;
t471 = Ifges(6,1) * t466 + Ifges(6,4) * t465 + Ifges(6,5) * t450;
t442 = Ifges(7,4) * t95;
t44 = Ifges(7,2) * t520 + Ifges(7,6) * t283 + t442;
t470 = -t44 / 0.2e1;
t469 = t44 / 0.2e1;
t89 = Ifges(7,4) * t520;
t45 = Ifges(7,1) * t95 + Ifges(7,5) * t283 + t89;
t468 = -t45 / 0.2e1;
t467 = t45 / 0.2e1;
t462 = -t520 / 0.2e1;
t461 = t520 / 0.2e1;
t460 = -t95 / 0.2e1;
t459 = t95 / 0.2e1;
t458 = Ifges(5,1) * t452 + Ifges(5,4) * t453 + Ifges(5,5) * t449;
t457 = -t358 / 0.2e1;
t456 = t358 / 0.2e1;
t455 = -t169 / 0.2e1;
t454 = t169 / 0.2e1;
t448 = -t283 / 0.2e1;
t447 = t283 / 0.2e1;
t446 = -t292 / 0.2e1;
t445 = t292 / 0.2e1;
t443 = mrSges(7,3) * t15;
t440 = pkin(5) * t169;
t438 = t16 * mrSges(7,3);
t303 = t313 * pkin(8);
t434 = mrSges(6,3) * t358;
t433 = mrSges(6,3) * t169;
t431 = Ifges(4,4) * t317;
t422 = cos(pkin(11));
t421 = sin(pkin(11));
t108 = -qJDD(3) * pkin(3) + qJDD(4) - t119;
t420 = t108 * t313;
t413 = t307 * t314;
t412 = t307 * t318;
t362 = t421 * t318;
t365 = t422 * t314;
t230 = t309 * t365 + t362;
t367 = t307 * t422;
t183 = t230 * t317 - t313 * t367;
t363 = t421 * t314;
t364 = t422 * t318;
t229 = -t309 * t364 + t363;
t408 = (-t183 * t294 + t229 * t295) * mrSges(7,1) + (-t183 * t295 - t229 * t294) * mrSges(7,2);
t232 = -t309 * t363 + t364;
t366 = t307 * t421;
t185 = t232 * t317 + t313 * t366;
t231 = t309 * t362 + t365;
t407 = (-t185 * t294 + t231 * t295) * mrSges(7,1) + (-t185 * t295 - t231 * t294) * mrSges(7,2);
t237 = t309 * t313 + t317 * t413;
t406 = (-t237 * t294 - t295 * t412) * mrSges(7,1) + (-t237 * t295 + t294 * t412) * mrSges(7,2);
t298 = pkin(8) * t398;
t377 = t306 * t398;
t244 = pkin(4) * t377 + t298;
t261 = pkin(4) * t415 + t303;
t400 = qJD(2) * t314;
t388 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t243;
t387 = Ifges(6,5) * t79 + Ifges(6,6) * t80 + Ifges(6,3) * t255;
t384 = mrSges(4,3) * t392;
t379 = t307 * t400;
t378 = t318 * t402;
t361 = -t391 / 0.2e1;
t360 = t391 / 0.2e1;
t348 = Ifges(4,5) * t317 - Ifges(4,6) * t313;
t347 = Ifges(5,5) * t308 - Ifges(5,6) * t306;
t180 = -t237 * t306 - t308 * t412;
t181 = t237 * t308 - t306 * t412;
t103 = t180 * t316 - t181 * t312;
t104 = t180 * t312 + t181 * t316;
t57 = t103 * t315 - t104 * t311;
t58 = t103 * t311 + t104 * t315;
t140 = -t218 * t315 + t219 * t311;
t141 = -t218 * t311 - t219 * t315;
t170 = -t257 * t311 - t315 * t343;
t171 = t257 * t315 - t311 * t343;
t341 = t388 - t483;
t236 = -t309 * t317 + t313 * t413;
t338 = t313 * (Ifges(4,1) * t317 - t432);
t182 = t230 * t313 + t317 * t367;
t184 = t232 * t313 - t317 * t366;
t333 = -g(1) * t184 - g(2) * t182 - g(3) * t236;
t147 = -pkin(4) * t250 + t192;
t323 = t317 * (Ifges(5,3) * t313 + t317 * t347);
t81 = -pkin(4) * t211 + t108;
t271 = -qJD(3) * mrSges(4,2) + t384;
t259 = t356 * qJD(2);
t225 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t263;
t224 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t262;
t213 = pkin(5) * t343 - t296;
t208 = -mrSges(5,1) * t392 - mrSges(5,3) * t251;
t207 = mrSges(5,2) * t392 + mrSges(5,3) * t250;
t205 = -pkin(8) * t414 + t249;
t191 = mrSges(4,1) * t262 + mrSges(4,2) * t263;
t187 = -qJD(3) * t236 + t317 * t378;
t186 = qJD(3) * t237 + t313 * t378;
t174 = pkin(5) * t218 + t261;
t161 = t251 * Ifges(5,4) + t250 * Ifges(5,2) - Ifges(5,6) * t392;
t159 = mrSges(5,1) * t262 - mrSges(5,3) * t212;
t158 = -mrSges(5,2) * t262 + mrSges(5,3) * t211;
t144 = t187 * t308 + t306 * t379;
t143 = -t187 * t306 + t308 * t379;
t139 = mrSges(6,1) * t292 - t433;
t138 = -mrSges(6,2) * t292 + t434;
t132 = -t209 * t311 - t210 * t315;
t131 = -t209 * t315 + t210 * t311;
t117 = -pkin(5) * t151 + t244;
t115 = t212 * Ifges(5,4) + t211 * Ifges(5,2) + t262 * Ifges(5,6);
t97 = -qJD(6) * t171 + t233 * t311 - t234 * t315;
t96 = qJD(6) * t170 - t233 * t315 - t234 * t311;
t87 = -pkin(5) * t358 + t147;
t75 = mrSges(7,1) * t283 - mrSges(7,3) * t95;
t74 = -mrSges(7,2) * t283 + mrSges(7,3) * t520;
t69 = -mrSges(6,2) * t255 + mrSges(6,3) * t80;
t68 = mrSges(6,1) * t255 - mrSges(6,3) * t79;
t61 = -qJD(6) * t141 - t150 * t311 + t151 * t315;
t60 = qJD(6) * t140 + t150 * t315 + t151 * t311;
t47 = -qJD(5) * t104 + t143 * t316 - t144 * t312;
t46 = qJD(5) * t103 + t143 * t312 + t144 * t316;
t40 = -pkin(5) * t80 + t81;
t22 = -mrSges(7,2) * t243 + mrSges(7,3) * t27;
t21 = mrSges(7,1) * t243 - mrSges(7,3) * t26;
t18 = t315 * t33 - t424;
t17 = -t311 * t33 - t423;
t14 = -qJD(6) * t58 - t311 * t46 + t315 * t47;
t13 = qJD(6) * t57 + t311 * t47 + t315 * t46;
t1 = [m(2) * qJDD(1) + t103 * t68 + t104 * t69 + t13 * t74 + t46 * t138 + t47 * t139 + t14 * t75 + t143 * t208 + t144 * t207 + t181 * t158 + t180 * t159 + t187 * t271 + t57 * t21 + t58 * t22 + t237 * t224 + (-t225 + t493) * t236 + (t501 + t523) * t186 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t319 - t191) * t318 + (-mrSges(3,1) * t319 - mrSges(3,2) * qJDD(2) - qJD(2) * t259) * t314) * t307 + (-m(2) - m(3) + t491) * g(3) + m(3) * (qJDD(1) * t309 ^ 2 + (t222 * t318 + t223 * t314) * t307) + m(4) * (t118 * t237 - t119 * t236 - t186 * t201 + t187 * t202 + (-t214 * t318 + t266 * t400) * t307) + m(7) * (t13 * t16 + t14 * t15 + t186 * t87 + t2 * t58 + t236 * t40 + t3 * t57) + m(6) * (t103 * t12 + t104 * t11 + t147 * t186 + t236 * t81 + t46 * t52 + t47 * t51) + m(5) * (t108 * t236 + t123 * t143 + t124 * t144 + t180 * t62 + t181 * t63 + t186 * t192); -(Ifges(5,5) * t212 + Ifges(5,6) * t211 + Ifges(5,3) * t262 + t387 + t388) * t317 / 0.2e1 + (t229 * t521 + t230 * t479) * g(2) + (t231 * t521 + t232 * t479) * g(1) + (-t202 * mrSges(4,3) + t492 / 0.2e1 + Ifges(6,5) * t454 + Ifges(7,5) * t459 + Ifges(6,6) * t456 + Ifges(7,6) * t461 + Ifges(6,3) * t445 + Ifges(7,3) * t447 + t478) * t399 + (-t411 * t62 - t415 * t63) * mrSges(5,3) + (t222 + t276) * mrSges(3,1) + (t140 * t2 - t141 * t3 - t15 * t60 + t16 * t61) * mrSges(7,3) + (-m(6) * t147 - m(7) * t87 + t482 - t523) * t313 * t381 + t431 * t516 + t350 * t517 + (-t11 * t218 + t12 * t219 - t150 * t51 + t151 * t52) * mrSges(6,3) + (-Ifges(6,4) * t219 - Ifges(6,2) * t218) * t465 + (-Ifges(6,5) * t219 - Ifges(6,6) * t218) * t450 + t81 * (mrSges(6,1) * t218 - mrSges(6,2) * t219) + (-Ifges(6,1) * t219 - Ifges(6,4) * t218) * t466 + t501 * t298 + (-pkin(2) * t214 + t495 * pkin(8) - (t266 * t314 + (-t201 * t313 + t202 * t317) * t318) * t404) * m(4) + t499 * t207 + t500 * t208 + (t205 * t62 + t206 * t63 + (t192 * t398 + t420) * pkin(8) + t499 * t124 + t500 * t123) * m(5) + t488 * t398 / 0.2e1 + (Ifges(6,4) * t150 + Ifges(6,2) * t151) * t456 + t338 * t360 + (-t201 * t398 + t495) * mrSges(4,3) + t151 * t464 + t60 * t467 + t61 * t469 - t219 * t471 - t218 * t472 + t141 * t476 + t140 * t477 + t411 * t458 + t150 * t463 + (-t223 + t277) * mrSges(3,2) + t354 * t420 + (t135 - t225) * t303 + (Ifges(6,5) * t150 + Ifges(6,6) * t151) * t445 + (t491 * (pkin(2) * t412 + pkin(8) * t413) + (t533 * t318 + (-m(6) * t441 - m(7) * t264 + t519) * t314) * t307) * g(3) - t214 * t356 + t261 * t39 + t244 * t98 + Ifges(3,3) * qJDD(2) + t205 * t159 + t206 * t158 - pkin(2) * t191 + t174 * t10 + t147 * (-mrSges(6,1) * t151 + mrSges(6,2) * t150) + t40 * (-mrSges(7,1) * t140 + mrSges(7,2) * t141) + t117 * t53 + t100 * t68 + t101 * t69 + t87 * (-mrSges(7,1) * t61 + mrSges(7,2) * t60) + t38 * t22 + t37 * t21 - t115 * t415 / 0.2e1 + t323 * t361 + qJDD(3) * (Ifges(4,5) * t313 + Ifges(4,6) * t317) + t507 * t138 + t508 * t139 + (t100 * t12 + t101 * t11 + t147 * t244 + t261 * t81 + t507 * t52 + t508 * t51) * m(6) + t496 * t513 + t514 * t74 + t515 * t75 + (t117 * t87 + t15 * t515 + t16 * t514 + t174 * t40 + t2 * t38 + t3 * t37) * m(7) + ((-Ifges(4,2) * t313 + t431) * t360 - t271 * t381 - Ifges(6,6) * t465 - Ifges(6,5) * t466 - Ifges(7,6) * t473 - Ifges(7,5) * t474 - Ifges(5,3) * t449 - Ifges(6,3) * t450 - Ifges(7,3) * t451 - Ifges(5,5) * t452 - Ifges(5,6) * t453 + Ifges(4,4) * t516 + Ifges(4,2) * t517 + pkin(8) * t224 + t63 * mrSges(5,2) - t62 * mrSges(5,1) + t481 + t483) * t317 + (Ifges(4,1) * t263 + Ifges(4,4) * t517 + t347 * t449 + t349 * t453 + t351 * t452) * t313 + (Ifges(6,1) * t150 + Ifges(6,4) * t151) * t454 + ((-t201 * t317 - t202 * t313) * t518 + t348 * t513 - t485) * qJD(3) + (Ifges(7,4) * t141 + Ifges(7,2) * t140) * t473 + (Ifges(7,4) * t60 + Ifges(7,2) * t61) * t461 + (Ifges(7,5) * t141 + Ifges(7,6) * t140) * t451 + (Ifges(7,5) * t60 + Ifges(7,6) * t61) * t447 - t161 * t377 / 0.2e1 + (Ifges(7,1) * t141 + Ifges(7,4) * t140) * t474 + (Ifges(7,1) * t60 + Ifges(7,4) * t61) * t459 + t259 * t382 - t271 * t386; t525 * t463 + t526 * t464 + (Ifges(7,1) * t459 + Ifges(7,4) * t461 + Ifges(7,5) * t447 - t443 + t467) * t96 + (Ifges(7,5) * t132 + Ifges(7,6) * t131) * t448 + (Ifges(7,4) * t132 + Ifges(7,2) * t131) * t462 + (t323 / 0.2e1 - t338 / 0.2e1) * t319 + t81 * (mrSges(6,1) * t343 + mrSges(6,2) * t257) - t343 * t472 + (Ifges(6,4) * t257 - Ifges(6,2) * t343) * t465 + (Ifges(6,1) * t257 - Ifges(6,4) * t343) * t466 + (Ifges(6,5) * t257 - Ifges(6,6) * t343) * t450 + (-Ifges(6,5) * t210 - Ifges(6,6) * t209) * t446 + (-Ifges(6,4) * t210 - Ifges(6,2) * t209) * t457 + (-Ifges(6,1) * t210 - Ifges(6,4) * t209) * t455 + (-t145 - t397) * t208 + (t158 * t308 - t159 * t306) * qJ(4) + (-t496 / 0.2e1 + t485) * qJD(2) + (Ifges(7,4) * t459 + Ifges(7,2) * t461 + Ifges(7,6) * t447 + t438 + t469) * t97 + ((-t132 + t96) * mrSges(7,2) + (t131 - t97) * mrSges(7,1)) * t87 + (-t146 + t396) * t207 + t505 * t139 + t502 * t53 + (-t123 * t145 - t124 * t146 - pkin(3) * t108 + (-t123 * t306 + t124 * t308) * qJD(4) + t494 * qJ(4)) * m(5) + t494 * mrSges(5,3) - t492 * t401 / 0.2e1 + (t182 * t487 + t183 * t486) * g(2) + (t236 * t487 + t237 * t486) * g(3) + (t184 * t487 + t185 * t486) * g(1) - (-Ifges(4,2) * t401 + t297 + t488) * t392 / 0.2e1 + (t385 + t482) * t202 + (-t11 * t343 - t12 * t257 - t51 * t525 + t52 * t526) * mrSges(6,3) + (-mrSges(6,1) * t526 + mrSges(6,2) * t525) * t147 + (Ifges(7,1) * t132 + Ifges(7,4) * t131) * t460 + t308 * t115 / 0.2e1 + t132 * t468 + t131 * t470 + t257 * t471 + (Ifges(7,4) * t171 + Ifges(7,2) * t170) * t473 + (Ifges(7,1) * t171 + Ifges(7,4) * t170) * t474 + t171 * t476 + t170 * t477 + (Ifges(5,5) * t306 + Ifges(5,6) * t308) * t449 + (Ifges(7,5) * t171 + Ifges(7,6) * t170) * t451 + (Ifges(5,1) * t306 + t429) * t452 + (Ifges(5,2) * t308 + t430) * t453 + t306 * t458 + (Ifges(6,5) * t455 + Ifges(7,5) * t460 + Ifges(6,6) * t457 + Ifges(7,6) * t462 + Ifges(6,3) * t446 + Ifges(7,3) * t448 - t478) * t401 + (t384 - t271) * t201 + (-Ifges(6,1) * t233 - Ifges(6,4) * t234) * t454 + (-Ifges(6,4) * t233 - Ifges(6,2) * t234) * t456 + (-Ifges(6,5) * t233 - Ifges(6,6) * t234) * t445 - t296 * t39 + Ifges(4,5) * t263 - Ifges(4,6) * t262 + Ifges(4,3) * qJDD(3) + t213 * t10 + t189 * t68 + t190 * t69 - t175 * t98 + t40 * (-mrSges(7,1) * t170 + mrSges(7,2) * t171) + (-t131 * t16 + t132 * t15 + t170 * t2 - t171 * t3) * mrSges(7,3) - pkin(3) * t135 - t118 * mrSges(4,2) + t119 * mrSges(4,1) + t71 * t21 + t72 * t22 + t348 * t361 + t506 * t138 + (t11 * t190 + t12 * t189 - t147 * t175 - t296 * t81 + t505 * t51 + t506 * t52) * m(6) + t510 * t74 + t511 * t75 + (t15 * t511 + t16 * t510 + t2 * t72 + t213 * t40 + t3 * t71 + t502 * t87) * m(7) + t108 * t355 + t161 * t380 / 0.2e1; -t358 * t138 + t169 * t139 - t250 * t207 + t251 * t208 - t520 * t74 + t95 * t75 + (t15 * t95 - t16 * t520 + t333 + t40) * m(7) + (t169 * t51 - t358 * t52 + t333 + t81) * m(6) + (t123 * t251 - t124 * t250 + t108 + t333) * m(5) + t493; (-t406 - (-t237 * t300 + t299 * t412) * mrSges(6,2) + t504 * (-t237 * t299 - t300 * t412)) * g(3) + (-(-t185 * t300 - t231 * t299) * mrSges(6,2) - t407 + t504 * (-t185 * t299 + t231 * t300)) * g(1) + (-(-t183 * t300 - t229 * t299) * mrSges(6,2) - t408 + t504 * (-t183 * t299 + t229 * t300)) * g(2) + (Ifges(6,1) * t358 - t428) * t455 + (Ifges(6,5) * t358 - Ifges(6,6) * t169) * t446 - t147 * (mrSges(6,1) * t169 + mrSges(6,2) * t358) - t481 + (t433 + t139) * t52 + (t434 - t138) * t51 + t387 - (mrSges(7,1) * t87 + Ifges(7,4) * t460 + Ifges(7,2) * t462 + Ifges(7,6) * t448 - t438 + t470) * t95 + (-mrSges(7,2) * t87 + Ifges(7,1) * t460 + Ifges(7,4) * t462 + Ifges(7,5) * t448 + t443 + t468) * t520 + (t2 * t311 + t3 * t315 + (-t15 * t311 + t16 * t315) * qJD(6)) * t475 + t84 * t454 + (-Ifges(6,2) * t169 + t165 + t85) * t457 + t341 - t18 * t74 - t17 * t75 - t53 * t440 - m(7) * (t15 * t17 + t16 * t18 + t440 * t87) + (t315 * t21 + t311 * t22 + (-t311 * t75 + t315 * t74) * qJD(6)) * pkin(5); -t87 * (mrSges(7,1) * t95 + mrSges(7,2) * t520) + t44 * t459 + (Ifges(7,5) * t520 - Ifges(7,6) * t95) * t448 + (Ifges(7,1) * t520 - t442) * t460 - t15 * t74 + t16 * t75 - g(1) * t407 - g(2) * t408 - g(3) * t406 + (t15 * t520 + t16 * t95) * mrSges(7,3) + t341 + (-Ifges(7,2) * t95 + t45 + t89) * t462;];
tau  = t1;
