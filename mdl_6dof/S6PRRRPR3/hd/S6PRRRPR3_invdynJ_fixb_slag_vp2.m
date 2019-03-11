% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:11:03
% EndTime: 2019-03-08 23:11:36
% DurationCPUTime: 17.88s
% Computational Cost: add. (8203->699), mult. (18195->910), div. (0->0), fcn. (13529->14), ass. (0->349)
t293 = sin(qJ(3));
t298 = -pkin(9) - pkin(8);
t255 = t298 * t293;
t296 = cos(qJ(3));
t256 = t298 * t296;
t292 = sin(qJ(4));
t464 = cos(qJ(4));
t180 = t292 * t255 - t256 * t464;
t407 = t292 * t296;
t238 = t293 * t464 + t407;
t368 = qJD(3) * t298;
t243 = t293 * t368;
t367 = t464 * t296;
t341 = qJD(3) * t367;
t297 = cos(qJ(2));
t289 = sin(pkin(6));
t397 = qJD(1) * t289;
t365 = t297 * t397;
t534 = qJD(4) * t180 - t238 * t365 + t292 * t243 - t298 * t341;
t533 = mrSges(5,1) - mrSges(6,2);
t532 = mrSges(5,2) - mrSges(6,3);
t286 = qJD(3) + qJD(4);
t358 = qJD(4) * t464;
t408 = t292 * t293;
t169 = t286 * t408 - t296 * t358 - t341;
t531 = -t169 * pkin(5) + t534;
t311 = t238 * qJD(4);
t170 = qJD(3) * t238 + t311;
t391 = qJD(3) * t293;
t279 = pkin(3) * t391;
t320 = qJ(5) * t169 - qJD(5) * t238 + t279;
t294 = sin(qJ(2));
t366 = t294 * t397;
t479 = pkin(4) + pkin(10);
t530 = t170 * t479 + t320 - t366;
t290 = cos(pkin(6));
t396 = qJD(1) * t290;
t267 = t296 * t396;
t246 = qJD(2) * pkin(8) + t366;
t350 = pkin(9) * qJD(2) + t246;
t171 = -t293 * t350 + t267;
t364 = t293 * t396;
t172 = t296 * t350 + t364;
t409 = t292 * t172;
t79 = t171 * t464 - t409;
t529 = pkin(3) * t358 - t79;
t291 = sin(qJ(6));
t295 = cos(qJ(6));
t520 = mrSges(7,1) * t291 + mrSges(7,2) * t295;
t528 = mrSges(7,3) + t533;
t527 = -t520 + t532;
t526 = -Ifges(5,5) + Ifges(6,4);
t525 = -Ifges(5,6) + Ifges(6,5);
t395 = qJD(2) * t293;
t227 = -qJD(2) * t367 + t292 * t395;
t178 = t227 * t291 + t286 * t295;
t215 = Ifges(5,4) * t227;
t228 = t238 * qJD(2);
t218 = qJD(6) + t228;
t509 = t218 * Ifges(7,3);
t177 = t227 * t295 - t286 * t291;
t510 = t177 * Ifges(7,6);
t524 = t228 * Ifges(5,1) + t286 * Ifges(5,5) + t178 * Ifges(7,5) - t215 + t509 + t510;
t507 = -qJD(5) - t529;
t386 = qJD(2) * qJD(3);
t244 = qJDD(2) * t293 + t296 * t386;
t317 = -qJDD(2) * t296 + t293 * t386;
t110 = qJD(4) * t227 - t464 * t244 + t292 * t317;
t108 = qJDD(6) - t110;
t111 = qJD(2) * t311 + t292 * t244 + t317 * t464;
t285 = qJDD(3) + qJDD(4);
t52 = qJD(6) * t177 + t111 * t291 + t285 * t295;
t30 = mrSges(7,1) * t108 - mrSges(7,3) * t52;
t357 = qJD(2) * t397;
t257 = t294 * t357;
t385 = qJDD(1) * t289;
t210 = t297 * t385 - t257;
t202 = -qJDD(2) * pkin(2) - t210;
t157 = pkin(3) * t317 + t202;
t302 = t110 * qJ(5) - t228 * qJD(5) + t157;
t18 = t111 * t479 + t302;
t456 = t228 * pkin(5);
t167 = qJD(3) * pkin(3) + t171;
t76 = -t464 * t167 + t409;
t324 = t76 + t456;
t519 = qJD(5) + t324;
t40 = -t286 * t479 + t519;
t276 = pkin(3) * t296 + pkin(2);
t216 = -qJD(2) * t276 - t365;
t308 = -qJ(5) * t228 + t216;
t75 = t227 * t479 + t308;
t24 = t291 * t40 + t295 * t75;
t431 = qJD(6) * t24;
t389 = qJD(4) * t292;
t185 = t246 * t296 + t364;
t258 = t297 * t357;
t211 = t294 * t385 + t258;
t203 = qJDD(2) * pkin(8) + t211;
t384 = qJDD(1) * t290;
t84 = -qJD(3) * t185 - t203 * t293 + t296 * t384;
t71 = qJDD(3) * pkin(3) - pkin(9) * t244 + t84;
t83 = qJD(3) * t267 + t296 * t203 - t246 * t391 + t293 * t384;
t74 = -pkin(9) * t317 + t83;
t17 = -t167 * t389 - t172 * t358 - t292 * t74 + t464 * t71;
t313 = qJDD(5) - t17;
t6 = -t110 * pkin(5) - t285 * t479 + t313;
t2 = -t18 * t291 + t295 * t6 - t431;
t23 = -t291 * t75 + t295 * t40;
t329 = t23 * t291 - t24 * t295;
t1 = qJD(6) * t23 + t18 * t295 + t291 * t6;
t459 = t1 * t291;
t306 = m(7) * (-qJD(6) * t329 + t2 * t295 + t459);
t53 = -qJD(6) * t178 + t111 * t295 - t285 * t291;
t31 = -mrSges(7,2) * t108 + mrSges(7,3) * t53;
t523 = t291 * t31 + t295 * t30 + t306;
t336 = mrSges(7,1) * t295 - mrSges(7,2) * t291;
t457 = t227 * pkin(5);
t166 = t464 * t172;
t77 = t292 * t167 + t166;
t73 = -t286 * qJ(5) - t77;
t49 = -t73 - t457;
t441 = t178 * Ifges(7,4);
t67 = t177 * Ifges(7,2) + t218 * Ifges(7,6) + t441;
t522 = -t295 * t67 / 0.2e1 + t336 * t49;
t415 = t289 * t294;
t223 = t290 * t296 - t293 * t415;
t433 = cos(pkin(11));
t352 = t433 * t297;
t288 = sin(pkin(11));
t417 = t288 * t294;
t222 = -t290 * t417 + t352;
t414 = t289 * t296;
t521 = -t222 * t293 + t288 * t414;
t103 = pkin(4) * t227 + t308;
t518 = -t216 * mrSges(5,1) + t103 * mrSges(6,2);
t517 = t23 * mrSges(7,1) + t216 * mrSges(5,2) - t24 * mrSges(7,2) - t103 * mrSges(6,3);
t481 = t52 / 0.2e1;
t480 = t53 / 0.2e1;
t516 = -m(7) - m(6);
t478 = t108 / 0.2e1;
t467 = t285 / 0.2e1;
t515 = -t317 / 0.2e1;
t514 = g(3) * t289;
t19 = -mrSges(7,1) * t53 + mrSges(7,2) * t52;
t93 = mrSges(6,1) * t111 - mrSges(6,3) * t285;
t513 = t19 - t93;
t237 = -t367 + t408;
t325 = -qJ(5) * t238 - t276;
t109 = t237 * t479 + t325;
t179 = -t464 * t255 - t256 * t292;
t126 = pkin(5) * t238 + t179;
t41 = -t109 * t291 + t126 * t295;
t512 = qJD(6) * t41 + t291 * t531 + t295 * t530;
t42 = t109 * t295 + t126 * t291;
t511 = -qJD(6) * t42 - t291 * t530 + t295 * t531;
t506 = t456 - t507;
t505 = -qJD(5) - t76;
t448 = mrSges(5,3) * t227;
t186 = -mrSges(5,2) * t286 - t448;
t453 = mrSges(6,1) * t227;
t188 = -mrSges(6,3) * t286 + t453;
t504 = t186 - t188;
t447 = mrSges(5,3) * t228;
t452 = mrSges(6,1) * t228;
t503 = t286 * t533 - t447 - t452;
t224 = t290 * t293 + t294 * t414;
t128 = -t223 * t464 + t224 * t292;
t410 = t291 * t297;
t260 = t289 * t410;
t99 = t128 * t295 + t260;
t287 = qJ(3) + qJ(4);
t283 = sin(t287);
t284 = cos(t287);
t502 = t532 * t283 - t284 * t533;
t501 = -t293 * t84 + t296 * t83;
t499 = mrSges(3,1) * t297 - mrSges(3,2) * t294;
t498 = m(5) - t516;
t497 = 0.2e1 * t467;
t151 = mrSges(5,1) * t227 + mrSges(5,2) * t228;
t152 = -mrSges(6,2) * t227 - mrSges(6,3) * t228;
t338 = -t296 * mrSges(4,1) + mrSges(4,2) * t293;
t496 = t338 * qJD(2) + t151 + t152;
t205 = t283 * t415 - t290 * t284;
t206 = t283 * t290 + t284 * t415;
t495 = t205 * t528 + t206 * t527;
t418 = t288 * t289;
t160 = t222 * t283 - t284 * t418;
t161 = t222 * t284 + t283 * t418;
t494 = t160 * t528 + t161 * t527;
t353 = t433 * t294;
t416 = t288 * t297;
t220 = t290 * t353 + t416;
t354 = t289 * t433;
t158 = t220 * t283 + t284 * t354;
t159 = t220 * t284 - t283 * t354;
t493 = t158 * t528 + t159 * t527;
t72 = -pkin(4) * t286 - t505;
t492 = m(6) * t72 - t503;
t491 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t16 = t167 * t358 - t172 * t389 + t292 * t71 + t464 * t74;
t10 = -qJ(5) * t285 - qJD(5) * t286 - t16;
t92 = -mrSges(5,2) * t285 - mrSges(5,3) * t111;
t490 = m(5) * t16 - m(6) * t10 + t92 - t93;
t14 = -t285 * pkin(4) + t313;
t91 = mrSges(5,1) * t285 + mrSges(5,3) * t110;
t94 = -t110 * mrSges(6,1) + t285 * mrSges(6,2);
t489 = -m(5) * t17 + m(6) * t14 - t91 + t94;
t488 = -m(5) * t77 + m(6) * t73 - t504;
t487 = m(5) * t76 + t492;
t486 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(5,3) - t336;
t316 = m(4) * pkin(2) - t338;
t485 = t283 * t520 + mrSges(3,1) + t316 - t502;
t90 = -mrSges(7,1) * t177 + mrSges(7,2) * t178;
t484 = -m(7) * t49 + t488 - t90;
t379 = m(4) * pkin(8) + mrSges(4,3);
t483 = mrSges(3,2) - t379 + t486;
t300 = qJD(2) ^ 2;
t482 = Ifges(7,1) * t481 + Ifges(7,4) * t480 + Ifges(7,5) * t478;
t476 = -t177 / 0.2e1;
t475 = -t178 / 0.2e1;
t474 = t178 / 0.2e1;
t473 = -t218 / 0.2e1;
t472 = -t227 / 0.2e1;
t471 = t227 / 0.2e1;
t470 = -t228 / 0.2e1;
t469 = t228 / 0.2e1;
t466 = -t286 / 0.2e1;
t465 = t286 / 0.2e1;
t463 = pkin(3) * t292;
t462 = pkin(10) * t160;
t461 = pkin(10) * t205;
t458 = t158 * pkin(10);
t446 = mrSges(7,3) * t295;
t445 = Ifges(4,4) * t293;
t444 = Ifges(4,4) * t296;
t443 = Ifges(7,4) * t291;
t442 = Ifges(7,4) * t295;
t440 = t228 * Ifges(5,4);
t439 = t228 * Ifges(6,6);
t434 = t188 - t90;
t432 = qJ(5) * t283;
t429 = t170 * t291;
t428 = t170 * t295;
t219 = -t290 * t352 + t417;
t427 = t219 * t284;
t221 = t290 * t416 + t353;
t426 = t221 * t284;
t424 = t228 * t291;
t423 = t237 * t291;
t422 = t237 * t295;
t419 = t284 * t297;
t413 = t289 * t297;
t406 = t295 * t297;
t401 = -t219 * t276 - t220 * t298;
t400 = -t221 * t276 - t222 * t298;
t394 = qJD(2) * t294;
t393 = qJD(2) * t296;
t392 = qJD(2) * t297;
t390 = qJD(3) * t296;
t388 = qJD(6) * t291;
t387 = qJD(6) * t295;
t382 = Ifges(7,5) * t52 + Ifges(7,6) * t53 + Ifges(7,3) * t108;
t381 = t464 * pkin(3);
t278 = pkin(3) * t395;
t378 = mrSges(4,3) * t395;
t377 = mrSges(4,3) * t393;
t374 = t289 * t406;
t363 = t289 * t394;
t361 = t293 * t392;
t360 = t296 * t392;
t355 = -t388 / 0.2e1;
t349 = -t158 * pkin(4) + t159 * qJ(5);
t348 = -t160 * pkin(4) + qJ(5) * t161;
t347 = -t205 * pkin(4) + qJ(5) * t206;
t78 = t171 * t292 + t166;
t345 = -pkin(4) * t427 - t219 * t432 + t401;
t344 = -pkin(4) * t426 - t221 * t432 + t400;
t343 = t293 * t365;
t342 = t296 * t365;
t275 = -t381 - pkin(4);
t339 = t521 * pkin(3);
t334 = Ifges(7,1) * t291 + t442;
t333 = t296 * Ifges(4,2) + t445;
t332 = Ifges(7,2) * t295 + t443;
t331 = Ifges(4,5) * t296 - Ifges(4,6) * t293;
t330 = Ifges(7,5) * t291 + Ifges(7,6) * t295;
t149 = pkin(4) * t228 + qJ(5) * t227;
t112 = -mrSges(7,2) * t218 + mrSges(7,3) * t177;
t113 = mrSges(7,1) * t218 - mrSges(7,3) * t178;
t328 = t112 * t295 - t113 * t291;
t184 = -t246 * t293 + t267;
t327 = -t184 * t293 + t185 * t296;
t326 = t223 * pkin(3);
t116 = t149 + t278;
t323 = -t128 * t291 + t374;
t129 = t292 * t223 + t224 * t464;
t322 = t237 * t387 + t429;
t321 = t237 * t388 - t428;
t247 = -qJD(2) * pkin(2) - t365;
t319 = t247 * (mrSges(4,1) * t293 + mrSges(4,2) * t296);
t318 = t293 * (Ifges(4,1) * t296 - t445);
t81 = -t464 * t243 - t255 * t358 - t256 * t389 - t368 * t407;
t315 = -g(1) * t160 - g(2) * t158 - g(3) * t205;
t314 = -t220 * t293 - t296 * t354;
t312 = t339 + t348;
t310 = t314 * pkin(3);
t309 = t326 + t347;
t305 = t310 + t349;
t304 = qJD(3) * t224 + t289 * t361;
t303 = qJD(3) * t223 + t289 * t360;
t11 = t52 * Ifges(7,4) + t53 * Ifges(7,2) + t108 * Ifges(7,6);
t121 = t286 * Ifges(6,5) + t227 * Ifges(6,3) - t439;
t214 = Ifges(6,6) * t227;
t122 = t286 * Ifges(6,4) - t228 * Ifges(6,2) + t214;
t123 = -t227 * Ifges(5,2) + t286 * Ifges(5,6) + t440;
t173 = Ifges(7,4) * t177;
t68 = t178 * Ifges(7,1) + t218 * Ifges(7,5) + t173;
t7 = -pkin(5) * t111 - t10;
t301 = (t214 + t122) * t472 + (-Ifges(5,1) * t470 - Ifges(7,5) * t475 + Ifges(6,2) * t469 - Ifges(7,6) * t476 - Ifges(7,3) * t473 + t466 * t526 + t517) * t227 + t526 * t110 + t7 * t520 + t522 * qJD(6) + (-t424 / 0.2e1 + t355) * t68 - (t177 * t332 + t178 * t334 + t218 * t330) * qJD(6) / 0.2e1 + (t388 + t424) * t23 * mrSges(7,3) + (-Ifges(5,2) * t228 - t215 + t524) * t471 + (Ifges(6,3) * t472 - t24 * t446 + t330 * t473 + t332 * t476 + t334 * t475 + t466 * t525 + t518 + t522) * t228 + t525 * t111 + t14 * mrSges(6,2) - t16 * mrSges(5,2) + t17 * mrSges(5,1) - t10 * mrSges(6,3) + (-t440 + t121) * t470 + (t439 + t123) * t469 + t76 * t448 + t72 * t453 + (Ifges(6,1) + Ifges(5,3)) * t285 - t291 * t11 / 0.2e1 - t73 * t452 + (Ifges(7,5) * t295 - Ifges(7,6) * t291) * t478 + (-Ifges(7,2) * t291 + t442) * t480 + (Ifges(7,1) * t295 - t443) * t481 + t295 * t482;
t277 = Ifges(4,4) * t393;
t273 = qJ(5) + t463;
t252 = -qJD(3) * mrSges(4,2) + t377;
t251 = qJD(3) * mrSges(4,1) - t378;
t242 = t276 * t413;
t226 = Ifges(4,1) * t395 + Ifges(4,5) * qJD(3) + t277;
t225 = Ifges(4,6) * qJD(3) + qJD(2) * t333;
t217 = t228 * pkin(10);
t213 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t244;
t212 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t317;
t174 = mrSges(4,1) * t317 + t244 * mrSges(4,2);
t156 = pkin(4) * t237 + t325;
t127 = -t237 * pkin(5) + t180;
t101 = t149 + t217;
t89 = t116 + t217;
t60 = t78 - t457;
t59 = t77 - t457;
t57 = pkin(4) * t170 + t320;
t47 = -pkin(5) * t170 - t81;
t37 = mrSges(5,1) * t111 - mrSges(5,2) * t110;
t36 = -mrSges(6,2) * t111 + mrSges(6,3) * t110;
t34 = qJD(4) * t129 + t292 * t303 + t304 * t464;
t29 = t101 * t295 + t291 * t59;
t28 = -t101 * t291 + t295 * t59;
t27 = t291 * t60 + t295 * t89;
t26 = -t291 * t89 + t295 * t60;
t25 = t111 * pkin(4) + t302;
t22 = qJD(6) * t99 + t291 * t34 + t295 * t363;
t21 = qJD(6) * t323 - t291 * t363 + t295 * t34;
t3 = [m(4) * (t327 * t290 * qJD(3) + t84 * t223 + t83 * t224) + t303 * t252 - t304 * t251 + t223 * t213 + t224 * t212 + m(7) * (-t1 * t323 + t2 * t99 + t21 * t23 + t22 * t24) + t99 * t30 - t323 * t31 + t22 * t112 + t21 * t113 + (m(3) * t290 ^ 2 + m(2)) * qJDD(1) + t487 * t34 + t489 * t128 + (-t174 - t37 - t36) * t413 + t496 * t363 + t484 * (-t223 * t358 + t224 * t389 + t292 * t304 - t303 * t464) + (m(7) * t7 + t19 + t490) * t129 + (-m(2) - m(3) - m(4) - t498) * g(3) + ((-mrSges(3,1) * t294 - mrSges(3,2) * t297) * t300 + t499 * qJDD(2) + m(4) * (t185 * (-t294 * t391 + t360) + t184 * (-t294 * t390 - t361) - t202 * t297 + t247 * t394) + m(5) * (-t157 * t297 + t216 * t394) + m(6) * (t103 * t394 - t25 * t297) + m(3) * (t210 * t297 + t211 * t294)) * t289; (t258 - t211) * mrSges(3,2) + (-m(5) * t242 + t516 * (t242 + (pkin(4) * t284 + t432) * t413) + (-m(7) * pkin(10) * t419 + t502 * t297 + (-mrSges(7,1) * t410 - mrSges(7,2) * t406) * t283 + (t298 * t498 + t486) * t294) * t289) * g(3) + (t318 + t296 * (-Ifges(4,2) * t293 + t444)) * t386 / 0.2e1 + (qJD(6) * t68 + t11) * t422 / 0.2e1 + t244 * t444 / 0.2e1 + (-mrSges(5,3) * t77 + mrSges(6,1) * t73 - Ifges(5,4) * t469 + Ifges(6,6) * t470 + Ifges(6,3) * t471 - Ifges(5,2) * t472 - t123 / 0.2e1 + t121 / 0.2e1 + t525 * t465 - t518) * t170 + (t157 * mrSges(5,1) + t10 * mrSges(6,1) - t25 * mrSges(6,2) - t16 * mrSges(5,3) + t330 * t478 + t332 * t480 + t334 * t481 - t7 * t336 + t67 * t355 + (Ifges(5,2) + Ifges(6,3)) * t111 + (Ifges(5,4) + Ifges(6,6)) * t110 + t525 * t497) * t237 + t218 * (Ifges(7,5) * t322 - Ifges(7,6) * t321) / 0.2e1 + t41 * t30 + t42 * t31 - t225 * t391 / 0.2e1 + t226 * t390 / 0.2e1 + (-t72 * mrSges(6,1) - t76 * mrSges(5,3) - t510 / 0.2e1 - t509 / 0.2e1 - Ifges(5,1) * t469 + Ifges(6,2) * t470 + Ifges(6,6) * t471 - Ifges(5,4) * t472 - Ifges(7,5) * t474 + t122 / 0.2e1 + t526 * t465 - t517 - t524 / 0.2e1) * t169 + (-pkin(2) * t202 - (t247 * t294 + t297 * t327) * t397) * m(4) + (t257 + t210) * mrSges(3,1) + t67 * t428 / 0.2e1 + t68 * t429 / 0.2e1 - t110 * Ifges(5,1) * t238 + t296 * (Ifges(4,4) * t244 - Ifges(4,2) * t317) / 0.2e1 + (-Ifges(5,4) * t111 + Ifges(5,5) * t285 + t382) * t238 / 0.2e1 + (t14 * mrSges(6,1) + t157 * mrSges(5,2) - t17 * mrSges(5,3) - t25 * mrSges(6,3) + Ifges(5,5) * t467 + Ifges(7,5) * t481 + Ifges(7,6) * t480 + Ifges(7,3) * t478 + (-Ifges(5,4) / 0.2e1 - Ifges(6,6)) * t111 - Ifges(6,2) * t110 - t497 * Ifges(6,4) + t491) * t238 + t251 * t343 + Ifges(3,3) * qJDD(2) - t252 * t342 + (Ifges(4,1) * t244 + Ifges(4,4) * t515) * t293 + (t1 * t422 - t2 * t423 - t23 * t322 - t24 * t321 - t419 * t514) * mrSges(7,3) + t511 * t113 + t512 * t112 + (t1 * t42 + t127 * t7 + t2 * t41 + t23 * t511 + t24 * t512 + t47 * t49) * m(7) + (Ifges(7,1) * t322 - Ifges(7,4) * t321) * t474 + (m(4) * ((-t184 * t296 - t185 * t293) * qJD(3) + t501) + t296 * t212 - t293 * t213 - t251 * t390 - t252 * t391) * pkin(8) + (-t184 * t390 - t185 * t391 + t501) * mrSges(4,3) + (-m(5) * t216 - m(6) * t103 - t496) * t366 + t202 * t338 + t177 * (Ifges(7,4) * t322 - Ifges(7,2) * t321) / 0.2e1 + t489 * t179 + t490 * t180 + t488 * t81 + t484 * (-t292 * t343 + t342 * t464) + (-m(7) * (-pkin(10) * t427 + t345) + mrSges(7,3) * t427 - m(5) * t401 - m(6) * t345 + t483 * t220 + t485 * t219) * g(2) + (-m(7) * (-pkin(10) * t426 + t344) + mrSges(7,3) * t426 - m(5) * t400 - m(6) * t344 + t483 * t222 + t485 * t221) * g(1) + t49 * (mrSges(7,1) * t321 + mrSges(7,2) * t322) + (t319 + t331 * qJD(3) / 0.2e1) * qJD(3) + t534 * t487 + t151 * t279 + m(5) * (-t157 * t276 + t216 * t279) + m(6) * (t103 * t57 + t156 * t25) - t276 * t37 + qJDD(3) * (Ifges(4,5) * t293 + Ifges(4,6) * t296) + t156 * t36 - pkin(2) * t174 + t47 * t90 + t423 * t482 + (-t294 * t379 - t297 * t316 - t499) * t514 + t333 * t515 + t127 * t19 + t57 * t152; (-t521 * mrSges(4,1) - (-t222 * t296 - t293 * t418) * mrSges(4,2) - m(7) * (t312 - t462) - m(5) * t339 - m(6) * t312 + t494) * g(1) + (t112 * t387 - t113 * t388 + t523) * (-pkin(10) + t275) + (-t216 * t278 - t76 * t78 - t77 * t79 + (t464 * t17 + t16 * t292 + (t292 * t76 + t464 * t77) * qJD(4)) * pkin(3)) * m(5) - (-Ifges(4,2) * t395 + t226 + t277) * t393 / 0.2e1 + (t378 + t251) * t185 + t301 + t225 * t395 / 0.2e1 - t300 * t318 / 0.2e1 - t2 * t446 - t331 * t386 / 0.2e1 - t151 * t278 + t77 * t447 + (t377 - t252) * t184 + (-t24 * t387 - t459) * mrSges(7,3) + Ifges(4,3) * qJDD(3) + (m(7) * (t23 * t295 + t24 * t291) + t295 * t113 + t291 * t112 + t492) * pkin(3) * t389 + t513 * t273 + (-t23 * t26 - t24 * t27 + t273 * t7 + t49 * t506) * m(7) + t506 * t90 + (-t10 * t273 - t103 * t116 + t14 * t275 + t507 * t73 - t72 * t78) * m(6) + t507 * t188 + t503 * t78 + (-m(7) * (t305 - t458) - m(5) * t310 - m(6) * t305 - t314 * mrSges(4,1) - (-t220 * t296 + t293 * t354) * mrSges(4,2) + t493) * g(2) + (-m(7) * (t309 - t461) - m(5) * t326 - m(6) * t309 - mrSges(4,1) * t223 + mrSges(4,2) * t224 + t495) * g(3) - qJD(2) * t319 - Ifges(4,6) * t317 + t92 * t463 + Ifges(4,5) * t244 + t275 * t94 + t529 * t186 - t83 * mrSges(4,2) + t84 * mrSges(4,1) + t91 * t381 - t27 * t112 - t26 * t113 - t116 * t152; t301 + (-(qJD(6) * t112 + t30) * t479 + (-t2 - t431) * mrSges(7,3)) * t295 - t434 * qJD(5) + (t447 + t503) * t77 + t504 * t76 + (-t1 * mrSges(7,3) - (-qJD(6) * t113 + t31) * t479) * t291 - t479 * t306 + t493 * g(2) + t494 * g(1) + t495 * g(3) + t513 * qJ(5) + t324 * t90 - pkin(4) * t94 - t29 * t112 - t28 * t113 - t149 * t152 + (qJ(5) * t7 - t23 * t28 - t24 * t29 + (-t349 + t458) * g(2) + (-t348 + t462) * g(1) + (-t347 + t461) * g(3) + t519 * t49) * m(7) + (-pkin(4) * t14 - t348 * g(1) - t349 * g(2) - t347 * g(3) - qJ(5) * t10 - t103 * t149 + t505 * t73 - t72 * t77) * m(6); t434 * t286 + t328 * qJD(6) + (t152 + t328) * t228 + t94 + (-t228 * t329 - t286 * t49 + t315) * m(7) + (t103 * t228 + t286 * t73 + t14 + t315) * m(6) + t523; -t49 * (mrSges(7,1) * t178 + mrSges(7,2) * t177) + (Ifges(7,1) * t177 - t441) * t475 + t67 * t474 + (Ifges(7,5) * t177 - Ifges(7,6) * t178) * t473 - t23 * t112 + t24 * t113 - g(1) * ((t160 * t295 - t221 * t291) * mrSges(7,1) + (-t160 * t291 - t221 * t295) * mrSges(7,2)) - g(2) * ((t158 * t295 - t219 * t291) * mrSges(7,1) + (-t158 * t291 - t219 * t295) * mrSges(7,2)) - g(3) * ((t205 * t295 + t260) * mrSges(7,1) + (-t205 * t291 + t374) * mrSges(7,2)) + (t177 * t23 + t178 * t24) * mrSges(7,3) + t382 + (-Ifges(7,2) * t178 + t173 + t68) * t476 + t491;];
tau  = t3;
