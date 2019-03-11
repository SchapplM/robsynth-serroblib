% Calculate vector of inverse dynamics joint torques for
% S6RRPRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:15
% EndTime: 2019-03-09 09:55:10
% DurationCPUTime: 37.71s
% Computational Cost: add. (7982->746), mult. (17778->895), div. (0->0), fcn. (12406->10), ass. (0->345)
t494 = -mrSges(7,2) - mrSges(6,3);
t464 = mrSges(5,2) + t494;
t495 = -mrSges(6,2) + mrSges(5,1);
t528 = -mrSges(7,1) - mrSges(4,3);
t527 = Ifges(5,4) + Ifges(6,6);
t490 = Ifges(5,5) + Ifges(7,5);
t526 = Ifges(6,4) - t490;
t491 = Ifges(7,4) + Ifges(6,5);
t525 = -Ifges(5,6) + t491;
t275 = cos(qJ(2));
t380 = qJD(1) * t275;
t251 = -qJD(4) + t380;
t426 = -t251 / 0.2e1;
t463 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t524 = t463 * t426;
t272 = sin(qJ(4));
t398 = cos(pkin(9));
t423 = cos(qJ(4));
t332 = t423 * t398;
t318 = t275 * t332;
t269 = sin(pkin(9));
t364 = t269 * t380;
t171 = qJD(1) * t318 - t272 * t364;
t375 = qJD(4) * t272;
t469 = qJD(4) * t332 - t269 * t375;
t474 = t171 - t469;
t343 = t398 * qJD(2);
t273 = sin(qJ(2));
t381 = qJD(1) * t273;
t361 = t269 * t381;
t298 = -t343 + t361;
t190 = t423 * t298;
t348 = t273 * t398;
t331 = qJD(1) * t348;
t374 = t269 * qJD(2);
t211 = t331 + t374;
t373 = qJD(1) * qJD(2);
t354 = t275 * t373;
t305 = qJDD(1) * t273 + t354;
t286 = t269 * qJDD(2) + t305 * t398;
t338 = t398 * qJDD(2);
t287 = t269 * t305 - t338;
t62 = qJD(4) * t190 + t211 * t375 + t272 * t287 - t286 * t423;
t449 = -t62 / 0.2e1;
t448 = t62 / 0.2e1;
t285 = t211 * t423 - t272 * t298;
t63 = qJD(4) * t285 + t272 * t286 + t287 * t423;
t446 = t63 / 0.2e1;
t142 = t211 * t272 + t190;
t439 = t142 / 0.2e1;
t262 = t275 * qJDD(1);
t355 = t273 * t373;
t220 = -t262 + t355;
t213 = qJDD(4) + t220;
t430 = t213 / 0.2e1;
t438 = -t285 / 0.2e1;
t352 = -t380 / 0.2e1;
t523 = -mrSges(5,3) - mrSges(6,1);
t493 = Ifges(5,1) + Ifges(7,3);
t492 = -Ifges(5,4) + Ifges(7,6);
t489 = Ifges(7,2) + Ifges(6,3);
t488 = Ifges(6,6) - Ifges(7,6);
t215 = t269 * t423 + t272 * t398;
t292 = t275 * t215;
t170 = qJD(1) * t292;
t193 = t215 * qJD(4);
t475 = t170 - t193;
t270 = -pkin(8) - qJ(3);
t522 = -m(7) * (pkin(5) - t270) - m(4) * qJ(3) + t528;
t138 = Ifges(7,6) * t285;
t403 = t285 * Ifges(6,6);
t484 = t489 * t142 - t491 * t251 + t138 - t403;
t404 = t285 * Ifges(5,4);
t73 = -t142 * Ifges(5,2) - t251 * Ifges(5,6) + t404;
t519 = -t484 / 0.2e1 + t73 / 0.2e1;
t268 = pkin(9) + qJ(4);
t260 = sin(t268);
t261 = cos(t268);
t518 = -t464 * t260 + t495 * t261;
t517 = t305 / 0.2e1;
t140 = Ifges(5,4) * t142;
t405 = Ifges(7,6) * t142;
t483 = -t490 * t251 + t493 * t285 - t140 + t405;
t516 = -t483 / 0.2e1;
t226 = t270 * t269;
t346 = t398 * qJ(3);
t227 = pkin(8) * t398 + t346;
t358 = qJD(4) * t423;
t377 = qJD(3) * t269;
t105 = t272 * (qJD(4) * t227 + t377) - qJD(3) * t332 - t226 * t358;
t321 = pkin(2) * t273 - qJ(3) * t275;
t217 = t321 * qJD(1);
t159 = pkin(7) * t361 + t398 * t217;
t347 = t275 * t398;
t304 = pkin(3) * t273 - pkin(8) * t347;
t124 = qJD(1) * t304 + t159;
t194 = t269 * t217;
t336 = pkin(7) * t348;
t389 = t269 * t275;
t296 = -pkin(8) * t389 - t336;
t145 = qJD(1) * t296 + t194;
t78 = t272 * t124 + t423 * t145;
t64 = -qJ(5) * t381 - t78;
t515 = t105 - t64;
t157 = t272 * t226 + t227 * t423;
t106 = qJD(3) * t215 + qJD(4) * t157;
t328 = -t124 * t423 + t272 * t145;
t514 = t106 - t328;
t258 = pkin(7) * t380;
t204 = pkin(3) * t364 + t258;
t513 = qJ(5) * t474 - qJD(5) * t215 - t204;
t274 = sin(qJ(1));
t276 = cos(qJ(1));
t512 = g(1) * t276 + g(2) * t274;
t511 = -t142 * pkin(5) + qJD(6);
t350 = -qJ(3) * t273 - pkin(1);
t222 = -pkin(2) * t275 + t350;
t198 = t222 * qJD(1);
t229 = qJD(2) * qJ(3) + t258;
t149 = t269 * t198 + t398 * t229;
t102 = -pkin(8) * t298 + t149;
t148 = t398 * t198 - t229 * t269;
t291 = -pkin(3) * t380 - pkin(8) * t211 + t148;
t97 = t423 * t291;
t36 = t102 * t272 - t97;
t316 = pkin(5) * t285 + t36;
t510 = t316 + qJD(5);
t326 = mrSges(3,1) * t273 + mrSges(3,2) * t275;
t408 = Ifges(3,4) * t273;
t509 = t273 * (Ifges(3,1) * t275 - t408) / 0.2e1 - pkin(1) * t326;
t290 = t272 * t291;
t376 = qJD(3) * t273;
t395 = qJDD(1) * pkin(1);
t132 = t220 * pkin(2) - qJ(3) * t305 - qJD(1) * t376 - t395;
t255 = pkin(7) * t262;
t257 = pkin(7) * t381;
t178 = qJDD(2) * qJ(3) + t255 + (qJD(3) - t257) * qJD(2);
t88 = t398 * t132 - t269 * t178;
t45 = t220 * pkin(3) - pkin(8) * t286 + t88;
t89 = t269 * t132 + t398 * t178;
t51 = -pkin(8) * t287 + t89;
t7 = -qJD(4) * t290 - t102 * t358 - t272 * t51 + t423 * t45;
t300 = qJDD(5) - t7;
t415 = pkin(4) + qJ(6);
t1 = -t62 * pkin(5) + t251 * qJD(6) - t213 * t415 + t300;
t6 = qJD(4) * t97 - t102 * t375 + t272 * t45 + t423 * t51;
t4 = -qJ(5) * t213 + qJD(5) * t251 - t6;
t2 = -pkin(5) * t63 + qJDD(6) - t4;
t5 = -t213 * pkin(4) + t300;
t508 = -t7 * mrSges(5,1) + t6 * mrSges(5,2) - t5 * mrSges(6,2) - t2 * mrSges(7,2) + t4 * mrSges(6,3) + t1 * mrSges(7,3);
t437 = t285 / 0.2e1;
t440 = -t142 / 0.2e1;
t507 = -Ifges(5,4) * t440 + Ifges(6,2) * t438 + t426 * t526 - t493 * t437 + t488 * t439 + t516;
t506 = -Ifges(5,2) * t440 + Ifges(6,6) * t438 + t426 * t525 + t492 * t437 + t489 * t439 - t519;
t23 = t251 * t415 + t510;
t37 = t102 * t423 + t290;
t32 = t251 * qJ(5) - t37;
t26 = -t32 + t511;
t31 = pkin(4) * t251 + qJD(5) + t36;
t323 = Ifges(3,2) * t275 + t408;
t505 = t26 * mrSges(7,2) + t31 * mrSges(6,2) - Ifges(3,6) * qJD(2) / 0.2e1 - qJD(1) * t323 / 0.2e1 + Ifges(4,5) * t211 / 0.2e1 - Ifges(4,6) * t298 / 0.2e1 + Ifges(4,3) * t352 + t526 * t438 + t524 + t525 * t439 - t23 * mrSges(7,3) - t32 * mrSges(6,3) - t36 * mrSges(5,1) - t37 * mrSges(5,2);
t210 = t305 * pkin(7);
t187 = -qJDD(2) * pkin(2) + qJDD(3) + t210;
t121 = pkin(3) * t287 + t187;
t277 = t62 * qJ(5) - qJD(5) * t285 + t121;
t3 = t142 * qJD(6) + t415 * t63 + t277;
t447 = -t63 / 0.2e1;
t500 = -t213 / 0.2e1;
t8 = t63 * pkin(4) + t277;
t504 = t121 * mrSges(5,1) + t4 * mrSges(6,1) - t2 * mrSges(7,1) - t8 * mrSges(6,2) - t6 * mrSges(5,3) + t3 * mrSges(7,3) + Ifges(5,6) * t500 + 0.2e1 * t489 * t446 + t492 * t449 + (t525 + t491) * t430 + (t446 - t447) * Ifges(5,2) + (t488 + t527) * t448;
t503 = t5 * mrSges(6,1) + t1 * mrSges(7,1) + t121 * mrSges(5,2) - t3 * mrSges(7,2) - t7 * mrSges(5,3) - t8 * mrSges(6,3) + Ifges(6,4) * t500 + 0.2e1 * t493 * t449 + t527 * t447 + (-t488 + t492) * t446 + (-t526 + t490) * t430 + (-t448 + t449) * Ifges(6,2);
t502 = -m(3) - m(4);
t501 = -m(6) - m(7);
t427 = t220 / 0.2e1;
t498 = t286 / 0.2e1;
t497 = -t287 / 0.2e1;
t496 = qJD(2) / 0.2e1;
t214 = t269 * t272 - t332;
t485 = qJD(6) * t214 - t475 * t415 + t513;
t482 = t475 * pkin(5) - t515;
t357 = t415 * t273;
t481 = -t474 * pkin(5) + qJD(1) * t357 + t514;
t327 = t275 * mrSges(3,1) - t273 * mrSges(3,2);
t480 = -mrSges(2,1) - t327;
t479 = -t475 * pkin(4) + t513;
t311 = -mrSges(4,1) * t398 + t269 * mrSges(4,2);
t295 = m(4) * pkin(2) - t311;
t478 = t275 * t295;
t410 = mrSges(5,3) * t142;
t107 = mrSges(5,2) * t251 - t410;
t414 = mrSges(6,1) * t142;
t110 = mrSges(6,3) * t251 + t414;
t477 = -t107 + t110;
t409 = mrSges(5,3) * t285;
t108 = -mrSges(5,1) * t251 - t409;
t413 = mrSges(6,1) * t285;
t112 = -mrSges(6,2) * t251 + t413;
t476 = t108 - t112;
t412 = mrSges(7,1) * t142;
t111 = -mrSges(7,2) * t251 - t412;
t384 = -t110 + t111;
t473 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t298 + t211 * mrSges(4,2) + mrSges(3,3) * t381;
t265 = t275 * pkin(4);
t396 = qJ(5) * t275;
t472 = t260 * t396 + t261 * t265;
t471 = -t37 - t511;
t470 = t522 * t273 - t478;
t468 = -t269 * t88 + t398 * t89;
t209 = -pkin(7) * t355 + t255;
t467 = t209 * t275 + t210 * t273;
t466 = t523 * t273;
t465 = -m(5) + t501;
t313 = mrSges(4,1) * t269 + mrSges(4,2) * t398;
t460 = -mrSges(3,3) - t313 + mrSges(2,2);
t221 = -qJD(2) * pkin(2) + qJD(3) + t257;
t459 = -t221 * t275 * t313 - t148 * (mrSges(4,1) * t273 - mrSges(4,3) * t347) - t149 * (-mrSges(4,2) * t273 - mrSges(4,3) * t389);
t406 = Ifges(4,4) * t269;
t310 = Ifges(4,1) * t398 - t406;
t359 = Ifges(4,4) * t398;
t312 = -Ifges(4,2) * t269 + t359;
t458 = (Ifges(4,1) * t211 - Ifges(4,4) * t298 - Ifges(4,5) * t380) * t347 + t211 * (Ifges(4,5) * t273 + t275 * t310) - (Ifges(4,6) * t273 + t275 * t312) * t298;
t457 = t463 * t213 + t525 * t63 + t526 * t62;
t456 = m(7) * qJ(6) + mrSges(7,3);
t454 = t456 + t495;
t139 = Ifges(6,6) * t142;
t72 = -Ifges(6,4) * t251 - Ifges(6,2) * t285 + t139;
t445 = t72 / 0.2e1;
t443 = Ifges(4,1) * t498 + Ifges(4,4) * t497 + Ifges(4,5) * t427;
t425 = t251 / 0.2e1;
t422 = pkin(3) * t269;
t419 = g(3) * t273;
t263 = t273 * pkin(7);
t411 = mrSges(7,1) * t285;
t407 = Ifges(3,4) * t275;
t397 = qJ(5) * t142;
t394 = t187 * t273;
t391 = t261 * t273;
t390 = t269 * t273;
t388 = t270 * t273;
t387 = t273 * t276;
t386 = t274 * t275;
t253 = pkin(3) * t398 + pkin(2);
t234 = t275 * t253;
t385 = t275 * t276;
t208 = t398 * t222;
t147 = -pkin(8) * t348 + t208 + (-pkin(7) * t269 - pkin(3)) * t275;
t167 = pkin(7) * t347 + t269 * t222;
t155 = -pkin(8) * t390 + t167;
t87 = t272 * t147 + t423 * t155;
t189 = qJD(2) * t321 - t376;
t379 = qJD(2) * t273;
t368 = pkin(7) * t379;
t153 = t398 * t189 + t269 * t368;
t378 = qJD(2) * t275;
t259 = pkin(7) * t378;
t363 = t275 * t374;
t205 = pkin(3) * t363 + t259;
t218 = pkin(3) * t390 + t263;
t382 = t276 * pkin(1) + t274 * pkin(7);
t367 = m(4) - t465;
t366 = t270 * t387;
t266 = t276 * pkin(7);
t365 = t274 * t388 + t276 * t422 + t266;
t20 = t62 * mrSges(7,2) + t63 * mrSges(7,3);
t44 = -t62 * mrSges(6,1) + t213 * mrSges(6,2);
t43 = -t63 * mrSges(7,1) + t213 * mrSges(7,2);
t41 = -t62 * mrSges(7,1) - t213 * mrSges(7,3);
t345 = -t373 / 0.2e1;
t340 = -qJ(5) * t260 - t253;
t156 = -t423 * t226 + t227 * t272;
t339 = t234 - t388;
t337 = t253 * t385 + t274 * t422 + t382;
t334 = -m(7) * t415 - mrSges(7,3);
t79 = -t87 + t396;
t186 = -t272 * t390 + t273 * t332;
t329 = -qJ(5) * t186 + t218;
t86 = t147 * t423 - t272 * t155;
t322 = Ifges(3,5) * t275 - Ifges(3,6) * t273;
t80 = t265 - t86;
t309 = Ifges(4,5) * t398 - Ifges(4,6) * t269;
t113 = qJD(2) * t304 + t153;
t176 = t269 * t189;
t125 = qJD(2) * t296 + t176;
t306 = -t113 * t423 + t272 * t125 + t147 * t375 + t155 * t358;
t24 = t272 * t113 + t423 * t125 + t147 * t358 - t155 * t375;
t303 = -t215 * qJ(5) - t253;
t117 = -qJD(2) * t318 + t193 * t273 + t272 * t363;
t302 = qJ(5) * t117 - qJD(5) * t186 + t205;
t181 = t260 * t385 - t274 * t261;
t182 = t274 * t260 + t261 * t385;
t299 = t182 * pkin(4) + t181 * qJ(5) + t337;
t179 = t260 * t386 + t261 * t276;
t297 = -g(1) * t181 - g(2) * t179 - t260 * t419;
t288 = t275 * (Ifges(4,3) * t273 + t275 * t309);
t21 = -qJ(5) * t379 + qJD(5) * t275 - t24;
t158 = pkin(3) * t298 + t221;
t281 = -qJ(5) * t285 + t158;
t256 = Ifges(3,4) * t380;
t232 = qJ(5) * t391;
t230 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t380;
t197 = Ifges(3,1) * t381 + Ifges(3,5) * qJD(2) + t256;
t185 = t215 * t273;
t180 = -t260 * t276 + t261 * t386;
t169 = -mrSges(4,1) * t380 - mrSges(4,3) * t211;
t168 = mrSges(4,2) * t380 - mrSges(4,3) * t298;
t166 = -pkin(7) * t389 + t208;
t165 = t286 * mrSges(4,2);
t160 = -pkin(7) * t331 + t194;
t154 = -qJD(2) * t336 + t176;
t137 = t214 * pkin(4) + t303;
t134 = Ifges(4,4) * t211 - Ifges(4,2) * t298 - Ifges(4,6) * t380;
t129 = t220 * mrSges(4,1) - mrSges(4,3) * t286;
t128 = -t220 * mrSges(4,2) - mrSges(4,3) * t287;
t118 = qJD(2) * t292 + t273 * t469;
t116 = -t214 * pkin(5) + t157;
t115 = pkin(5) * t215 + t156;
t109 = mrSges(7,3) * t251 + t411;
t104 = mrSges(4,1) * t287 + t165;
t101 = t214 * t415 + t303;
t98 = pkin(4) * t185 + t329;
t93 = Ifges(4,4) * t286 - t287 * Ifges(4,2) + Ifges(4,6) * t220;
t85 = -mrSges(6,2) * t142 - mrSges(6,3) * t285;
t84 = mrSges(5,1) * t142 + mrSges(5,2) * t285;
t83 = pkin(4) * t285 + t397;
t82 = -mrSges(7,2) * t285 + mrSges(7,3) * t142;
t81 = t185 * t415 + t329;
t65 = -pkin(4) * t381 + t328;
t53 = t62 * mrSges(5,2);
t52 = t62 * mrSges(6,3);
t48 = -pkin(5) * t185 - t79;
t47 = t142 * pkin(4) + t281;
t46 = t186 * pkin(5) + t275 * qJ(6) + t80;
t42 = mrSges(6,1) * t63 - mrSges(6,3) * t213;
t40 = -mrSges(5,2) * t213 - mrSges(5,3) * t63;
t39 = mrSges(5,1) * t213 + mrSges(5,3) * t62;
t35 = t285 * t415 + t397;
t30 = pkin(4) * t118 + t302;
t29 = t142 * t415 + t281;
t22 = -pkin(4) * t379 + t306;
t19 = t63 * mrSges(5,1) - t53;
t18 = -t63 * mrSges(6,2) + t52;
t17 = qJD(6) * t185 + t118 * t415 + t302;
t10 = -pkin(5) * t118 - t21;
t9 = -t117 * pkin(5) - qJD(2) * t357 + t275 * qJD(6) + t306;
t11 = [(-t348 * t88 - t390 * t89) * mrSges(4,3) + (-qJDD(2) * mrSges(3,1) + mrSges(3,3) * t305 + t104) * t263 + t407 * t517 + t467 * mrSges(3,3) + t503 * t186 + t504 * t185 + (mrSges(5,1) * t158 + t32 * mrSges(6,1) - t26 * mrSges(7,1) - t47 * mrSges(6,2) - t37 * mrSges(5,3) + mrSges(7,3) * t29 + t506) * t118 + (-t31 * mrSges(6,1) - t23 * mrSges(7,1) - mrSges(5,2) * t158 + mrSges(7,2) * t29 - t36 * mrSges(5,3) + mrSges(6,3) * t47 + t445 + t507) * t117 + t458 * t496 - (t408 + t323) * t220 / 0.2e1 - t306 * t108 + m(5) * (t121 * t218 + t158 * t205 + t24 * t37 + t306 * t36 + t6 * t87 + t7 * t86) + t327 * t395 + t509 * t373 + (Ifges(6,4) * t438 + Ifges(5,6) * t440 + t490 * t437 + t491 * t439 + t505 + t524) * t379 - t230 * t368 - pkin(1) * (t220 * mrSges(3,1) + mrSges(3,2) * t305) - t93 * t390 / 0.2e1 + m(4) * (t148 * t153 + t149 * t154 + t166 * t88 + t167 * t89 + (t221 * t378 + t394) * pkin(7)) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t467) + t166 * t129 + t167 * t128 + t154 * t168 + t153 * t169 + t24 * t107 + t9 * t109 + t21 * t110 + t10 * t111 + t22 * t112 + t98 * t18 + t30 * t85 + t86 * t39 + t87 * t40 + t79 * t42 + t80 * t44 + t81 * t20 + t17 * t82 + t46 * t41 + t48 * t43 - t134 * t363 / 0.2e1 + (Ifges(3,5) * qJDD(2) + t310 * t498 + t312 * t497 + t305 * Ifges(3,1) + t309 * t427 + (m(7) * pkin(5) - t528) * t274 * g(1)) * t273 + (-m(5) * t365 + t501 * (-t180 * pkin(4) - qJ(5) * t179 + t365) + t502 * t266 + t460 * t276 + t454 * t180 - t464 * t179 + (m(3) * pkin(1) - m(4) * t350 + t478 + t465 * (-pkin(1) - t234) - t466 - t480) * t274) * g(1) + m(7) * (t1 * t46 + t10 * t26 + t17 * t29 + t2 * t48 + t23 * t9 + t3 * t81) + m(6) * (t21 * t32 + t22 * t31 + t30 * t47 + t4 * t79 + t5 * t80 + t8 * t98) + t473 * t259 + Ifges(2,3) * qJDD(1) + t348 * t443 + t197 * t378 / 0.2e1 + ((-t497 + t287 / 0.2e1) * Ifges(4,6) + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * qJDD(2) - t463 * t430 + t89 * mrSges(4,2) - Ifges(5,6) * t447 - Ifges(6,4) * t448 - t88 * mrSges(4,1) + (-Ifges(3,2) * t273 + t407) * t373 / 0.2e1 - t491 * t446 - t490 * t449 - Ifges(4,3) * t427 + t508 + (-t498 - t286 / 0.2e1) * Ifges(4,5) - t457 / 0.2e1 + (-pkin(7) * mrSges(3,3) - Ifges(4,3) / 0.2e1 - Ifges(3,2) / 0.2e1) * t220 + Ifges(3,4) * t517) * t275 + (-m(5) * (t337 - t366) - m(6) * (t299 - t366) - m(7) * t299 + t523 * t387 + t502 * t382 + t460 * t274 - t454 * t182 + t464 * t181 + (t470 + t480) * t276) * g(2) + (t322 * t496 - t459) * qJD(2) + t205 * t84 + t218 * t19 + t288 * t345 + t313 * t394; (-Ifges(5,2) * t439 + Ifges(6,6) * t437 + t425 * t525 + t492 * t438 + t489 * t440 + t519) * t170 + (Ifges(5,4) * t439 - Ifges(6,2) * t437 - t425 * t526 + t493 * t438 - t488 * t440 + t516) * t171 + (t44 - t39) * t156 + (-t159 - t377) * t169 + (-t36 * t474 + t37 * t475) * mrSges(5,3) + (-t31 * t474 - t32 * t475) * mrSges(6,1) + (t256 + t197) * t352 + t503 * t215 + t504 * t214 + (Ifges(6,4) * t437 - Ifges(3,2) * t352 + Ifges(5,6) * t439 + t463 * t425 + t490 * t438 + t491 * t440 - t505) * t381 + t506 * t193 + (Ifges(4,2) * t398 + t406) * t497 + (Ifges(4,1) * t269 + t359) * t498 + t328 * t108 + (-t23 * t474 + t26 * t475) * mrSges(7,1) + (-t458 / 0.2e1 + t459 + (t288 / 0.2e1 - t509) * qJD(1)) * qJD(1) + (-t121 * t253 - t156 * t7 + t157 * t6 - t158 * t204 + (-t105 - t78) * t37 + t514 * t36) * m(5) - t507 * t469 + Ifges(3,5) * t305 + (-t148 * t159 - t149 * t160 - t221 * t258 - t187 * pkin(2) + (-t148 * t269 + t149 * t398) * qJD(3) + t468 * qJ(3)) * m(4) + t468 * mrSges(4,3) + t137 * t18 + t115 * t41 + t116 * t43 - pkin(2) * t104 - t78 * t107 - t64 * t110 - t65 * t112 + t101 * t20 + t479 * t85 + t269 * t443 + (-m(5) * t339 - m(6) * (t339 + t472) - m(7) * (t234 + t472) - t327 + (-t261 * t456 - t518) * t275 + t466 + t470) * g(3) + (-t42 + t40) * t157 + (t137 * t8 + t156 * t5 - t157 * t4 + t479 * t47 + t515 * t32 + (t106 - t65) * t31) * m(6) - t473 * t258 + t474 * t445 + t477 * t105 + (mrSges(6,2) * t475 + mrSges(6,3) * t474) * t47 + (mrSges(7,2) * t474 - mrSges(7,3) * t475) * t29 + (-mrSges(5,1) * t475 - mrSges(5,2) * t474) * t158 - t476 * t106 + t398 * t93 / 0.2e1 + Ifges(3,3) * qJDD(2) - t269 * qJ(3) * t129 + t481 * t109 + t482 * t111 + t485 * t82 + (t1 * t115 + t101 * t3 + t116 * t2 + t23 * t481 + t26 * t482 + t29 * t485) * m(7) + t230 * t257 + t512 * (t326 + ((m(5) + m(6)) * t270 + t522 + t523) * t275 + (m(5) * t253 + t295 - m(6) * (-pkin(4) * t261 + t340) - m(7) * t340 - t261 * t334 + t518) * t273) - t204 * t84 - t209 * mrSges(3,2) - t210 * mrSges(3,1) - Ifges(3,6) * t220 + t322 * t345 + t128 * t346 + (Ifges(4,5) * t269 + Ifges(4,6) * t398) * t427 - t253 * t19 + t187 * t311 + t134 * t364 / 0.2e1 + (qJD(3) * t398 - t160) * t168; -t168 * t343 - (-t269 * t354 + t338) * mrSges(4,1) - t53 + t52 + t165 + t211 * t169 + t495 * t63 - (-t107 - t384) * t142 + (-t109 + t476) * t285 + t367 * t275 * g(3) + ((mrSges(4,1) * qJDD(1) + qJD(1) * t168) * t269 - t512 * t367) * t273 + t20 + (t142 * t26 - t23 * t285 + t3) * m(7) + (-t142 * t32 - t285 * t31 + t8) * m(6) + (t142 * t37 - t285 * t36 + t121) * m(5) + (t148 * t211 + t149 * t298 + t187) * m(4); (t142 * t526 + t285 * t525) * t425 + (-g(3) * t232 + qJ(5) * t2 - t1 * t415 + t471 * t23 + t26 * t510 - t29 * t35) * m(7) + (Ifges(6,2) * t142 + t403 + t73) * t437 + t316 * t111 - t415 * t41 - t508 - t83 * t85 - t35 * t82 - t32 * t413 + (-pkin(4) * t5 - qJ(5) * t4 - qJD(5) * t32 - t47 * t83) * m(6) + t471 * t109 + (-t42 + t43) * qJ(5) + (-m(6) * t32 + t410 - t477) * t36 + t384 * qJD(5) + (-m(6) * t31 + t409 + t476) * t37 + t457 - (-mrSges(5,1) * t260 - mrSges(5,2) * t261) * t419 - t158 * (mrSges(5,1) * t285 - mrSges(5,2) * t142) - t47 * (-mrSges(6,2) * t285 + mrSges(6,3) * t142) - t29 * (mrSges(7,2) * t142 + mrSges(7,3) * t285) + (-Ifges(5,2) * t285 - t140 + t483) * t439 + (t285 * t489 + t139 - t405 + t72) * t440 + (-t142 * t493 + t138 - t404 + t484) * t438 + (-m(6) * t232 + (t494 * t261 + (m(6) * pkin(4) - mrSges(6,2) - t334) * t260) * t273) * g(3) + (t501 * (-t181 * pkin(4) + qJ(5) * t182) + t464 * t182 + t454 * t181) * g(1) + (t501 * (-t179 * pkin(4) + qJ(5) * t180) + t464 * t180 + t454 * t179) * g(2) + t26 * t411 + t23 * t412 + t31 * t414 - pkin(4) * t44; t384 * t251 + (t82 + t85) * t285 + t41 + t44 + (t251 * t26 + t285 * t29 + t1 + t297) * m(7) + (-t251 * t32 + t285 * t47 + t297 + t5) * m(6); -t251 * t109 - t142 * t82 + (-g(1) * t182 - g(2) * t180 - g(3) * t391 - t29 * t142 - t23 * t251 + t2) * m(7) + t43;];
tau  = t11;
