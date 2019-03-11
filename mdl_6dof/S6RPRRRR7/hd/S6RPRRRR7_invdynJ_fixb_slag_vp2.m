% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:57
% EndTime: 2019-03-09 07:17:22
% DurationCPUTime: 13.75s
% Computational Cost: add. (16146->670), mult. (31814->886), div. (0->0), fcn. (21986->14), ass. (0->321)
t281 = qJD(3) + qJD(4);
t271 = qJD(5) + t281;
t287 = sin(qJ(6));
t292 = cos(qJ(6));
t289 = sin(qJ(4));
t290 = sin(qJ(3));
t294 = cos(qJ(4));
t295 = cos(qJ(3));
t219 = -t289 * t295 - t294 * t290;
t206 = t219 * qJD(1);
t375 = qJD(1) * t295;
t376 = qJD(1) * t290;
t207 = -t289 * t376 + t294 * t375;
t288 = sin(qJ(5));
t293 = cos(qJ(5));
t323 = t206 * t288 + t293 * t207;
t125 = t271 * t292 - t287 * t323;
t126 = t271 * t287 + t292 * t323;
t421 = mrSges(6,3) * t323;
t489 = -mrSges(6,1) * t271 - mrSges(7,1) * t125 + mrSges(7,2) * t126 + t421;
t298 = -pkin(1) - pkin(7);
t245 = qJD(1) * t298 + qJD(2);
t392 = t245 * t290;
t189 = -pkin(8) * t376 + t392;
t178 = t289 * t189;
t226 = t295 * t245;
t190 = -pkin(8) * t375 + t226;
t181 = qJD(3) * pkin(3) + t190;
t134 = t294 * t181 - t178;
t196 = t207 * pkin(9);
t114 = t134 - t196;
t106 = pkin(4) * t281 + t114;
t179 = t294 * t189;
t135 = t181 * t289 + t179;
t437 = pkin(9) * t206;
t115 = t135 + t437;
t403 = t115 * t288;
t67 = t106 * t293 - t403;
t65 = -pkin(5) * t271 - t67;
t470 = -m(7) * t65 - t489;
t505 = -m(6) * t67 - t470;
t239 = pkin(3) * t376 + qJD(1) * qJ(2);
t169 = -pkin(4) * t206 + t239;
t448 = -t271 / 0.2e1;
t343 = t293 * t206 - t207 * t288;
t453 = -t343 / 0.2e1;
t148 = qJD(6) - t343;
t454 = -t148 / 0.2e1;
t456 = -t126 / 0.2e1;
t457 = -t125 / 0.2e1;
t381 = t293 * t115;
t68 = t106 * t288 + t381;
t66 = pkin(10) * t271 + t68;
t80 = -pkin(5) * t343 - pkin(10) * t323 + t169;
t23 = t287 * t80 + t292 * t66;
t490 = t23 * mrSges(7,2);
t22 = -t287 * t66 + t292 * t80;
t491 = t22 * mrSges(7,1);
t504 = -t169 * mrSges(6,1) + Ifges(7,5) * t456 - Ifges(6,2) * t453 - Ifges(6,6) * t448 + Ifges(7,6) * t457 + Ifges(7,3) * t454 + t490 - t491;
t501 = t295 / 0.2e1;
t278 = t290 * pkin(3);
t286 = qJ(3) + qJ(4);
t273 = sin(t286);
t274 = cos(t286);
t339 = -mrSges(5,1) * t273 - mrSges(5,2) * t274;
t340 = mrSges(4,1) * t290 + mrSges(4,2) * t295;
t500 = m(5) * t278 - t339 + t340;
t277 = qJ(5) + t286;
t259 = sin(t277);
t260 = cos(t277);
t499 = mrSges(6,1) * t259 + (mrSges(6,2) - mrSges(7,3)) * t260;
t425 = mrSges(7,2) * t287;
t428 = mrSges(7,1) * t292;
t498 = t425 - t428;
t364 = qJD(1) * qJD(3);
t224 = qJDD(1) * t295 - t290 * t364;
t225 = -qJDD(1) * t290 - t295 * t364;
t307 = t219 * qJD(4);
t129 = qJD(1) * t307 + t224 * t294 + t225 * t289;
t280 = qJDD(3) + qJDD(4);
t244 = qJDD(1) * t298 + qJDD(2);
t374 = qJD(3) * t290;
t170 = t295 * t244 - t245 * t374;
t143 = qJDD(3) * pkin(3) - pkin(8) * t224 + t170;
t373 = qJD(3) * t295;
t171 = t290 * t244 + t245 * t373;
t150 = pkin(8) * t225 + t171;
t71 = -qJD(4) * t135 + t294 * t143 - t150 * t289;
t38 = pkin(4) * t280 - pkin(9) * t129 + t71;
t380 = t294 * t295;
t321 = t289 * t290 - t380;
t130 = qJD(1) * qJD(4) * t321 - t224 * t289 + t225 * t294;
t370 = qJD(4) * t294;
t371 = qJD(4) * t289;
t70 = t289 * t143 + t294 * t150 + t181 * t370 - t189 * t371;
t45 = pkin(9) * t130 + t70;
t15 = -qJD(5) * t68 - t288 * t45 + t293 * t38;
t270 = qJDD(5) + t280;
t12 = -pkin(5) * t270 - t15;
t368 = qJD(5) * t293;
t369 = qJD(5) * t288;
t14 = t106 * t368 - t115 * t369 + t288 * t38 + t293 * t45;
t147 = Ifges(6,4) * t343;
t11 = pkin(10) * t270 + t14;
t283 = qJD(1) * qJD(2);
t246 = qJDD(1) * qJ(2) + t283;
t177 = -pkin(3) * t225 + t246;
t105 = -pkin(4) * t130 + t177;
t63 = qJD(5) * t343 + t129 * t293 + t130 * t288;
t64 = -qJD(5) * t323 - t129 * t288 + t130 * t293;
t17 = -pkin(5) * t64 - pkin(10) * t63 + t105;
t2 = qJD(6) * t22 + t11 * t292 + t17 * t287;
t3 = -qJD(6) * t23 - t11 * t287 + t17 * t292;
t337 = mrSges(7,1) * t287 + mrSges(7,2) * t292;
t312 = t65 * t337;
t331 = Ifges(7,5) * t292 - Ifges(7,6) * t287;
t414 = Ifges(7,4) * t292;
t333 = -Ifges(7,2) * t287 + t414;
t415 = Ifges(7,4) * t287;
t335 = Ifges(7,1) * t292 - t415;
t367 = qJD(6) * t287;
t350 = -t367 / 0.2e1;
t124 = Ifges(7,4) * t125;
t59 = t126 * Ifges(7,1) + t148 * Ifges(7,5) + t124;
t405 = t292 * t59;
t352 = t405 / 0.2e1;
t366 = qJD(6) * t292;
t416 = Ifges(6,4) * t323;
t419 = mrSges(7,3) * t292;
t420 = mrSges(7,3) * t287;
t447 = t287 / 0.2e1;
t451 = t323 / 0.2e1;
t452 = -t323 / 0.2e1;
t62 = qJDD(6) - t64;
t459 = t62 / 0.2e1;
t36 = -qJD(6) * t126 + t270 * t292 - t287 * t63;
t461 = t36 / 0.2e1;
t35 = qJD(6) * t125 + t270 * t287 + t292 * t63;
t462 = t35 / 0.2e1;
t57 = t126 * Ifges(7,5) + t125 * Ifges(7,6) + t148 * Ifges(7,3);
t413 = t126 * Ifges(7,4);
t58 = t125 * Ifges(7,2) + t148 * Ifges(7,6) + t413;
t8 = t35 * Ifges(7,4) + t36 * Ifges(7,2) + t62 * Ifges(7,6);
t9 = t35 * Ifges(7,1) + t36 * Ifges(7,4) + t62 * Ifges(7,5);
t96 = Ifges(6,2) * t343 + t271 * Ifges(6,6) + t416;
t97 = Ifges(6,1) * t323 + t271 * Ifges(6,5) + t147;
t497 = -t14 * mrSges(6,2) - t3 * t420 + t2 * t419 + t12 * t498 + t15 * mrSges(6,1) + Ifges(6,3) * t270 + (Ifges(7,1) * t287 + t414) * t462 + (Ifges(7,2) * t292 + t415) * t461 + (Ifges(7,5) * t287 + Ifges(7,6) * t292) * t459 + t58 * t350 + t9 * t447 + Ifges(6,6) * t64 + Ifges(6,5) * t63 + t292 * t8 / 0.2e1 + (-t22 * t366 - t23 * t367) * mrSges(7,3) + (t312 + t352) * qJD(6) + (t125 * t333 + t126 * t335 + t148 * t331) * qJD(6) / 0.2e1 + t96 * t451 + (t147 + t97) * t453 + (t57 - t416) * t452;
t495 = -t169 * mrSges(6,2) + Ifges(6,1) * t452 + Ifges(6,5) * t448 + t22 * t419 + t23 * t420 + t331 * t454 + t333 * t457 + t335 * t456 - t312 - t405 / 0.2e1 + t58 * t447;
t494 = -m(4) - m(5);
t493 = -m(7) - m(6);
t348 = -pkin(5) * t259 + t260 * pkin(10);
t492 = m(7) * t348;
t164 = t219 * t288 - t293 * t321;
t431 = pkin(8) - t298;
t229 = t431 * t290;
t230 = t431 * t295;
t168 = -t294 * t229 - t289 * t230;
t296 = cos(qJ(1));
t360 = t260 * t425;
t426 = mrSges(6,2) * t259;
t486 = (-t360 - t426) * t296;
t291 = sin(qJ(1));
t357 = t260 * t428;
t390 = t260 * t291;
t391 = t259 * t291;
t485 = -mrSges(6,1) * t390 - mrSges(7,3) * t391 - (t357 - t360) * t291;
t484 = -t259 * t498 + t499;
t165 = t281 * t380 - t289 * t374 - t290 * t371;
t166 = qJD(3) * t219 + t307;
t322 = t293 * t219 + t288 * t321;
t84 = qJD(5) * t322 - t165 * t288 + t293 * t166;
t314 = t164 * t366 + t287 * t84;
t441 = pkin(4) * t274;
t444 = pkin(3) * t295;
t228 = t441 + t444;
t483 = m(5) * t444 + m(6) * t228;
t482 = -t170 * t295 - t171 * t290;
t18 = mrSges(7,1) * t62 - mrSges(7,3) * t35;
t19 = -mrSges(7,2) * t62 + mrSges(7,3) * t36;
t481 = -t287 * t18 + t292 * t19;
t427 = mrSges(5,2) * t273;
t480 = t426 + t427;
t479 = mrSges(6,1) * t260 + t259 * mrSges(7,3) + t357;
t478 = -g(1) * t291 + g(2) * t296;
t341 = mrSges(4,1) * t295 - mrSges(4,2) * t290;
t417 = Ifges(4,4) * t295;
t477 = (-Ifges(4,1) * t290 - t417) * t501 + qJ(2) * t341;
t422 = mrSges(6,3) * t343;
t131 = -mrSges(6,2) * t271 + t422;
t90 = -mrSges(7,2) * t148 + mrSges(7,3) * t125;
t91 = mrSges(7,1) * t148 - mrSges(7,3) * t126;
t476 = -t287 * t91 + t292 * t90 + t131;
t475 = -t296 * t427 + t486;
t474 = mrSges(5,1) * t274 + t479;
t16 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t473 = -m(7) * t12 - t16;
t237 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t376;
t238 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t375;
t472 = (t237 * t295 - t238 * t290) * qJD(3);
t471 = -t134 * t166 - t135 * t165 + t219 * t70 + t321 * t71;
t469 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t468 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t467 = -m(4) * t482 + t295 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t224) + t290 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t225);
t378 = pkin(5) * t390 + pkin(10) * t391;
t466 = -m(7) * t378 + t485;
t103 = pkin(5) * t323 - pkin(10) * t343;
t465 = -mrSges(3,3) + t492 + mrSges(2,2) - t499 - t500;
t330 = t22 * t292 + t23 * t287;
t304 = -qJD(6) * t330 + t2 * t292 - t287 * t3;
t464 = m(7) * t304 - t91 * t366 - t90 * t367 + t481;
t463 = qJD(1) ^ 2;
t458 = -m(3) - m(4);
t297 = -pkin(8) - pkin(7);
t455 = t126 / 0.2e1;
t449 = t207 / 0.2e1;
t445 = pkin(3) * t294;
t443 = pkin(4) * t207;
t442 = pkin(4) * t273;
t440 = pkin(4) * t288;
t439 = pkin(4) * t293;
t434 = t67 * mrSges(6,3);
t433 = t68 * mrSges(6,3);
t424 = mrSges(5,3) * t206;
t423 = mrSges(5,3) * t207;
t418 = Ifges(4,4) * t290;
t412 = t207 * Ifges(5,4);
t400 = t164 * t287;
t399 = t164 * t292;
t389 = t274 * t291;
t388 = t287 * t291;
t387 = t287 * t296;
t386 = t288 * t289;
t385 = t289 * t293;
t383 = t291 * t292;
t382 = t292 * t296;
t256 = qJ(2) + t278;
t141 = t294 * t190 - t178;
t263 = pkin(4) + t445;
t201 = pkin(3) * t385 + t288 * t263;
t377 = t296 * pkin(1) + t291 * qJ(2);
t365 = qJDD(1) * mrSges(3,2);
t247 = pkin(3) * t373 + qJD(2);
t362 = m(6) * t441;
t361 = Ifges(7,5) * t35 + Ifges(7,6) * t36 + Ifges(7,3) * t62;
t266 = pkin(3) * t375;
t276 = t296 * qJ(2);
t349 = -pkin(1) * t291 + t276;
t347 = -t364 / 0.2e1;
t345 = (t246 + t283) * qJ(2);
t140 = -t190 * t289 - t179;
t167 = t229 * t289 - t294 * t230;
t182 = -pkin(4) * t219 + t256;
t149 = pkin(4) * t165 + t247;
t342 = -pkin(5) * t260 - pkin(10) * t259;
t336 = t295 * Ifges(4,1) - t418;
t334 = -t290 * Ifges(4,2) + t417;
t332 = -Ifges(4,5) * t290 - Ifges(4,6) * t295;
t329 = -t22 * t287 + t23 * t292;
t328 = -t287 * t90 - t292 * t91;
t138 = pkin(9) * t321 + t167;
t139 = pkin(9) * t219 + t168;
t93 = t138 * t288 + t139 * t293;
t94 = -pkin(5) * t322 - pkin(10) * t164 + t182;
t42 = t287 * t94 + t292 * t93;
t41 = -t287 * t93 + t292 * t94;
t325 = t293 * t138 - t139 * t288;
t318 = t348 - t442;
t200 = -pkin(3) * t386 + t263 * t293;
t315 = t140 - t437;
t313 = t164 * t367 - t292 * t84;
t310 = t290 * (-Ifges(4,2) * t295 - t418);
t214 = t431 * t374;
t215 = qJD(3) * t230;
t108 = t289 * t214 - t294 * t215 + t229 * t371 - t230 * t370;
t88 = t103 + t443;
t109 = -qJD(4) * t168 + t294 * t214 + t215 * t289;
t303 = qJD(5) * t164 + t293 * t165 + t166 * t288;
t302 = -pkin(9) * t166 + t109;
t145 = t206 * Ifges(5,2) + t281 * Ifges(5,6) + t412;
t195 = Ifges(5,4) * t206;
t146 = t207 * Ifges(5,1) + t281 * Ifges(5,5) + t195;
t300 = -t239 * (mrSges(5,1) * t207 + mrSges(5,2) * t206) + t145 * t449 - t207 * (Ifges(5,1) * t206 - t412) / 0.2e1 + t135 * t423 + t134 * t424 - t281 * (Ifges(5,5) * t206 - Ifges(5,6) * t207) / 0.2e1 + Ifges(5,3) * t280 + Ifges(5,5) * t129 + Ifges(5,6) * t130 - t70 * mrSges(5,2) + t71 * mrSges(5,1) - (-Ifges(5,2) * t207 + t146 + t195) * t206 / 0.2e1 + (t433 + t504) * t323 + (t434 + t495) * t343 + t497;
t284 = -pkin(9) + t297;
t269 = -pkin(1) * qJDD(1) + qJDD(2);
t262 = -pkin(5) - t439;
t242 = mrSges(5,1) * t389;
t227 = t278 + t442;
t221 = t340 * qJD(1);
t205 = Ifges(4,5) * qJD(3) + qJD(1) * t336;
t204 = Ifges(4,6) * qJD(3) + qJD(1) * t334;
t191 = -pkin(5) - t200;
t188 = t259 * t382 - t388;
t187 = t259 * t387 + t383;
t186 = t259 * t383 + t387;
t185 = -t259 * t388 + t382;
t176 = mrSges(5,1) * t281 - t423;
t175 = -mrSges(5,2) * t281 + t424;
t174 = t266 + t443;
t158 = -mrSges(5,1) * t206 + mrSges(5,2) * t207;
t120 = -mrSges(5,2) * t280 + mrSges(5,3) * t130;
t119 = mrSges(5,1) * t280 - mrSges(5,3) * t129;
t116 = -t196 + t141;
t102 = -mrSges(6,1) * t343 + mrSges(6,2) * t323;
t89 = -pkin(9) * t165 + t108;
t82 = t266 + t88;
t77 = t293 * t116 + t288 * t315;
t73 = t114 * t293 - t403;
t72 = t114 * t288 + t381;
t53 = -mrSges(6,2) * t270 + mrSges(6,3) * t64;
t52 = mrSges(6,1) * t270 - mrSges(6,3) * t63;
t30 = t103 * t287 + t292 * t67;
t29 = t103 * t292 - t287 * t67;
t28 = t287 * t82 + t292 * t77;
t27 = -t287 * t77 + t292 * t82;
t26 = t287 * t88 + t292 * t73;
t25 = -t287 * t73 + t292 * t88;
t24 = pkin(5) * t303 - pkin(10) * t84 + t149;
t20 = qJD(5) * t325 + t288 * t302 + t293 * t89;
t5 = -qJD(6) * t42 - t20 * t287 + t24 * t292;
t4 = qJD(6) * t41 + t20 * t292 + t24 * t287;
t1 = [t467 * t298 - (Ifges(7,3) * t459 + Ifges(7,6) * t461 + Ifges(7,5) * t462 + t361 / 0.2e1 - Ifges(6,4) * t63 - Ifges(6,2) * t64 - Ifges(6,6) * t270 + t105 * mrSges(6,1) - t14 * mrSges(6,3) + t468) * t322 - (mrSges(5,2) * t177 + Ifges(5,1) * t129 + Ifges(5,4) * t130 + Ifges(5,5) * t280) * t321 + (t340 + 0.2e1 * mrSges(3,3)) * t246 + m(4) * t345 + t471 * mrSges(5,3) + (-mrSges(5,1) * t177 + Ifges(5,4) * t129 + Ifges(5,2) * t130 + Ifges(5,6) * t280) * t219 + (Ifges(4,1) * t224 + Ifges(4,4) * t225) * t501 + t310 * t347 + t84 * t352 + (-Ifges(7,1) * t313 - Ifges(7,4) * t314 + Ifges(7,5) * t303) * t455 + (Ifges(5,1) * t166 - Ifges(5,4) * t165) * t449 + m(7) * (t2 * t42 + t22 * t5 + t23 * t4 + t3 * t41) + m(6) * (t105 * t182 + t14 * t93 + t149 * t169 + t20 * t68) + t298 * t472 + (-t2 * t400 + t22 * t313 - t23 * t314 - t3 * t399) * mrSges(7,3) + t482 * mrSges(4,3) - t314 * t58 / 0.2e1 + m(5) * (t108 * t135 + t109 * t134 + t167 * t71 + t168 * t70 + t177 * t256 + t239 * t247) + m(3) * (-pkin(1) * t269 + t345) - t204 * t373 / 0.2e1 - t205 * t374 / 0.2e1 + t65 * (mrSges(7,1) * t314 - mrSges(7,2) * t313) + qJD(3) ^ 2 * t332 / 0.2e1 + t225 * t334 / 0.2e1 + t224 * t336 / 0.2e1 - pkin(1) * t365 + t9 * t399 / 0.2e1 - t8 * t400 / 0.2e1 - (-m(6) * t15 - t473 - t52) * t325 - t84 * t434 + (t105 * mrSges(6,2) - t15 * mrSges(6,3) + Ifges(6,1) * t63 + Ifges(6,4) * t64 + Ifges(6,5) * t270 + t12 * t337 + t331 * t459 + t333 * t461 + t335 * t462 + t350 * t59) * t164 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + qJDD(3) * (Ifges(4,5) * t295 - Ifges(4,6) * t290) + t281 * (Ifges(5,5) * t166 - Ifges(5,6) * t165) / 0.2e1 + t269 * mrSges(3,2) + t247 * t158 + t256 * (-mrSges(5,1) * t130 + mrSges(5,2) * t129) + t239 * (mrSges(5,1) * t165 + mrSges(5,2) * t166) + qJ(2) * (-mrSges(4,1) * t225 + mrSges(4,2) * t224) + qJD(2) * t221 + t206 * (Ifges(5,4) * t166 - Ifges(5,2) * t165) / 0.2e1 + t108 * t175 + t109 * t176 + t182 * (-mrSges(6,1) * t64 + mrSges(6,2) * t63) + t166 * t146 / 0.2e1 + t167 * t119 + t168 * t120 - t165 * t145 / 0.2e1 + t149 * t102 + t20 * t131 + t93 * t53 + t84 * t97 / 0.2e1 + t4 * t90 + t5 * t91 + t41 * t18 + t42 * t19 + (Ifges(6,1) * t84 - Ifges(6,4) * t303) * t451 + t148 * (-Ifges(7,5) * t313 - Ifges(7,6) * t314 + Ifges(7,3) * t303) / 0.2e1 + t125 * (-Ifges(7,4) * t313 - Ifges(7,2) * t314 + Ifges(7,6) * t303) / 0.2e1 + t271 * (Ifges(6,5) * t84 - Ifges(6,6) * t303) / 0.2e1 + t169 * (mrSges(6,1) * t303 + mrSges(6,2) * t84) - t303 * t96 / 0.2e1 + t303 * t57 / 0.2e1 - t303 * t433 + t343 * (Ifges(6,4) * t84 - Ifges(6,2) * t303) / 0.2e1 + t477 * t364 - t303 * t490 + (-t186 * mrSges(7,1) - t185 * mrSges(7,2) + t493 * (t291 * t227 - t284 * t296 + t377) + (-m(3) + t494) * t377 + (-m(4) * pkin(7) + m(5) * t297 - t469) * t296 + t465 * t291) * g(2) + (-m(3) * t349 - t188 * mrSges(7,1) + t187 * mrSges(7,2) + t493 * (t296 * t227 + t291 * t284 + t349) + t494 * t276 + (-m(4) * t298 - m(5) * (-pkin(1) + t297) + t469) * t291 + t465 * t296) * g(1) + t303 * t491 - t290 * (Ifges(4,4) * t224 + Ifges(4,2) * t225) / 0.2e1 + t505 * (qJD(5) * t93 + t288 * t89 - t293 * t302); t365 - t321 * t119 - t219 * t120 + t165 * t175 + t166 * t176 - t489 * t84 - (t16 - t52) * t164 + t472 + t476 * t303 + (qJ(2) * t458 - mrSges(3,3)) * t463 - (qJD(6) * t328 + t481 + t53) * t322 + m(6) * (-t14 * t322 + t15 * t164 + t303 * t68 + t67 * t84) - m(5) * t471 + m(3) * t269 + m(7) * (-t12 * t164 + t303 * t329 - t304 * t322 - t65 * t84) + (-m(5) * t239 - m(6) * t169 - m(7) * t330 - t102 - t158 - t221 + t328) * qJD(1) + t478 * (m(5) - t458 - t493) + t467; t464 * (pkin(10) + t201) + (t12 * t191 - t22 * t27 - t23 * t28) * m(7) + (t14 * t201 + t15 * t200 - t169 * t174 - t68 * t77) * m(6) + t332 * t347 + t119 * t445 - t237 * t226 - t158 * t266 - m(5) * (t134 * t140 + t135 * t141 + t239 * t266) + (t175 * t370 + m(5) * (t289 * t70 + t294 * t71 + (-t134 * t289 + t135 * t294) * qJD(4)) - t176 * t371 + t289 * t120) * pkin(3) + t478 * t341 + ((-m(7) * (-t228 + t342) + t474 + t483) * t296 + t475) * g(2) + (-t242 + (-m(7) * t228 + t480 - t483) * t291 + t466) * g(1) + t300 + t204 * t375 / 0.2e1 + t205 * t376 / 0.2e1 + Ifges(4,6) * t225 + Ifges(4,5) * t224 + t200 * t52 + t201 * t53 + t191 * t16 - t141 * t175 - t140 * t176 + t170 * mrSges(4,1) - t171 * mrSges(4,2) - t174 * t102 + Ifges(4,3) * qJDD(3) - t77 * t131 - t28 * t90 - t27 * t91 + (-m(7) * (-t278 + t318) + m(6) * t227 + t484 + t500) * g(3) + (m(6) * t68 + m(7) * t329 + t476) * (t263 * t368 + (-t289 * t369 + (t293 * t294 - t386) * qJD(4)) * pkin(3)) + (t310 / 0.2e1 - t477) * t463 + t238 * t392 + t505 * (-t116 * t288 + t293 * t315 + t263 * t369 + (t289 * t368 + (t288 * t294 + t385) * qJD(4)) * pkin(3)); -t102 * t443 - t73 * t131 - t134 * t175 + t135 * t176 + t262 * t16 - t25 * t91 - t26 * t90 + t52 * t439 + t53 * t440 + t300 + (-t22 * t25 - t23 * t26 - t65 * t72 + t12 * t262 + (t288 * t65 + t293 * t329) * qJD(5) * pkin(4)) * m(7) + ((t14 * t288 + t15 * t293 + (-t288 * t67 + t293 * t68) * qJD(5)) * pkin(4) - t169 * t443 + t67 * t72 - t68 * t73) * m(6) + t476 * pkin(4) * t368 + (m(6) * t442 - m(7) * t318 - t339 + t484) * g(3) + ((-m(7) * (t342 - t441) + t362 + t474) * t296 + t475) * g(2) + (-m(7) * (pkin(4) * t389 + t378) - t242 + (-t362 + t480) * t291 + t485) * g(1) + t464 * (pkin(10) + t440) + t489 * (pkin(4) * t369 - t72); -m(7) * (t22 * t29 + t23 * t30) - t30 * t90 - t29 * t91 + (t422 - t131) * t67 + (t484 - t492) * g(3) + ((-m(7) * t342 + t479) * t296 + t486) * g(2) + (mrSges(6,2) * t391 + t466) * g(1) + t473 * pkin(5) + (t421 + t470) * t68 + t504 * t323 + t495 * t343 + t464 * pkin(10) + t497; -t65 * (mrSges(7,1) * t126 + mrSges(7,2) * t125) + (Ifges(7,1) * t125 - t413) * t456 + t58 * t455 + (Ifges(7,5) * t125 - Ifges(7,6) * t126) * t454 - t22 * t90 + t23 * t91 - g(1) * (mrSges(7,1) * t185 - mrSges(7,2) * t186) - g(2) * (mrSges(7,1) * t187 + mrSges(7,2) * t188) + g(3) * t337 * t260 + (t125 * t22 + t126 * t23) * mrSges(7,3) + t361 + (-Ifges(7,2) * t126 + t124 + t59) * t457 + t468;];
tau  = t1;
