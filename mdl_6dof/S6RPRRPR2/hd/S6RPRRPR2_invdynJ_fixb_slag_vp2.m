% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:25
% EndTime: 2019-03-09 05:01:05
% DurationCPUTime: 28.75s
% Computational Cost: add. (12458->829), mult. (26632->1115), div. (0->0), fcn. (18323->18), ass. (0->358)
t301 = sin(pkin(10));
t281 = pkin(1) * t301 + pkin(7);
t261 = t281 * qJD(1);
t307 = sin(qJ(3));
t311 = cos(qJ(3));
t209 = qJD(2) * t311 - t307 * t261;
t351 = pkin(3) * t307 - pkin(8) * t311;
t252 = t351 * qJD(1);
t306 = sin(qJ(4));
t310 = cos(qJ(4));
t151 = -t209 * t306 + t310 * t252;
t392 = t310 * t311;
t331 = pkin(4) * t307 - qJ(5) * t392;
t304 = -qJ(5) - pkin(8);
t355 = qJD(4) * t304;
t522 = -qJD(1) * t331 - qJD(5) * t306 + t310 * t355 - t151;
t152 = t310 * t209 + t306 * t252;
t386 = qJD(1) * t311;
t368 = t306 * t386;
t379 = qJD(5) * t310;
t521 = -qJ(5) * t368 - t306 * t355 + t152 - t379;
t300 = sin(pkin(11));
t302 = cos(pkin(11));
t240 = t300 * t310 + t302 * t306;
t327 = t240 * t311;
t193 = qJD(1) * t327;
t222 = t240 * qJD(4);
t520 = t193 - t222;
t334 = t300 * t306 - t302 * t310;
t326 = t334 * t311;
t194 = qJD(1) * t326;
t223 = t334 * qJD(4);
t519 = t194 - t223;
t495 = t521 * t300 + t302 * t522;
t494 = t300 * t522 - t521 * t302;
t518 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t298 = qJ(4) + pkin(11);
t291 = cos(t298);
t425 = pkin(4) * t310;
t257 = pkin(5) * t291 + t425;
t250 = pkin(3) + t257;
t294 = qJ(6) + t298;
t279 = sin(t294);
t280 = cos(t294);
t287 = pkin(3) + t425;
t289 = sin(t298);
t347 = -mrSges(5,1) * t310 + mrSges(5,2) * t306;
t517 = m(5) * pkin(3) + m(6) * t287 + m(7) * t250 + mrSges(6,1) * t291 + mrSges(7,1) * t280 - mrSges(6,2) * t289 - mrSges(7,2) * t279 - t347;
t516 = Ifges(5,3) + Ifges(6,3);
t387 = qJD(1) * t307;
t515 = -pkin(5) * t387 - pkin(9) * t519 + t495;
t514 = pkin(9) * t520 + t494;
t297 = -pkin(9) + t304;
t512 = -m(5) * pkin(8) + m(6) * t304 + m(7) * t297 - t518;
t378 = qJD(1) * qJD(3);
t254 = qJDD(1) * t311 - t307 * t378;
t511 = t254 / 0.2e1;
t510 = t378 / 0.2e1;
t384 = qJD(3) * t310;
t248 = -t306 * t387 + t384;
t249 = qJD(3) * t306 + t310 * t387;
t165 = t248 * t300 + t249 * t302;
t509 = pkin(9) * t165;
t210 = t307 * qJD(2) + t311 * t261;
t382 = qJD(4) * t306;
t486 = -t210 + (-t368 + t382) * pkin(4);
t508 = qJD(2) * qJD(3) + t281 * qJDD(1);
t380 = qJD(4) * t310;
t383 = qJD(3) * t311;
t321 = t306 * t383 + t307 * t380;
t299 = qJ(1) + pkin(10);
t290 = sin(t299);
t292 = cos(t299);
t507 = g(1) * t292 + g(2) * t290;
t451 = m(6) + m(7);
t477 = -t451 - m(4) - m(5);
t385 = qJD(3) * t307;
t305 = sin(qJ(6));
t309 = cos(qJ(6));
t353 = t302 * t248 - t249 * t300;
t506 = -t165 * t305 + t309 * t353;
t98 = t165 * t309 + t305 * t353;
t255 = qJDD(1) * t307 + t311 * t378;
t158 = qJD(4) * t248 + qJDD(3) * t306 + t255 * t310;
t238 = qJDD(4) - t254;
t192 = qJD(3) * pkin(8) + t210;
t352 = pkin(3) * t311 + pkin(8) * t307;
t333 = -pkin(2) - t352;
t303 = cos(pkin(10));
t431 = pkin(1) * t303;
t234 = t333 - t431;
t195 = t234 * qJD(1);
t124 = t192 * t310 + t195 * t306;
t145 = t307 * qJDD(2) - t261 * t385 + t311 * t508;
t133 = qJDD(3) * pkin(8) + t145;
t283 = -pkin(2) - t431;
t260 = t283 * qJDD(1);
t157 = -pkin(3) * t254 - pkin(8) * t255 + t260;
t54 = -qJD(4) * t124 - t133 * t306 + t310 * t157;
t34 = pkin(4) * t238 - qJ(5) * t158 - qJD(5) * t249 + t54;
t159 = -qJD(4) * t249 + qJDD(3) * t310 - t255 * t306;
t53 = t310 * t133 + t306 * t157 - t192 * t382 + t195 * t380;
t36 = qJ(5) * t159 + qJD(5) * t248 + t53;
t12 = t300 * t34 + t302 * t36;
t87 = -t158 * t300 + t159 * t302;
t10 = pkin(9) * t87 + t12;
t273 = qJD(4) - t386;
t108 = qJ(5) * t248 + t124;
t103 = t300 * t108;
t123 = -t192 * t306 + t310 * t195;
t107 = -qJ(5) * t249 + t123;
t99 = pkin(4) * t273 + t107;
t51 = t302 * t99 - t103;
t37 = pkin(5) * t273 - t509 + t51;
t499 = pkin(9) * t353;
t399 = t302 * t108;
t52 = t300 * t99 + t399;
t40 = t52 + t499;
t13 = -t305 * t40 + t309 * t37;
t11 = -t300 * t36 + t302 * t34;
t88 = t158 * t302 + t159 * t300;
t8 = pkin(5) * t238 - pkin(9) * t88 + t11;
t2 = qJD(6) * t13 + t10 * t309 + t305 * t8;
t14 = t305 * t37 + t309 * t40;
t3 = -qJD(6) * t14 - t10 * t305 + t309 * t8;
t505 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t504 = -t54 * mrSges(5,1) - t11 * mrSges(6,1) + t53 * mrSges(5,2) + t12 * mrSges(6,2);
t418 = Ifges(4,4) * t307;
t342 = t311 * Ifges(4,2) + t418;
t503 = t14 * mrSges(7,2) + t52 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t342 / 0.2e1 - t13 * mrSges(7,1) - t51 * mrSges(6,1);
t471 = m(6) * pkin(4);
t24 = qJD(6) * t506 + t305 * t87 + t309 * t88;
t470 = t24 / 0.2e1;
t25 = -qJD(6) * t98 - t305 * t88 + t309 * t87;
t469 = t25 / 0.2e1;
t457 = t87 / 0.2e1;
t456 = t88 / 0.2e1;
t450 = t158 / 0.2e1;
t449 = t159 / 0.2e1;
t232 = qJDD(6) + t238;
t444 = t232 / 0.2e1;
t443 = t238 / 0.2e1;
t48 = -t87 * mrSges(6,1) + t88 * mrSges(6,2);
t9 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t500 = -t48 - t9;
t498 = -qJD(1) / 0.2e1;
t265 = t304 * t306;
t267 = t304 * t310;
t174 = t302 * t265 + t267 * t300;
t140 = -pkin(9) * t240 + t174;
t175 = t300 * t265 - t302 * t267;
t141 = -pkin(9) * t334 + t175;
t70 = t140 * t309 - t141 * t305;
t497 = qJD(6) * t70 + t305 * t515 + t309 * t514;
t71 = t140 * t305 + t141 * t309;
t496 = -qJD(6) * t71 - t305 * t514 + t309 * t515;
t427 = pkin(4) * t302;
t282 = pkin(5) + t427;
t428 = pkin(4) * t300;
t215 = t282 * t309 - t305 * t428;
t57 = -t107 * t300 - t399;
t43 = t57 - t499;
t58 = t302 * t107 - t103;
t44 = t58 - t509;
t493 = t215 * qJD(6) - t305 * t43 - t309 * t44;
t216 = t282 * t305 + t309 * t428;
t492 = -t216 * qJD(6) + t305 * t44 - t309 * t43;
t491 = t471 + mrSges(5,1);
t94 = -mrSges(5,1) * t159 + mrSges(5,2) * t158;
t490 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t255 + t94;
t488 = t307 * t507;
t485 = -pkin(5) * t520 + t486;
t251 = t281 * t392;
t167 = t306 * t234 + t251;
t372 = mrSges(4,3) * t387;
t484 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t248 + mrSges(5,2) * t249 + t372;
t483 = Ifges(5,5) * t158 + Ifges(6,5) * t88 + Ifges(5,6) * t159 + Ifges(6,6) * t87 + t238 * t516;
t235 = Ifges(5,4) * t248;
t149 = t249 * Ifges(5,1) + t273 * Ifges(5,5) + t235;
t288 = Ifges(4,4) * t386;
t482 = Ifges(4,1) * t387 + Ifges(4,5) * qJD(3) + t310 * t149 + t288;
t115 = mrSges(5,1) * t238 - mrSges(5,3) * t158;
t116 = -mrSges(5,2) * t238 + mrSges(5,3) * t159;
t481 = -t306 * t115 + t310 * t116;
t146 = qJDD(2) * t311 - t261 * t383 - t307 * t508;
t480 = t145 * t311 - t146 * t307;
t479 = -t306 * t54 + t310 * t53;
t269 = qJD(6) + t273;
t478 = t249 * Ifges(5,5) + t165 * Ifges(6,5) + t98 * Ifges(7,5) + t248 * Ifges(5,6) + Ifges(6,6) * t353 + t506 * Ifges(7,6) + t269 * Ifges(7,3) + t273 * t516;
t426 = pkin(4) * t306;
t256 = pkin(5) * t289 + t426;
t475 = -m(7) * t256 + mrSges(3,2) - mrSges(4,3);
t349 = t311 * mrSges(4,1) - mrSges(4,2) * t307;
t474 = t307 * t518 + mrSges(3,1) + t349;
t473 = Ifges(7,4) * t470 + Ifges(7,2) * t469 + Ifges(7,6) * t444;
t472 = Ifges(7,1) * t470 + Ifges(7,4) * t469 + Ifges(7,5) * t444;
t468 = Ifges(6,4) * t456 + Ifges(6,2) * t457 + Ifges(6,6) * t443;
t467 = Ifges(6,1) * t456 + Ifges(6,4) * t457 + Ifges(6,5) * t443;
t432 = Ifges(7,4) * t98;
t46 = Ifges(7,2) * t506 + Ifges(7,6) * t269 + t432;
t466 = -t46 / 0.2e1;
t465 = t46 / 0.2e1;
t91 = Ifges(7,4) * t506;
t47 = Ifges(7,1) * t98 + Ifges(7,5) * t269 + t91;
t464 = -t47 / 0.2e1;
t463 = t47 / 0.2e1;
t462 = Ifges(5,1) * t450 + Ifges(5,4) * t449 + Ifges(5,5) * t443;
t85 = Ifges(6,4) * t165 + Ifges(6,2) * t353 + Ifges(6,6) * t273;
t461 = -t85 / 0.2e1;
t460 = t85 / 0.2e1;
t86 = Ifges(6,1) * t165 + Ifges(6,4) * t353 + Ifges(6,5) * t273;
t459 = -t86 / 0.2e1;
t458 = t86 / 0.2e1;
t455 = -t506 / 0.2e1;
t454 = t506 / 0.2e1;
t453 = -t98 / 0.2e1;
t452 = t98 / 0.2e1;
t448 = -t353 / 0.2e1;
t447 = t353 / 0.2e1;
t446 = -t165 / 0.2e1;
t445 = t165 / 0.2e1;
t441 = t249 / 0.2e1;
t440 = -t269 / 0.2e1;
t439 = t269 / 0.2e1;
t438 = -t273 / 0.2e1;
t437 = t273 / 0.2e1;
t435 = mrSges(6,3) * t51;
t434 = mrSges(6,3) * t52;
t433 = mrSges(7,3) * t14;
t308 = sin(qJ(1));
t430 = pkin(1) * t308;
t429 = pkin(4) * t249;
t422 = g(3) * t307;
t421 = t13 * mrSges(7,3);
t312 = cos(qJ(1));
t296 = t312 * pkin(1);
t253 = t351 * qJD(3);
t367 = t281 * t385;
t388 = t310 * t253 + t306 * t367;
t69 = -t307 * t379 + t331 * qJD(3) + (-t251 + (qJ(5) * t307 - t234) * t306) * qJD(4) + t388;
t389 = t234 * t380 + t306 * t253;
t395 = t307 * t310;
t79 = (-qJ(5) * qJD(4) - qJD(3) * t281) * t395 + (-qJD(5) * t307 + (-qJ(5) * qJD(3) - qJD(4) * t281) * t311) * t306 + t389;
t39 = t300 * t69 + t302 * t79;
t417 = Ifges(4,4) * t311;
t416 = Ifges(5,4) * t306;
t415 = Ifges(5,4) * t310;
t414 = t123 * mrSges(5,3);
t413 = t124 * mrSges(5,3);
t412 = t249 * Ifges(5,4);
t403 = t290 * t306;
t402 = t290 * t311;
t401 = t292 * t306;
t400 = t292 * t311;
t397 = t306 * t307;
t396 = t306 * t311;
t268 = t307 * t281;
t212 = t310 * t234;
t137 = -qJ(5) * t395 + t212 + (-t281 * t306 - pkin(4)) * t311;
t150 = -qJ(5) * t397 + t167;
t76 = t300 * t137 + t302 * t150;
t176 = t279 * t402 + t280 * t292;
t177 = t279 * t292 - t280 * t402;
t391 = -t176 * mrSges(7,1) + t177 * mrSges(7,2);
t178 = -t279 * t400 + t280 * t290;
t179 = t279 * t290 + t280 * t400;
t390 = t178 * mrSges(7,1) - t179 * mrSges(7,2);
t219 = pkin(4) * t397 + t268;
t381 = qJD(4) * t307;
t376 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t232;
t371 = mrSges(4,3) * t386;
t258 = t281 * t383;
t173 = pkin(4) * t321 + t258;
t148 = t248 * Ifges(5,2) + t273 * Ifges(5,6) + t412;
t364 = -t306 * t148 / 0.2e1;
t38 = -t300 * t79 + t302 * t69;
t75 = t302 * t137 - t150 * t300;
t348 = mrSges(4,1) * t307 + mrSges(4,2) * t311;
t346 = mrSges(5,1) * t306 + mrSges(5,2) * t310;
t345 = -mrSges(7,1) * t279 - mrSges(7,2) * t280;
t344 = Ifges(5,1) * t310 - t416;
t343 = Ifges(5,1) * t306 + t415;
t341 = -Ifges(5,2) * t306 + t415;
t340 = Ifges(5,2) * t310 + t416;
t339 = Ifges(4,5) * t311 - Ifges(4,6) * t307;
t338 = Ifges(5,5) * t310 - Ifges(5,6) * t306;
t337 = Ifges(5,5) * t306 + Ifges(5,6) * t310;
t208 = t334 * t307;
t61 = -pkin(5) * t311 + pkin(9) * t208 + t75;
t207 = t240 * t307;
t62 = -pkin(9) * t207 + t76;
t28 = -t305 * t62 + t309 * t61;
t29 = t305 * t61 + t309 * t62;
t131 = -t207 * t309 + t208 * t305;
t132 = -t207 * t305 - t208 * t309;
t160 = -t240 * t305 - t309 * t334;
t161 = t240 * t309 - t305 * t334;
t336 = t250 * t311 - t297 * t307;
t335 = t287 * t311 - t304 * t307;
t332 = t376 - t505;
t204 = t290 * t310 - t292 * t396;
t202 = t290 * t396 + t292 * t310;
t191 = -qJD(3) * pkin(3) - t209;
t330 = t191 * t346;
t329 = t283 * qJD(1) * t348;
t328 = t307 * (Ifges(4,1) * t311 - t418);
t322 = -t306 * t381 + t310 * t383;
t134 = -qJDD(3) * pkin(3) - t146;
t320 = Ifges(5,5) * t307 + t311 * t344;
t319 = Ifges(5,6) * t307 + t311 * t341;
t318 = Ifges(5,3) * t307 + t311 * t338;
t155 = -pkin(4) * t248 + qJD(5) + t191;
t80 = -pkin(4) * t159 + qJDD(5) + t134;
t317 = (-t123 * t310 - t124 * t306) * qJD(4) + t479;
t264 = -qJD(3) * mrSges(4,2) + t371;
t233 = t346 * t307;
t213 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t254;
t205 = t292 * t392 + t403;
t203 = -t290 * t392 + t401;
t196 = pkin(5) * t334 - t287;
t189 = mrSges(5,1) * t273 - mrSges(5,3) * t249;
t188 = -mrSges(5,2) * t273 + mrSges(5,3) * t248;
t185 = t289 * t290 + t291 * t400;
t184 = -t289 * t400 + t290 * t291;
t183 = t289 * t292 - t291 * t402;
t182 = t289 * t402 + t291 * t292;
t166 = -t281 * t396 + t212;
t156 = pkin(5) * t207 + t219;
t136 = -qJD(3) * t326 - t222 * t307;
t135 = -qJD(3) * t327 + t334 * t381;
t130 = mrSges(6,1) * t273 - mrSges(6,3) * t165;
t129 = -mrSges(6,2) * t273 + mrSges(6,3) * t353;
t120 = -t193 * t305 - t194 * t309;
t119 = -t193 * t309 + t194 * t305;
t118 = pkin(5) * t165 + t429;
t110 = -qJD(4) * t167 + t388;
t109 = (-t307 * t384 - t311 * t382) * t281 + t389;
t102 = -mrSges(6,1) * t353 + mrSges(6,2) * t165;
t101 = -pkin(5) * t135 + t173;
t100 = -pkin(5) * t353 + t155;
t93 = -qJD(6) * t161 - t222 * t309 + t223 * t305;
t92 = qJD(6) * t160 - t222 * t305 - t223 * t309;
t78 = mrSges(7,1) * t269 - mrSges(7,3) * t98;
t77 = -mrSges(7,2) * t269 + mrSges(7,3) * t506;
t72 = t158 * Ifges(5,4) + t159 * Ifges(5,2) + t238 * Ifges(5,6);
t66 = mrSges(6,1) * t238 - mrSges(6,3) * t88;
t65 = -mrSges(6,2) * t238 + mrSges(6,3) * t87;
t56 = -qJD(6) * t132 + t135 * t309 - t136 * t305;
t55 = qJD(6) * t131 + t135 * t305 + t136 * t309;
t50 = -mrSges(7,1) * t506 + mrSges(7,2) * t98;
t49 = -pkin(5) * t87 + t80;
t27 = pkin(9) * t135 + t39;
t26 = pkin(5) * t385 - pkin(9) * t136 + t38;
t20 = -mrSges(7,2) * t232 + mrSges(7,3) * t25;
t19 = mrSges(7,1) * t232 - mrSges(7,3) * t24;
t5 = -qJD(6) * t29 + t26 * t309 - t27 * t305;
t4 = qJD(6) * t28 + t26 * t305 + t27 * t309;
t1 = [(Ifges(7,4) * t132 + Ifges(7,2) * t131) * t469 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t454 + (-Ifges(6,5) * t208 - Ifges(6,6) * t207 + t307 * t338) * t443 + (t307 * Ifges(4,1) + t417 / 0.2e1 + t283 * mrSges(4,2)) * t255 + m(4) * (t260 * t283 + (-t210 * t385 + t480) * t281) + (Ifges(7,1) * t132 + Ifges(7,4) * t131) * t470 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t452 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t439 + (Ifges(7,5) * t132 + Ifges(7,6) * t131) * t444 + (Ifges(6,1) * t136 + Ifges(6,4) * t135) * t445 + (-Ifges(6,1) * t208 - Ifges(6,4) * t207) * t456 + (Ifges(6,5) * t136 + Ifges(6,6) * t135 + qJD(3) * t318 - t337 * t381) * t437 + (-t13 * t55 + t131 * t2 - t132 * t3 + t14 * t56) * mrSges(7,3) + (-t123 * t322 - t124 * t321 - t395 * t54 - t397 * t53) * mrSges(5,3) + (t11 * t208 - t12 * t207 + t135 * t52 - t136 * t51) * mrSges(6,3) + t248 * (qJD(3) * t319 - t340 * t381) / 0.2e1 + t490 * t268 - (t310 * t148 + t306 * t149) * t381 / 0.2e1 + (-t209 * t383 + t480) * mrSges(4,3) + m(5) * (t109 * t124 + t110 * t123 + t166 * t54 + t167 * t53 + (t134 * t307 + t191 * t383) * t281) + (Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(3,1) * t303 - 0.2e1 * mrSges(3,2) * t301 + m(3) * (t301 ^ 2 + t303 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (t478 / 0.2e1 - t124 * mrSges(5,2) + t123 * mrSges(5,1) + Ifges(7,3) * t439 + Ifges(6,5) * t445 + Ifges(6,6) * t447 + Ifges(7,5) * t452 + Ifges(7,6) * t454 - t210 * mrSges(4,3) + Ifges(6,3) * t437 - t503) * t385 + t56 * t465 - t208 * t467 - t207 * t468 + t132 * t472 + t131 * t473 + t364 * t383 - t264 * t367 + (qJD(3) * t320 - t343 * t381) * t441 - t72 * t397 / 0.2e1 + (Ifges(6,4) * t136 + Ifges(6,2) * t135) * t447 + (-Ifges(6,4) * t208 - Ifges(6,2) * t207) * t457 + (-t403 * t471 - m(3) * t296 - mrSges(2,1) * t312 - t205 * mrSges(5,1) - t185 * mrSges(6,1) - t179 * mrSges(7,1) + mrSges(2,2) * t308 - t204 * mrSges(5,2) - t184 * mrSges(6,2) - t178 * mrSges(7,2) + t477 * (t292 * pkin(2) + t290 * pkin(7) + t296) + t475 * t290 + (-m(5) * t352 - m(6) * t335 - m(7) * t336 - t474) * t292) * g(2) + (-t401 * t471 + m(3) * t430 + mrSges(2,1) * t308 - t203 * mrSges(5,1) - t183 * mrSges(6,1) - t177 * mrSges(7,1) + mrSges(2,2) * t312 - t202 * mrSges(5,2) - t182 * mrSges(6,2) - t176 * mrSges(7,2) + t477 * (t292 * pkin(7) - t430) + t475 * t292 + (-m(7) * (-pkin(2) - t336) - m(6) * (-pkin(2) - t335) - m(5) * t333 + m(4) * pkin(2) + t474) * t290) * g(1) + t80 * (mrSges(6,1) * t207 - mrSges(6,2) * t208) - t283 * mrSges(4,1) * t254 + t307 * t341 * t449 + t307 * t344 * t450 + qJDD(3) * Ifges(4,5) * t307 + t136 * t458 + t135 * t460 + t395 * t462 + t55 * t463 + ((-Ifges(4,2) * t307 + t417) * t510 - Ifges(7,6) * t469 - Ifges(7,5) * t470 - Ifges(7,3) * t444 - Ifges(6,5) * t456 - Ifges(6,6) * t457 - t516 * t443 - t483 / 0.2e1 + Ifges(4,4) * t255 / 0.2e1 + Ifges(4,2) * t511 - t376 / 0.2e1 - Ifges(5,6) * t449 - Ifges(5,5) * t450 + (-m(4) * qJD(3) * t209 + t213) * t281 + Ifges(4,6) * qJDD(3) + t504 + t505) * t311 + t482 * t383 / 0.2e1 + t484 * t258 + t328 * t510 + (t418 + t342) * t511 + qJD(3) * t329 + t134 * t233 + t219 * t48 + t109 * t188 + t110 * t189 + t166 * t115 + t167 * t116 + t173 * t102 + t155 * (-mrSges(6,1) * t135 + mrSges(6,2) * t136) + t156 * t9 + t49 * (-mrSges(7,1) * t131 + mrSges(7,2) * t132) + t39 * t129 + t38 * t130 + t100 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t101 * t50 + t4 * t77 + t5 * t78 + t75 * t66 + t76 * t65 + qJD(3) ^ 2 * t339 / 0.2e1 + t28 * t19 + t29 * t20 + m(7) * (t100 * t101 + t13 * t5 + t14 * t4 + t156 * t49 + t2 * t29 + t28 * t3) + m(6) * (t11 * t75 + t12 * t76 + t155 * t173 + t219 * t80 + t38 * t51 + t39 * t52) - t260 * t349 + t191 * (mrSges(5,1) * t321 + mrSges(5,2) * t322); m(3) * qJDD(2) + t136 * t129 + t135 * t130 + t131 * t19 + t132 * t20 - t207 * t66 - t208 * t65 + t55 * t77 + t56 * t78 + (-m(3) + t477) * g(3) + ((t188 * t310 - t189 * t306 + t264) * qJD(3) - t490 + t500) * t311 + (t213 + (-t306 * t188 - t310 * t189) * qJD(4) + (t102 + t50 + t484) * qJD(3) + t481) * t307 + m(7) * (t100 * t385 + t13 * t56 + t131 * t3 + t132 * t2 + t14 * t55 - t311 * t49) + m(6) * (-t11 * t207 - t12 * t208 + t135 * t51 + t136 * t52 + t155 * t385 - t311 * t80) + m(4) * (t145 * t307 + t146 * t311 + (-t209 * t307 + t210 * t311) * qJD(3)) + m(5) * ((-t134 + (-t123 * t306 + t124 * t310) * qJD(3)) * t311 + (qJD(3) * t191 + t317) * t307); t507 * (t307 * t517 + t348) + (-g(3) * t517 + t507 * t512) * t311 + (t307 * t512 - t349) * g(3) + (-mrSges(6,1) * t520 + mrSges(6,2) * t519) * t155 + (t248 * t319 + t249 * t320 + t273 * t318) * t498 + t494 * t129 + (t11 * t174 + t12 * t175 + t155 * t486 - t287 * t80 + t494 * t52 + t495 * t51) * m(6) + t495 * t130 - t478 * t387 / 0.2e1 + (Ifges(7,1) * t120 + Ifges(7,4) * t119) * t453 - t330 * t386 + (Ifges(6,5) * t446 + Ifges(7,5) * t453 + Ifges(6,6) * t448 + Ifges(7,6) * t455 + Ifges(6,3) * t438 + Ifges(7,3) * t440 + t503) * t387 + (-pkin(3) * t134 - t123 * t151 - t124 * t152) * m(5) - t339 * t378 / 0.2e1 + t496 * t78 + (t100 * t485 + t13 * t496 + t14 * t497 + t196 * t49 + t2 * t71 + t3 * t70) * m(7) + t497 * t77 + (-t123 * (mrSges(5,1) * t307 - mrSges(5,3) * t392) - t124 * (-mrSges(5,2) * t307 - mrSges(5,3) * t396) - t329 + t328 * t498) * qJD(1) + t119 * t466 + t240 * t467 + (Ifges(7,4) * t161 + Ifges(7,2) * t160) * t469 + (Ifges(7,1) * t161 + Ifges(7,4) * t160) * t470 + t161 * t472 + t160 * t473 - (Ifges(6,1) * t445 + Ifges(6,4) * t447 + Ifges(6,5) * t437 - t435 + t458) * t223 + (Ifges(7,4) * t120 + Ifges(7,2) * t119) * t455 + (Ifges(7,5) * t161 + Ifges(7,6) * t160) * t444 + (Ifges(7,4) * t452 + Ifges(7,2) * t454 + Ifges(7,6) * t439 + t433 + t465) * t93 - t382 * t413 + (Ifges(7,5) * t120 + Ifges(7,6) * t119) * t440 - (Ifges(6,4) * t445 + Ifges(6,2) * t447 + Ifges(6,6) * t437 + t434 + t460) * t222 + (-Ifges(6,4) * t194 - Ifges(6,2) * t193) * t448 + (-Ifges(6,1) * t194 - Ifges(6,4) * t193) * t446 + (-t11 * t240 - t12 * t334 + t193 * t52 - t194 * t51) * mrSges(6,3) + t485 * t50 + t486 * t102 + ((-t120 + t92) * mrSges(7,2) + (t119 - t93) * mrSges(7,1)) * t100 + (t248 * t341 + t249 * t344 + t273 * t338) * qJD(4) / 0.2e1 + t340 * t449 + t343 * t450 - t194 * t459 - t193 * t461 + t306 * t462 + t120 * t464 + (Ifges(7,1) * t452 + Ifges(7,4) * t454 + Ifges(7,5) * t439 - t421 + t463) * t92 + t479 * mrSges(5,3) + (m(5) * t317 - t188 * t382 - t189 * t380 + t481) * pkin(8) - (-Ifges(4,2) * t387 + t288 + t482) * t386 / 0.2e1 + (-m(5) * t191 + t372 - t484) * t210 + t148 * t368 / 0.2e1 + (t330 + t364) * qJD(4) + t310 * t72 / 0.2e1 - t287 * t48 + (-t119 * t14 + t120 * t13 + t160 * t2 - t161 * t3) * mrSges(7,3) + Ifges(4,5) * t255 + Ifges(4,6) * t254 + t196 * t9 - t152 * t188 - t151 * t189 + t174 * t66 + t175 * t65 + t49 * (-mrSges(7,1) * t160 + mrSges(7,2) * t161) - t145 * mrSges(4,2) + t146 * mrSges(4,1) - pkin(3) * t94 + t134 * t347 + t70 * t19 + t71 * t20 + Ifges(4,3) * qJDD(3) - t334 * t468 + (Ifges(6,5) * t240 - Ifges(6,6) * t334 + t337) * t443 + (Ifges(6,1) * t240 - Ifges(6,4) * t334) * t456 + (Ifges(6,4) * t240 - Ifges(6,2) * t334) * t457 + t80 * (mrSges(6,1) * t334 + mrSges(6,2) * t240) + (-Ifges(6,5) * t194 - Ifges(6,6) * t193) * t438 + (-t264 + t371) * t209 + (t149 / 0.2e1 - t414) * t380; -(mrSges(7,1) * t100 + Ifges(7,4) * t453 + Ifges(7,2) * t455 + Ifges(7,6) * t440 - t433 + t466) * t98 + t483 + (-mrSges(7,2) * t100 + Ifges(7,1) * t453 + Ifges(7,4) * t455 + Ifges(7,5) * t440 + t421 + t464) * t506 + (t182 * mrSges(6,1) - t183 * mrSges(6,2) - m(7) * (-t256 * t402 - t257 * t292) - t391 - mrSges(5,2) * t203 + t491 * t202) * g(2) + (-t184 * mrSges(6,1) + t185 * mrSges(6,2) - m(7) * (-t256 * t400 + t257 * t290) - t390 + mrSges(5,2) * t205 - t491 * t204) * g(1) + t492 * t78 + (-t100 * t118 + t13 * t492 + t14 * t493 + t2 * t216 + t215 * t3 + t256 * t422) * m(7) + t493 * t77 + t332 + (t11 * t302 + t12 * t300) * t471 + t249 * t413 + t248 * t414 + t66 * t427 + t65 * t428 + (Ifges(5,5) * t248 - Ifges(5,6) * t249) * t438 + t148 * t441 - t249 * (Ifges(5,1) * t248 - t412) / 0.2e1 + (m(6) * t426 + mrSges(6,1) * t289 + mrSges(6,2) * t291 - t345) * t422 - t102 * t429 - m(6) * (t155 * t429 + t51 * t57 + t52 * t58) - (-Ifges(5,2) * t249 + t149 + t235) * t248 / 0.2e1 + (-mrSges(6,2) * t155 + Ifges(6,1) * t446 + Ifges(6,4) * t448 + Ifges(6,5) * t438 + t435 + t459) * t353 - t191 * (mrSges(5,1) * t249 + mrSges(5,2) * t248) + g(3) * t233 + t215 * t19 + t216 * t20 - t123 * t188 + t124 * t189 - t58 * t129 - t57 * t130 - t118 * t50 - t504 - (mrSges(6,1) * t155 + Ifges(6,4) * t446 + Ifges(6,2) * t448 + Ifges(6,6) * t438 - t434 + t461) * t165; t451 * t311 * g(3) - t353 * t129 + t165 * t130 - t506 * t77 + t98 * t78 + (t13 * t98 - t14 * t506 - t488 + t49) * m(7) + (t165 * t51 - t353 * t52 - t488 + t80) * m(6) - t500; -t100 * (mrSges(7,1) * t98 + mrSges(7,2) * t506) + (Ifges(7,1) * t506 - t432) * t453 + t46 * t452 + (Ifges(7,5) * t506 - Ifges(7,6) * t98) * t440 - t13 * t77 + t14 * t78 - g(1) * t390 - g(2) * t391 - t345 * t422 + (t13 * t506 + t14 * t98) * mrSges(7,3) + t332 + (-Ifges(7,2) * t98 + t47 + t91) * t455;];
tau  = t1;
