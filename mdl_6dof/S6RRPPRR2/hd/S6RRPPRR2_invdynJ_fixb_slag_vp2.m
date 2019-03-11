% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:51
% EndTime: 2019-03-09 08:50:37
% DurationCPUTime: 29.18s
% Computational Cost: add. (19462->856), mult. (44846->1121), div. (0->0), fcn. (34707->18), ass. (0->376)
t310 = sin(pkin(11));
t312 = cos(pkin(11));
t348 = -mrSges(5,1) * t312 + mrSges(5,2) * t310;
t526 = -m(5) * pkin(3) - mrSges(4,1) + t348;
t311 = sin(pkin(10));
t321 = cos(qJ(2));
t404 = cos(pkin(10));
t356 = t404 * t321;
t317 = sin(qJ(2));
t383 = qJD(1) * t317;
t252 = -qJD(1) * t356 + t311 * t383;
t316 = sin(qJ(5));
t320 = cos(qJ(5));
t337 = t310 * t316 - t312 * t320;
t173 = t337 * t252;
t256 = t337 * qJD(5);
t495 = -t256 - t173;
t272 = t310 * t320 + t312 * t316;
t172 = t272 * t252;
t257 = t272 * qJD(5);
t494 = -t257 - t172;
t511 = m(6) + m(7);
t374 = m(5) + t511;
t525 = m(4) + t374;
t271 = t311 * t321 + t317 * t404;
t254 = t271 * qJD(1);
t368 = pkin(2) * t383;
t177 = pkin(3) * t254 + qJ(4) * t252 + t368;
t313 = -qJ(3) - pkin(7);
t286 = t313 * t321;
t276 = qJD(1) * t286;
t258 = t311 * t276;
t284 = t313 * t317;
t275 = qJD(1) * t284;
t212 = t275 * t404 + t258;
t130 = t310 * t177 + t312 * t212;
t400 = t252 * t310;
t113 = pkin(8) * t400 + t130;
t435 = pkin(2) * t311;
t290 = qJ(4) + t435;
t424 = pkin(8) + t290;
t263 = t424 * t310;
t264 = t424 * t312;
t375 = qJD(5) * t320;
t377 = qJD(4) * t312;
t378 = qJD(4) * t310;
t129 = t312 * t177 - t212 * t310;
t399 = t252 * t312;
t97 = pkin(4) * t254 + pkin(8) * t399 + t129;
t504 = -t263 * t375 + (-t113 + t377) * t320 + (-qJD(5) * t264 - t378 - t97) * t316;
t192 = -t316 * t263 + t320 * t264;
t503 = -t272 * qJD(4) - qJD(5) * t192 + t113 * t316 - t320 * t97;
t524 = mrSges(4,2) - mrSges(6,3) - mrSges(7,3);
t225 = qJD(2) * t310 + t254 * t312;
t352 = t312 * qJD(2) - t254 * t310;
t155 = t225 * t320 + t316 * t352;
t315 = sin(qJ(6));
t319 = cos(qJ(6));
t513 = -t225 * t316 + t320 * t352;
t523 = -t155 * t315 + t319 * t513;
t92 = t155 * t319 + t315 * t513;
t522 = pkin(9) * t494 + t504;
t521 = -pkin(5) * t254 - pkin(9) * t495 + t503;
t305 = t312 * pkin(4);
t294 = t305 + pkin(3);
t308 = pkin(11) + qJ(5);
t302 = cos(t308);
t274 = pkin(5) * t302 + t294;
t304 = qJ(6) + t308;
t292 = sin(t304);
t293 = cos(t304);
t300 = sin(t308);
t518 = -m(6) * t294 - m(7) * t274 - mrSges(6,1) * t302 - mrSges(7,1) * t293 + mrSges(6,2) * t300 + mrSges(7,2) * t292;
t285 = -mrSges(3,1) * t321 + mrSges(3,2) * t317;
t309 = qJ(2) + pkin(10);
t301 = sin(t309);
t303 = cos(t309);
t517 = t301 * t524 + t303 * t526 + t285;
t516 = m(5) * qJ(4) + mrSges(5,3);
t318 = sin(qJ(1));
t322 = cos(qJ(1));
t515 = g(1) * t322 + g(2) * t318;
t514 = -m(3) * pkin(1) - mrSges(2,1) + t517;
t373 = qJD(1) * qJD(2);
t360 = t317 * t373;
t372 = qJDD(1) * t321;
t277 = -t360 + t372;
t278 = qJDD(1) * t317 + t321 * t373;
t214 = t311 * t277 + t278 * t404;
t180 = qJDD(2) * t312 - t214 * t310;
t181 = qJDD(2) * t310 + t214 * t312;
t70 = qJD(5) * t513 + t180 * t316 + t181 * t320;
t71 = -qJD(5) * t155 + t180 * t320 - t181 * t316;
t30 = qJD(6) * t523 + t315 * t71 + t319 * t70;
t475 = t30 / 0.2e1;
t31 = -qJD(6) * t92 - t315 * t70 + t319 * t71;
t474 = t31 / 0.2e1;
t467 = t70 / 0.2e1;
t466 = t71 / 0.2e1;
t454 = t180 / 0.2e1;
t453 = t181 / 0.2e1;
t213 = -t404 * t277 + t278 * t311;
t210 = qJDD(5) + t213;
t205 = qJDD(6) + t210;
t452 = t205 / 0.2e1;
t451 = t210 / 0.2e1;
t450 = t213 / 0.2e1;
t510 = t277 / 0.2e1;
t425 = qJD(2) / 0.2e1;
t191 = -t320 * t263 - t264 * t316;
t162 = -pkin(9) * t272 + t191;
t163 = -pkin(9) * t337 + t192;
t102 = t162 * t315 + t163 * t319;
t509 = -qJD(6) * t102 - t315 * t522 + t319 * t521;
t101 = t162 * t319 - t163 * t315;
t508 = qJD(6) * t101 + t315 * t521 + t319 * t522;
t306 = t321 * pkin(2);
t296 = t306 + pkin(1);
t507 = t225 * Ifges(5,5);
t506 = t352 * Ifges(5,6);
t405 = qJDD(2) / 0.2e1;
t476 = m(7) * pkin(5);
t505 = mrSges(6,1) + t476;
t502 = Ifges(4,5) * qJD(2);
t501 = Ifges(4,6) * qJD(2);
t115 = t172 * t319 - t173 * t315;
t207 = t272 * t319 - t315 * t337;
t139 = -qJD(6) * t207 + t256 * t315 - t257 * t319;
t498 = -t115 + t139;
t116 = t172 * t315 + t173 * t319;
t206 = -t272 * t315 - t319 * t337;
t138 = qJD(6) * t206 - t256 * t319 - t257 * t315;
t497 = -t116 + t138;
t333 = -t311 * t317 + t356;
t201 = -pkin(3) * t333 - qJ(4) * t271 - t296;
t217 = t311 * t284 - t286 * t404;
t145 = t312 * t201 - t217 * t310;
t393 = t271 * t312;
t111 = -pkin(4) * t333 - pkin(8) * t393 + t145;
t146 = t310 * t201 + t312 * t217;
t394 = t271 * t310;
t120 = -pkin(8) * t394 + t146;
t61 = t316 * t111 + t320 * t120;
t496 = qJD(2) * mrSges(4,1) + mrSges(5,1) * t352 - mrSges(5,2) * t225 - mrSges(4,3) * t254;
t297 = pkin(7) * t372;
t267 = -pkin(7) * t360 + t297;
t268 = t278 * pkin(7);
t493 = t267 * t321 + t268 * t317;
t402 = qJDD(1) * pkin(1);
t243 = -pkin(2) * t277 + qJDD(3) - t402;
t112 = pkin(3) * t213 - qJ(4) * t214 - qJD(4) * t254 + t243;
t380 = qJD(3) * t317;
t197 = qJDD(2) * pkin(2) - qJ(3) * t278 - qJD(1) * t380 - t268;
t381 = qJD(2) * t317;
t366 = pkin(7) * t381;
t379 = qJD(3) * t321;
t208 = qJ(3) * t277 + t297 + (-t366 + t379) * qJD(1);
t135 = t311 * t197 + t404 * t208;
t128 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t135;
t62 = t312 * t112 - t128 * t310;
t63 = t310 * t112 + t312 * t128;
t492 = -t310 * t62 + t312 * t63;
t357 = t404 * t276;
t211 = t275 * t311 - t357;
t161 = -pkin(4) * t400 + t211;
t491 = -pkin(5) * t494 - t161;
t249 = qJD(5) + t252;
t241 = qJD(6) + t249;
t490 = t155 * Ifges(6,5) + t92 * Ifges(7,5) + Ifges(6,6) * t513 + Ifges(7,6) * t523 + t252 * Ifges(5,3) + t249 * Ifges(6,3) + t241 * Ifges(7,3) + t506 + t507;
t489 = 0.2e1 * t405;
t281 = -qJD(1) * t296 + qJD(3);
t166 = pkin(3) * t252 - qJ(4) * t254 + t281;
t266 = qJD(2) * pkin(2) + t275;
t203 = t311 * t266 - t357;
t190 = qJD(2) * qJ(4) + t203;
t117 = t312 * t166 - t190 * t310;
t82 = pkin(4) * t252 - pkin(8) * t225 + t117;
t118 = t310 * t166 + t312 * t190;
t98 = pkin(8) * t352 + t118;
t51 = -t316 * t98 + t320 * t82;
t35 = -pkin(9) * t155 + t51;
t32 = pkin(5) * t249 + t35;
t52 = t316 * t82 + t320 * t98;
t36 = pkin(9) * t513 + t52;
t406 = t319 * t36;
t14 = t315 * t32 + t406;
t202 = t266 * t404 + t258;
t183 = -qJD(2) * pkin(3) + qJD(4) - t202;
t147 = -pkin(4) * t352 + t183;
t85 = -pkin(5) * t513 + t147;
t488 = -t85 * mrSges(7,1) + mrSges(7,3) * t14;
t407 = t315 * t36;
t13 = t319 * t32 - t407;
t487 = t85 * mrSges(7,2) - mrSges(7,3) * t13;
t347 = mrSges(5,1) * t310 + mrSges(5,2) * t312;
t486 = -m(3) * pkin(7) + m(5) * t313 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - t347;
t46 = pkin(4) * t213 - pkin(8) * t181 + t62;
t57 = pkin(8) * t180 + t63;
t12 = -qJD(5) * t52 - t316 * t57 + t320 * t46;
t6 = pkin(5) * t210 - pkin(9) * t70 + t12;
t376 = qJD(5) * t316;
t11 = t316 * t46 + t320 * t57 + t82 * t375 - t376 * t98;
t7 = pkin(9) * t71 + t11;
t2 = qJD(6) * t13 + t315 * t6 + t319 * t7;
t3 = -qJD(6) * t14 - t315 * t7 + t319 * t6;
t485 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t484 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t418 = Ifges(5,4) * t312;
t343 = -Ifges(5,2) * t310 + t418;
t419 = Ifges(5,4) * t310;
t345 = Ifges(5,1) * t312 - t419;
t439 = t312 / 0.2e1;
t481 = (t225 * Ifges(5,1) + Ifges(5,4) * t352 + Ifges(5,5) * t252) * t439 - t310 * (t225 * Ifges(5,4) + Ifges(5,2) * t352 + Ifges(5,6) * t252) / 0.2e1 + t183 * t347 + t281 * mrSges(4,2) - t202 * mrSges(4,3) + t225 * t345 / 0.2e1 + t352 * t343 / 0.2e1;
t480 = t203 * mrSges(4,3) + t118 * mrSges(5,2) + t14 * mrSges(7,2) + t52 * mrSges(6,2) - t117 * mrSges(5,1) - t13 * mrSges(7,1) - t507 / 0.2e1 - t281 * mrSges(4,1) - t506 / 0.2e1 - t51 * mrSges(6,1);
t478 = Ifges(7,4) * t475 + Ifges(7,2) * t474 + Ifges(7,6) * t452;
t477 = Ifges(7,1) * t475 + Ifges(7,4) * t474 + Ifges(7,5) * t452;
t473 = Ifges(6,4) * t467 + Ifges(6,2) * t466 + Ifges(6,6) * t451;
t472 = Ifges(6,1) * t467 + Ifges(6,4) * t466 + Ifges(6,5) * t451;
t436 = Ifges(7,4) * t92;
t44 = Ifges(7,2) * t523 + Ifges(7,6) * t241 + t436;
t471 = -t44 / 0.2e1;
t470 = t44 / 0.2e1;
t86 = Ifges(7,4) * t523;
t45 = Ifges(7,1) * t92 + Ifges(7,5) * t241 + t86;
t469 = -t45 / 0.2e1;
t468 = t45 / 0.2e1;
t417 = Ifges(6,4) * t155;
t80 = Ifges(6,2) * t513 + Ifges(6,6) * t249 + t417;
t465 = t80 / 0.2e1;
t152 = Ifges(6,4) * t513;
t81 = Ifges(6,1) * t155 + Ifges(6,5) * t249 + t152;
t464 = t81 / 0.2e1;
t463 = -t523 / 0.2e1;
t462 = t523 / 0.2e1;
t461 = -t92 / 0.2e1;
t460 = t92 / 0.2e1;
t459 = Ifges(5,1) * t453 + Ifges(5,4) * t454 + Ifges(5,5) * t450;
t458 = -t513 / 0.2e1;
t457 = t513 / 0.2e1;
t456 = -t155 / 0.2e1;
t455 = t155 / 0.2e1;
t449 = -t241 / 0.2e1;
t448 = t241 / 0.2e1;
t447 = -t249 / 0.2e1;
t446 = t249 / 0.2e1;
t445 = -t252 / 0.2e1;
t444 = t252 / 0.2e1;
t441 = t254 / 0.2e1;
t433 = pkin(4) * t310;
t432 = pkin(5) * t155;
t430 = pkin(5) * t300;
t429 = pkin(7) * t321;
t426 = g(3) * t301;
t314 = -pkin(8) - qJ(4);
t423 = mrSges(6,3) * t513;
t422 = mrSges(6,3) * t155;
t421 = Ifges(3,4) * t317;
t420 = Ifges(3,4) * t321;
t413 = t254 * Ifges(4,4);
t412 = t301 * mrSges(5,3);
t403 = qJ(4) * t301;
t255 = t333 * qJD(2);
t398 = t255 * t310;
t397 = t255 * t312;
t307 = -pkin(9) + t314;
t392 = t301 * t307;
t391 = t301 * t314;
t390 = t303 * t322;
t389 = t318 * t292;
t388 = t318 * t293;
t387 = t318 * t300;
t386 = t318 * t302;
t253 = t271 * qJD(2);
t367 = pkin(2) * t381;
t157 = pkin(3) * t253 - qJ(4) * t255 - qJD(4) * t271 + t367;
t358 = qJD(2) * t313;
t250 = t317 * t358 + t379;
t251 = t321 * t358 - t380;
t179 = t250 * t404 + t311 * t251;
t105 = t310 * t157 + t312 * t179;
t227 = t293 * t322 + t303 * t389;
t228 = t292 * t322 - t303 * t388;
t385 = -t227 * mrSges(7,1) + t228 * mrSges(7,2);
t229 = -t292 * t390 + t388;
t230 = t293 * t390 + t389;
t384 = t229 * mrSges(7,1) - t230 * mrSges(7,2);
t382 = qJD(1) * t321;
t371 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t383) * t429;
t370 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t205;
t369 = Ifges(6,5) * t70 + Ifges(6,6) * t71 + Ifges(6,3) * t210;
t363 = t404 * pkin(2);
t37 = -t71 * mrSges(6,1) + t70 * mrSges(6,2);
t10 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t354 = t213 * mrSges(4,1) + t214 * mrSges(4,2);
t126 = -t180 * mrSges(5,1) + t181 * mrSges(5,2);
t60 = t320 * t111 - t120 * t316;
t104 = t312 * t157 - t179 * t310;
t178 = t250 * t311 - t404 * t251;
t216 = -t404 * t284 - t286 * t311;
t295 = -t363 - pkin(3);
t148 = pkin(4) * t398 + t178;
t174 = pkin(4) * t394 + t216;
t350 = mrSges(3,1) * t317 + mrSges(3,2) * t321;
t346 = -mrSges(7,1) * t292 - mrSges(7,2) * t293;
t344 = t321 * Ifges(3,2) + t421;
t342 = Ifges(3,5) * t321 - Ifges(3,6) * t317;
t341 = Ifges(5,5) * t312 - Ifges(5,6) * t310;
t185 = t337 * t271;
t53 = -pkin(5) * t333 + pkin(9) * t185 + t60;
t184 = t272 * t271;
t56 = -pkin(9) * t184 + t61;
t21 = -t315 * t56 + t319 * t53;
t22 = t315 * t53 + t319 * t56;
t134 = t197 * t404 - t311 * t208;
t340 = t117 * t310 - t118 * t312;
t124 = -t184 * t319 + t185 * t315;
t125 = -t184 * t315 - t185 * t319;
t339 = t274 * t303 - t392;
t338 = t294 * t303 - t391;
t280 = -t305 + t295;
t336 = t370 + t485;
t335 = pkin(1) * t350;
t233 = -t300 * t390 + t386;
t231 = t302 * t322 + t303 * t387;
t334 = t317 * (Ifges(3,1) * t321 - t421);
t76 = pkin(4) * t253 - pkin(8) * t397 + t104;
t84 = -pkin(8) * t398 + t105;
t23 = t111 * t375 - t120 * t376 + t316 * t76 + t320 * t84;
t131 = -qJDD(2) * pkin(3) + qJDD(4) - t134;
t24 = -qJD(5) * t61 - t316 * t84 + t320 * t76;
t87 = -t180 * pkin(4) + t131;
t298 = Ifges(3,4) * t382;
t287 = t322 * t296;
t283 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t382;
t279 = t430 + t433;
t262 = Ifges(3,1) * t383 + Ifges(3,5) * qJD(2) + t298;
t261 = Ifges(3,6) * qJD(2) + qJD(1) * t344;
t248 = Ifges(4,4) * t252;
t235 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t252;
t234 = t302 * t390 + t387;
t232 = t300 * t322 - t303 * t386;
t222 = pkin(5) * t337 + t280;
t198 = mrSges(4,1) * t252 + mrSges(4,2) * t254;
t194 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t214;
t193 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t213;
t187 = t254 * Ifges(4,1) - t248 + t502;
t186 = -t252 * Ifges(4,2) + t413 + t501;
t168 = mrSges(5,1) * t252 - mrSges(5,3) * t225;
t167 = -mrSges(5,2) * t252 + mrSges(5,3) * t352;
t137 = mrSges(5,1) * t213 - mrSges(5,3) * t181;
t136 = -mrSges(5,2) * t213 + mrSges(5,3) * t180;
t133 = mrSges(6,1) * t249 - t422;
t132 = -mrSges(6,2) * t249 + t423;
t127 = pkin(5) * t184 + t174;
t123 = -t255 * t272 + t256 * t271;
t122 = -t255 * t337 - t257 * t271;
t99 = t181 * Ifges(5,4) + t180 * Ifges(5,2) + t213 * Ifges(5,6);
t96 = -mrSges(6,1) * t513 + mrSges(6,2) * t155;
t75 = mrSges(7,1) * t241 - mrSges(7,3) * t92;
t74 = -mrSges(7,2) * t241 + mrSges(7,3) * t523;
t69 = -pkin(5) * t123 + t148;
t65 = -mrSges(6,2) * t210 + mrSges(6,3) * t71;
t64 = mrSges(6,1) * t210 - mrSges(6,3) * t70;
t54 = -mrSges(7,1) * t523 + mrSges(7,2) * t92;
t48 = -qJD(6) * t125 - t122 * t315 + t123 * t319;
t47 = qJD(6) * t124 + t122 * t319 + t123 * t315;
t40 = -t71 * pkin(5) + t87;
t26 = -mrSges(7,2) * t205 + mrSges(7,3) * t31;
t25 = mrSges(7,1) * t205 - mrSges(7,3) * t30;
t20 = pkin(9) * t123 + t23;
t17 = pkin(5) * t253 - pkin(9) * t122 + t24;
t16 = t319 * t35 - t407;
t15 = -t315 * t35 - t406;
t5 = -qJD(6) * t22 + t17 * t319 - t20 * t315;
t4 = qJD(6) * t21 + t17 * t315 + t20 * t319;
t1 = [-(Ifges(5,5) * t181 + Ifges(5,6) * t180 + Ifges(5,3) * t213 + t369 + t370) * t333 / 0.2e1 - (Ifges(6,3) * t451 + Ifges(7,3) * t452 + Ifges(5,5) * t453 + Ifges(5,6) * t454 - t63 * mrSges(5,2) + t62 * mrSges(5,1) - Ifges(4,4) * t214 + Ifges(4,2) * t213 + t243 * mrSges(4,1) + Ifges(6,6) * t466 + Ifges(6,5) * t467 + Ifges(7,6) * t474 + Ifges(7,5) * t475 + Ifges(5,3) * t450 - t135 * mrSges(4,3) - t489 * Ifges(4,6) + t484 + t485) * t333 - qJDD(2) * mrSges(3,2) * t429 + (-t117 * t397 - t118 * t398 - t393 * t62 - t394 * t63) * mrSges(5,3) + (t277 * t429 + t493) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t493) + t344 * t510 + (-t480 + Ifges(6,5) * t455 - t186 / 0.2e1 + Ifges(6,6) * t457 + Ifges(7,5) * t460 + Ifges(7,6) * t462 - Ifges(4,6) * t425 - Ifges(4,4) * t441 + Ifges(5,3) * t444 - Ifges(4,2) * t445 + Ifges(6,3) * t446 + Ifges(7,3) * t448 + t490 / 0.2e1) * t253 + (Ifges(6,5) * t122 + Ifges(6,6) * t123) * t446 + (-m(4) * t202 + m(5) * t183 - t496) * t178 + (Ifges(6,1) * t122 + Ifges(6,4) * t123) * t455 + Ifges(3,6) * t321 * t405 + (Ifges(7,1) * t47 + Ifges(7,4) * t48) * t460 + (Ifges(7,1) * t125 + Ifges(7,4) * t124) * t475 + (-m(4) * t134 + m(5) * t131 + t126 - t194) * t216 + (t321 * (-Ifges(3,2) * t317 + t420) + t334) * t373 / 0.2e1 + (t187 / 0.2e1 + Ifges(4,5) * t425 + Ifges(4,1) * t441 + t341 * t444 + Ifges(4,4) * t445 + t481) * t255 + m(7) * (t127 * t40 + t13 * t5 + t14 * t4 + t2 * t22 + t21 * t3 + t69 * t85) + m(6) * (t11 * t61 + t12 * t60 + t147 * t148 + t174 * t87 + t23 * t52 + t24 * t51) + (Ifges(7,4) * t47 + Ifges(7,2) * t48) * t462 + (Ifges(7,4) * t125 + Ifges(7,2) * t124) * t474 + t321 * t262 * t425 + (t243 * mrSges(4,2) - t134 * mrSges(4,3) + Ifges(4,1) * t214 - Ifges(4,4) * t213 + Ifges(4,5) * t489 + t131 * t347 + t341 * t450 + t343 * t454 + t345 * t453) * t271 + t21 * t25 + t22 * t26 + (-Ifges(6,5) * t185 - Ifges(6,6) * t184) * t451 + (-Ifges(6,1) * t185 - Ifges(6,4) * t184) * t467 + (-t11 * t184 + t12 * t185 - t122 * t51 + t123 * t52) * mrSges(6,3) + (-Ifges(6,4) * t185 - Ifges(6,2) * t184) * t466 + t87 * (mrSges(6,1) * t184 - mrSges(6,2) * t185) + (Ifges(7,5) * t125 + Ifges(7,6) * t124) * t452 + (Ifges(7,5) * t47 + Ifges(7,6) * t48) * t448 + (Ifges(6,4) * t122 + Ifges(6,2) * t123) * t457 + Ifges(2,3) * qJDD(1) + t321 * (Ifges(3,4) * t278 + Ifges(3,2) * t277 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t285 * t402 - t99 * t394 / 0.2e1 - pkin(1) * (-mrSges(3,1) * t277 + mrSges(3,2) * t278) - t261 * t381 / 0.2e1 - t335 * t373 - t283 * t366 + t179 * t235 + (-t232 * mrSges(6,1) - t228 * mrSges(7,1) - t231 * mrSges(6,2) - t227 * mrSges(7,2) + (m(4) * t313 - m(6) * (-t313 + t433) - m(7) * (t279 - t313) + t486) * t322 + (m(4) * t296 - m(6) * (-t296 - t338) - m(5) * (-t296 - t403) + t412 - m(7) * (-t296 - t339) - t514) * t318) * g(1) + t217 * t193 + t278 * t420 / 0.2e1 - t296 * t354 + t174 * t37 + t105 * t167 + t104 * t168 + t393 * t459 + t122 * t464 + t123 * t465 + t47 * t468 + t48 * t470 + m(5) * (t104 * t117 + t105 * t118 + t145 * t62 + t146 * t63) + m(4) * (t135 * t217 + t179 * t203 - t243 * t296 + t281 * t367) + (-m(5) * t287 - t234 * mrSges(6,1) - t230 * mrSges(7,1) - t233 * mrSges(6,2) - t229 * mrSges(7,2) + (-m(4) - t511) * (-t318 * t313 + t287) + (-m(6) * t433 - m(7) * t279 + t486) * t318 + (-m(6) * t338 - m(7) * t339 - t301 * t516 + t514) * t322) * g(2) + t198 * t367 + (t342 * t425 - t371) * qJD(2) + t60 * t64 + t61 * t65 + (Ifges(3,1) * t278 + Ifges(3,4) * t510 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t278) + t489 * Ifges(3,5)) * t317 + t69 * t54 + t4 * t74 + t5 * t75 + (t124 * t2 - t125 * t3 - t13 * t47 + t14 * t48) * mrSges(7,3) - t185 * t472 - t184 * t473 + t125 * t477 + t124 * t478 + t85 * (-mrSges(7,1) * t48 + mrSges(7,2) * t47) + t40 * (-mrSges(7,1) * t124 + mrSges(7,2) * t125) + t127 * t10 + t23 * t132 + t24 * t133 + t145 * t137 + t146 * t136 + t147 * (-mrSges(6,1) * t123 + mrSges(6,2) * t122) + t148 * t96; t515 * (t350 + t525 * pkin(2) * t317 + (m(6) * t314 + m(7) * t307 - t516 + t524) * t303 + (-t518 - t526) * t301) + (-t11 * t337 - t12 * t272 + t494 * t52 - t495 * t51) * mrSges(6,3) + (Ifges(6,5) * t272 - Ifges(6,6) * t337) * t451 + t87 * (mrSges(6,1) * t337 + mrSges(6,2) * t272) + (Ifges(6,4) * t272 - Ifges(6,2) * t337) * t466 + (Ifges(6,1) * t272 - Ifges(6,4) * t337) * t467 - t337 * t473 + (Ifges(7,1) * t116 + Ifges(7,4) * t115) * t461 + (-t115 * t14 + t116 * t13 + t2 * t206 - t207 * t3) * mrSges(7,3) + t491 * t54 + (-t117 * t399 - t118 * t400 + t492) * mrSges(5,3) + (-t340 * qJD(4) - t117 * t129 - t118 * t130 + t131 * t295 - t183 * t211 + t290 * t492) * m(5) + (Ifges(7,4) * t116 + Ifges(7,2) * t115) * t463 + t131 * t348 + (pkin(7) * t283 + t261 / 0.2e1) * t383 + (-t248 + t187) * t444 + (-t378 - t129) * t168 + (-t130 + t377) * t167 + (t202 * t211 - t203 * t212 - t281 * t368 + (t134 * t404 + t135 * t311) * pkin(2)) * m(4) + (-mrSges(6,1) * t494 + mrSges(6,2) * t495) * t147 + t496 * t211 + (t371 + (-t334 / 0.2e1 + t335) * qJD(1)) * qJD(1) + (t136 * t312 - t137 * t310) * t290 + t310 * t459 + (Ifges(6,5) * t173 + Ifges(6,6) * t172) * t447 + (Ifges(7,5) * t116 + Ifges(7,6) * t115) * t449 + (Ifges(7,5) * t207 + Ifges(7,6) * t206) * t452 - (-Ifges(3,2) * t383 + t262 + t298) * t382 / 0.2e1 + (Ifges(7,1) * t460 + Ifges(7,4) * t462 + Ifges(7,5) * t448 + t468 + t487) * t138 + (Ifges(7,4) * t460 + Ifges(7,2) * t462 + Ifges(7,6) * t448 + t470 + t488) * t139 + (Ifges(6,1) * t173 + Ifges(6,4) * t172) * t456 - (-Ifges(4,1) * t252 - t413 + t490) * t254 / 0.2e1 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (Ifges(5,2) * t312 + t419) * t454 + (Ifges(5,1) * t310 + t418) * t453 + (-Ifges(6,1) * t256 - Ifges(6,4) * t257) * t455 + (-Ifges(6,4) * t256 - Ifges(6,2) * t257) * t457 + (-Ifges(6,5) * t256 - Ifges(6,6) * t257) * t446 + t295 * t126 + Ifges(3,6) * t277 + Ifges(3,5) * t278 + t280 * t37 - t267 * mrSges(3,2) - t268 * mrSges(3,1) - t342 * t373 / 0.2e1 - t198 * t368 - t212 * t235 + t222 * t10 - Ifges(4,6) * t213 + Ifges(4,5) * t214 + t40 * (-mrSges(7,1) * t206 + mrSges(7,2) * t207) + t191 * t64 + t192 * t65 - t172 * t80 / 0.2e1 - t173 * t81 / 0.2e1 - t161 * t96 - t256 * t464 - t257 * t465 + t116 * t469 + (t480 + Ifges(6,5) * t456 + Ifges(6,6) * t458 + Ifges(7,5) * t461 + Ifges(7,6) * t463 - Ifges(4,2) * t444 + Ifges(5,3) * t445 + Ifges(6,3) * t447 + Ifges(7,3) * t449 + t501 / 0.2e1) * t254 + (t502 / 0.2e1 - t341 * t445 + t481) * t252 + t503 * t133 + t504 * t132 + (t11 * t192 + t12 * t191 - t147 * t161 + t280 * t87 + t503 * t51 + t504 * t52) * m(6) + (-m(4) * t306 - m(6) * (t306 - t391) - m(5) * (t306 + t403) - t412 - m(7) * (t306 - t392) + t518 * t303 + t517) * g(3) + t508 * t74 + t509 * t75 + (t101 * t3 + t102 * t2 + t13 * t509 + t14 * t508 + t222 * t40 + t491 * t85) * m(7) + (Ifges(6,4) * t173 + Ifges(6,2) * t172) * t458 + t115 * t471 + t272 * t472 + (Ifges(7,4) * t207 + Ifges(7,2) * t206) * t474 + (Ifges(7,1) * t207 + Ifges(7,4) * t206) * t475 + t207 * t477 + t206 * t478 + t194 * t363 + t193 * t435 + t99 * t439 + t186 * t441 + (Ifges(5,5) * t310 + Ifges(5,6) * t312) * t450 + t101 * t25 + t102 * t26 - t85 * (-mrSges(7,1) * t115 + mrSges(7,2) * t116) + t134 * mrSges(4,1) - t135 * mrSges(4,2); t498 * t75 + t497 * t74 + t494 * t133 + t495 * t132 + (-t54 - t96 + t496) * t254 + (t167 * t312 - t168 * t310 + t235) * t252 + t312 * t137 + t310 * t136 - t337 * t64 + t272 * t65 + t206 * t25 + t207 * t26 + t354 + (-g(1) * t318 + g(2) * t322) * t525 + (t13 * t498 + t14 * t497 + t2 * t207 + t206 * t3 - t254 * t85) * m(7) + (t11 * t272 - t12 * t337 - t147 * t254 + t494 * t51 + t495 * t52) * m(6) + (-t183 * t254 - t252 * t340 + t310 * t63 + t312 * t62) * m(5) + (t202 * t254 + t203 * t252 + t243) * m(4); -t513 * t132 + t155 * t133 - t352 * t167 + t225 * t168 - t523 * t74 + t92 * t75 + t10 + t126 + t37 + (t13 * t92 - t14 * t523 + t40) * m(7) + (t155 * t51 - t513 * t52 + t87) * m(6) + (t117 * t225 - t118 * t352 + t131) * m(5) + (t303 * g(3) - t301 * t515) * t374; (-mrSges(6,2) * t232 + t231 * t505 - t385) * g(2) + (m(7) * t430 + mrSges(6,1) * t300 + mrSges(6,2) * t302 - t346) * t426 + (mrSges(6,2) * t234 - t233 * t505 - t384) * g(1) + t336 + t369 - (Ifges(7,4) * t461 + Ifges(7,2) * t463 + Ifges(7,6) * t449 + t471 - t488) * t92 + (Ifges(7,1) * t461 + Ifges(7,4) * t463 + Ifges(7,5) * t449 + t469 - t487) * t523 + t80 * t455 + (Ifges(6,1) * t513 - t417) * t456 + (t422 + t133) * t52 + (-Ifges(6,2) * t155 + t152 + t81) * t458 + (t423 - t132) * t51 - t54 * t432 - m(7) * (t13 * t15 + t14 * t16 + t432 * t85) - t147 * (mrSges(6,1) * t155 + mrSges(6,2) * t513) + t484 - t16 * t74 - t15 * t75 + (t2 * t315 + t3 * t319 + (-t13 * t315 + t14 * t319) * qJD(6)) * t476 + (Ifges(6,5) * t513 - Ifges(6,6) * t155) * t447 + ((-t315 * t75 + t319 * t74) * qJD(6) + t319 * t25 + t315 * t26) * pkin(5); -t85 * (mrSges(7,1) * t92 + mrSges(7,2) * t523) + (Ifges(7,1) * t523 - t436) * t461 + t44 * t460 + (Ifges(7,5) * t523 - Ifges(7,6) * t92) * t449 - t13 * t74 + t14 * t75 - g(1) * t384 - g(2) * t385 - t346 * t426 + (t13 * t523 + t14 * t92) * mrSges(7,3) + t336 + (-Ifges(7,2) * t92 + t45 + t86) * t463;];
tau  = t1;
