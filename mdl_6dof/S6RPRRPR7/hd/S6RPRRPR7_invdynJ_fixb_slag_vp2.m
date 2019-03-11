% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:20:04
% EndTime: 2019-03-09 05:20:21
% DurationCPUTime: 12.64s
% Computational Cost: add. (12330->650), mult. (24676->853), div. (0->0), fcn. (16998->14), ass. (0->308)
t277 = qJD(3) + qJD(4);
t284 = sin(qJ(6));
t288 = cos(qJ(6));
t285 = sin(qJ(4));
t286 = sin(qJ(3));
t289 = cos(qJ(4));
t290 = cos(qJ(3));
t216 = -t285 * t290 - t289 * t286;
t203 = t216 * qJD(1);
t362 = qJD(1) * t290;
t363 = qJD(1) * t286;
t204 = -t285 * t363 + t289 * t362;
t282 = sin(pkin(10));
t283 = cos(pkin(10));
t317 = t203 * t282 + t283 * t204;
t126 = t277 * t288 - t284 * t317;
t127 = t277 * t284 + t288 * t317;
t392 = mrSges(6,1) * t277 + mrSges(7,1) * t126 - mrSges(7,2) * t127 - mrSges(6,3) * t317;
t293 = -pkin(1) - pkin(7);
t242 = qJD(1) * t293 + qJD(2);
t378 = t242 * t286;
t188 = -pkin(8) * t363 + t378;
t173 = t285 * t188;
t223 = t290 * t242;
t189 = -pkin(8) * t362 + t223;
t176 = qJD(3) * pkin(3) + t189;
t133 = t289 * t176 - t173;
t194 = t204 * qJ(5);
t112 = t133 - t194;
t104 = pkin(4) * t277 + t112;
t174 = t289 * t188;
t134 = t176 * t285 + t174;
t391 = qJ(5) * t203;
t113 = t134 + t391;
t389 = t113 * t282;
t56 = t104 * t283 - t389;
t54 = -pkin(5) * t277 - t56;
t486 = -m(6) * t56 + m(7) * t54 - t392;
t281 = qJ(3) + qJ(4);
t267 = pkin(10) + t281;
t252 = sin(t267);
t269 = sin(t281);
t461 = mrSges(5,2) * t269 + mrSges(6,2) * t252;
t485 = t290 / 0.2e1;
t105 = t283 * t113;
t57 = t282 * t104 + t105;
t55 = pkin(9) * t277 + t57;
t236 = pkin(3) * t363 + qJD(1) * qJ(2);
t164 = -pkin(4) * t203 + qJD(5) + t236;
t337 = t283 * t203 - t204 * t282;
t74 = -pkin(5) * t337 - pkin(9) * t317 + t164;
t20 = -t284 * t55 + t288 * t74;
t484 = t20 * mrSges(7,1);
t21 = t284 * t74 + t288 * t55;
t483 = t21 * mrSges(7,2);
t273 = t286 * pkin(3);
t270 = cos(t281);
t332 = -mrSges(5,1) * t269 - mrSges(5,2) * t270;
t333 = mrSges(4,1) * t286 + mrSges(4,2) * t290;
t482 = m(5) * t273 - t332 + t333;
t253 = cos(t267);
t481 = t252 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t253;
t408 = mrSges(7,2) * t284;
t411 = mrSges(7,1) * t288;
t480 = t408 - t411;
t353 = qJD(1) * qJD(3);
t221 = qJDD(1) * t290 - t286 * t353;
t222 = -qJDD(1) * t286 - t290 * t353;
t300 = t216 * qJD(4);
t128 = qJD(1) * t300 + t221 * t289 + t222 * t285;
t275 = qJDD(3) + qJDD(4);
t241 = qJDD(1) * t293 + qJDD(2);
t361 = qJD(3) * t286;
t165 = t290 * t241 - t242 * t361;
t141 = qJDD(3) * pkin(3) - pkin(8) * t221 + t165;
t360 = qJD(3) * t290;
t166 = t286 * t241 + t242 * t360;
t151 = pkin(8) * t222 + t166;
t60 = -qJD(4) * t134 + t289 * t141 - t151 * t285;
t29 = pkin(4) * t275 - qJ(5) * t128 - qJD(5) * t204 + t60;
t367 = t289 * t290;
t315 = t285 * t286 - t367;
t129 = qJD(1) * qJD(4) * t315 - t221 * t285 + t222 * t289;
t357 = qJD(4) * t289;
t358 = qJD(4) * t285;
t59 = t285 * t141 + t289 * t151 + t176 * t357 - t188 * t358;
t31 = qJ(5) * t129 + qJD(5) * t203 + t59;
t15 = t282 * t29 + t283 * t31;
t12 = pkin(9) * t275 + t15;
t80 = -t128 * t282 + t129 * t283;
t81 = t128 * t283 + t129 * t282;
t279 = qJD(1) * qJD(2);
t244 = qJDD(1) * qJ(2) + t279;
t172 = -pkin(3) * t222 + t244;
t98 = -pkin(4) * t129 + qJDD(5) + t172;
t17 = -pkin(5) * t80 - pkin(9) * t81 + t98;
t2 = qJD(6) * t20 + t12 * t288 + t17 * t284;
t3 = -qJD(6) * t21 - t12 * t284 + t17 * t288;
t479 = t2 * t288 - t3 * t284;
t46 = qJD(6) * t126 + t275 * t284 + t288 * t81;
t445 = t46 / 0.2e1;
t47 = -qJD(6) * t127 + t275 * t288 - t284 * t81;
t444 = t47 / 0.2e1;
t79 = qJDD(6) - t80;
t442 = t79 / 0.2e1;
t477 = -m(5) - m(4);
t476 = -m(7) - m(6);
t434 = t275 / 0.2e1;
t475 = -t317 / 0.2e1;
t474 = t317 / 0.2e1;
t473 = -t337 / 0.2e1;
t16 = -mrSges(7,1) * t47 + mrSges(7,2) * t46;
t72 = mrSges(6,1) * t275 - mrSges(6,3) * t81;
t472 = t16 - t72;
t471 = Ifges(6,4) * t317;
t470 = Ifges(6,4) * t337;
t469 = Ifges(5,5) * t315;
t468 = Ifges(5,6) * t216;
t330 = mrSges(7,1) * t284 + mrSges(7,2) * t288;
t467 = t330 * t54;
t131 = -mrSges(6,2) * t277 + mrSges(6,3) * t337;
t143 = qJD(6) - t337;
t87 = -mrSges(7,2) * t143 + mrSges(7,3) * t126;
t88 = mrSges(7,1) * t143 - mrSges(7,3) * t127;
t321 = -t284 * t88 + t288 * t87;
t309 = -t131 - t321;
t414 = pkin(8) - t293;
t226 = t414 * t286;
t227 = t414 * t290;
t163 = -t289 * t226 - t285 * t227;
t465 = pkin(5) * t252 - t253 * pkin(9);
t464 = -t480 * t252 + t481;
t160 = t277 * t367 - t285 * t361 - t286 * t358;
t161 = qJD(3) * t216 + t300;
t100 = t160 * t282 - t283 * t161;
t159 = t216 * t282 - t283 * t315;
t355 = qJD(6) * t288;
t306 = -t100 * t284 + t159 * t355;
t425 = pkin(4) * t270;
t428 = pkin(3) * t290;
t225 = t425 + t428;
t463 = m(5) * t428 + m(6) * t225;
t318 = t165 * t290 + t166 * t286;
t18 = mrSges(7,1) * t79 - mrSges(7,3) * t46;
t19 = -mrSges(7,2) * t79 + mrSges(7,3) * t47;
t462 = -t284 * t18 + t288 * t19;
t287 = sin(qJ(1));
t291 = cos(qJ(1));
t460 = -g(1) * t287 + g(2) * t291;
t334 = mrSges(4,1) * t290 - mrSges(4,2) * t286;
t404 = Ifges(4,4) * t290;
t459 = qJ(2) * t334 + (-Ifges(4,1) * t286 - t404) * t485;
t458 = 0.2e1 * t434;
t349 = t253 * t408;
t457 = (-t349 - t461) * t291;
t348 = t253 * t411;
t375 = t270 * t287;
t376 = t253 * t287;
t377 = t252 * t287;
t456 = -mrSges(5,1) * t375 - mrSges(6,1) * t376 - mrSges(7,3) * t377 - (t348 - t349) * t287;
t455 = mrSges(5,1) * t270 + mrSges(6,1) * t253 + t252 * mrSges(7,3) + t348;
t234 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t363;
t235 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t362;
t454 = (t234 * t290 - t235 * t286) * qJD(3);
t453 = -t133 * t161 - t134 * t160 + t216 * t59 + t315 * t60;
t452 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t451 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t450 = m(4) * t318 + t290 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t221) + t286 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t222);
t449 = -m(7) * t465 + mrSges(2,2) - mrSges(3,3) - t481 - t482;
t323 = t20 * t288 + t21 * t284;
t297 = -qJD(6) * t323 + t479;
t356 = qJD(6) * t284;
t448 = m(7) * t297 - t88 * t355 - t87 * t356 + t462;
t447 = qJD(1) ^ 2;
t446 = Ifges(7,1) * t445 + Ifges(7,4) * t444 + Ifges(7,5) * t442;
t441 = -m(3) - m(4);
t292 = -pkin(8) - pkin(7);
t440 = -t126 / 0.2e1;
t439 = -t127 / 0.2e1;
t438 = t127 / 0.2e1;
t437 = -t143 / 0.2e1;
t435 = t204 / 0.2e1;
t433 = -t277 / 0.2e1;
t431 = t288 / 0.2e1;
t429 = pkin(3) * t289;
t427 = pkin(4) * t204;
t426 = pkin(4) * t269;
t424 = pkin(4) * t282;
t423 = pkin(4) * t283;
t417 = t56 * mrSges(6,3);
t416 = t57 * mrSges(6,3);
t407 = mrSges(5,3) * t203;
t406 = mrSges(5,3) * t204;
t405 = Ifges(4,4) * t286;
t403 = Ifges(7,4) * t284;
t402 = Ifges(7,4) * t288;
t401 = pkin(3) * qJD(4);
t400 = t127 * Ifges(7,4);
t399 = t204 * Ifges(5,4);
t386 = t337 * t284;
t385 = t337 * t288;
t384 = t159 * t284;
t383 = t159 * t288;
t374 = t282 * t285;
t373 = t283 * t285;
t372 = t284 * t287;
t371 = t284 * t291;
t369 = t287 * t288;
t368 = t288 * t291;
t256 = qJ(2) + t273;
t139 = t289 * t189 - t173;
t365 = pkin(5) * t376 + pkin(9) * t377;
t259 = pkin(4) + t429;
t193 = pkin(3) * t373 + t282 * t259;
t364 = t291 * pkin(1) + t287 * qJ(2);
t354 = qJDD(1) * mrSges(3,2);
t245 = pkin(3) * t360 + qJD(2);
t351 = m(6) * t425;
t350 = Ifges(7,5) * t46 + Ifges(7,6) * t47 + Ifges(7,3) * t79;
t262 = pkin(3) * t362;
t123 = Ifges(7,4) * t126;
t53 = t127 * Ifges(7,1) + t143 * Ifges(7,5) + t123;
t345 = t53 * t431;
t343 = -t80 * mrSges(6,1) + t81 * mrSges(6,2);
t342 = -t356 / 0.2e1;
t272 = t291 * qJ(2);
t341 = -pkin(1) * t287 + t272;
t340 = -t353 / 0.2e1;
t338 = (t244 + t279) * qJ(2);
t138 = -t189 * t285 - t174;
t162 = t226 * t285 - t289 * t227;
t177 = -pkin(4) * t216 + t256;
t146 = pkin(4) * t160 + t245;
t336 = -pkin(5) * t253 - pkin(9) * t252;
t335 = -g(1) * t291 - g(2) * t287;
t329 = t290 * Ifges(4,1) - t405;
t328 = Ifges(7,1) * t288 - t403;
t327 = -t286 * Ifges(4,2) + t404;
t326 = -Ifges(7,2) * t284 + t402;
t325 = -Ifges(4,5) * t286 - Ifges(4,6) * t290;
t324 = Ifges(7,5) * t288 - Ifges(7,6) * t284;
t322 = -t20 * t284 + t21 * t288;
t14 = -t282 * t31 + t283 * t29;
t137 = qJ(5) * t216 + t163;
t307 = qJ(5) * t315 + t162;
t86 = t283 * t137 + t282 * t307;
t316 = -t283 * t216 - t282 * t315;
t89 = pkin(5) * t316 - pkin(9) * t159 + t177;
t33 = t284 * t89 + t288 * t86;
t32 = -t284 * t86 + t288 * t89;
t320 = -t284 * t87 - t288 * t88;
t319 = t283 * t160 + t161 * t282;
t312 = -t426 - t465;
t192 = -pkin(3) * t374 + t259 * t283;
t308 = t138 - t391;
t305 = t100 * t288 + t159 * t356;
t303 = t286 * (-Ifges(4,2) * t290 - t405);
t211 = t414 * t361;
t212 = qJD(3) * t227;
t107 = t285 * t211 - t289 * t212 + t226 * t358 - t227 * t357;
t84 = pkin(5) * t317 - pkin(9) * t337 + t427;
t108 = -qJD(4) * t163 + t289 * t211 + t212 * t285;
t296 = -qJ(5) * t161 + qJD(5) * t315 + t108;
t11 = -pkin(5) * t275 - t14;
t144 = t203 * Ifges(5,2) + t277 * Ifges(5,6) + t399;
t195 = Ifges(5,4) * t203;
t145 = t204 * Ifges(5,1) + t277 * Ifges(5,5) + t195;
t51 = t127 * Ifges(7,5) + t126 * Ifges(7,6) + t143 * Ifges(7,3);
t52 = t126 * Ifges(7,2) + t143 * Ifges(7,6) + t400;
t8 = t46 * Ifges(7,4) + t47 * Ifges(7,2) + t79 * Ifges(7,6);
t91 = Ifges(6,2) * t337 + t277 * Ifges(6,6) + t471;
t92 = Ifges(6,1) * t317 + t277 * Ifges(6,5) + t470;
t295 = (t342 + t386 / 0.2e1) * t52 + (-t164 * mrSges(6,2) + Ifges(6,1) * t475 + Ifges(6,5) * t433 + t324 * t437 + t326 * t440 + t328 * t439 + t417 - t467) * t337 + (Ifges(5,3) + Ifges(6,3)) * t275 + (-t164 * mrSges(6,1) + Ifges(7,5) * t439 - Ifges(6,2) * t473 - Ifges(6,6) * t433 + Ifges(7,6) * t440 + Ifges(7,3) * t437 + t416 + t483 - t484) * t317 + (t126 * t326 + t127 * t328 + t143 * t324) * qJD(6) / 0.2e1 - (-Ifges(5,2) * t204 + t145 + t195) * t203 / 0.2e1 + (t467 + t345) * qJD(6) - t236 * (mrSges(5,1) * t204 + mrSges(5,2) * t203) + t91 * t474 + t284 * t446 + t8 * t431 + (Ifges(5,5) * t203 - Ifges(5,6) * t204) * t433 + t144 * t435 + (Ifges(7,5) * t284 + Ifges(7,6) * t288) * t442 + (Ifges(7,2) * t288 + t403) * t444 + (Ifges(7,1) * t284 + t402) * t445 + ((-t356 + t386) * t21 + (-t355 + t385) * t20 + t479) * mrSges(7,3) + t11 * t480 + t134 * t406 + t133 * t407 + Ifges(5,5) * t128 + Ifges(5,6) * t129 + Ifges(6,6) * t80 + Ifges(6,5) * t81 + t60 * mrSges(5,1) - t59 * mrSges(5,2) + (-t471 + t51) * t475 + t14 * mrSges(6,1) - t15 * mrSges(6,2) - t204 * (Ifges(5,1) * t203 - t399) / 0.2e1 - t53 * t385 / 0.2e1 + (t470 + t92) * t473;
t276 = -qJ(5) + t292;
t266 = -pkin(1) * qJDD(1) + qJDD(2);
t255 = -pkin(5) - t423;
t224 = t273 + t426;
t218 = t333 * qJD(1);
t202 = Ifges(4,5) * qJD(3) + qJD(1) * t329;
t201 = Ifges(4,6) * qJD(3) + qJD(1) * t327;
t184 = -pkin(5) - t192;
t183 = t252 * t368 - t372;
t182 = t252 * t371 + t369;
t181 = t252 * t369 + t371;
t180 = -t252 * t372 + t368;
t171 = mrSges(5,1) * t277 - t406;
t170 = -mrSges(5,2) * t277 + t407;
t169 = t262 + t427;
t153 = -mrSges(5,1) * t203 + mrSges(5,2) * t204;
t119 = -mrSges(5,2) * t275 + mrSges(5,3) * t129;
t118 = mrSges(5,1) * t275 - mrSges(5,3) * t128;
t115 = -t194 + t139;
t97 = -mrSges(6,1) * t337 + mrSges(6,2) * t317;
t82 = t262 + t84;
t73 = -qJ(5) * t160 + qJD(5) * t216 + t107;
t71 = -mrSges(6,2) * t275 + mrSges(6,3) * t80;
t64 = t283 * t115 + t282 * t308;
t62 = t112 * t283 - t389;
t61 = t112 * t282 + t105;
t42 = pkin(5) * t319 + pkin(9) * t100 + t146;
t27 = t282 * t296 + t283 * t73;
t25 = t284 * t82 + t288 * t64;
t24 = -t284 * t64 + t288 * t82;
t23 = t284 * t84 + t288 * t62;
t22 = -t284 * t62 + t288 * t84;
t5 = -qJD(6) * t33 - t27 * t284 + t288 * t42;
t4 = qJD(6) * t32 + t27 * t288 + t284 * t42;
t1 = [t486 * (t282 * t73 - t283 * t296) - t319 * t483 - t100 * t345 + (-Ifges(6,1) * t100 - Ifges(6,4) * t319) * t474 + (Ifges(5,5) * t161 - Ifges(6,5) * t100 - Ifges(5,6) * t160 - Ifges(6,6) * t319) * t277 / 0.2e1 + t164 * (mrSges(6,1) * t319 - mrSges(6,2) * t100) + t337 * (-Ifges(6,4) * t100 - Ifges(6,2) * t319) / 0.2e1 - t100 * t92 / 0.2e1 + t100 * t417 + (mrSges(5,2) * t256 - Ifges(5,1) * t315 + Ifges(5,4) * t216) * t128 + (-mrSges(5,1) * t256 - Ifges(5,4) * t315 + Ifges(5,2) * t216) * t129 + t172 * (-mrSges(5,1) * t216 - mrSges(5,2) * t315) + (t333 + 0.2e1 * mrSges(3,3)) * t244 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t319 * t484 + (t98 * mrSges(6,2) - t14 * mrSges(6,3) + Ifges(6,1) * t81 + Ifges(6,4) * t80 + Ifges(6,5) * t458 + t11 * t330 + t324 * t442 + t326 * t444 + t328 * t445 + t53 * t342) * t159 + (Ifges(7,3) * t442 + Ifges(7,6) * t444 + Ifges(7,5) * t445 + t98 * mrSges(6,1) - Ifges(6,2) * t80 - Ifges(6,4) * t81 + t350 / 0.2e1 - t15 * mrSges(6,3) - t458 * Ifges(6,6) + t451) * t316 + t459 * t353 - t318 * mrSges(4,3) - t306 * t52 / 0.2e1 + qJDD(3) * (Ifges(4,5) * t290 - Ifges(4,6) * t286) + t303 * t340 + t126 * (-Ifges(7,4) * t305 - Ifges(7,2) * t306 + Ifges(7,6) * t319) / 0.2e1 + t143 * (-Ifges(7,5) * t305 - Ifges(7,6) * t306 + Ifges(7,3) * t319) / 0.2e1 + (-Ifges(7,1) * t305 - Ifges(7,4) * t306 + Ifges(7,5) * t319) * t438 - t319 * t91 / 0.2e1 + t319 * t51 / 0.2e1 - t319 * t416 + (Ifges(4,1) * t221 + Ifges(4,4) * t222) * t485 + (-t2 * t384 + t20 * t305 - t21 * t306 - t3 * t383) * mrSges(7,3) + t54 * (mrSges(7,1) * t306 - mrSges(7,2) * t305) + t293 * t454 + (Ifges(5,1) * t161 - Ifges(5,4) * t160) * t435 + t383 * t446 + m(5) * (t107 * t134 + t108 * t133 + t162 * t60 + t163 * t59 + t172 * t256 + t236 * t245) + m(4) * t338 + t450 * t293 + t453 * mrSges(5,3) + t236 * (mrSges(5,1) * t160 + mrSges(5,2) * t161) + t245 * t153 + qJ(2) * (-mrSges(4,1) * t222 + mrSges(4,2) * t221) + qJD(2) * t218 + t203 * (Ifges(5,4) * t161 - Ifges(5,2) * t160) / 0.2e1 + t107 * t170 + t108 * t171 - t160 * t144 / 0.2e1 + t161 * t145 / 0.2e1 + t162 * t118 + t163 * t119 + t146 * t97 + t27 * t131 + t86 * t71 + t4 * t87 + t5 * t88 + t32 * t18 + t33 * t19 - t8 * t384 / 0.2e1 - t201 * t360 / 0.2e1 - t202 * t361 / 0.2e1 - pkin(1) * t354 + t177 * t343 + m(3) * (-pkin(1) * t266 + t338) + (-t469 / 0.2e1 + t468 / 0.2e1) * t275 + (t468 - t469) * t434 + qJD(3) ^ 2 * t325 / 0.2e1 + t222 * t327 / 0.2e1 + t221 * t329 / 0.2e1 + (-m(6) * t14 + m(7) * t11 + t472) * (t137 * t282 - t283 * t307) + (-t181 * mrSges(7,1) - t180 * mrSges(7,2) + t476 * (t287 * t224 - t276 * t291 + t364) + (-m(3) + t477) * t364 + (-m(4) * pkin(7) + m(5) * t292 - t452) * t291 + t449 * t287) * g(2) + (-m(3) * t341 - t183 * mrSges(7,1) + t182 * mrSges(7,2) + t476 * (t291 * t224 + t287 * t276 + t341) + t477 * t272 + (-m(5) * (-pkin(1) + t292) - m(4) * t293 + t452) * t287 + t449 * t291) * g(1) + t266 * mrSges(3,2) + m(7) * (t2 * t33 + t20 * t5 + t21 * t4 + t3 * t32) + m(6) * (t146 * t164 + t15 * t86 + t177 * t98 + t27 * t57) - t286 * (Ifges(4,4) * t221 + Ifges(4,2) * t222) / 0.2e1; t354 - t315 * t118 - t216 * t119 + t160 * t170 + t161 * t171 - t472 * t159 - t392 * t100 + t454 + (qJ(2) * t441 - mrSges(3,3)) * t447 - t309 * t319 + (qJD(6) * t320 + t462 + t71) * t316 + m(6) * (-t100 * t56 + t14 * t159 + t15 * t316 + t319 * t57) - m(5) * t453 + m(3) * t266 + m(7) * (t100 * t54 - t11 * t159 + t297 * t316 + t319 * t322) + (-m(5) * t236 - m(6) * t164 - m(7) * t323 - t153 - t218 + t320 - t97) * qJD(1) + t460 * (m(5) - t441 - t476) + t450; t486 * (-t115 * t282 + t283 * t308 + (t282 * t289 + t373) * t401) + (t303 / 0.2e1 - t459) * t447 + t460 * t334 + ((-m(7) * (-t225 + t336) + t455 + t463) * t291 + t457) * g(2) + (-m(7) * t365 + (-m(7) * t225 + t461 - t463) * t287 + t456) * g(1) + t448 * (pkin(9) + t193) + t325 * t340 + t295 + (-m(7) * (-t273 + t312) + m(6) * t224 + t464 + t482) * g(3) + t118 * t429 + t235 * t378 + (t11 * t184 - t20 * t24 - t21 * t25) * m(7) + (t14 * t192 + t15 * t193 - t164 * t169 - t57 * t64) * m(6) + Ifges(4,5) * t221 + Ifges(4,6) * t222 + t192 * t72 + t193 * t71 + t184 * t16 + t165 * mrSges(4,1) - t166 * mrSges(4,2) - t169 * t97 - t139 * t170 - t138 * t171 - t64 * t131 - t25 * t87 - t24 * t88 + Ifges(4,3) * qJDD(3) - t234 * t223 + t201 * t362 / 0.2e1 + t202 * t363 / 0.2e1 - t153 * t262 - m(5) * (t133 * t138 + t134 * t139 + t236 * t262) + (m(6) * t57 + m(7) * t322 - t309) * (t283 * t289 - t374) * t401 + (m(5) * (t285 * t59 + t289 * t60 + (-t133 * t285 + t134 * t289) * qJD(4)) + t170 * t357 - t171 * t358 + t285 * t119) * pkin(3); -t62 * t131 - t133 * t170 + t134 * t171 + t255 * t16 - t22 * t88 - t23 * t87 + t72 * t423 + t71 * t424 - t97 * t427 + t295 + t392 * t61 + (t11 * t255 - t20 * t22 - t21 * t23 - t54 * t61) * m(7) + ((t14 * t283 + t15 * t282) * pkin(4) - t164 * t427 + t56 * t61 - t57 * t62) * m(6) + (m(6) * t426 - m(7) * t312 - t332 + t464) * g(3) + ((-m(7) * (t336 - t425) + t351 + t455) * t291 + t457) * g(2) + (-m(7) * (pkin(4) * t375 + t365) + (-t351 + t461) * t287 + t456) * g(1) + t448 * (pkin(9) + t424); t321 * qJD(6) + t392 * t317 + t309 * t337 + t288 * t18 + t284 * t19 + t343 + (t143 * t322 + t2 * t284 + t3 * t288 - t317 * t54 + t335) * m(7) + (t317 * t56 - t337 * t57 + t335 + t98) * m(6); -t54 * (mrSges(7,1) * t127 + mrSges(7,2) * t126) + (Ifges(7,1) * t126 - t400) * t439 + t52 * t438 + (Ifges(7,5) * t126 - Ifges(7,6) * t127) * t437 - t20 * t87 + t21 * t88 - g(1) * (mrSges(7,1) * t180 - mrSges(7,2) * t181) - g(2) * (mrSges(7,1) * t182 + mrSges(7,2) * t183) + g(3) * t330 * t253 + (t126 * t20 + t127 * t21) * mrSges(7,3) + t350 + (-Ifges(7,2) * t127 + t123 + t53) * t440 + t451;];
tau  = t1;
