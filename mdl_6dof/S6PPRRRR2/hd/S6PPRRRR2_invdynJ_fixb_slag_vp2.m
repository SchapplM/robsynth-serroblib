% Calculate vector of inverse dynamics joint torques for
% S6PPRRRR2
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:17
% EndTime: 2019-03-08 19:03:48
% DurationCPUTime: 16.44s
% Computational Cost: add. (9977->678), mult. (24773->972), div. (0->0), fcn. (21676->18), ass. (0->323)
t365 = cos(pkin(6));
t227 = qJD(1) * t365 + qJD(2);
t238 = sin(pkin(13));
t240 = sin(pkin(6));
t245 = sin(qJ(3));
t241 = cos(pkin(13));
t364 = cos(pkin(7));
t308 = t241 * t364;
t300 = t245 * t308;
t392 = cos(qJ(3));
t258 = (t238 * t392 + t300) * t240;
t239 = sin(pkin(7));
t354 = t239 * t245;
t123 = qJD(1) * t258 + t227 * t354;
t244 = sin(qJ(4));
t248 = cos(qJ(4));
t298 = pkin(4) * t244 - pkin(10) * t248;
t462 = t298 * qJD(4) - t123;
t218 = -pkin(4) * t248 - pkin(10) * t244 - pkin(3);
t243 = sin(qJ(5));
t279 = t392 * t308;
t273 = t240 * t279;
t264 = qJD(1) * t273;
t356 = t238 * t240;
t322 = qJD(1) * t356;
t323 = t239 * t392;
t252 = -t227 * t323 + t245 * t322 - t264;
t247 = cos(qJ(5));
t337 = qJD(5) * t247;
t339 = qJD(5) * t243;
t341 = qJD(4) * t247;
t349 = t247 * t248;
t445 = t252 * t349 + t218 * t337 + (-t244 * t341 - t248 * t339) * pkin(9) + t462 * t243;
t342 = qJD(4) * t244;
t333 = pkin(9) * t342;
t352 = t243 * t248;
t461 = t243 * t333 + t247 * t462 - t252 * t352;
t460 = -m(5) - m(6);
t229 = pkin(9) * t349;
t278 = pkin(5) * t244 - pkin(11) * t349;
t459 = t278 * qJD(4) + (-t229 + (pkin(11) * t244 - t218) * t243) * qJD(5) + t461;
t340 = qJD(4) * t248;
t265 = t243 * t340 + t244 * t337;
t458 = pkin(11) * t265 - t445;
t343 = qJD(3) * t248;
t320 = t243 * t343;
t249 = -pkin(11) - pkin(10);
t324 = qJD(5) * t249;
t208 = t298 * qJD(3);
t121 = qJD(3) * pkin(9) + t123;
t328 = t240 * t241 * t239;
t164 = -qJD(1) * t328 + t227 * t364;
t88 = -t244 * t121 + t164 * t248;
t66 = t243 * t208 + t247 * t88;
t457 = pkin(11) * t320 + t243 * t324 - t66;
t65 = t247 * t208 - t243 * t88;
t456 = -qJD(3) * t278 + t247 * t324 - t65;
t232 = pkin(5) * t247 + pkin(4);
t237 = qJ(5) + qJ(6);
t234 = sin(t237);
t235 = cos(t237);
t294 = -mrSges(6,1) * t247 + mrSges(6,2) * t243;
t455 = m(6) * pkin(4) + m(7) * t232 + mrSges(7,1) * t235 - mrSges(7,2) * t234 - t294;
t454 = m(6) * pkin(10) - m(7) * t249 + mrSges(6,3) + mrSges(7,3);
t336 = qJD(3) * qJD(4);
t212 = qJDD(3) * t248 - t244 * t336;
t196 = qJDD(5) - t212;
t189 = qJDD(6) + t196;
t399 = t189 / 0.2e1;
t345 = qJD(3) * t244;
t199 = -t243 * t345 + t341;
t213 = qJDD(3) * t244 + t248 * t336;
t131 = qJD(5) * t199 + qJDD(4) * t243 + t213 * t247;
t200 = qJD(4) * t243 + t247 * t345;
t132 = -qJD(5) * t200 + qJDD(4) * t247 - t213 * t243;
t242 = sin(qJ(6));
t246 = cos(qJ(6));
t138 = t199 * t242 + t200 * t246;
t51 = -qJD(6) * t138 - t131 * t242 + t132 * t246;
t411 = t51 / 0.2e1;
t306 = t246 * t199 - t200 * t242;
t50 = qJD(6) * t306 + t131 * t246 + t132 * t242;
t412 = t50 / 0.2e1;
t413 = Ifges(7,1) * t412 + Ifges(7,4) * t411 + Ifges(7,5) * t399;
t414 = Ifges(7,4) * t412 + Ifges(7,2) * t411 + Ifges(7,6) * t399;
t295 = -mrSges(5,1) * t248 + mrSges(5,2) * t244;
t432 = m(7) - t460;
t453 = pkin(3) * t432 + t244 * t454 + t248 * t455 + mrSges(4,1) - t295;
t405 = t131 / 0.2e1;
t404 = t132 / 0.2e1;
t398 = t196 / 0.2e1;
t452 = t212 / 0.2e1;
t451 = t213 / 0.2e1;
t198 = t247 * t218;
t351 = t244 * t247;
t139 = -pkin(11) * t351 + t198 + (-pkin(9) * t243 - pkin(5)) * t248;
t166 = t243 * t218 + t229;
t353 = t243 * t244;
t151 = -pkin(11) * t353 + t166;
t93 = t139 * t242 + t151 * t246;
t450 = -qJD(6) * t93 + t242 * t458 + t246 * t459;
t92 = t139 * t246 - t151 * t242;
t449 = qJD(6) * t92 + t242 * t459 - t246 * t458;
t222 = t249 * t243;
t223 = t249 * t247;
t156 = t222 * t246 + t223 * t242;
t448 = qJD(6) * t156 + t242 * t456 + t246 * t457;
t415 = m(7) * pkin(5);
t447 = -mrSges(6,1) - t415;
t157 = t222 * t242 - t223 * t246;
t446 = -qJD(6) * t157 - t242 * t457 + t246 * t456;
t444 = -qJD(5) * t166 + t461;
t82 = -mrSges(6,1) * t132 + mrSges(6,2) * t131;
t443 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t213 + t82;
t358 = t164 * t244;
t388 = pkin(5) * t243;
t442 = pkin(5) * t339 - t358 - (qJD(3) * t388 + t121) * t248;
t228 = qJD(5) - t343;
t224 = qJD(6) + t228;
t441 = t200 * Ifges(6,5) + t138 * Ifges(7,5) + t199 * Ifges(6,6) + Ifges(7,6) * t306 + t228 * Ifges(6,3) + t224 * Ifges(7,3);
t282 = t242 * t243 - t246 * t247;
t171 = t282 * t244;
t203 = t295 * qJD(3);
t440 = mrSges(4,1) * qJD(3) - t203;
t330 = mrSges(5,3) * t345;
t439 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t199 + mrSges(6,2) * t200 + t330;
t363 = cos(pkin(12));
t297 = t365 * t363;
t362 = sin(pkin(12));
t256 = t238 * t362 - t241 * t297;
t310 = t240 * t363;
t438 = t239 * t310 + t256 * t364;
t296 = t365 * t362;
t257 = t238 * t363 + t241 * t296;
t309 = t240 * t362;
t437 = -t239 * t309 + t257 * t364;
t152 = -mrSges(5,1) * t212 + mrSges(5,2) * t213;
t436 = mrSges(4,1) * qJDD(3) - t152;
t193 = Ifges(6,4) * t199;
t126 = t200 * Ifges(6,1) + t228 * Ifges(6,5) + t193;
t233 = Ifges(5,4) * t343;
t435 = Ifges(5,1) * t345 + Ifges(5,5) * qJD(4) + t247 * t126 + t233;
t226 = qJDD(1) * t365 + qJDD(2);
t162 = -qJDD(1) * t328 + t226 * t364;
t281 = t240 * t300;
t302 = qJD(3) * t322;
t303 = qJD(3) * t323;
t317 = qJDD(1) * t356;
t76 = qJD(3) * t264 + qJDD(1) * t281 + t226 * t354 + t227 * t303 - t245 * t302 + t392 * t317;
t74 = qJDD(3) * pkin(9) + t76;
t31 = -t121 * t342 + t244 * t162 + t164 * t340 + t248 * t74;
t32 = -t121 * t340 + t162 * t248 - t164 * t342 - t244 * t74;
t434 = -t244 * t32 + t248 * t31;
t113 = qJD(3) * t218 + t252;
t29 = qJDD(4) * pkin(10) + t31;
t344 = qJD(3) * t245;
t321 = t239 * t344;
t77 = -qJD(3) * qJD(1) * t281 + qJDD(1) * t273 + t226 * t323 - t227 * t321 - t245 * t317 - t392 * t302;
t75 = -qJDD(3) * pkin(3) - t77;
t52 = -t212 * pkin(4) - t213 * pkin(10) + t75;
t89 = t121 * t248 + t358;
t87 = qJD(4) * pkin(10) + t89;
t8 = t113 * t337 + t243 * t52 + t247 * t29 - t339 * t87;
t43 = t113 * t243 + t247 * t87;
t9 = -qJD(5) * t43 - t243 * t29 + t247 * t52;
t433 = -t243 * t9 + t247 * t8;
t431 = qJD(5) + qJD(6);
t430 = mrSges(5,1) + t455;
t429 = mrSges(5,2) - t454;
t41 = pkin(11) * t199 + t43;
t367 = t246 * t41;
t42 = t247 * t113 - t243 * t87;
t40 = -pkin(11) * t200 + t42;
t39 = pkin(5) * t228 + t40;
t11 = t242 * t39 + t367;
t86 = -qJD(4) * pkin(4) - t88;
t69 = -pkin(5) * t199 + t86;
t428 = -mrSges(7,1) * t69 + mrSges(7,3) * t11;
t370 = t242 * t41;
t10 = t246 * t39 - t370;
t427 = mrSges(7,2) * t69 - mrSges(7,3) * t10;
t425 = m(4) + m(3) + t432;
t423 = -m(6) * t86 - t439;
t422 = -t9 * mrSges(6,1) + t8 * mrSges(6,2);
t6 = pkin(5) * t196 - pkin(11) * t131 + t9;
t7 = pkin(11) * t132 + t8;
t2 = qJD(6) * t10 + t242 * t6 + t246 * t7;
t3 = -qJD(6) * t11 - t242 * t7 + t246 * t6;
t421 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t293 = t243 * mrSges(6,1) + t247 * mrSges(6,2);
t326 = pkin(9) + t388;
t420 = -m(7) * t326 - t234 * mrSges(7,1) - t235 * mrSges(7,2) + pkin(9) * t460 + mrSges(4,2) - t293;
t377 = Ifges(5,4) * t244;
t290 = t248 * Ifges(5,2) + t377;
t419 = t11 * mrSges(7,2) + Ifges(5,6) * qJD(4) / 0.2e1 + qJD(3) * t290 / 0.2e1 - t10 * mrSges(7,1);
t90 = -mrSges(7,1) * t306 + mrSges(7,2) * t138;
t418 = -m(7) * t69 + t423 - t90;
t18 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t30 = -qJDD(4) * pkin(4) - t32;
t19 = -pkin(5) * t132 + t30;
t417 = -m(5) * t32 + m(6) * t30 + m(7) * t19 + t18 + t443;
t416 = -m(5) * t88 - t418;
t250 = qJD(3) ^ 2;
t410 = Ifges(6,1) * t405 + Ifges(6,4) * t404 + Ifges(6,5) * t398;
t372 = Ifges(7,4) * t138;
t79 = Ifges(7,2) * t306 + Ifges(7,6) * t224 + t372;
t409 = -t79 / 0.2e1;
t408 = t79 / 0.2e1;
t133 = Ifges(7,4) * t306;
t80 = Ifges(7,1) * t138 + Ifges(7,5) * t224 + t133;
t407 = -t80 / 0.2e1;
t406 = t80 / 0.2e1;
t403 = -t306 / 0.2e1;
t402 = t306 / 0.2e1;
t401 = -t138 / 0.2e1;
t400 = t138 / 0.2e1;
t396 = t200 / 0.2e1;
t395 = -t224 / 0.2e1;
t394 = t224 / 0.2e1;
t389 = pkin(5) * t200;
t174 = t238 * t297 + t241 * t362;
t104 = t174 * t245 + t392 * t438;
t105 = t174 * t392 - t245 * t438;
t142 = t239 * t256 - t310 * t364;
t60 = t105 * t248 + t142 * t244;
t382 = (t104 * t235 - t234 * t60) * mrSges(7,1) + (-t104 * t234 - t235 * t60) * mrSges(7,2);
t175 = -t238 * t296 + t241 * t363;
t106 = t175 * t245 + t392 * t437;
t107 = t175 * t392 - t245 * t437;
t143 = t239 * t257 + t309 * t364;
t62 = t107 * t248 + t143 * t244;
t381 = (t106 * t235 - t234 * t62) * mrSges(7,1) + (-t106 * t234 - t235 * t62) * mrSges(7,2);
t311 = t239 * t365;
t141 = t245 * t311 + t258;
t263 = t364 * t365 - t328;
t110 = t141 * t248 + t244 * t263;
t280 = t392 * t311;
t355 = t238 * t245;
t140 = t240 * t355 - t273 - t280;
t380 = (-t110 * t234 + t140 * t235) * mrSges(7,1) + (-t110 * t235 - t140 * t234) * mrSges(7,2);
t379 = mrSges(6,3) * t199;
t378 = mrSges(6,3) * t200;
t376 = Ifges(5,4) * t248;
t375 = Ifges(6,4) * t200;
t374 = Ifges(6,4) * t243;
t373 = Ifges(6,4) * t247;
t371 = pkin(5) * qJD(6);
t369 = t244 * t30;
t360 = mrSges(4,2) * qJD(3);
t346 = mrSges(4,2) * qJDD(3);
t338 = qJD(5) * t244;
t335 = Ifges(7,5) * t50 + Ifges(7,6) * t51 + Ifges(7,3) * t189;
t332 = pkin(9) * t340;
t329 = mrSges(5,3) * t343;
t327 = Ifges(6,5) * t131 + Ifges(6,6) * t132 + Ifges(6,3) * t196;
t125 = t199 * Ifges(6,2) + t228 * Ifges(6,6) + t375;
t319 = -t243 * t125 / 0.2e1;
t307 = t336 / 0.2e1;
t292 = Ifges(6,1) * t247 - t374;
t291 = Ifges(6,1) * t243 + t373;
t289 = -Ifges(6,2) * t243 + t373;
t288 = Ifges(6,2) * t247 + t374;
t287 = Ifges(5,5) * t248 - Ifges(5,6) * t244;
t286 = Ifges(6,5) * t247 - Ifges(6,6) * t243;
t285 = Ifges(6,5) * t243 + Ifges(6,6) * t247;
t63 = -t110 * t243 + t140 * t247;
t64 = t110 * t247 + t140 * t243;
t26 = -t242 * t64 + t246 * t63;
t27 = t242 * t63 + t246 * t64;
t179 = t244 * t364 + t248 * t354;
t149 = -t243 * t179 - t247 * t323;
t267 = -t247 * t179 + t243 * t323;
t94 = t149 * t246 + t242 * t267;
t95 = t149 * t242 - t246 * t267;
t202 = t242 * t247 + t243 * t246;
t277 = t335 - t421;
t276 = t86 * t293;
t120 = -qJD(3) * pkin(3) + t252;
t275 = t120 * (mrSges(5,1) * t244 + mrSges(5,2) * t248);
t274 = t244 * (Ifges(5,1) * t248 - t377);
t272 = t202 * t248;
t271 = t282 * t248;
t266 = -t243 * t338 + t247 * t340;
t178 = t244 * t354 - t248 * t364;
t262 = Ifges(6,5) * t244 + t248 * t292;
t261 = Ifges(6,6) * t244 + t248 * t289;
t260 = Ifges(6,3) * t244 + t248 * t286;
t109 = t141 * t244 - t248 * t263;
t145 = t431 * t202;
t220 = -qJD(4) * mrSges(5,2) + t329;
t214 = t326 * t244;
t172 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t212;
t170 = t202 * t244;
t168 = qJD(3) * t271;
t167 = qJD(3) * t272;
t165 = -pkin(9) * t352 + t198;
t163 = pkin(5) * t265 + t332;
t159 = mrSges(6,1) * t228 - t378;
t158 = -mrSges(6,2) * t228 + t379;
t147 = -qJD(4) * t178 + t248 * t303;
t130 = t141 * qJD(3);
t129 = (t280 + (t279 - t355) * t240) * qJD(3);
t118 = mrSges(7,1) * t224 - mrSges(7,3) * t138;
t117 = -mrSges(7,2) * t224 + mrSges(7,3) * t306;
t112 = -mrSges(6,2) * t196 + mrSges(6,3) * t132;
t111 = mrSges(6,1) * t196 - mrSges(6,3) * t131;
t98 = -qJD(4) * t272 + t171 * t431;
t97 = -qJD(4) * t271 - t145 * t244;
t85 = qJD(5) * t267 - t243 * t147 + t247 * t321;
t84 = qJD(5) * t149 + t247 * t147 + t243 * t321;
t67 = t131 * Ifges(6,4) + t132 * Ifges(6,2) + t196 * Ifges(6,6);
t58 = -qJD(4) * t109 + t129 * t248;
t45 = -mrSges(7,2) * t189 + mrSges(7,3) * t51;
t44 = mrSges(7,1) * t189 - mrSges(7,3) * t50;
t21 = -qJD(6) * t95 - t242 * t84 + t246 * t85;
t20 = qJD(6) * t94 + t242 * t85 + t246 * t84;
t17 = qJD(5) * t63 + t130 * t243 + t247 * t58;
t16 = -qJD(5) * t64 + t130 * t247 - t243 * t58;
t13 = t246 * t40 - t370;
t12 = -t242 * t40 - t367;
t5 = -qJD(6) * t27 + t16 * t246 - t17 * t242;
t4 = qJD(6) * t26 + t16 * t242 + t17 * t246;
t1 = [m(7) * (t10 * t5 + t11 * t4 + t2 * t27 + t26 * t3) + m(6) * (t16 * t42 + t17 * t43 + t63 * t9 + t64 * t8) + m(5) * (t110 * t31 + t58 * t89) - t141 * t346 - t129 * t360 + m(3) * (t226 * t365 + (t238 ^ 2 + t241 ^ 2) * t240 ^ 2 * qJDD(1)) + t58 * t220 + t110 * t172 + t17 * t158 + t16 * t159 + t4 * t117 + t5 * t118 + t63 * t111 + t64 * t112 + t26 * t44 + t27 * t45 + m(4) * (t123 * t129 + t76 * t141 + t162 * t263) + m(2) * qJDD(1) + (-m(4) * t77 + m(5) * t75 - t436) * t140 + (m(4) * t252 + m(5) * t120 - t440) * t130 + t416 * (qJD(4) * t110 + t129 * t244) + t417 * t109 + (-m(2) - t425) * g(3); m(5) * (t89 * t147 + t31 * t179 + (t120 * t344 - t392 * t75) * t239) + m(4) * (t162 * t364 + (t392 * t77 + t245 * t76 + (t123 * t392 + t245 * t252) * qJD(3)) * t239) + t203 * t321 + m(7) * (t10 * t21 + t11 * t20 + t2 * t95 + t3 * t94) + m(6) * (t149 * t9 - t267 * t8 + t42 * t85 + t43 * t84) + m(3) * t226 + t147 * t220 + t179 * t172 + t84 * t158 + t85 * t159 + t149 * t111 - t267 * t112 + t20 * t117 + t21 * t118 + t94 * t44 + t95 * t45 + (-mrSges(4,1) * t250 - t346) * t354 + (-mrSges(4,2) * t250 + t436) * t323 + t417 * t178 + t416 * (qJD(4) * t179 + t244 * t303) + (-g(1) * t309 + g(2) * t310 - g(3) * t365) * t425; (Ifges(7,5) * t97 + Ifges(7,6) * t98) * t394 + (-t418 * t244 + t248 * t220 + m(5) * (-t244 * t88 + t248 * t89) - t360) * t252 + (Ifges(5,4) * t451 + Ifges(5,2) * t452 + pkin(9) * t172 - Ifges(6,6) * t404 - Ifges(6,5) * t405 - Ifges(7,6) * t411 - Ifges(7,5) * t412 - Ifges(6,3) * t398 - Ifges(7,3) * t399 + (-Ifges(5,2) * t244 + t376) * t307 + t421 + t422) * t248 + (-t120 * t123 - pkin(3) * t75 + ((-t244 * t89 - t248 * t88) * qJD(4) + t434) * pkin(9)) * m(5) + (-g(1) * t107 - g(2) * t105 - g(3) * t141 - t340 * t88 + t434) * mrSges(5,3) + (t42 * mrSges(6,1) - t43 * mrSges(6,2) + Ifges(7,5) * t400 + Ifges(7,6) * t402 + Ifges(7,3) * t394 - t419 - t89 * mrSges(5,3) + t441 / 0.2e1) * t342 + (Ifges(7,4) * t97 + Ifges(7,2) * t98) * t402 + (t106 * t453 + t107 * t420) * g(1) + (t140 * t453 + t141 * t420) * g(3) + (t104 * t453 + t105 * t420) * g(2) + (-Ifges(7,4) * t171 - Ifges(7,2) * t170) * t411 + (-t10 * t97 + t11 * t98 - t170 * t2 + t171 * t3) * mrSges(7,3) + t19 * (mrSges(7,1) * t170 - mrSges(7,2) * t171) + (-Ifges(7,1) * t171 - Ifges(7,4) * t170) * t412 + (-Ifges(7,5) * t171 - Ifges(7,6) * t170) * t399 - (t247 * t125 + t243 * t126) * t338 / 0.2e1 - (t335 + t327) * t248 / 0.2e1 + t449 * t117 + t450 * t118 - t220 * t333 + (t10 * t450 + t11 * t449 + t163 * t69 + t19 * t214 + t2 * t93 + t3 * t92) * m(7) + t439 * t332 + t440 * t123 + t445 * t158 + (t165 * t9 + t166 * t8 + (t340 * t86 + t369) * pkin(9) + t445 * t43 + t444 * t42) * m(6) + t199 * (qJD(4) * t261 - t288 * t338) / 0.2e1 + t228 * (qJD(4) * t260 - t285 * t338) / 0.2e1 + t443 * pkin(9) * t244 + (-t265 * t43 - t266 * t42 - t351 * t9 - t353 * t8) * mrSges(6,3) + (Ifges(5,1) * t213 + Ifges(5,4) * t452 + t286 * t398 + t289 * t404 + t292 * t405) * t244 + t435 * t340 / 0.2e1 - t67 * t353 / 0.2e1 + t376 * t451 + t290 * t452 + t75 * t295 + qJD(4) ^ 2 * t287 / 0.2e1 + Ifges(4,3) * qJDD(3) + qJDD(4) * (Ifges(5,5) * t244 + Ifges(5,6) * t248) + t214 * t18 + t165 * t111 + t166 * t112 - pkin(3) * t152 + t163 * t90 + t86 * (mrSges(6,1) * t265 + mrSges(6,2) * t266) + t69 * (-mrSges(7,1) * t98 + mrSges(7,2) * t97) + t92 * t44 + t93 * t45 - t76 * mrSges(4,2) + t77 * mrSges(4,1) + t97 * t406 + t98 * t408 + t351 * t410 - t171 * t413 - t170 * t414 + (qJD(4) * t262 - t291 * t338) * t396 + t293 * t369 + t319 * t340 + t274 * t307 + qJD(4) * t275 + (Ifges(7,1) * t97 + Ifges(7,4) * t98) * t400 + t444 * t159; (-Ifges(7,5) * t168 - Ifges(7,6) * t167) * t395 + (-Ifges(7,4) * t168 - Ifges(7,2) * t167) * t403 + (-Ifges(7,1) * t168 - Ifges(7,4) * t167) * t401 - t69 * (mrSges(7,1) * t167 - mrSges(7,2) * t168) + (-pkin(4) * t30 - t42 * t65 - t43 * t66) * m(6) + (t19 * mrSges(7,1) - 0.2e1 * t414 - t2 * mrSges(7,3) - (Ifges(7,1) * t400 + Ifges(7,4) * t402 + Ifges(7,5) * t394 + t406 + t427) * t431) * t282 + (mrSges(7,2) * t19 - mrSges(7,3) * t3 + 0.2e1 * t413) * t202 + (t199 * t289 + t200 * t292 + t228 * t286) * qJD(5) / 0.2e1 - (t199 * t261 + t200 * t262 + t228 * t260) * qJD(3) / 0.2e1 + (t10 * t446 + t11 * t448 + t156 * t3 + t157 * t2 - t19 * t232 + t442 * t69) * m(7) + t448 * t117 + (t330 + t423) * t89 + t446 * t118 + (-t10 * t168 + t11 * t167) * mrSges(7,3) + t125 * t320 / 0.2e1 + (t329 - t220) * t88 + (t319 + t276) * qJD(5) + (-t275 - t42 * (mrSges(6,1) * t244 - mrSges(6,3) * t349) - t43 * (-mrSges(6,2) * t244 - mrSges(6,3) * t352)) * qJD(3) + (Ifges(7,5) * t401 + Ifges(7,6) * t403 + Ifges(7,3) * t395 + t419) * t345 - t287 * t336 / 0.2e1 + t126 * t337 / 0.2e1 + (m(6) * ((-t43 * t243 - t42 * t247) * qJD(5) + t433) + t247 * t112 - t243 * t111 - t159 * t337 - t158 * t339) * pkin(10) + (-t337 * t42 - t339 * t43 + t433) * mrSges(6,3) - (-Ifges(5,2) * t345 + t233 + t435) * t343 / 0.2e1 - t250 * t274 / 0.2e1 + (t429 * t62 - t430 * (-t107 * t244 + t143 * t248)) * g(1) + (t109 * t430 + t110 * t429) * g(3) + (t429 * t60 - t430 * (-t105 * t244 + t142 * t248)) * g(2) - (Ifges(7,4) * t400 + Ifges(7,2) * t402 + Ifges(7,6) * t394 + t408 + t428) * t145 + t30 * t294 + Ifges(5,3) * qJDD(4) + t247 * t67 / 0.2e1 - t232 * t18 + Ifges(5,5) * t213 + Ifges(5,6) * t212 + t156 * t44 + t157 * t45 - t66 * t158 - t65 * t159 - pkin(4) * t82 + t288 * t404 + t291 * t405 - t168 * t407 - t167 * t409 + t243 * t410 + t285 * t398 - t31 * mrSges(5,2) + t32 * mrSges(5,1) - t276 * t343 - t441 * t345 / 0.2e1 + t442 * t90; -t422 + (-t381 - (-t106 * t243 - t247 * t62) * mrSges(6,2) + t447 * (t106 * t247 - t243 * t62)) * g(1) + (-t382 - (-t104 * t243 - t247 * t60) * mrSges(6,2) + t447 * (t104 * t247 - t243 * t60)) * g(2) - (-Ifges(6,2) * t200 + t126 + t193) * t199 / 0.2e1 - (Ifges(7,4) * t401 + Ifges(7,2) * t403 + Ifges(7,6) * t395 + t409 - t428) * t138 + (Ifges(7,1) * t401 + Ifges(7,4) * t403 + Ifges(7,5) * t395 + t407 - t427) * t306 + (t159 + t378) * t43 - m(7) * (t10 * t12 + t11 * t13 + t389 * t69) - t90 * t389 + (-t158 + t379) * t42 + t277 + (mrSges(6,2) * t64 + t447 * t63 - t380) * g(3) + t327 - t200 * (Ifges(6,1) * t199 - t375) / 0.2e1 + (-t242 * t371 - t12) * t118 + (t246 * t371 - t13) * t117 - t228 * (Ifges(6,5) * t199 - Ifges(6,6) * t200) / 0.2e1 - t86 * (mrSges(6,1) * t200 + mrSges(6,2) * t199) + (t242 * t45 + t246 * t44) * pkin(5) + (t2 * t242 + t246 * t3 + (-t10 * t242 + t11 * t246) * qJD(6)) * t415 + t125 * t396; -t69 * (mrSges(7,1) * t138 + mrSges(7,2) * t306) + (Ifges(7,1) * t306 - t372) * t401 + t79 * t400 + (Ifges(7,5) * t306 - Ifges(7,6) * t138) * t395 - t10 * t117 + t11 * t118 - g(1) * t381 - g(2) * t382 - g(3) * t380 + (t10 * t306 + t11 * t138) * mrSges(7,3) + t277 + (-Ifges(7,2) * t138 + t133 + t80) * t403;];
tau  = t1;
