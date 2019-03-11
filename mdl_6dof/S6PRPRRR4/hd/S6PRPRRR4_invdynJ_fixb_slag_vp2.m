% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:36:03
% EndTime: 2019-03-08 20:36:38
% DurationCPUTime: 20.87s
% Computational Cost: add. (11283->698), mult. (26769->943), div. (0->0), fcn. (21446->18), ass. (0->318)
t262 = cos(qJ(5));
t242 = pkin(5) * t262 + pkin(4);
t250 = qJ(5) + qJ(6);
t245 = sin(t250);
t246 = cos(t250);
t258 = sin(qJ(5));
t296 = -mrSges(6,1) * t262 + mrSges(6,2) * t258;
t427 = -m(6) * pkin(4) - m(7) * t242 - mrSges(7,1) * t246 + mrSges(7,2) * t245 - mrSges(5,1) + t296;
t265 = -pkin(10) - pkin(9);
t426 = -m(6) * pkin(9) + m(7) * t265 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t254 = cos(pkin(12));
t263 = cos(qJ(4));
t251 = sin(pkin(12));
t259 = sin(qJ(4));
t351 = t251 * t259;
t211 = -t263 * t254 + t351;
t206 = t211 * qJD(4);
t212 = t251 * t263 + t254 * t259;
t207 = t212 * qJD(4);
t260 = sin(qJ(2));
t253 = sin(pkin(6));
t338 = qJD(1) * t253;
t310 = t260 * t338;
t481 = pkin(4) * t207 + pkin(9) * t206 - t310;
t264 = cos(qJ(2));
t345 = t253 * t264;
t270 = t211 * t345;
t380 = pkin(8) + qJ(3);
t222 = t380 * t251;
t223 = t380 * t254;
t435 = -t263 * t222 - t223 * t259;
t437 = qJD(1) * t270 - t211 * qJD(3) + qJD(4) * t435;
t249 = pkin(12) + qJ(4);
t243 = sin(t249);
t244 = cos(t249);
t480 = t426 * t243 + t427 * t244;
t241 = pkin(3) * t254 + pkin(2);
t151 = pkin(4) * t211 - pkin(9) * t212 - t241;
t171 = -t222 * t259 + t223 * t263;
t330 = qJD(5) * t262;
t331 = qJD(5) * t258;
t447 = t151 * t330 - t171 * t331 + t481 * t258 + t437 * t262;
t479 = -t437 * t258 + t481 * t262;
t162 = t262 * t171;
t355 = t206 * t262;
t477 = pkin(10) * t355 + pkin(5) * t207 + (-t162 + (pkin(10) * t212 - t151) * t258) * qJD(5) + t479;
t276 = -t258 * t206 + t212 * t330;
t476 = pkin(10) * t276 - t447;
t313 = qJD(5) * t265;
t336 = qJD(2) * t254;
t204 = -qJD(2) * t351 + t263 * t336;
t357 = t204 * t258;
t205 = t212 * qJD(2);
t149 = pkin(4) * t205 - pkin(9) * t204;
t220 = qJD(2) * qJ(3) + t310;
t255 = cos(pkin(6));
t337 = qJD(1) * t255;
t235 = t254 * t337;
t165 = t235 + (-pkin(8) * qJD(2) - t220) * t251;
t181 = t254 * t220 + t251 * t337;
t166 = pkin(8) * t336 + t181;
t90 = t165 * t263 - t259 * t166;
t60 = t258 * t149 + t262 * t90;
t475 = pkin(10) * t357 + t258 * t313 - t60;
t356 = t204 * t262;
t59 = t262 * t149 - t258 * t90;
t474 = -pkin(5) * t205 + pkin(10) * t356 + t262 * t313 - t59;
t471 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t470 = t331 - t357;
t290 = -t264 * t338 + qJD(3);
t193 = -qJD(2) * t241 + t290;
t442 = Ifges(5,5) * qJD(4);
t469 = -t193 * mrSges(5,2) - t442 / 0.2e1;
t261 = cos(qJ(6));
t197 = qJD(5) - t204;
t178 = qJD(4) * t258 + t205 * t262;
t105 = -pkin(4) * t204 - pkin(9) * t205 + t193;
t91 = t165 * t259 + t166 * t263;
t87 = qJD(4) * pkin(9) + t91;
t51 = t262 * t105 - t258 * t87;
t42 = -pkin(10) * t178 + t51;
t34 = pkin(5) * t197 + t42;
t257 = sin(qJ(6));
t177 = qJD(4) * t262 - t205 * t258;
t52 = t105 * t258 + t262 * t87;
t43 = pkin(10) * t177 + t52;
t365 = t257 * t43;
t15 = t261 * t34 - t365;
t364 = t261 * t43;
t16 = t257 * t34 + t364;
t441 = Ifges(5,6) * qJD(4);
t468 = -t193 * mrSges(5,1) - t51 * mrSges(6,1) - t15 * mrSges(7,1) + t52 * mrSges(6,2) + t16 * mrSges(7,2) + t441 / 0.2e1;
t415 = m(7) * pkin(5);
t194 = Ifges(5,4) * t204;
t466 = Ifges(5,2) * t204;
t465 = Ifges(6,6) * t177;
t464 = Ifges(6,3) * t197;
t462 = t251 ^ 2 + t254 ^ 2;
t295 = mrSges(6,1) * t258 + mrSges(6,2) * t262;
t461 = -t245 * mrSges(7,1) - t246 * mrSges(7,2) - t258 * t415 - mrSges(5,3) - t295;
t299 = t261 * t177 - t178 * t257;
t155 = -qJD(2) * t206 + qJDD(2) * t212;
t84 = qJD(5) * t177 + qJDD(4) * t258 + t155 * t262;
t85 = -qJD(5) * t178 + qJDD(4) * t262 - t155 * t258;
t27 = qJD(6) * t299 + t257 * t85 + t261 * t84;
t414 = t27 / 0.2e1;
t104 = t177 * t257 + t178 * t261;
t28 = -qJD(6) * t104 - t257 * t84 + t261 * t85;
t413 = t28 / 0.2e1;
t406 = t84 / 0.2e1;
t405 = t85 / 0.2e1;
t323 = qJDD(2) * t254;
t324 = qJDD(2) * t251;
t156 = -qJD(2) * t207 - t259 * t324 + t263 * t323;
t148 = qJDD(5) - t156;
t143 = qJDD(6) + t148;
t400 = t143 / 0.2e1;
t399 = t148 / 0.2e1;
t353 = t212 * t262;
t81 = t262 * t151 - t171 * t258;
t65 = pkin(5) * t211 - pkin(10) * t353 + t81;
t354 = t212 * t258;
t82 = t258 * t151 + t162;
t70 = -pkin(10) * t354 + t82;
t31 = t257 * t65 + t261 * t70;
t459 = -qJD(6) * t31 + t476 * t257 + t261 * t477;
t30 = -t257 * t70 + t261 * t65;
t458 = qJD(6) * t30 + t257 * t477 - t476 * t261;
t192 = qJD(6) + t197;
t453 = Ifges(6,5) * t178 + Ifges(7,5) * t104 + Ifges(7,6) * t299 + Ifges(7,3) * t192 + t464 + t465;
t226 = t265 * t258;
t227 = t265 * t262;
t175 = t226 * t261 + t227 * t257;
t452 = qJD(6) * t175 + t257 * t474 + t261 * t475;
t176 = t226 * t257 - t227 * t261;
t451 = -qJD(6) * t176 - t257 * t475 + t261 * t474;
t450 = mrSges(4,3) * t462;
t448 = -qJD(5) * t82 + t479;
t210 = -mrSges(4,1) * t323 + mrSges(4,2) * t324;
t75 = -t156 * mrSges(5,1) + t155 * mrSges(5,2);
t446 = t210 + t75;
t445 = -t415 - mrSges(6,1);
t45 = -mrSges(6,1) * t85 + mrSges(6,2) * t84;
t444 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t155 + t45;
t443 = pkin(5) * t470 - t91;
t284 = t181 * t254 - (-t220 * t251 + t235) * t251;
t440 = t264 * t284;
t282 = t257 * t258 - t261 * t262;
t137 = t282 * t212;
t119 = t282 * t204;
t428 = qJD(5) + qJD(6);
t163 = t428 * t282;
t439 = -t163 + t119;
t216 = t257 * t262 + t258 * t261;
t118 = t216 * t204;
t164 = t428 * t216;
t438 = -t164 + t118;
t376 = mrSges(5,3) * t205;
t436 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t177 + mrSges(6,2) * t178 + t376;
t298 = -mrSges(4,1) * t254 + mrSges(4,2) * t251;
t434 = mrSges(5,1) * t204 - mrSges(5,2) * t205 - t298 * qJD(2);
t305 = qJD(2) * t338;
t229 = t264 * t305;
t327 = qJDD(1) * t253;
t196 = t260 * t327 + t229;
t179 = t196 + t471;
t326 = qJDD(1) * t255;
t233 = t254 * t326;
t144 = -t179 * t251 + t233;
t145 = t254 * t179 + t251 * t326;
t431 = -t144 * t251 + t145 * t254;
t120 = t233 + (-pkin(8) * qJDD(2) - t179) * t251;
t121 = pkin(8) * t323 + t145;
t332 = qJD(4) * t263;
t333 = qJD(4) * t259;
t38 = t259 * t120 + t263 * t121 + t165 * t332 - t166 * t333;
t36 = qJDD(4) * pkin(9) + t38;
t228 = t260 * t305;
t195 = t264 * t327 - t228;
t278 = qJDD(3) - t195;
t167 = -qJDD(2) * t241 + t278;
t66 = -pkin(4) * t156 - pkin(9) * t155 + t167;
t11 = t105 * t330 + t258 * t66 + t262 * t36 - t331 * t87;
t12 = -qJD(5) * t52 - t258 * t36 + t262 * t66;
t430 = t11 * t262 - t12 * t258;
t429 = m(7) + m(6) + m(5);
t86 = -qJD(4) * pkin(4) - t90;
t67 = -pkin(5) * t177 + t86;
t424 = -mrSges(7,1) * t67 + t16 * mrSges(7,3);
t423 = mrSges(7,2) * t67 - t15 * mrSges(7,3);
t422 = -m(6) * t86 - t436;
t6 = pkin(5) * t148 - pkin(10) * t84 + t12;
t7 = pkin(10) * t85 + t11;
t2 = qJD(6) * t15 + t257 * t6 + t261 * t7;
t3 = -qJD(6) * t16 - t257 * t7 + t261 * t6;
t421 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t420 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t307 = m(4) * qJ(3) + mrSges(4,3);
t419 = mrSges(3,2) - t307 + t461;
t274 = m(4) * pkin(2) - t298;
t418 = mrSges(3,1) + t274 - t480;
t417 = Ifges(7,4) * t414 + Ifges(7,2) * t413 + Ifges(7,6) * t400;
t416 = Ifges(7,1) * t414 + Ifges(7,4) * t413 + Ifges(7,5) * t400;
t412 = Ifges(6,1) * t406 + Ifges(6,4) * t405 + Ifges(6,5) * t399;
t369 = Ifges(7,4) * t104;
t49 = Ifges(7,2) * t299 + Ifges(7,6) * t192 + t369;
t411 = -t49 / 0.2e1;
t410 = t49 / 0.2e1;
t99 = Ifges(7,4) * t299;
t50 = Ifges(7,1) * t104 + Ifges(7,5) * t192 + t99;
t409 = -t50 / 0.2e1;
t408 = t50 / 0.2e1;
t404 = -t299 / 0.2e1;
t403 = t299 / 0.2e1;
t402 = -t104 / 0.2e1;
t401 = t104 / 0.2e1;
t398 = -t177 / 0.2e1;
t397 = -t178 / 0.2e1;
t396 = t178 / 0.2e1;
t395 = -t192 / 0.2e1;
t394 = t192 / 0.2e1;
t393 = -t197 / 0.2e1;
t392 = -t204 / 0.2e1;
t390 = t205 / 0.2e1;
t387 = t262 / 0.2e1;
t386 = pkin(5) * t178;
t384 = g(3) * t253;
t363 = cos(pkin(11));
t302 = t363 * t260;
t252 = sin(pkin(11));
t347 = t252 * t264;
t201 = t255 * t302 + t347;
t303 = t253 * t363;
t158 = t201 * t244 - t243 * t303;
t301 = t363 * t264;
t348 = t252 * t260;
t200 = -t255 * t301 + t348;
t379 = (-t158 * t245 + t200 * t246) * mrSges(7,1) + (-t158 * t246 - t200 * t245) * mrSges(7,2);
t203 = -t255 * t348 + t301;
t349 = t252 * t253;
t160 = t203 * t244 + t243 * t349;
t202 = t255 * t347 + t302;
t378 = (-t160 * t245 + t202 * t246) * mrSges(7,1) + (-t160 * t246 - t202 * t245) * mrSges(7,2);
t377 = mrSges(5,3) * t204;
t375 = mrSges(6,3) * t177;
t374 = mrSges(6,3) * t178;
t373 = Ifges(5,4) * t205;
t372 = Ifges(6,4) * t178;
t371 = Ifges(6,4) * t258;
t370 = Ifges(6,4) * t262;
t346 = t253 * t260;
t188 = t243 * t255 + t244 * t346;
t342 = (-t188 * t245 - t246 * t345) * mrSges(7,1) + (-t188 * t246 + t245 * t345) * mrSges(7,2);
t335 = qJD(2) * t260;
t322 = Ifges(7,5) * t27 + Ifges(7,6) * t28 + Ifges(7,3) * t143;
t321 = Ifges(6,5) * t84 + Ifges(6,6) * t85 + Ifges(6,3) * t148;
t56 = -mrSges(7,1) * t299 + mrSges(7,2) * t104;
t319 = t56 + t436;
t318 = m(4) + t429;
t316 = t258 * t345;
t315 = t262 * t345;
t174 = Ifges(6,4) * t177;
t78 = t178 * Ifges(6,1) + t197 * Ifges(6,5) + t174;
t314 = t78 * t387;
t309 = t253 * t335;
t304 = -t331 / 0.2e1;
t293 = Ifges(6,1) * t262 - t371;
t292 = -Ifges(6,2) * t258 + t370;
t291 = Ifges(6,5) * t262 - Ifges(6,6) * t258;
t198 = -t251 * t346 + t254 * t255;
t199 = t251 * t255 + t254 * t346;
t125 = t198 * t259 + t199 * t263;
t107 = -t125 * t258 - t315;
t280 = -t125 * t262 + t316;
t57 = t107 * t261 + t257 * t280;
t58 = t107 * t257 - t261 * t280;
t113 = -mrSges(6,2) * t197 + t375;
t114 = mrSges(6,1) * t197 - t374;
t287 = t113 * t262 - t114 * t258;
t283 = t263 * t198 - t199 * t259;
t39 = t120 * t263 - t259 * t121 - t165 * t333 - t166 * t332;
t281 = t322 + t421;
t277 = t86 * t295;
t275 = t212 * t331 + t355;
t271 = t212 * t345;
t37 = -qJDD(4) * pkin(4) - t39;
t116 = qJD(3) * t212 + qJD(4) * t171;
t266 = qJD(2) ^ 2;
t214 = -qJD(2) * pkin(2) + t290;
t185 = -qJD(4) * mrSges(5,2) + t377;
t184 = -qJDD(2) * pkin(2) + t278;
t136 = t216 * t212;
t130 = t205 * Ifges(5,1) + t194 + t442;
t129 = t373 + t441 + t466;
t128 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t156;
t117 = pkin(5) * t354 - t435;
t89 = qJD(2) * t271 + qJD(4) * t125;
t88 = -qJD(2) * t270 + qJD(4) * t283;
t77 = t177 * Ifges(6,2) + t197 * Ifges(6,6) + t372;
t74 = mrSges(7,1) * t192 - mrSges(7,3) * t104;
t73 = -mrSges(7,2) * t192 + mrSges(7,3) * t299;
t72 = pkin(5) * t276 + t116;
t62 = -mrSges(6,2) * t148 + mrSges(6,3) * t85;
t61 = mrSges(6,1) * t148 - mrSges(6,3) * t84;
t55 = t137 * t428 + t216 * t206;
t54 = -t164 * t212 + t206 * t282;
t47 = qJD(5) * t280 - t258 * t88 + t262 * t309;
t46 = qJD(5) * t107 + t258 * t309 + t262 * t88;
t32 = t84 * Ifges(6,4) + t85 * Ifges(6,2) + t148 * Ifges(6,6);
t23 = -pkin(5) * t85 + t37;
t22 = -mrSges(7,2) * t143 + mrSges(7,3) * t28;
t21 = mrSges(7,1) * t143 - mrSges(7,3) * t27;
t18 = t261 * t42 - t365;
t17 = -t257 * t42 - t364;
t14 = -qJD(6) * t58 - t257 * t46 + t261 * t47;
t13 = qJD(6) * t57 + t257 * t47 + t261 * t46;
t10 = -mrSges(7,1) * t28 + mrSges(7,2) * t27;
t1 = [t107 * t61 - t280 * t62 + t46 * t113 + t47 * t114 + t125 * t128 + t13 * t73 + t14 * t74 + t88 * t185 + t57 * t21 + t58 * t22 + (-t198 * t251 + t199 * t254) * qJDD(2) * mrSges(4,3) + t319 * t89 - (t10 + t444) * t283 + (-m(2) - m(3) - t318) * g(3) + m(4) * (t144 * t198 + t145 * t199) + m(5) * (t125 * t38 + t283 * t39 + t88 * t91 - t89 * t90) + m(7) * (t13 * t16 + t14 * t15 + t2 * t58 - t23 * t283 + t3 * t57 + t67 * t89) + m(6) * (t107 * t12 - t11 * t280 - t283 * t37 + t46 * t52 + t47 * t51 + t86 * t89) + ((mrSges(3,1) * qJDD(2) - t446) * t264 + (-qJDD(2) * mrSges(3,2) - qJD(2) * t434) * t260 + m(4) * (qJD(2) * t440 - t184 * t264 + t214 * t335) + m(3) * (t195 * t264 + t196 * t260) + m(5) * (-t167 * t264 + t193 * t335) + (-t260 * mrSges(3,1) + (-mrSges(3,2) + t450) * t264) * t266) * t253 + (m(3) * t255 ^ 2 + m(2)) * qJDD(1); (mrSges(5,1) * t167 - mrSges(5,3) * t38 - Ifges(5,4) * t155 + Ifges(6,5) * t406 + Ifges(7,5) * t414 - Ifges(5,2) * t156 - Ifges(5,6) * qJDD(4) + Ifges(6,6) * t405 + Ifges(7,6) * t413 + Ifges(6,3) * t399 + Ifges(7,3) * t400 + t420 + t421) * t211 + t197 * (-Ifges(6,5) * t275 - Ifges(6,6) * t276) / 0.2e1 + (t462 * (-t229 + t471) + t431) * mrSges(4,3) + (-mrSges(5,3) * t91 - Ifges(5,4) * t390 - t466 / 0.2e1 + Ifges(7,3) * t394 + Ifges(6,5) * t396 - t129 / 0.2e1 + Ifges(7,5) * t401 + Ifges(7,6) * t403 + t453 / 0.2e1 + t464 / 0.2e1 + t465 / 0.2e1 - t468) * t207 + (t90 * mrSges(5,3) - Ifges(5,1) * t390 - t194 / 0.2e1 - t130 / 0.2e1 - t314 + t469) * t206 + (Ifges(7,4) * t54 + Ifges(7,2) * t55) * t403 + (t167 * mrSges(5,2) - t39 * mrSges(5,3) + Ifges(5,1) * t155 + Ifges(5,4) * t156 + Ifges(5,5) * qJDD(4) + t291 * t399 + t292 * t405 + t293 * t406 + t295 * t37 + t304 * t78) * t212 + (-t116 * t90 - t167 * t241 + t171 * t38 - t193 * t310 + t39 * t435 + t437 * t91) * m(5) + (t11 * t82 + t116 * t86 + t12 * t81 - t37 * t435 + t447 * t52 + t448 * t51) * m(6) - t444 * t435 + (-Ifges(7,4) * t137 - Ifges(7,2) * t136) * t413 + (-Ifges(7,1) * t137 - Ifges(7,4) * t136) * t414 + (-t136 * t2 + t137 * t3 - t15 * t54 + t16 * t55) * mrSges(7,3) + (-Ifges(7,5) * t137 - Ifges(7,6) * t136) * t400 + t23 * (mrSges(7,1) * t136 - mrSges(7,2) * t137) + (t322 + t321) * t211 / 0.2e1 + t177 * (-Ifges(6,4) * t275 - Ifges(6,2) * t276) / 0.2e1 + (Ifges(7,1) * t54 + Ifges(7,4) * t55) * t401 + ((t480 * t264 + (-t380 * t429 + t461) * t260) * t253 - t429 * t241 * t345) * g(3) + (Ifges(4,4) * t251 + Ifges(4,2) * t254) * t323 + (-t264 * t384 + t195 + t228) * mrSges(3,1) + (-t11 * t354 - t12 * t353 + t275 * t51 - t276 * t52) * mrSges(6,3) + (-Ifges(6,1) * t275 - Ifges(6,4) * t276) * t396 + (Ifges(4,1) * t251 + Ifges(4,4) * t254) * t324 + (Ifges(7,5) * t54 + Ifges(7,6) * t55) * t394 + (t260 * t384 - t196 + t229) * mrSges(3,2) + Ifges(3,3) * qJDD(2) - (t260 * t307 + t264 * t274) * t384 - t32 * t354 / 0.2e1 + (-t429 * (-t202 * t241 + t203 * t380) + t419 * t203 + t418 * t202) * g(1) + (-t429 * (-t200 * t241 + t201 * t380) + t419 * t201 + t418 * t200) * g(2) - t241 * t75 - pkin(2) * t210 + t171 * t128 + t117 * t10 + t81 * t61 + t82 * t62 + t72 * t56 + t67 * (-mrSges(7,1) * t55 + mrSges(7,2) * t54) + t30 * t21 + t31 * t22 - t137 * t416 - t136 * t417 + t54 * t408 + t55 * t410 + t353 * t412 + (m(5) * t90 - m(7) * t67 + t422 - t56) * qJD(1) * t271 - t276 * t77 / 0.2e1 + t434 * t310 + t436 * t116 + t437 * t185 + (-pkin(2) * t184 + t284 * qJD(3) + t431 * qJ(3) - (t214 * t260 + t440) * t338) * m(4) + t447 * t113 + t448 * t114 + t458 * t73 + t459 * t74 + (t117 * t23 + t15 * t459 + t16 * t458 + t2 * t31 + t3 * t30 + t67 * t72) * m(7) + t86 * (mrSges(6,1) * t276 - mrSges(6,2) * t275) + t184 * t298; t438 * t74 + t439 * t73 - t266 * t450 - t319 * t205 + t262 * t61 + t258 * t62 - t282 * t21 + t216 * t22 + (-t185 - t287) * t204 + t287 * qJD(5) + (t15 * t438 + t16 * t439 + t2 * t216 - t205 * t67 - t282 * t3) * m(7) + (t11 * t258 + t12 * t262 - t205 * t86 + t197 * (-t258 * t51 + t262 * t52)) * m(6) + (-t204 * t91 + t205 * t90 + t167) * m(5) + (-qJD(2) * t284 + t184) * m(4) + (-g(1) * t202 - g(2) * t200 + g(3) * t345) * t318 + t446; (t130 + t194) * t392 + (t376 + t422) * t91 - (Ifges(7,1) * t401 + Ifges(7,4) * t403 + Ifges(7,5) * t394 + t408 + t423) * t163 - (Ifges(7,4) * t401 + Ifges(7,2) * t403 + Ifges(7,6) * t394 + t410 + t424) * t164 + (t426 * t158 + t427 * (-t201 * t243 - t244 * t303)) * g(2) + (t426 * t160 + t427 * (-t203 * t243 + t244 * t349)) * g(1) + (t426 * t188 + t427 * (-t243 * t346 + t244 * t255)) * g(3) + (-t470 * t52 + (-t330 + t356) * t51 + t430) * mrSges(6,3) + (Ifges(6,5) * t397 + Ifges(7,5) * t402 - Ifges(5,2) * t392 + Ifges(6,6) * t398 + Ifges(7,6) * t404 + Ifges(6,3) * t393 + Ifges(7,3) * t395 + t468) * t205 + (t291 * t393 + t292 * t398 + t293 * t397 - t277 + t469) * t204 + (t277 + t314) * qJD(5) + (-pkin(4) * t37 - t51 * t59 - t52 * t60) * m(6) + (-Ifges(7,1) * t119 - Ifges(7,4) * t118) * t402 + (-Ifges(7,4) * t119 - Ifges(7,2) * t118) * t404 - t67 * (mrSges(7,1) * t118 - mrSges(7,2) * t119) + (-Ifges(7,5) * t119 - Ifges(7,6) * t118) * t395 + (t177 * t292 + t178 * t293 + t197 * t291) * qJD(5) / 0.2e1 + (t377 - t185) * t90 + t32 * t387 + t129 * t390 + (t357 / 0.2e1 + t304) * t77 + Ifges(5,3) * qJDD(4) + (t118 * t16 - t119 * t15 - t2 * t282 - t216 * t3) * mrSges(7,3) + t23 * (mrSges(7,1) * t282 + mrSges(7,2) * t216) + (Ifges(7,1) * t216 - Ifges(7,4) * t282) * t414 + (Ifges(7,5) * t216 - Ifges(7,6) * t282) * t400 + (Ifges(7,4) * t216 - Ifges(7,2) * t282) * t413 - t282 * t417 - t78 * t356 / 0.2e1 - t242 * t10 + t175 * t21 + t176 * t22 + Ifges(5,5) * t155 + Ifges(5,6) * t156 - t60 * t113 - t59 * t114 - pkin(4) * t45 - t38 * mrSges(5,2) + t39 * mrSges(5,1) + t216 * t416 + (Ifges(6,5) * t258 + Ifges(6,6) * t262) * t399 + (Ifges(6,2) * t262 + t371) * t405 + (Ifges(6,1) * t258 + t370) * t406 - t119 * t409 - t118 * t411 + t258 * t412 + (m(6) * ((-t258 * t52 - t262 * t51) * qJD(5) + t430) + t262 * t62 - t258 * t61 - t114 * t330 - t113 * t331) * pkin(9) + t443 * t56 + t451 * t74 + t452 * t73 + (t15 * t451 + t16 * t452 + t175 * t3 + t176 * t2 - t23 * t242 + t443 * t67) * m(7) - (Ifges(5,1) * t204 - t373 + t453) * t205 / 0.2e1 + t37 * t296; -(Ifges(7,4) * t402 + Ifges(7,2) * t404 + Ifges(7,6) * t395 + t411 - t424) * t104 + (Ifges(7,1) * t402 + Ifges(7,4) * t404 + Ifges(7,5) * t395 + t409 - t423) * t299 + (-t342 - (-t188 * t262 + t316) * mrSges(6,2) + t445 * (-t188 * t258 - t315)) * g(3) + (-t378 - (-t160 * t262 - t202 * t258) * mrSges(6,2) + t445 * (-t160 * t258 + t202 * t262)) * g(1) + (-t379 - (-t158 * t262 - t200 * t258) * mrSges(6,2) + t445 * (-t158 * t258 + t200 * t262)) * g(2) + t281 + t420 + (Ifges(6,5) * t177 - Ifges(6,6) * t178) * t393 + t77 * t396 + (Ifges(6,1) * t177 - t372) * t397 - t56 * t386 - m(7) * (t15 * t17 + t16 * t18 + t386 * t67) + (t374 + t114) * t52 + (t375 - t113) * t51 + (-Ifges(6,2) * t178 + t174 + t78) * t398 - t86 * (mrSges(6,1) * t178 + mrSges(6,2) * t177) - t18 * t73 - t17 * t74 + (t2 * t257 + t261 * t3 + (-t15 * t257 + t16 * t261) * qJD(6)) * t415 + t321 + ((-t257 * t74 + t261 * t73) * qJD(6) + t21 * t261 + t22 * t257) * pkin(5); -t67 * (mrSges(7,1) * t104 + mrSges(7,2) * t299) + (Ifges(7,1) * t299 - t369) * t402 + t49 * t401 + (Ifges(7,5) * t299 - Ifges(7,6) * t104) * t395 - t15 * t73 + t16 * t74 - g(1) * t378 - g(2) * t379 - g(3) * t342 + (t104 * t16 + t15 * t299) * mrSges(7,3) + t281 + (-Ifges(7,2) * t104 + t50 + t99) * t404;];
tau  = t1;
