% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:46
% EndTime: 2019-03-09 03:58:17
% DurationCPUTime: 18.80s
% Computational Cost: add. (10989->732), mult. (21949->962), div. (0->0), fcn. (15159->14), ass. (0->340)
t246 = sin(qJ(5));
t243 = sin(pkin(10));
t385 = pkin(3) * t243;
t222 = pkin(8) + t385;
t375 = pkin(9) + t222;
t299 = qJD(5) * t375;
t247 = sin(qJ(3));
t350 = cos(pkin(10));
t297 = t350 * t247;
t387 = cos(qJ(3));
t304 = qJD(1) * t387;
t181 = -qJD(1) * t297 - t243 * t304;
t347 = t181 * t246;
t253 = -pkin(1) - pkin(7);
t211 = qJD(1) * t253 + qJD(2);
t326 = qJD(1) * t247;
t168 = -qJ(4) * t326 + t211 * t247;
t157 = t243 * t168;
t202 = t387 * t211;
t169 = -qJ(4) * t304 + t202;
t112 = t169 * t350 - t157;
t289 = t350 * t387;
t180 = -qJD(1) * t289 + t243 * t326;
t290 = pkin(3) * t304;
t114 = -t180 * pkin(4) - t181 * pkin(8) + t290;
t250 = cos(qJ(5));
t59 = t250 * t112 + t246 * t114;
t466 = -pkin(9) * t347 + t246 * t299 + t59;
t346 = t181 * t250;
t58 = -t112 * t246 + t250 * t114;
t465 = pkin(5) * t180 + pkin(9) * t346 - t250 * t299 - t58;
t324 = qJD(5) * t246;
t464 = t324 - t347;
t203 = pkin(3) * t326 + qJD(1) * qJ(2) + qJD(4);
t388 = t246 / 0.2e1;
t439 = Ifges(5,5) * qJD(3);
t146 = qJD(3) * t250 + t180 * t246;
t171 = qJD(5) - t181;
t147 = qJD(3) * t246 - t180 * t250;
t365 = t147 * Ifges(6,4);
t73 = t146 * Ifges(6,2) + t171 * Ifges(6,6) + t365;
t143 = Ifges(6,4) * t146;
t74 = t147 * Ifges(6,1) + t171 * Ifges(6,5) + t143;
t463 = -t203 * mrSges(5,2) - t250 * t74 / 0.2e1 + t73 * t388 - t439 / 0.2e1;
t249 = cos(qJ(6));
t101 = -pkin(4) * t181 + pkin(8) * t180 + t203;
t164 = qJD(3) * pkin(3) + t169;
t298 = t350 * t168;
t98 = t243 * t164 + t298;
t89 = qJD(3) * pkin(8) + t98;
t54 = t250 * t101 - t246 * t89;
t40 = -pkin(9) * t147 + t54;
t34 = pkin(5) * t171 + t40;
t245 = sin(qJ(6));
t55 = t101 * t246 + t250 * t89;
t41 = pkin(9) * t146 + t55;
t355 = t245 * t41;
t13 = t249 * t34 - t355;
t353 = t249 * t41;
t14 = t245 * t34 + t353;
t438 = Ifges(5,6) * qJD(3);
t462 = t203 * mrSges(5,1) + t54 * mrSges(6,1) + t13 * mrSges(7,1) - t55 * mrSges(6,2) - t14 * mrSges(7,2) - t438 / 0.2e1;
t461 = t387 / 0.2e1;
t251 = cos(qJ(1));
t381 = g(2) * t251;
t170 = Ifges(5,4) * t181;
t459 = t146 * Ifges(6,6);
t458 = t171 * Ifges(6,3);
t457 = t181 * Ifges(5,2);
t240 = qJ(3) + pkin(10);
t231 = sin(t240);
t232 = cos(t240);
t270 = t247 * mrSges(4,1) + mrSges(4,2) * t387;
t456 = t231 * mrSges(5,1) + t270 + (mrSges(5,2) - mrSges(6,3) - mrSges(7,3)) * t232;
t415 = m(7) * pkin(5);
t291 = t249 * t146 - t147 * t245;
t318 = qJD(1) * qJD(3);
t200 = qJDD(1) * t387 - t247 * t318;
t303 = qJD(3) * t387;
t201 = -qJD(1) * t303 - t247 * qJDD(1);
t137 = t200 * t350 + t243 * t201;
t77 = qJD(5) * t146 + qJDD(3) * t246 + t137 * t250;
t78 = -qJD(5) * t147 + qJDD(3) * t250 - t137 * t246;
t27 = qJD(6) * t291 + t245 * t78 + t249 * t77;
t414 = t27 / 0.2e1;
t84 = t146 * t245 + t147 * t249;
t28 = -qJD(6) * t84 - t245 * t77 + t249 * t78;
t413 = t28 / 0.2e1;
t402 = -m(3) - m(4);
t136 = -t243 * t200 + t201 * t350;
t135 = qJDD(5) - t136;
t132 = qJDD(6) + t135;
t401 = t132 / 0.2e1;
t453 = t231 * t381;
t190 = t375 * t246;
t191 = t375 * t250;
t129 = -t190 * t245 + t191 * t249;
t448 = -qJD(6) * t129 + t245 * t466 + t465 * t249;
t128 = -t190 * t249 - t191 * t245;
t447 = qJD(6) * t128 + t465 * t245 - t249 * t466;
t165 = qJD(6) + t171;
t446 = t147 * Ifges(6,5) + t84 * Ifges(7,5) + Ifges(7,6) * t291 + t165 * Ifges(7,3) + t458 + t459;
t445 = -mrSges(6,1) - t415;
t193 = t243 * t387 + t297;
t277 = t245 * t246 - t249 * t250;
t120 = t277 * t193;
t325 = qJD(3) * t247;
t182 = -qJD(3) * t289 + t243 * t325;
t196 = t245 * t250 + t246 * t249;
t429 = qJD(5) + qJD(6);
t444 = t277 * qJD(1) + t120 * t429 + t196 * t182;
t139 = t429 * t196;
t443 = -t196 * qJD(1) - t139 * t193 + t182 * t277;
t125 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t137;
t38 = -mrSges(6,1) * t78 + mrSges(6,2) * t77;
t442 = t38 - t125;
t111 = t169 * t243 + t298;
t441 = pkin(5) * t464 - t111;
t364 = t180 * mrSges(5,3);
t440 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t146 + mrSges(6,2) * t147 - t364;
t192 = t243 * t247 - t289;
t119 = t196 * t192;
t248 = sin(qJ(1));
t437 = t231 * t248;
t238 = t247 * pkin(3);
t224 = qJ(2) + t238;
t130 = pkin(4) * t193 + pkin(8) * t192 + t224;
t330 = qJ(4) - t253;
t204 = t330 * t247;
t307 = t387 * t253;
t205 = -qJ(4) * t387 + t307;
t141 = -t204 * t350 + t243 * t205;
t133 = t250 * t141;
t71 = t246 * t130 + t133;
t110 = t277 * t181;
t138 = t429 * t277;
t436 = -t138 + t110;
t109 = t196 * t181;
t435 = -t139 + t109;
t377 = t250 * pkin(5);
t227 = pkin(4) + t377;
t242 = qJ(5) + qJ(6);
t234 = sin(t242);
t235 = cos(t242);
t287 = -mrSges(6,1) * t250 + mrSges(6,2) * t246;
t434 = m(6) * pkin(4) + m(7) * t227 + t235 * mrSges(7,1) - t234 * mrSges(7,2) - t287;
t319 = qJD(1) * qJD(2);
t212 = qJDD(1) * qJ(2) + t319;
t210 = qJDD(1) * t253 + qJDD(2);
t144 = t387 * t210 - t211 * t325;
t145 = t247 * t210 + t211 * t303;
t433 = -t144 * t387 - t145 * t247;
t52 = mrSges(6,1) * t135 - mrSges(6,3) * t77;
t53 = -mrSges(6,2) * t135 + mrSges(6,3) * t78;
t432 = -t246 * t52 + t250 * t53;
t323 = qJD(5) * t250;
t321 = t247 * qJD(4);
t104 = qJ(4) * t201 - qJD(1) * t321 + t145;
t302 = t387 * qJD(4);
t94 = qJDD(3) * pkin(3) - t200 * qJ(4) - qJD(1) * t302 + t144;
t57 = t350 * t104 + t243 * t94;
t51 = qJDD(3) * pkin(8) + t57;
t153 = -pkin(3) * t201 + qJDD(4) + t212;
t63 = -pkin(4) * t136 - pkin(8) * t137 + t153;
t11 = t101 * t323 + t246 * t63 + t250 * t51 - t324 * t89;
t12 = -qJD(5) * t55 - t246 * t51 + t250 * t63;
t431 = t11 * t250 - t12 * t246;
t430 = -g(1) * t248 + t381;
t320 = -m(5) - m(6) - m(7);
t273 = mrSges(4,1) * t387 - mrSges(4,2) * t247;
t311 = Ifges(4,4) * t387;
t428 = qJ(2) * t273 + (-Ifges(4,1) * t247 - t311) * t461;
t351 = qJDD(3) / 0.2e1;
t427 = -qJDD(3) / 0.2e1 - t351;
t426 = -m(4) * t433 + t247 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t201);
t97 = t164 * t350 - t157;
t88 = -qJD(3) * pkin(4) - t97;
t69 = -t146 * pkin(5) + t88;
t425 = -t69 * mrSges(7,1) + t14 * mrSges(7,3);
t424 = t69 * mrSges(7,2) - t13 * mrSges(7,3);
t423 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t183 = t193 * qJD(3);
t56 = -t243 * t104 + t350 * t94;
t422 = t182 * t98 + t183 * t97 + t192 * t56 - t193 * t57;
t10 = pkin(9) * t78 + t11;
t8 = pkin(5) * t135 - pkin(9) * t77 + t12;
t2 = qJD(6) * t13 + t10 * t249 + t245 * t8;
t3 = -qJD(6) * t14 - t10 * t245 + t249 * t8;
t421 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t420 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t252 = -pkin(9) - pkin(8);
t341 = t232 * t252;
t383 = pkin(8) * t232;
t419 = -m(6) * (pkin(4) * t231 - t383) + mrSges(2,2) - m(7) * (t227 * t231 + t341) - mrSges(3,3) - t456;
t418 = qJD(1) ^ 2;
t417 = Ifges(7,4) * t414 + Ifges(7,2) * t413 + Ifges(7,6) * t401;
t416 = Ifges(7,1) * t414 + Ifges(7,4) * t413 + Ifges(7,5) * t401;
t386 = Ifges(7,4) * t84;
t36 = Ifges(7,2) * t291 + Ifges(7,6) * t165 + t386;
t412 = -t36 / 0.2e1;
t411 = t36 / 0.2e1;
t80 = Ifges(7,4) * t291;
t37 = Ifges(7,1) * t84 + Ifges(7,5) * t165 + t80;
t410 = -t37 / 0.2e1;
t409 = t37 / 0.2e1;
t408 = t77 / 0.2e1;
t407 = t78 / 0.2e1;
t406 = -t291 / 0.2e1;
t405 = t291 / 0.2e1;
t404 = -t84 / 0.2e1;
t403 = t84 / 0.2e1;
t400 = t135 / 0.2e1;
t399 = -t146 / 0.2e1;
t398 = -t147 / 0.2e1;
t397 = t147 / 0.2e1;
t396 = -t165 / 0.2e1;
t395 = t165 / 0.2e1;
t394 = -t171 / 0.2e1;
t393 = -t180 / 0.2e1;
t391 = -t181 / 0.2e1;
t384 = pkin(5) * t147;
t380 = g(3) * t232;
t374 = mrSges(5,3) * t181;
t373 = mrSges(6,3) * t146;
t372 = mrSges(6,3) * t147;
t371 = Ifges(4,4) * t247;
t370 = Ifges(6,4) * t246;
t369 = Ifges(6,4) * t250;
t363 = t180 * Ifges(5,4);
t50 = -qJDD(3) * pkin(4) - t56;
t360 = t192 * t50;
t345 = t183 * t250;
t344 = t192 * t246;
t343 = t192 * t250;
t340 = t234 * t248;
t339 = t234 * t251;
t338 = t235 * t248;
t337 = t235 * t251;
t336 = t246 * t248;
t335 = t246 * t251;
t208 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t304;
t333 = t247 * t208;
t332 = t248 * t250;
t331 = t250 * t251;
t160 = -t231 * t340 + t337;
t161 = t231 * t338 + t339;
t329 = t160 * mrSges(7,1) - t161 * mrSges(7,2);
t162 = t231 * t339 + t338;
t163 = t231 * t337 - t340;
t328 = t162 * mrSges(7,1) + t163 * mrSges(7,2);
t327 = t251 * pkin(1) + t248 * qJ(2);
t322 = qJDD(1) * mrSges(3,2);
t213 = pkin(3) * t303 + qJD(2);
t46 = -mrSges(7,1) * t291 + mrSges(7,2) * t84;
t316 = t46 + t440;
t315 = Ifges(7,5) * t27 + Ifges(7,6) * t28 + Ifges(7,3) * t132;
t314 = Ifges(6,5) * t77 + Ifges(6,6) * t78 + Ifges(6,3) * t135;
t313 = t387 * pkin(3);
t306 = t350 * pkin(3);
t301 = t323 / 0.2e1;
t237 = t251 * qJ(2);
t300 = -pkin(1) * t248 + t237;
t296 = -t318 / 0.2e1;
t295 = -t136 * mrSges(5,1) + t137 * mrSges(5,2);
t293 = (t212 + t319) * qJ(2);
t166 = t325 * t330 - t302;
t167 = qJD(3) * t205 - t321;
t106 = t243 * t166 + t167 * t350;
t113 = -pkin(4) * t182 + pkin(8) * t183 + t213;
t292 = -t106 * t246 + t250 * t113;
t70 = t250 * t130 - t141 * t246;
t105 = -t350 * t166 + t167 * t243;
t140 = -t204 * t243 - t350 * t205;
t223 = -t306 - pkin(4);
t286 = mrSges(6,1) * t246 + mrSges(6,2) * t250;
t285 = -mrSges(7,1) * t234 - mrSges(7,2) * t235;
t284 = Ifges(6,1) * t250 - t370;
t283 = -Ifges(6,2) * t246 + t369;
t282 = Ifges(6,5) * t250 - Ifges(6,6) * t246;
t49 = pkin(5) * t193 + pkin(9) * t343 + t70;
t62 = pkin(9) * t344 + t71;
t23 = -t245 * t62 + t249 * t49;
t24 = t245 * t49 + t249 * t62;
t281 = t55 * t246 + t54 * t250;
t280 = t246 * t54 - t250 * t55;
t100 = mrSges(6,1) * t171 - t372;
t99 = -mrSges(6,2) * t171 + t373;
t279 = -t250 * t100 - t246 * t99;
t278 = -t100 * t246 + t250 * t99;
t274 = t315 + t421;
t272 = Ifges(4,1) * t387 - t371;
t271 = -Ifges(4,2) * t247 + t311;
t269 = -Ifges(4,5) * t247 - Ifges(4,6) * t387;
t155 = -qJD(3) * mrSges(5,2) + t374;
t268 = -t155 - t278;
t174 = t231 * t335 + t332;
t172 = -t231 * t336 + t331;
t267 = t88 * t286;
t265 = t183 * t246 + t192 * t323;
t264 = t192 * t324 - t345;
t32 = t250 * t106 + t246 * t113 + t130 * t323 - t141 * t324;
t260 = t247 * (-Ifges(4,2) * t387 - t371);
t255 = -qJD(5) * t281 + t431;
t244 = -qJ(4) - pkin(7);
t230 = -pkin(1) * qJDD(1) + qJDD(2);
t219 = t248 * t313;
t207 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t326;
t206 = t223 - t377;
t197 = t270 * qJD(1);
t185 = Ifges(4,5) * qJD(3) + qJD(1) * t272;
t184 = Ifges(4,6) * qJD(3) + qJD(1) * t271;
t176 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t200;
t175 = t231 * t331 - t336;
t173 = t231 * t332 + t335;
t126 = -mrSges(5,1) * t181 - mrSges(5,2) * t180;
t124 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t136;
t121 = t277 * t192;
t118 = t196 * t193;
t117 = -t180 * Ifges(5,1) + t170 + t439;
t116 = -t363 + t438 + t457;
t95 = -pkin(5) * t344 + t140;
t66 = mrSges(7,1) * t165 - mrSges(7,3) * t84;
t65 = -mrSges(7,2) * t165 + mrSges(7,3) * t291;
t64 = -pkin(5) * t265 + t105;
t45 = -t138 * t192 + t183 * t196;
t43 = t119 * t429 + t277 * t183;
t33 = -qJD(5) * t71 + t292;
t31 = -t78 * pkin(5) + t50;
t30 = t77 * Ifges(6,1) + t78 * Ifges(6,4) + t135 * Ifges(6,5);
t29 = t77 * Ifges(6,4) + t78 * Ifges(6,2) + t135 * Ifges(6,6);
t22 = pkin(9) * t265 + t32;
t21 = pkin(9) * t345 - pkin(5) * t182 + (-t133 + (-pkin(9) * t192 - t130) * t246) * qJD(5) + t292;
t20 = -mrSges(7,2) * t132 + mrSges(7,3) * t28;
t19 = mrSges(7,1) * t132 - mrSges(7,3) * t27;
t16 = t249 * t40 - t355;
t15 = -t245 * t40 - t353;
t9 = -mrSges(7,1) * t28 + mrSges(7,2) * t27;
t5 = -qJD(6) * t24 + t21 * t249 - t22 * t245;
t4 = qJD(6) * t23 + t21 * t245 + t22 * t249;
t1 = [(Ifges(7,5) * t121 + Ifges(7,6) * t119) * t401 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (Ifges(6,1) * t264 + Ifges(6,4) * t265) * t397 + t121 * t416 + t119 * t417 + t43 * t409 + t45 * t411 + (t207 * t303 - t208 * t325 + t426) * t253 + t201 * t271 / 0.2e1 + t200 * t272 / 0.2e1 + qJD(3) ^ 2 * t269 / 0.2e1 + t88 * (-mrSges(6,1) * t265 + mrSges(6,2) * t264) + m(6) * (t11 * t71 + t12 * t70 + t32 * t55 + t33 * t54) + m(5) * (t106 * t98 + t141 * t57 + t153 * t224 + t203 * t213) + (Ifges(7,5) * t43 + Ifges(7,6) * t45) * t395 + t433 * mrSges(4,3) + (-m(5) * t97 + m(6) * t88 + t440) * t105 + (-m(5) * t56 + m(6) * t50 + t442) * t140 + (-t153 * mrSges(5,2) + Ifges(5,5) * t427 - t282 * t400 - t283 * t407 - t284 * t408 + t73 * t301) * t192 + (t153 * mrSges(5,1) + Ifges(6,5) * t408 + Ifges(7,5) * t414 + Ifges(5,6) * t427 + Ifges(6,6) * t407 + Ifges(7,6) * t413 + Ifges(6,3) * t400 + Ifges(7,3) * t401 + t420 + t421) * t193 + t428 * t318 + t146 * (Ifges(6,4) * t264 + Ifges(6,2) * t265) / 0.2e1 + (t270 + 0.2e1 * mrSges(3,3)) * t212 + (t11 * t344 + t12 * t343 - t264 * t54 + t265 * t55) * mrSges(6,3) + t422 * mrSges(5,3) + (Ifges(7,4) * t121 + Ifges(7,2) * t119) * t413 + (Ifges(7,1) * t121 + Ifges(7,4) * t119) * t414 + m(7) * (t13 * t5 + t14 * t4 + t2 * t24 + t23 * t3 + t31 * t95 + t64 * t69) + (Ifges(7,1) * t43 + Ifges(7,4) * t45) * t403 + t171 * (Ifges(6,5) * t264 + Ifges(6,6) * t265) / 0.2e1 + (t119 * t2 - t121 * t3 - t13 * t43 + t14 * t45) * mrSges(7,3) + (Ifges(7,4) * t43 + Ifges(7,2) * t45) * t405 + m(4) * t293 + (-t335 * t415 - t173 * mrSges(6,1) - t161 * mrSges(7,1) - t172 * mrSges(6,2) - t160 * mrSges(7,2) + t402 * t327 + t320 * (t248 * t238 - t244 * t251 + t327) + (-m(4) * pkin(7) + t423) * t251 + t419 * t248) * g(2) + (t336 * t415 - m(3) * t300 - m(4) * t237 - t175 * mrSges(6,1) - t163 * mrSges(7,1) + t174 * mrSges(6,2) + t162 * mrSges(7,2) + t320 * (t251 * t238 + t248 * t244 + t300) + (-m(4) * t253 - t423) * t248 + t419 * t251) * g(1) + (qJD(5) * t74 + t29) * t344 / 0.2e1 + (t315 + t314) * t193 / 0.2e1 + (-Ifges(5,1) * t393 - t117 / 0.2e1 - t170 / 0.2e1 + t463) * t183 + (-t446 / 0.2e1 - Ifges(7,5) * t403 - Ifges(7,6) * t405 - t458 / 0.2e1 - t459 / 0.2e1 + Ifges(5,4) * t393 - Ifges(7,3) * t395 - Ifges(6,5) * t397 + t457 / 0.2e1 + t116 / 0.2e1 - t462) * t182 + t176 * t307 - t247 * (Ifges(4,4) * t200 + Ifges(4,2) * t201 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t230 * mrSges(3,2) + t213 * t126 + qJ(2) * (-mrSges(4,1) * t201 + mrSges(4,2) * t200) + t137 * (-Ifges(5,1) * t192 - Ifges(5,4) * t193) + qJD(2) * t197 + t136 * (-Ifges(5,4) * t192 - Ifges(5,2) * t193) + (Ifges(4,5) * t387 - Ifges(4,6) * t247) * t351 + t106 * t155 + t141 * t124 + t31 * (-mrSges(7,1) * t119 + mrSges(7,2) * t121) + t95 * t9 - t286 * t360 + t32 * t99 + t33 * t100 + t5 * t66 + t69 * (-mrSges(7,1) * t45 + mrSges(7,2) * t43) + t70 * t52 + t71 * t53 - t30 * t343 / 0.2e1 + t64 * t46 + t4 * t65 - t185 * t325 / 0.2e1 - pkin(1) * t322 + t23 * t19 + t24 * t20 - t184 * t303 / 0.2e1 + t260 * t296 + (Ifges(4,1) * t200 + Ifges(4,4) * t201 + Ifges(4,5) * qJDD(3)) * t461 + m(3) * (-pkin(1) * t230 + t293) + t224 * t295; t322 + t387 * t176 - t118 * t19 - t120 * t20 + t444 * t66 + t443 * t65 + (t207 * t387 - t333) * qJD(3) + (qJ(2) * t402 - mrSges(3,3)) * t418 + (t9 + t442) * t192 + t316 * t183 + t268 * t182 + (qJD(5) * t279 + t124 + t432) * t193 + (-m(5) * t203 - t126 - t197 + t279) * qJD(1) - m(5) * t422 + m(3) * t230 + t430 * (-t320 - t402) + (-t118 * t3 - t120 * t2 + t13 * t444 + t14 * t443 + t183 * t69 + t192 * t31) * m(7) + (-qJD(1) * t281 + t182 * t280 + t183 * t88 + t193 * t255 + t360) * m(6) + t426; (t170 + t117) * t391 + (-(-t227 * t232 + t231 * t252 - t313) * t381 + t128 * t3 + t129 * t2 + t206 * t31 + t441 * t69 + t447 * t14 + t448 * t13) * m(7) + (t223 * t50 - g(1) * (pkin(8) * t437 + t219) - t111 * t88 - t54 * t58 - t55 * t59 - (-pkin(4) * t232 - pkin(8) * t231 - t313) * t381) * m(6) + (Ifges(5,1) * t181 + t363 + t446) * t180 / 0.2e1 + t447 * t65 + t448 * t66 + t50 * t287 + t196 * t416 + (Ifges(6,5) * t246 + Ifges(6,6) * t250) * t400 + (Ifges(6,2) * t250 + t370) * t407 + (Ifges(6,1) * t246 + t369) * t408 - t110 * t410 - t109 * t412 + (mrSges(6,1) * t331 + mrSges(7,1) * t337 - mrSges(6,2) * t335 - mrSges(7,2) * t339) * g(2) * t232 + (t109 * t14 - t110 * t13 - t196 * t3 - t2 * t277 + t453) * mrSges(7,3) - t277 * t417 + (Ifges(7,1) * t196 - Ifges(7,4) * t277) * t414 + t31 * (mrSges(7,1) * t277 + mrSges(7,2) * t196) + (Ifges(7,5) * t196 - Ifges(7,6) * t277) * t401 + (-m(7) * t219 + (-(-m(7) * t252 + mrSges(7,3)) * t231 - t434 * t232) * t248) * g(1) + (m(6) * t255 - t100 * t323 - t324 * t99 + t432) * t222 - t440 * t111 + t441 * t46 + (t260 / 0.2e1 - t428) * t418 + (-m(6) * (-t238 + t383) + m(5) * t238 - m(7) * (-t238 - t341) + t434 * t231 + t456) * g(3) + (Ifges(7,4) * t196 - Ifges(7,2) * t277) * t413 - (Ifges(7,1) * t403 + Ifges(7,4) * t405 + Ifges(7,5) * t395 + t409 + t424) * t138 - (Ifges(7,4) * t403 + Ifges(7,2) * t405 + Ifges(7,6) * t395 + t411 + t425) * t139 + (m(5) * t313 + mrSges(5,1) * t232 - mrSges(5,2) * t231 + t273) * t430 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + ((t243 * t57 + t350 * t56) * pkin(3) + t97 * t111 - t98 * t112 - t203 * t290) * m(5) + (-Ifges(7,1) * t110 - Ifges(7,4) * t109) * t404 + (-Ifges(7,5) * t110 - Ifges(7,6) * t109) * t396 + (-Ifges(7,4) * t110 - Ifges(7,2) * t109) * t406 - t69 * (mrSges(7,1) * t109 - mrSges(7,2) * t110) + (t146 * t283 + t147 * t284 + t171 * t282) * qJD(5) / 0.2e1 + (-g(1) * t437 + t453 - t464 * t55 + (-t323 + t346) * t54 + t431) * mrSges(6,3) + (-Ifges(6,5) * t398 - Ifges(7,5) * t404 + Ifges(5,2) * t391 - Ifges(6,6) * t399 - Ifges(7,6) * t406 - Ifges(6,3) * t394 - Ifges(7,3) * t396 + t462) * t180 + (t282 * t394 + t283 * t399 + t284 * t398 - t267 + t463) * t181 + t124 * t385 + t30 * t388 + t116 * t393 + t97 * t374 + t211 * t333 + t125 * t306 + qJD(5) * t267 + t250 * t29 / 0.2e1 + t223 * t38 + t206 * t9 + Ifges(4,6) * t201 + Ifges(4,5) * t200 - t112 * t155 + t144 * mrSges(4,1) - t145 * mrSges(4,2) + Ifges(5,5) * t137 + Ifges(5,6) * t136 + t128 * t19 + t129 * t20 - t98 * t364 - t59 * t99 - t58 * t100 + t56 * mrSges(5,1) - t57 * mrSges(5,2) - t73 * t324 / 0.2e1 + t185 * t326 / 0.2e1 - t207 * t202 + t184 * t304 / 0.2e1 - t126 * t290 + t269 * t296 + t74 * t301; -t277 * t19 + t196 * t20 + t246 * t53 + t250 * t52 + t435 * t66 + t436 * t65 + t278 * qJD(5) + t268 * t181 + t316 * t180 + t295 + (g(1) * t251 + g(2) * t248) * t320 + (t13 * t435 + t14 * t436 + t180 * t69 + t196 * t2 - t277 * t3) * m(7) + (t11 * t246 + t12 * t250 - t171 * t280 + t180 * t88) * m(6) + (-t180 * t97 - t181 * t98 + t153) * m(5); (t372 + t100) * t55 + (-mrSges(6,2) * t175 + t174 * t445 - t328) * g(2) + (t2 * t245 + t249 * t3 + (-t13 * t245 + t14 * t249) * qJD(6)) * t415 + (-Ifges(6,2) * t147 + t143 + t74) * t399 + t420 + t274 + (mrSges(6,2) * t173 + t172 * t445 - t329) * g(1) + (t246 * t415 - t285 + t286) * t380 + (t373 - t99) * t54 - (Ifges(7,4) * t404 + Ifges(7,2) * t406 + Ifges(7,6) * t396 + t412 - t425) * t84 + (Ifges(7,1) * t404 + Ifges(7,4) * t406 + Ifges(7,5) * t396 + t410 - t424) * t291 + t314 + (Ifges(6,5) * t146 - Ifges(6,6) * t147) * t394 + t73 * t397 + (Ifges(6,1) * t146 - t365) * t398 - t46 * t384 - m(7) * (t13 * t15 + t14 * t16 + t384 * t69) - t88 * (mrSges(6,1) * t147 + mrSges(6,2) * t146) - t15 * t66 - t16 * t65 + (t19 * t249 + t20 * t245 + (-t245 * t66 + t249 * t65) * qJD(6)) * pkin(5); -t69 * (mrSges(7,1) * t84 + mrSges(7,2) * t291) + (Ifges(7,1) * t291 - t386) * t404 + t36 * t403 + (Ifges(7,5) * t291 - Ifges(7,6) * t84) * t396 - t13 * t65 + t14 * t66 - g(1) * t329 - g(2) * t328 - t285 * t380 + (t13 * t291 + t14 * t84) * mrSges(7,3) + t274 + (-Ifges(7,2) * t84 + t37 + t80) * t406;];
tau  = t1;
