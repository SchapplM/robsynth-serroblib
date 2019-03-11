% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:15
% EndTime: 2019-03-08 20:18:39
% DurationCPUTime: 14.59s
% Computational Cost: add. (4164->586), mult. (8779->769), div. (0->0), fcn. (5990->10), ass. (0->271)
t420 = Ifges(6,1) + Ifges(7,1);
t403 = Ifges(6,5) + Ifges(7,4);
t418 = Ifges(6,6) - Ifges(7,6);
t176 = sin(qJ(4));
t179 = cos(qJ(4));
t181 = -pkin(2) - pkin(8);
t180 = cos(qJ(2));
t172 = sin(pkin(6));
t314 = qJD(1) * t172;
t272 = t180 * t314;
t210 = qJD(3) - t272;
t123 = qJD(2) * t181 + t210;
t174 = cos(pkin(6));
t296 = qJDD(1) * t174;
t306 = qJD(4) * t179;
t177 = sin(qJ(2));
t311 = qJD(2) * t177;
t270 = t172 * t311;
t155 = qJD(1) * t270;
t325 = t172 * t180;
t110 = qJDD(1) * t325 - t155;
t201 = qJDD(3) - t110;
t313 = qJD(1) * t174;
t415 = -qJD(4) * t313 + qJDD(2) * t181 + t201;
t19 = t123 * t306 + t176 * t415 + t179 * t296;
t307 = qJD(4) * t176;
t20 = -t123 * t307 - t176 * t296 + t179 * t415;
t80 = t123 * t179 - t176 * t313;
t425 = t176 * t19 + t179 * t20 - t307 * t80;
t297 = qJD(2) * qJD(4);
t141 = qJDD(2) * t179 - t176 * t297;
t175 = sin(qJ(5));
t178 = cos(qJ(5));
t298 = t178 * qJD(4);
t310 = qJD(2) * t179;
t136 = t175 * t310 - t298;
t304 = qJD(5) * t136;
t59 = qJDD(4) * t175 + t141 * t178 - t304;
t373 = t59 / 0.2e1;
t308 = qJD(4) * t175;
t137 = t178 * t310 + t308;
t60 = qJD(5) * t137 - t178 * qJDD(4) + t141 * t175;
t371 = t60 / 0.2e1;
t370 = -m(5) - m(4);
t408 = m(6) + m(7);
t142 = -qJDD(2) * t176 - t179 * t297;
t133 = qJDD(5) - t142;
t369 = t133 / 0.2e1;
t367 = t136 / 0.2e1;
t424 = -t137 / 0.2e1;
t299 = t176 * qJD(2);
t167 = qJD(5) + t299;
t423 = -t167 / 0.2e1;
t422 = qJD(4) / 0.2e1;
t421 = mrSges(6,3) + mrSges(7,2);
t419 = Ifges(7,2) + Ifges(6,3);
t273 = t177 * t314;
t391 = -t418 * t175 + t403 * t178;
t337 = Ifges(7,5) * t175;
t340 = Ifges(6,4) * t175;
t389 = t420 * t178 + t337 - t340;
t305 = qJD(4) * t181;
t372 = -t60 / 0.2e1;
t416 = (-Ifges(6,4) + Ifges(7,5)) * t371 + t420 * t373 + t403 * t369;
t232 = mrSges(5,1) * t176 + mrSges(5,2) * t179;
t404 = mrSges(3,2) - mrSges(4,3);
t414 = -t232 + t404;
t271 = t179 * t313;
t81 = t123 * t176 + t271;
t65 = qJD(4) * pkin(9) + t81;
t248 = -pkin(9) * t179 + qJ(3);
t360 = pkin(4) * t176;
t143 = t248 + t360;
t92 = qJD(2) * t143 + t273;
t22 = -t175 * t65 + t178 * t92;
t23 = t175 * t92 + t178 * t65;
t16 = qJDD(4) * pkin(9) + t19;
t302 = qJD(5) * t178;
t303 = qJD(5) * t175;
t294 = qJDD(2) * qJ(3);
t295 = qJDD(1) * t177;
t85 = t172 * t295 + t294 + (qJD(3) + t272) * qJD(2);
t39 = -pkin(4) * t142 - pkin(9) * t141 + t85;
t3 = t178 * t16 + t175 * t39 + t92 * t302 - t303 * t65;
t4 = -qJD(5) * t23 - t16 * t175 + t178 * t39;
t233 = -t175 * t4 + t178 * t3;
t413 = -t22 * t302 - t23 * t303 + t233;
t396 = qJD(6) - t22;
t15 = -pkin(5) * t167 + t396;
t18 = qJ(6) * t167 + t23;
t1 = qJ(6) * t133 + qJD(6) * t167 + t3;
t2 = -pkin(5) * t133 + qJDD(6) - t4;
t234 = t1 * t178 + t175 * t2;
t412 = t15 * t302 - t18 * t303 + t234;
t343 = Ifges(5,4) * t176;
t227 = t179 * Ifges(5,1) - t343;
t131 = Ifges(6,4) * t136;
t338 = Ifges(7,5) * t136;
t401 = t137 * t420 + t403 * t167 - t131 + t338;
t130 = Ifges(7,5) * t137;
t48 = Ifges(7,6) * t167 + Ifges(7,3) * t136 + t130;
t411 = Ifges(5,5) * qJD(4) + qJD(2) * t227 + t175 * t48 + t178 * t401;
t410 = t408 - t370;
t342 = Ifges(5,4) * t179;
t222 = -Ifges(5,2) * t176 + t342;
t409 = Ifges(5,6) * t422 + qJD(2) * t222 / 0.2e1 + t419 * t423 + t403 * t424 + t418 * t367;
t406 = t141 / 0.2e1;
t405 = t142 / 0.2e1;
t351 = mrSges(3,1) - mrSges(4,2);
t31 = mrSges(6,1) * t133 - mrSges(6,3) * t59;
t32 = -t133 * mrSges(7,1) + t59 * mrSges(7,2);
t349 = t32 - t31;
t33 = -mrSges(6,2) * t133 - mrSges(6,3) * t60;
t34 = -mrSges(7,2) * t60 + mrSges(7,3) * t133;
t348 = t33 + t34;
t345 = mrSges(6,3) * t136;
t88 = -mrSges(6,2) * t167 - t345;
t91 = -mrSges(7,2) * t136 + mrSges(7,3) * t167;
t347 = t88 + t91;
t344 = mrSges(6,3) * t137;
t89 = mrSges(6,1) * t167 - t344;
t90 = -mrSges(7,1) * t167 + mrSges(7,2) * t137;
t400 = t89 - t90;
t11 = mrSges(6,1) * t60 + mrSges(6,2) * t59;
t399 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t141 - t11;
t286 = mrSges(5,3) * t310;
t398 = -qJD(4) * mrSges(5,1) + mrSges(6,1) * t136 + mrSges(6,2) * t137 + t286;
t211 = pkin(5) * t175 - qJ(6) * t178;
t397 = -t271 - (-qJD(2) * t211 + t123) * t176 + qJD(5) * t211 - qJD(6) * t175;
t319 = t176 * t181;
t395 = t175 * t143 + t178 * t319;
t394 = -t176 * t391 + t179 * t419;
t393 = -t176 * t389 + t179 * t403;
t392 = t175 * t403 + t178 * t418;
t336 = Ifges(7,5) * t178;
t339 = Ifges(6,4) * t178;
t390 = t175 * t420 - t336 + t339;
t386 = t133 * t419 + t403 * t59 - t418 * t60;
t383 = -t59 * Ifges(7,5) / 0.2e1 - t133 * Ifges(7,6) / 0.2e1 + Ifges(6,4) * t373 + Ifges(6,6) * t369 + (Ifges(7,3) + Ifges(6,2)) * t372;
t382 = mrSges(5,3) + t351;
t212 = pkin(5) * t178 + qJ(6) * t175;
t229 = -t178 * mrSges(7,1) - t175 * mrSges(7,3);
t231 = mrSges(6,1) * t178 - mrSges(6,2) * t175;
t381 = -m(7) * t212 - mrSges(5,1) + t229 - t231;
t235 = pkin(4) * t179 + pkin(9) * t176;
t134 = qJD(4) * t235 + qJD(3);
t379 = -qJD(5) * t395 + t134 * t178;
t245 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t64 = -qJD(4) * pkin(4) - t80;
t378 = -m(6) * t64 - t398;
t244 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t279 = t176 * t325;
t122 = t174 * t179 - t279;
t171 = sin(pkin(10));
t173 = cos(pkin(10));
t323 = t174 * t180;
t119 = t171 * t323 + t173 * t177;
t326 = t172 * t179;
t70 = t119 * t176 + t171 * t326;
t117 = t171 * t177 - t173 * t323;
t72 = -t117 * t176 + t173 * t326;
t376 = -g(1) * t70 + g(2) * t72 - g(3) * t122;
t375 = qJ(3) * t370 + t179 * t421 - t248 * t408 + t414;
t374 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t182 = qJD(2) ^ 2;
t368 = -t136 / 0.2e1;
t365 = t137 / 0.2e1;
t363 = t167 / 0.2e1;
t353 = -qJD(2) / 0.2e1;
t352 = qJD(5) / 0.2e1;
t341 = Ifges(6,4) * t137;
t17 = -qJDD(4) * pkin(4) - t20;
t335 = t17 * t179;
t139 = t235 * qJD(2);
t42 = t175 * t139 + t178 * t80;
t140 = qJD(2) * qJ(3) + t273;
t330 = t140 * t180;
t329 = t143 * t178;
t328 = t172 * t176;
t327 = t172 * t177;
t324 = t174 * t177;
t322 = t175 * t176;
t320 = t176 * t178;
t318 = t177 * t178;
t315 = pkin(2) * t325 + qJ(3) * t327;
t309 = qJD(2) * t180;
t301 = qJD(5) * t179;
t300 = qJDD(2) * mrSges(4,2);
t10 = mrSges(7,1) * t60 - mrSges(7,3) * t59;
t293 = t10 - t399;
t67 = mrSges(7,1) * t136 - mrSges(7,3) * t137;
t290 = t67 + t398;
t287 = mrSges(5,3) * t299;
t281 = t176 * t327;
t280 = t177 * t326;
t266 = t179 * t305;
t277 = t175 * t134 + t143 * t302 + t178 * t266;
t276 = pkin(8) * t325 + t315;
t269 = t172 * t309;
t268 = t175 * t307;
t265 = t181 * t303;
t263 = qJD(1) * t309;
t254 = t302 / 0.2e1;
t253 = -t301 / 0.2e1;
t251 = t175 * t181 - pkin(5);
t108 = t117 * pkin(2);
t250 = -pkin(8) * t117 - t108;
t109 = t119 * pkin(2);
t249 = -pkin(8) * t119 - t109;
t246 = -t297 / 0.2e1;
t243 = t176 * t273;
t240 = t172 * t263;
t230 = mrSges(6,1) * t175 + mrSges(6,2) * t178;
t228 = mrSges(7,1) * t175 - mrSges(7,3) * t178;
t221 = -Ifges(6,2) * t175 + t339;
t220 = Ifges(6,2) * t178 + t340;
t217 = -Ifges(5,5) * t176 - Ifges(5,6) * t179;
t214 = Ifges(7,3) * t175 + t336;
t213 = -Ifges(7,3) * t178 + t337;
t209 = t15 * t178 - t18 * t175;
t207 = t23 * t175 + t22 * t178;
t41 = t139 * t178 - t175 * t80;
t206 = qJ(3) * t85 + qJD(3) * t140;
t205 = -t181 + t211;
t204 = t140 * t309 + t177 * t85;
t203 = -t122 * t175 + t172 * t318;
t76 = t122 * t178 + t175 * t327;
t121 = t174 * t176 + t179 * t325;
t200 = t140 * (mrSges(5,1) * t179 - mrSges(5,2) * t176);
t199 = t176 * (-Ifges(5,2) * t179 - t343);
t198 = t179 * (-Ifges(5,1) * t176 - t342);
t97 = (t175 * t180 + t176 * t318) * t172;
t192 = -t175 * t347 - t178 * t400;
t186 = Ifges(6,6) * t179 - t176 * t221;
t185 = Ifges(7,6) * t179 - t176 * t214;
t151 = -qJD(4) * mrSges(5,2) - t287;
t144 = -pkin(4) - t212;
t138 = t232 * qJD(2);
t135 = -qJD(2) * pkin(2) + t210;
t120 = -t171 * t324 + t173 * t180;
t118 = t171 * t180 + t173 * t324;
t113 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t142;
t111 = (t263 + t295) * t172;
t98 = t205 * t179;
t95 = -qJDD(2) * pkin(2) + t201;
t93 = -t175 * t319 + t329;
t87 = qJD(1) * t97;
t86 = t175 * t243 - t178 * t272;
t83 = t176 * t251 - t329;
t82 = qJ(6) * t176 + t395;
t78 = -mrSges(5,1) * t142 + mrSges(5,2) * t141;
t74 = -qJD(4) * t279 + t174 * t306 - t179 * t270;
t73 = -qJD(4) * t121 + t176 * t270;
t71 = t117 * t179 + t173 * t328;
t69 = t119 * t179 - t171 * t328;
t66 = pkin(5) * t137 + qJ(6) * t136;
t51 = -Ifges(6,2) * t136 + Ifges(6,6) * t167 + t341;
t40 = (qJD(5) * t212 - qJD(6) * t178) * t179 - t205 * t307;
t38 = -t175 * t266 + t379;
t37 = -t176 * t265 + t277;
t36 = -pkin(5) * t310 - t41;
t35 = qJ(6) * t310 + t42;
t29 = -t118 * t178 - t175 * t72;
t27 = -t120 * t178 + t175 * t70;
t25 = t251 * t306 - t379;
t24 = pkin(5) * t136 - qJ(6) * t137 + t64;
t21 = qJ(6) * t306 + (qJD(6) - t265) * t176 + t277;
t13 = qJD(5) * t203 + t175 * t269 + t178 * t73;
t12 = qJD(5) * t76 + t175 * t73 - t178 * t269;
t5 = pkin(5) * t60 - qJ(6) * t59 - qJD(6) * t137 + t17;
t6 = [t122 * t113 + t73 * t151 + t348 * t76 - t349 * t203 + t347 * t13 - t400 * t12 + t290 * t74 + t293 * t121 + (-m(2) - m(3) - t410) * g(3) + m(5) * (-t121 * t20 + t122 * t19 + t73 * t81 - t74 * t80) + m(7) * (t1 * t76 + t12 * t15 + t121 * t5 + t13 * t18 - t2 * t203 + t24 * t74) + m(6) * (-t12 * t22 + t121 * t17 + t13 * t23 + t203 * t4 + t3 * t76 + t64 * t74) + (t138 * t309 + t177 * t78 + m(5) * t204 + m(3) * (t110 * t180 + t111 * t177) + m(4) * (t135 * t311 - t180 * t95 + t204) + (-t177 * t351 - t180 * t404) * t182 + (-t177 * t404 + t180 * t351) * qJDD(2)) * t172 + (m(2) + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t174 ^ 2) * qJDD(1); (t425 * t181 + t206 - (t176 * t81 * t177 + t330) * t314) * m(5) - t425 * mrSges(5,3) + (t95 - t155) * mrSges(4,2) + (-t243 + t266) * t151 - qJDD(4) * Ifges(5,6) * t176 + t386 * t176 / 0.2e1 + (-t213 * t367 - t220 * t368 - t363 * t392 - t365 * t390) * t301 + (qJD(2) * qJD(3) - t240 + t294 + t85) * mrSges(4,3) + (-t111 + t240) * mrSges(3,2) + (-mrSges(6,1) * t64 - mrSges(7,1) * t24 + mrSges(7,2) * t18 + mrSges(6,3) * t23) * (-t178 * t301 + t268) + (-mrSges(6,2) * t64 - mrSges(7,2) * t15 + mrSges(6,3) * t22 + mrSges(7,3) * t24) * (t175 * t301 + t176 * t298) + (t3 * t395 + t4 * t93 + (t307 * t64 - t335) * t181 + (t37 - t87) * t23 + (t38 + t86) * t22) * m(6) + t395 * t33 + t198 * t297 / 0.2e1 - pkin(2) * t300 + t210 * t138 + t398 * t176 * t305 + t401 * t175 * t253 + (t110 + t155) * mrSges(3,1) + (t268 / 0.2e1 + t178 * t253) * t51 + ((-t273 * t80 + t305 * t81) * m(5) + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) - t383) * t175 + t399 * t181 + (mrSges(7,2) * t2 - mrSges(6,3) * t4 + t416) * t178 + (m(7) * t24 - t378 + t67) * t273 + Ifges(5,5) * qJDD(4) + Ifges(5,1) * t406 + Ifges(5,4) * t405 + t214 * t371 + t221 * t372 + t5 * t228 + t48 * t254 + t391 * t369 + t389 * t373) * t179 + t85 * t232 + t222 * t405 + t227 * t406 + t230 * t335 + t113 * t319 + t400 * t86 - t347 * t87 + t40 * t67 + (t22 * mrSges(6,1) - t15 * mrSges(7,1) - t23 * mrSges(6,2) - t81 * mrSges(5,3) + t18 * mrSges(7,3) - t409) * t306 - t411 * t307 / 0.2e1 + (m(4) * t109 - m(5) * t249 - t408 * (t120 * t360 + t249) - t245 * (-t119 * t175 + t120 * t320) + t244 * (t119 * t178 + t120 * t322) + t382 * t119 + t375 * t120) * g(1) + (m(4) * t108 - m(5) * t250 - t408 * (t118 * t360 + t250) - t245 * (-t117 * t175 + t118 * t320) + t244 * (t117 * t178 + t118 * t322) + t382 * t117 + t375 * t118) * g(2) + qJ(3) * t78 + t82 * t34 + t83 * t32 + t37 * t88 + t38 * t89 + t25 * t90 + t21 * t91 + t93 * t31 + t98 * t10 + t199 * t246 + (-Ifges(5,4) * t141 / 0.2e1 - Ifges(5,2) * t142 / 0.2e1 + Ifges(7,6) * t371 + Ifges(6,6) * t372 + t403 * t373 + t419 * t369 + t374) * t176 + (-m(4) * t315 - m(5) * t276 - t408 * (pkin(4) * t281 - pkin(9) * t280 + t276) - t245 * t97 + t244 * (t175 * t281 - t178 * t325) + t421 * t280 + (t177 * t414 - t382 * t180) * t172) * g(3) + (t185 * t367 + t186 * t368 + t217 * t422 + t394 * t363 + t393 * t365 + t200) * qJD(4) + (Ifges(3,3) + Ifges(4,1)) * qJDD(2) + (t1 * t82 + t2 * t83 + t24 * t40 + t5 * t98 + (t21 - t87) * t18 + (t25 - t86) * t15) * m(7) + (-pkin(2) * t95 + t206 - (t135 * t177 + t330) * t314) * m(4); t300 - t182 * mrSges(4,3) + m(4) * t95 + (-m(6) * t207 + m(7) * t209 + t140 * t370 - t138 + t192) * qJD(2) + ((-t175 * t400 + t178 * t347 + t151) * qJD(4) + m(5) * (qJD(4) * t81 + t20) + m(6) * (-t22 * t308 + t23 * t298 - t17) + m(7) * (t15 * t308 + t18 * t298 - t5) - t293) * t179 + (t113 + t348 * t178 + t349 * t175 + t290 * qJD(4) + t192 * qJD(5) + m(5) * (-qJD(4) * t80 + t19) + m(6) * (qJD(4) * t64 + t413) + m(7) * (qJD(4) * t24 + t412)) * t176 + t410 * (-g(1) * t119 - g(2) * t117 + g(3) * t325); (t352 * t389 + t353 * t393) * t137 + (pkin(9) * t348 + t383) * t178 + (t286 + t378) * t81 - t19 * mrSges(5,2) + t20 * mrSges(5,1) - pkin(4) * t11 + (-pkin(4) * t17 + (-qJD(5) * t207 + t233) * pkin(9) - t22 * t41 - t23 * t42) * m(6) + ((t186 / 0.2e1 - t185 / 0.2e1) * t136 - t200 - t15 * (-mrSges(7,1) * t179 - mrSges(7,2) * t320) - t22 * (mrSges(6,1) * t179 + mrSges(6,3) * t320) - t23 * (-mrSges(6,2) * t179 + mrSges(6,3) * t322) - t18 * (mrSges(7,2) * t322 + mrSges(7,3) * t179)) * qJD(2) + t390 * t373 + t392 * t369 + (t199 / 0.2e1 - t198 / 0.2e1) * t182 + (t352 * t391 + t353 * t394 + t228 * t24 + t230 * t64 - t175 * t51 / 0.2e1) * t167 + (-t221 / 0.2e1 + t214 / 0.2e1) * t304 - t400 * pkin(9) * t302 + (t48 / 0.2e1 - t347 * pkin(9)) * t303 + (-t287 - t151) * t80 + t5 * t229 - t17 * t231 + Ifges(5,3) * qJDD(4) + t213 * t371 + t220 * t372 + (-t15 * t36 - t18 * t35 + t144 * t5 + (qJD(5) * t209 + t234) * pkin(9) + t397 * t24) * m(7) + t397 * t67 + t401 * t254 + t409 * t310 + t411 * t299 / 0.2e1 + (t376 + t412) * mrSges(7,2) + (t376 + t413) * mrSges(6,3) + (mrSges(5,2) * t122 - t408 * (-t121 * pkin(4) + pkin(9) * t122) - t381 * t121) * g(3) + (-mrSges(5,2) * t72 - t408 * (t71 * pkin(4) - pkin(9) * t72) + t381 * t71) * g(2) + (mrSges(5,2) * t70 - t408 * (t69 * pkin(4) + pkin(9) * t70) + t381 * t69) * g(1) - t42 * t88 - t41 * t89 - t36 * t90 - t35 * t91 + Ifges(5,5) * t141 + Ifges(5,6) * t142 + t144 * t10 + t217 * t246 + (t349 * pkin(9) + t416) * t175; -pkin(5) * t32 + qJ(6) * t34 + (-Ifges(6,2) * t137 - t131 + t401) * t367 + (t136 * t15 + t137 * t18) * mrSges(7,2) + (-pkin(5) * t2 + qJ(6) * t1 - t15 * t23 + t18 * t396 - t24 * t66) * m(7) + t374 + (t344 + t400) * t23 + (-t345 - t347) * t22 + (-t403 * t136 - t137 * t418) * t423 + (-t136 * t420 + t130 - t341 + t48) * t424 + (-t203 * t245 + t244 * t76) * g(3) + (t244 * (t120 * t175 + t178 * t70) + t245 * t27) * g(1) + (t244 * (t118 * t175 - t178 * t72) + t245 * t29) * g(2) + t51 * t365 + (Ifges(7,3) * t137 - t338) * t368 - t66 * t67 + qJD(6) * t91 - t24 * (mrSges(7,1) * t137 + mrSges(7,3) * t136) - t64 * (mrSges(6,1) * t137 - mrSges(6,2) * t136) + t386; t137 * t67 - t167 * t91 + (-g(1) * t27 - g(2) * t29 + g(3) * t203 + t24 * t137 - t18 * t167 + t2) * m(7) + t32;];
tau  = t6;
