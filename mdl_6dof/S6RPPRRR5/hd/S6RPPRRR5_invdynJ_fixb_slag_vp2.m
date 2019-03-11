% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:16
% EndTime: 2019-03-09 02:28:29
% DurationCPUTime: 8.11s
% Computational Cost: add. (5539->534), mult. (10261->713), div. (0->0), fcn. (6208->10), ass. (0->242)
t186 = qJDD(4) + qJDD(5);
t196 = sin(qJ(6));
t200 = cos(qJ(6));
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t268 = qJD(1) * qJD(4);
t132 = qJDD(1) * t201 - t198 * t268;
t133 = -qJDD(1) * t198 - t201 * t268;
t197 = sin(qJ(5));
t339 = cos(qJ(5));
t213 = -t197 * t201 - t198 * t339;
t267 = qJD(1) * qJD(5);
t65 = t132 * t339 + t197 * t133 + t213 * t267;
t251 = qJD(1) * t339;
t278 = qJD(1) * t198;
t123 = -t197 * t278 + t201 * t251;
t187 = qJD(4) + qJD(5);
t98 = -t123 * t196 + t187 * t200;
t33 = qJD(6) * t98 + t186 * t196 + t200 * t65;
t350 = t33 / 0.2e1;
t99 = t123 * t200 + t187 * t196;
t34 = -qJD(6) * t99 + t186 * t200 - t196 * t65;
t349 = t34 / 0.2e1;
t212 = -t197 * t198 + t201 * t339;
t66 = -t197 * t132 + t133 * t339 - t212 * t267;
t63 = qJDD(6) - t66;
t348 = t63 / 0.2e1;
t318 = t123 * mrSges(6,3);
t303 = mrSges(6,1) * t187 + mrSges(7,1) * t98 - mrSges(7,2) * t99 - t318;
t17 = mrSges(7,1) * t63 - mrSges(7,3) * t33;
t18 = -mrSges(7,2) * t63 + mrSges(7,3) * t34;
t229 = -t196 * t17 + t200 * t18;
t270 = qJD(6) * t200;
t271 = qJD(6) * t196;
t277 = qJD(1) * t201;
t122 = -t197 * t277 - t198 * t251;
t117 = qJD(6) - t122;
t68 = -mrSges(7,2) * t117 + mrSges(7,3) * t98;
t69 = mrSges(7,1) * t117 - mrSges(7,3) * t99;
t375 = -t69 * t270 - t68 * t271 + t229;
t193 = qJ(4) + qJ(5);
t179 = sin(t193);
t238 = t198 * mrSges(5,1) + t201 * mrSges(5,2);
t180 = cos(t193);
t313 = t180 * mrSges(6,2);
t374 = t179 * mrSges(6,1) + t238 + t313;
t324 = mrSges(7,2) * t196;
t220 = mrSges(6,2) * t179 + t180 * t324;
t373 = Ifges(7,1) * t350 + Ifges(7,4) * t349 + Ifges(7,5) * t348;
t344 = -m(5) - m(4);
t372 = -m(6) - m(7);
t371 = -mrSges(4,2) - mrSges(3,3);
t370 = -mrSges(4,3) + mrSges(3,2);
t10 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t54 = mrSges(6,1) * t186 - mrSges(6,3) * t65;
t369 = t10 - t54;
t236 = mrSges(7,1) * t196 + mrSges(7,2) * t200;
t177 = qJD(1) * qJ(2) + qJD(3);
t152 = -qJD(1) * pkin(7) + t177;
t134 = t201 * t152;
t109 = -pkin(8) * t277 + t134;
t106 = qJD(4) * pkin(4) + t109;
t108 = (-pkin(8) * qJD(1) + t152) * t198;
t290 = t197 * t108;
t70 = t106 * t339 - t290;
t59 = -t187 * pkin(5) - t70;
t368 = t236 * t59;
t323 = mrSges(6,3) * t122;
t102 = -mrSges(6,2) * t187 + t323;
t226 = t196 * t69 - t200 * t68;
t219 = t102 - t226;
t199 = sin(qJ(1));
t326 = mrSges(7,1) * t200;
t264 = t180 * t326;
t294 = t180 * t199;
t296 = t179 * t199;
t366 = -mrSges(6,1) * t294 - mrSges(7,3) * t296 - t199 * t264;
t202 = cos(qJ(1));
t293 = t180 * t202;
t295 = t179 * t202;
t365 = -mrSges(6,1) * t293 - mrSges(7,3) * t295 - t202 * t264;
t168 = t180 * pkin(9);
t241 = pkin(5) * t179 - t168;
t189 = qJD(1) * qJD(2);
t364 = qJDD(1) * qJ(2) + t189;
t184 = t198 * pkin(4);
t195 = pkin(1) + qJ(3);
t243 = -t195 - t184;
t337 = pkin(4) * t201;
t363 = -m(6) * t337 + t220;
t249 = t339 * qJD(5);
t362 = t339 * qJD(4) + t249;
t253 = t339 * t108;
t71 = t197 * t106 + t253;
t60 = t187 * pkin(9) + t71;
t153 = qJD(1) * t195 - qJD(2);
t125 = pkin(4) * t278 + t153;
t67 = -pkin(5) * t122 - pkin(9) * t123 + t125;
t20 = -t196 * t60 + t200 * t67;
t21 = t196 * t67 + t200 * t60;
t360 = -t20 * t196 + t21 * t200;
t142 = qJDD(3) + t364;
t129 = -qJDD(1) * pkin(7) + t142;
t274 = qJD(4) * t198;
t93 = t201 * t129 - t152 * t274;
t273 = qJD(4) * t201;
t94 = t198 * t129 + t152 * t273;
t359 = -t198 * t94 - t201 * t93;
t358 = g(1) * t202 + g(2) * t199;
t272 = qJD(5) * t197;
t73 = qJDD(4) * pkin(4) - pkin(8) * t132 + t93;
t76 = pkin(8) * t133 + t94;
t15 = t106 * t249 - t108 * t272 + t197 * t73 + t339 * t76;
t12 = pkin(9) * t186 + t15;
t269 = qJD(1) * qJD(3);
t130 = qJDD(1) * t195 - qJDD(2) + t269;
t100 = -pkin(4) * t133 + t130;
t19 = -pkin(5) * t66 - pkin(9) * t65 + t100;
t3 = -qJD(6) * t21 - t12 * t196 + t19 * t200;
t333 = t196 * t3;
t357 = -t20 * t270 - t21 * t271 - t333;
t16 = -qJD(5) * t71 - t197 * t76 + t339 * t73;
t91 = -t197 * t273 - t198 * t362 - t201 * t272;
t92 = -t197 * t274 - t198 * t272 + t201 * t362;
t356 = t15 * t213 - t16 * t212 - t70 * t91 - t71 * t92;
t2 = qJD(6) * t20 + t12 * t200 + t19 * t196;
t355 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t354 = -m(5) * t359 + t201 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t132) + t198 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t133);
t167 = t180 * mrSges(7,3);
t353 = mrSges(2,1) - t167 - t370 + t374;
t352 = m(5) * pkin(7) + mrSges(2,2) + mrSges(5,3) + mrSges(6,3) + t371;
t204 = qJD(1) ^ 2;
t351 = m(7) * pkin(5);
t347 = -t98 / 0.2e1;
t346 = -t99 / 0.2e1;
t345 = t99 / 0.2e1;
t343 = -t117 / 0.2e1;
t340 = t123 / 0.2e1;
t338 = pkin(4) * t197;
t332 = t2 * t200;
t329 = t99 * Ifges(7,4);
t194 = qJ(2) - pkin(7);
t327 = pkin(8) - t194;
t322 = Ifges(5,4) * t198;
t321 = Ifges(5,4) * t201;
t320 = Ifges(7,4) * t196;
t319 = Ifges(7,4) * t200;
t317 = t123 * Ifges(6,4);
t13 = -t186 * pkin(5) - t16;
t316 = t212 * t13;
t306 = t200 * t91;
t302 = qJDD(1) * pkin(1);
t301 = t122 * t196;
t300 = t122 * t200;
t299 = t212 * t196;
t298 = t212 * t200;
t140 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t277;
t297 = t140 * t198;
t292 = t196 * t199;
t291 = t196 * t202;
t288 = t199 * t200;
t287 = t200 * t202;
t283 = t179 * t324 + t167;
t282 = pkin(5) * t294 + pkin(9) * t296;
t281 = pkin(5) * t293 + pkin(9) * t295;
t279 = t202 * pkin(1) + t199 * qJ(2);
t276 = qJD(2) * t198;
t275 = qJD(2) * t201;
t155 = pkin(4) * t273 + qJD(3);
t266 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t63;
t265 = t339 * pkin(4);
t262 = pkin(4) * t277;
t260 = t344 + t372;
t43 = t98 * Ifges(7,2) + t117 * Ifges(7,6) + t329;
t255 = -t196 * t43 / 0.2e1;
t254 = t202 * qJ(3) + t279;
t138 = t327 * t201;
t248 = t270 / 0.2e1;
t247 = -t268 / 0.2e1;
t246 = (t198 ^ 2 + t201 ^ 2) * t152;
t242 = -t66 * mrSges(6,1) + t65 * mrSges(6,2);
t89 = pkin(5) * t123 - pkin(9) * t122;
t240 = -t133 * mrSges(5,1) + t132 * mrSges(5,2);
t239 = mrSges(5,1) * t201 - mrSges(5,2) * t198;
t235 = t201 * Ifges(5,1) - t322;
t234 = Ifges(7,1) * t200 - t320;
t233 = -t198 * Ifges(5,2) + t321;
t232 = -Ifges(7,2) * t196 + t319;
t231 = -Ifges(5,5) * t198 - Ifges(5,6) * t201;
t230 = Ifges(7,5) * t200 - Ifges(7,6) * t196;
t228 = t196 * t21 + t20 * t200;
t225 = -t196 * t68 - t200 * t69;
t86 = -pkin(5) * t213 - pkin(9) * t212 - t243;
t137 = t327 * t198;
t96 = -t137 * t339 - t197 * t138;
t45 = -t196 * t96 + t200 * t86;
t46 = t196 * t86 + t200 * t96;
t222 = qJD(3) * t153 + t130 * t195;
t217 = -t196 * t91 - t212 * t270;
t216 = -t212 * t271 + t306;
t214 = t197 * t137 - t138 * t339;
t211 = t153 * t239;
t210 = t198 * (-Ifges(5,2) * t201 - t322);
t209 = t201 * (-Ifges(5,1) * t198 - t321);
t208 = t274 * t327 + t275;
t207 = -qJD(6) * t228 - t333;
t206 = t207 + t332;
t116 = Ifges(6,4) * t122;
t42 = t99 * Ifges(7,5) + t98 * Ifges(7,6) + t117 * Ifges(7,3);
t97 = Ifges(7,4) * t98;
t44 = t99 * Ifges(7,1) + t117 * Ifges(7,5) + t97;
t6 = t33 * Ifges(7,4) + t34 * Ifges(7,2) + t63 * Ifges(7,6);
t81 = t122 * Ifges(6,2) + t187 * Ifges(6,6) + t317;
t82 = t123 * Ifges(6,1) + t187 * Ifges(6,5) + t116;
t205 = (t117 * t230 + t232 * t98 + t234 * t99) * qJD(6) / 0.2e1 - (Ifges(6,1) * t122 - t317 + t42) * t123 / 0.2e1 - (-Ifges(6,2) * t123 + t116 + t82) * t122 / 0.2e1 + (Ifges(7,6) * t123 + t122 * t232) * t347 + (Ifges(7,5) * t196 + Ifges(7,6) * t200) * t348 + (Ifges(7,2) * t200 + t320) * t349 + (Ifges(7,1) * t196 + t319) * t350 + mrSges(7,3) * t332 + t81 * t340 + (Ifges(7,3) * t123 + t122 * t230) * t343 + (Ifges(7,5) * t123 + t122 * t234) * t346 + t70 * t323 - t21 * (-mrSges(7,2) * t123 - mrSges(7,3) * t301) - t20 * (mrSges(7,1) * t123 - mrSges(7,3) * t300) - t125 * (mrSges(6,1) * t123 + mrSges(6,2) * t122) - t122 * t368 + (-t300 / 0.2e1 + t248) * t44 + (t368 + t255) * qJD(6) + t196 * t373 + t13 * (t324 - t326) + t43 * t301 / 0.2e1 - t15 * mrSges(6,2) + t16 * mrSges(6,1) + Ifges(6,5) * t65 + Ifges(6,6) * t66 + Ifges(6,3) * t186 - t187 * (Ifges(6,5) * t122 - Ifges(6,6) * t123) / 0.2e1 + t200 * t6 / 0.2e1;
t203 = -pkin(8) - pkin(7);
t183 = t202 * qJ(2);
t178 = qJDD(2) - t302;
t171 = -t265 - pkin(5);
t139 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t278;
t131 = t238 * qJD(1);
t121 = Ifges(5,5) * qJD(4) + qJD(1) * t235;
t120 = Ifges(5,6) * qJD(4) + qJD(1) * t233;
t115 = t179 * t287 - t292;
t114 = -t179 * t291 - t288;
t113 = -t179 * t288 - t291;
t112 = t179 * t292 - t287;
t107 = -qJD(4) * t138 + t276;
t88 = -mrSges(6,1) * t122 + mrSges(6,2) * t123;
t80 = t89 + t262;
t75 = t109 * t339 - t290;
t74 = t109 * t197 + t253;
t55 = -mrSges(6,2) * t186 + mrSges(6,3) * t66;
t47 = pkin(5) * t92 - pkin(9) * t91 + t155;
t40 = qJD(5) * t214 + t107 * t339 + t197 * t208;
t32 = t196 * t89 + t200 * t70;
t31 = -t196 * t70 + t200 * t89;
t28 = t196 * t80 + t200 * t75;
t27 = -t196 * t75 + t200 * t80;
t9 = -qJD(6) * t46 - t196 * t40 + t200 * t47;
t8 = qJD(6) * t45 + t196 * t47 + t200 * t40;
t1 = [-(-t100 * mrSges(6,2) - t65 * Ifges(6,1) - t66 * Ifges(6,4) - t186 * Ifges(6,5) - t230 * t348 - t232 * t349 - t234 * t350 + t248 * t43) * t212 - (Ifges(7,3) * t348 + Ifges(7,6) * t349 + Ifges(7,5) * t350 + t266 / 0.2e1 - Ifges(6,4) * t65 - Ifges(6,2) * t66 + t100 * mrSges(6,1) - Ifges(6,6) * t186 + t355) * t213 + t130 * t238 + (-m(6) * t70 + m(7) * t59 - t303) * (qJD(5) * t96 + t197 * t107 - t208 * t339) + m(3) * (-pkin(1) * t178 + (t364 + t189) * qJ(2)) + m(6) * (-t100 * t243 + t125 * t155 + t15 * t96 + t40 * t71) - t243 * t242 + (mrSges(4,3) * t195 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + m(5) * (qJD(2) * t246 + t222) + (t139 * t273 - t140 * t274 + t354) * t194 + t356 * mrSges(6,3) + t359 * mrSges(5,3) + (t142 + t364) * mrSges(4,2) + (-t302 + t178) * mrSges(3,2) + (Ifges(6,1) * t91 - Ifges(6,4) * t92) * t340 + (Ifges(7,1) * t216 + Ifges(7,4) * t217 + Ifges(7,5) * t92) * t345 + t140 * t275 + t139 * t276 - t198 * (Ifges(5,4) * t132 + Ifges(5,2) * t133) / 0.2e1 - (qJD(6) * t44 + t6) * t299 / 0.2e1 + t195 * t240 + t133 * t233 / 0.2e1 + t132 * t235 / 0.2e1 + (-t2 * t299 - t20 * t216 + t21 * t217 - t298 * t3) * mrSges(7,3) + m(4) * (qJ(2) * t142 + qJD(2) * t177 + t222) + t98 * (Ifges(7,4) * t216 + Ifges(7,2) * t217 + Ifges(7,6) * t92) / 0.2e1 + t117 * (Ifges(7,5) * t216 + Ifges(7,6) * t217 + Ifges(7,3) * t92) / 0.2e1 + t59 * (-mrSges(7,1) * t217 + mrSges(7,2) * t216) + t20 * mrSges(7,1) * t92 - t21 * mrSges(7,2) * t92 + m(7) * (t2 * t46 + t20 * t9 + t21 * t8 + t3 * t45) + t201 * (Ifges(5,1) * t132 + Ifges(5,4) * t133) / 0.2e1 + (t269 + t130) * mrSges(4,3) + t298 * t373 + (t211 + t231 * qJD(4) / 0.2e1) * qJD(4) + (-m(3) * t279 - t115 * mrSges(7,1) - t114 * mrSges(7,2) + t344 * t254 + t372 * (t202 * t184 + t199 * t203 + t254) + t352 * t199 + (-m(7) * t241 - t353) * t202) * g(2) + (-t113 * mrSges(7,1) - t112 * mrSges(7,2) + t372 * (t202 * t203 + t183) + (-m(3) + t344) * t183 + t352 * t202 + (-m(7) * (-t241 + t243) - m(6) * t243 + m(3) * pkin(1) - t344 * t195 + t353) * t199) * g(1) + t236 * t316 + t44 * t306 / 0.2e1 - t120 * t273 / 0.2e1 - t121 * t274 / 0.2e1 + t209 * t268 / 0.2e1 + 0.2e1 * t364 * mrSges(3,3) + t45 * t17 + t46 * t18 + t8 * t68 + t9 * t69 + t91 * t82 / 0.2e1 + t92 * t42 / 0.2e1 - t92 * t81 / 0.2e1 + t96 * t55 + t210 * t247 + t91 * t255 + t40 * t102 + t122 * (Ifges(6,4) * t91 - Ifges(6,2) * t92) / 0.2e1 + t125 * (mrSges(6,1) * t92 + mrSges(6,2) * t91) + qJD(3) * t131 + t155 * t88 + t187 * (Ifges(6,5) * t91 - Ifges(6,6) * t92) / 0.2e1 - (-m(6) * t16 + m(7) * t13 + t369) * t214 + qJDD(4) * (Ifges(5,5) * t201 - Ifges(5,6) * t198); -t200 * t17 - t196 * t18 - t303 * t123 + t370 * qJDD(1) + t226 * qJD(6) + (-m(3) * qJ(2) + t371) * t204 + t219 * t122 + (-m(4) * t177 - t198 * t139 - t201 * t140) * qJD(1) + m(3) * t178 - m(4) * t130 - t240 - t242 + (-g(1) * t199 + g(2) * t202) * (m(3) - t260) + (-t117 * t360 + t123 * t59 - t2 * t196 - t3 * t200) * m(7) + (t122 * t71 - t123 * t70 - t100) * m(6) + (-qJD(1) * t246 - t130) * m(5); qJDD(1) * mrSges(4,2) - t204 * mrSges(4,3) + t303 * t91 - t369 * t212 + (t139 * t201 - t297) * qJD(4) + t219 * t92 - (qJD(6) * t225 + t229 + t55) * t213 - m(6) * t356 + m(4) * t142 + m(7) * (-t206 * t213 + t360 * t92 - t59 * t91 - t316) + (-m(6) * t125 - m(7) * t228 + t153 * t344 - t131 + t225 - t88) * qJD(1) + t358 * t260 + t354; (t363 * t202 + t365) * g(1) + (t363 * t199 + t366) * g(2) + (-t209 / 0.2e1 + t210 / 0.2e1) * t204 + (-m(7) * (-t184 - t241) + t179 * t326 - t283 + m(6) * t184 + t374) * g(3) + (m(7) * t206 + t375) * (pkin(9) + t338) - t303 * (pkin(4) * t272 - t74) + ((t339 * t16 + t15 * t197 + (-t197 * t70 + t339 * t71) * qJD(5)) * pkin(4) - t125 * t262 + t70 * t74 - t71 * t75) * m(6) + (-t20 * t27 - t21 * t28 - t59 * t74 + t13 * t171 + (t197 * t59 + t339 * t360) * qJD(5) * pkin(4) - g(1) * (t202 * t337 + t281) - g(2) * (t199 * t337 + t282)) * m(7) + t55 * t338 + t71 * t318 + t152 * t297 + t219 * pkin(4) * t249 + t54 * t265 - qJD(1) * t211 - t358 * t239 + t357 * mrSges(7,3) + Ifges(5,3) * qJDD(4) - t139 * t134 + t120 * t277 / 0.2e1 + t121 * t278 / 0.2e1 - t88 * t262 - t28 * t68 - t27 * t69 + t93 * mrSges(5,1) - t94 * mrSges(5,2) + t231 * t247 - t75 * t102 + Ifges(5,5) * t132 + Ifges(5,6) * t133 + t171 * t10 + t205; -pkin(5) * t10 + t207 * mrSges(7,3) + (-m(7) * t281 + t202 * t220 + t365) * g(1) + (-m(7) * t282 + t199 * t220 + t366) * g(2) + (-m(7) * t168 + t313 + (mrSges(6,1) + t326 + t351) * t179 - t283) * g(3) + (m(7) * (t332 + t357) + t375) * pkin(9) - t13 * t351 - m(7) * (t20 * t31 + t21 * t32 + t59 * t71) + (t303 + t318) * t71 - t32 * t68 - t31 * t69 - t70 * t102 + t205; -t59 * (mrSges(7,1) * t99 + mrSges(7,2) * t98) + (Ifges(7,1) * t98 - t329) * t346 + t43 * t345 + (Ifges(7,5) * t98 - Ifges(7,6) * t99) * t343 - t20 * t68 + t21 * t69 - g(1) * (mrSges(7,1) * t114 - mrSges(7,2) * t115) - g(2) * (-mrSges(7,1) * t112 + mrSges(7,2) * t113) + g(3) * t236 * t180 + (t20 * t98 + t21 * t99) * mrSges(7,3) + t266 + (-Ifges(7,2) * t99 + t44 + t97) * t347 + t355;];
tau  = t1;
