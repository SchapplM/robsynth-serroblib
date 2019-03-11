% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:12:05
% EndTime: 2019-03-09 10:12:29
% DurationCPUTime: 11.93s
% Computational Cost: add. (10527->556), mult. (27546->724), div. (0->0), fcn. (20689->8), ass. (0->276)
t399 = -Ifges(6,6) - Ifges(5,4);
t217 = sin(qJ(2));
t303 = sin(pkin(10));
t304 = cos(pkin(10));
t335 = cos(qJ(2));
t248 = t217 * t303 - t304 * t335;
t334 = cos(qJ(4));
t238 = t334 * t248;
t169 = qJD(1) * t238;
t247 = -t304 * t217 - t303 * t335;
t185 = t247 * qJD(1);
t216 = sin(qJ(4));
t134 = t185 * t216 - t169;
t214 = qJD(2) + qJD(4);
t215 = sin(qJ(6));
t218 = cos(qJ(6));
t114 = -t134 * t215 + t214 * t218;
t347 = t114 / 0.2e1;
t403 = Ifges(7,5) * t347;
t112 = -t134 * t218 - t214 * t215;
t402 = t112 / 0.2e1;
t244 = qJD(1) * t248;
t237 = t216 * t244;
t225 = -t185 * t334 - t237;
t129 = qJD(6) + t225;
t401 = t129 / 0.2e1;
t345 = t134 / 0.2e1;
t342 = -t225 / 0.2e1;
t400 = Ifges(5,1) * t342;
t338 = -t214 / 0.2e1;
t398 = Ifges(6,3) + Ifges(5,2);
t397 = Ifges(6,2) * t342;
t387 = Ifges(7,6) * t402 + Ifges(7,3) * t401;
t350 = pkin(4) + pkin(9);
t382 = pkin(5) * t225;
t286 = t335 * pkin(7);
t204 = qJ(3) * t335 + t286;
t198 = t204 * qJD(1);
t186 = t303 * t198;
t203 = (-qJ(3) - pkin(7)) * t217;
t197 = qJD(1) * t203;
t192 = qJD(2) * pkin(2) + t197;
t137 = t304 * t192 - t186;
t331 = t185 * pkin(8);
t106 = qJD(2) * pkin(3) + t137 + t331;
t271 = t304 * t198;
t138 = t303 * t192 + t271;
t241 = pkin(8) * t244;
t109 = -t241 + t138;
t63 = -t334 * t106 + t109 * t216;
t256 = t63 + t382;
t389 = qJD(5) + t256;
t38 = -t214 * t350 + t389;
t277 = qJD(1) * t335;
t143 = -qJD(1) * pkin(1) - pkin(2) * t277 + pkin(3) * t244 + qJD(3);
t65 = -pkin(4) * t134 - qJ(5) * t225 + t143;
t46 = -pkin(9) * t134 + t65;
t14 = -t215 * t46 + t218 * t38;
t15 = t215 * t38 + t218 * t46;
t337 = t214 / 0.2e1;
t395 = -t345 * t399 + t403;
t374 = -qJD(5) - t63;
t61 = -pkin(4) * t214 - t374;
t388 = -t61 * mrSges(6,1) - t14 * mrSges(7,1) - t143 * mrSges(5,2) + t15 * mrSges(7,2) - t63 * mrSges(5,3) + t65 * mrSges(6,3) + Ifges(6,4) * t337 + Ifges(5,5) * t338 - t387 - t395 + t397 + t400;
t396 = t387 - t388;
t344 = -t134 / 0.2e1;
t394 = mrSges(6,2) - mrSges(5,1);
t393 = mrSges(5,3) + mrSges(6,1);
t141 = -t197 * t303 - t271;
t115 = t241 + t141;
t142 = t304 * t197 - t186;
t116 = t142 + t331;
t282 = pkin(2) * t303;
t207 = t216 * t282;
t281 = t304 * pkin(2);
t209 = t281 + pkin(3);
t275 = qJD(4) * t334;
t376 = -qJD(4) * t207 - t216 * t115 - t116 * t334 + t209 * t275;
t332 = pkin(5) * t134;
t381 = -Ifges(5,5) + Ifges(6,4);
t380 = Ifges(6,5) - Ifges(5,6);
t392 = Ifges(4,4) * t247;
t377 = -qJD(5) - t376;
t302 = qJ(5) * t134;
t259 = t14 * t215 - t15 * t218;
t182 = qJD(2) * t203 + qJD(3) * t335;
t158 = t182 * qJD(1);
t183 = -qJD(2) * t204 - t217 * qJD(3);
t159 = t183 * qJD(1);
t111 = -t158 * t303 + t304 * t159;
t243 = qJD(2) * t248;
t236 = qJD(1) * t243;
t222 = pkin(8) * t236 + t111;
t64 = t216 * t106 + t109 * t334;
t113 = t304 * t158 + t303 * t159;
t242 = qJD(2) * t247;
t235 = qJD(1) * t242;
t97 = pkin(8) * t235 + t113;
t28 = qJD(4) * t64 + t216 * t97 - t334 * t222;
t228 = qJD(2) * t238;
t293 = qJD(4) * t216;
t87 = qJD(1) * t228 + qJD(4) * t169 - t185 * t293 - t216 * t235;
t11 = -t87 * pkin(5) + t28;
t290 = qJD(1) * qJD(2);
t274 = t217 * t290;
t208 = pkin(2) * t274;
t144 = -pkin(3) * t235 + t208;
t227 = t334 * t242;
t88 = -qJD(1) * t227 - qJD(4) * t237 - t185 * t275 - t216 * t236;
t31 = t88 * pkin(4) + t87 * qJ(5) - qJD(5) * t225 + t144;
t12 = t88 * pkin(9) + t31;
t1 = qJD(6) * t14 + t11 * t215 + t12 * t218;
t301 = qJD(6) * t15;
t2 = t11 * t218 - t12 * t215 - t301;
t368 = -t1 * t215 - t2 * t218;
t229 = m(7) * (-qJD(6) * t259 - t368);
t55 = qJD(6) * t112 + t215 * t88;
t32 = -mrSges(7,1) * t87 - mrSges(7,3) * t55;
t56 = -qJD(6) * t114 + t218 * t88;
t33 = mrSges(7,2) * t87 + mrSges(7,3) * t56;
t391 = t215 * t33 + t218 * t32 + t229;
t264 = mrSges(7,1) * t218 - mrSges(7,2) * t215;
t62 = -t214 * qJ(5) - t64;
t39 = -t62 + t332;
t317 = Ifges(7,4) * t114;
t51 = t112 * Ifges(7,2) + t129 * Ifges(7,6) + t317;
t390 = t39 * t264 - t218 * t51 / 0.2e1;
t27 = t106 * t275 - t109 * t293 + t216 * t222 + t334 * t97;
t386 = m(5) * t27 - t88 * mrSges(5,3);
t385 = t143 * mrSges(5,1) + t62 * mrSges(6,1) - t65 * mrSges(6,2) - t64 * mrSges(5,3) + Ifges(6,5) * t337 + Ifges(5,6) * t338 - t342 * t399 + t344 * t398;
t357 = t55 / 0.2e1;
t356 = t56 / 0.2e1;
t352 = -t87 / 0.2e1;
t383 = pkin(4) * t225;
t378 = t382 - t377;
t375 = t225 * t350;
t117 = -mrSges(5,2) * t214 + mrSges(5,3) * t134;
t119 = -mrSges(6,1) * t134 - mrSges(6,3) * t214;
t373 = -t119 + t117;
t372 = -t214 * t394 - t225 * t393;
t295 = qJD(1) * t217;
t371 = t217 * pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t277) + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t295) * t286;
t291 = qJD(6) * t218;
t370 = -t225 * t218 - t291;
t367 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t366 = -m(5) * t64 + m(6) * t62 - t373;
t365 = m(5) * t63 + m(6) * t61 - t372;
t260 = t14 * t218 + t15 * t215;
t363 = -m(7) * t260 - t365;
t261 = Ifges(7,5) * t215 + Ifges(7,6) * t218;
t316 = Ifges(7,4) * t215;
t262 = Ifges(7,2) * t218 + t316;
t315 = Ifges(7,4) * t218;
t263 = Ifges(7,1) * t215 + t315;
t336 = -t215 / 0.2e1;
t346 = -t129 / 0.2e1;
t348 = -t114 / 0.2e1;
t349 = -t112 / 0.2e1;
t110 = Ifges(7,4) * t112;
t52 = Ifges(7,1) * t114 + Ifges(7,5) * t129 + t110;
t360 = t261 * t346 + t262 * t349 + t263 * t348 + t52 * t336 - t385 + t390;
t358 = Ifges(7,1) * t357 + Ifges(7,4) * t356 + Ifges(7,5) * t352;
t341 = t225 / 0.2e1;
t339 = -t185 / 0.2e1;
t145 = t304 * t203 - t204 * t303;
t125 = pkin(8) * t247 + t145;
t146 = t303 * t203 + t304 * t204;
t126 = -pkin(8) * t248 + t146;
t74 = -t334 * t125 + t126 * t216;
t329 = t28 * t74;
t324 = t87 * mrSges(6,1);
t323 = t87 * mrSges(5,3);
t322 = t88 * mrSges(6,1);
t319 = Ifges(3,4) * t217;
t318 = Ifges(4,4) * t185;
t310 = t185 * mrSges(4,3);
t245 = t216 * t248;
t93 = -qJD(4) * t245 - t216 * t243 - t247 * t275 - t227;
t308 = t215 * t93;
t306 = t218 * t93;
t69 = -mrSges(7,1) * t112 + mrSges(7,2) * t114;
t305 = t119 - t69;
t300 = t225 * t215;
t139 = -t216 * t247 + t238;
t298 = t139 * t215;
t297 = t139 * t218;
t124 = t304 * t182 + t303 * t183;
t210 = -pkin(2) * t335 - pkin(1);
t296 = qJD(1) * t210;
t294 = qJD(2) * t217;
t292 = qJD(6) * t215;
t288 = Ifges(7,5) * t55 + Ifges(7,6) * t56 - Ifges(7,3) * t87;
t213 = pkin(2) * t294;
t212 = pkin(2) * t295;
t285 = Ifges(3,4) * t335;
t279 = t88 * mrSges(5,1) - t87 * mrSges(5,2);
t278 = -t88 * mrSges(6,2) + t87 * mrSges(6,3);
t276 = qJD(2) * t335;
t272 = -t292 / 0.2e1;
t150 = -pkin(3) * t185 + t212;
t66 = -t334 * t115 + t116 * t216;
t268 = qJD(1) * t276;
t266 = t276 / 0.2e1;
t180 = t209 * t334 - t207;
t140 = -t247 * t334 - t245;
t161 = pkin(3) * t248 + t210;
t223 = -t140 * qJ(5) + t161;
t57 = t139 * t350 + t223;
t58 = pkin(5) * t140 + t74;
t30 = t215 * t58 + t218 * t57;
t29 = -t215 * t57 + t218 * t58;
t71 = -mrSges(7,2) * t129 + mrSges(7,3) * t112;
t72 = mrSges(7,1) * t129 - mrSges(7,3) * t114;
t258 = -t215 * t72 + t218 * t71;
t123 = -t182 * t303 + t304 * t183;
t177 = -pkin(4) - t180;
t257 = t150 - t302;
t255 = Ifges(3,5) * t335 - Ifges(3,6) * t217;
t254 = t139 * t291 + t308;
t253 = t139 * t292 - t306;
t75 = t216 * t125 + t126 * t334;
t104 = pkin(8) * t242 + t124;
t224 = pkin(8) * t243 + t123;
t34 = -t334 * t104 - t125 * t275 + t126 * t293 - t216 * t224;
t251 = pkin(1) * (mrSges(3,1) * t217 + mrSges(3,2) * t335);
t250 = t217 * (Ifges(3,1) * t335 - t319);
t249 = (Ifges(3,2) * t335 + t319) * qJD(1);
t181 = t216 * t209 + t282 * t334;
t24 = -qJD(5) * t214 - t27;
t246 = Ifges(4,4) * t248;
t240 = mrSges(4,3) * t244;
t233 = -t243 / 0.2e1;
t232 = t242 / 0.2e1;
t35 = qJD(4) * t75 + t216 * t104 - t334 * t224;
t231 = mrSges(4,3) * t236;
t230 = mrSges(4,3) * t235;
t151 = -pkin(3) * t242 + t213;
t226 = -mrSges(4,1) * t235 - mrSges(4,2) * t236;
t92 = qJD(4) * t238 - t216 * t242 - t247 * t293 + t228;
t36 = t93 * pkin(4) + t92 * qJ(5) - t140 * qJD(5) + t151;
t10 = -pkin(5) * t88 - t24;
t7 = t55 * Ifges(7,4) + t56 * Ifges(7,2) - t87 * Ifges(7,6);
t221 = -t27 * mrSges(5,2) - t24 * mrSges(6,3) + t218 * t358 + t7 * t336 + t14 * mrSges(7,3) * t292 + (Ifges(7,1) * t218 - t316) * t357 + (-Ifges(7,2) * t215 + t315) * t356 + t52 * t272 + (Ifges(7,5) * t218 - Ifges(7,6) * t215) * t352 + t10 * (mrSges(7,1) * t215 + mrSges(7,2) * t218) + t380 * t88 + t381 * t87 + t394 * t28 + t390 * qJD(6) - (t112 * t262 + t114 * t263 + t129 * t261) * qJD(6) / 0.2e1;
t220 = qJD(2) * (-Ifges(4,1) * t248 + t392);
t211 = Ifges(3,4) * t277;
t200 = qJD(3) + t296;
t191 = Ifges(3,1) * t295 + Ifges(3,5) * qJD(2) + t211;
t190 = Ifges(3,6) * qJD(2) + t249;
t179 = Ifges(4,4) * t244;
t176 = qJ(5) + t181;
t156 = qJD(2) * mrSges(4,1) + t310;
t155 = -qJD(2) * mrSges(4,2) - t240;
t136 = mrSges(4,1) * t244 - t185 * mrSges(4,2);
t131 = -Ifges(4,1) * t185 + Ifges(4,5) * qJD(2) - t179;
t130 = -Ifges(4,2) * t244 + Ifges(4,6) * qJD(2) - t318;
t91 = mrSges(6,2) * t134 - mrSges(6,3) * t225;
t90 = -mrSges(5,1) * t134 + mrSges(5,2) * t225;
t89 = -t302 + t383;
t73 = t139 * pkin(4) + t223;
t68 = t257 + t383;
t60 = -t302 + t375;
t59 = -t139 * pkin(5) + t75;
t49 = t257 + t375;
t47 = t66 + t332;
t41 = t64 + t332;
t23 = t215 * t41 + t218 * t60;
t22 = -t215 * t60 + t218 * t41;
t21 = -mrSges(7,1) * t56 + mrSges(7,2) * t55;
t20 = t93 * pkin(9) + t36;
t19 = -t92 * pkin(5) + t35;
t18 = -pkin(5) * t93 - t34;
t17 = t215 * t47 + t218 * t49;
t16 = -t215 * t49 + t218 * t47;
t4 = -qJD(6) * t30 + t19 * t218 - t20 * t215;
t3 = qJD(6) * t29 + t19 * t215 + t20 * t218;
t5 = [(t139 * t87 + t342 * t93) * Ifges(6,6) + (Ifges(7,1) * t254 - Ifges(7,4) * t253) * t347 + (t1 * t297 - t14 * t254 - t15 * t253 - t2 * t298) * mrSges(7,3) + t298 * t358 + t220 * t339 + (-m(6) * t24 - t322 + t386) * t75 + (t200 * (-mrSges(4,1) * t247 - mrSges(4,2) * t248) - t371) * qJD(2) + t52 * t308 / 0.2e1 + t51 * t306 / 0.2e1 + t191 * t266 + (t250 - 0.2e1 * t251) * t290 + t136 * t213 + m(6) * (t31 * t73 + t36 * t65 + t329) + m(5) * (t143 * t151 + t144 * t161 + t329) + (-t323 - t324) * t74 + t365 * t35 + t366 * t34 + t210 * t226 + (-Ifges(3,2) * t217 + t285) * t268 + (-Ifges(5,1) * t341 - Ifges(5,4) * t345 + t344 * Ifges(6,6) + t381 * t337 - t396 + t397 - t403) * t92 + ((-Ifges(4,1) * t247 - t246) * t233 + (-Ifges(4,2) * t248 - t392) * t232 + (Ifges(3,1) * t217 + t285) * t266 - t247 * t220 / 0.2e1) * qJD(1) - (t249 + t190) * t294 / 0.2e1 + (qJD(6) * t52 + t7) * t297 / 0.2e1 + (-Ifges(4,5) * t248 + Ifges(4,6) * t247 + t255) * qJD(2) ^ 2 / 0.2e1 - (Ifges(4,2) * t247 - t246) * t236 + t73 * t278 + t161 * t279 + (mrSges(4,1) * t248 - mrSges(4,2) * t247) * t208 + (t111 * t247 - t113 * t248 + t137 * t243 + t138 * t242) * mrSges(4,3) + m(7) * (t1 * t30 + t10 * t59 + t14 * t4 + t15 * t3 + t18 * t39 + t2 * t29) + (-Ifges(5,4) * t341 - Ifges(5,2) * t345 + Ifges(6,3) * t344 + t380 * t337 + t385) * t93 + t39 * (mrSges(7,1) * t253 + mrSges(7,2) * t254) + t151 * t90 + t124 * t155 + t123 * t156 + t36 * t91 + t18 * t69 + t3 * t71 + t4 * t72 + t59 * t21 + t29 * t32 + t30 * t33 + m(4) * (t111 * t145 + t113 * t146 + t123 * t137 + t124 * t138 + (t200 + t296) * t213) + (t144 * mrSges(5,1) + t24 * mrSges(6,1) - t31 * mrSges(6,2) - t27 * mrSges(5,3) - t10 * t264 + t261 * t352 + t262 * t356 + t263 * t357 + t51 * t272 + t398 * t88 + (-t352 + t87 / 0.2e1) * Ifges(5,4)) * t139 + (t393 * t28 + Ifges(7,6) * t356 + Ifges(7,5) * t357 - t31 * mrSges(6,3) + t144 * mrSges(5,2) + (Ifges(7,3) + Ifges(5,1)) * t352 + t367 + t288 / 0.2e1 + (-Ifges(6,2) - Ifges(5,1) / 0.2e1) * t87 + t399 * t88) * t140 + t130 * t232 + t131 * t233 + (Ifges(7,5) * t254 - Ifges(7,6) * t253) * t401 + (Ifges(7,4) * t254 - Ifges(7,2) * t253) * t402 + t145 * t231 + t146 * t230; (-t143 * t150 - t180 * t28 + t376 * t64 - t63 * t66) * m(5) + (t371 + (t251 - t250 / 0.2e1) * qJD(1)) * qJD(1) - t177 * t324 + t130 * t339 + ((t215 * t71 + t218 * t72 - t363) * qJD(4) + t386) * t181 + (t14 * t300 + t15 * t370 + t368) * mrSges(7,3) + t372 * t66 - t138 * t310 + (t291 * t71 - t292 * t72 + t391) * (-pkin(9) + t177) - Ifges(3,6) * t274 + Ifges(3,5) * t268 + (-t322 + t21) * t176 + t230 * t282 - t136 * t212 + t180 * t323 + (Ifges(5,4) * t344 + Ifges(7,5) * t348 - Ifges(6,2) * t341 - Ifges(6,6) * t345 + Ifges(7,6) * t349 + Ifges(7,3) * t346 - t381 * t338 + t388 + t400) * t134 - (-Ifges(3,2) * t295 + t191 + t211) * t277 / 0.2e1 + t185 * (-Ifges(4,1) * t244 + t318) / 0.2e1 - t200 * (-t185 * mrSges(4,1) - mrSges(4,2) * t244) - qJD(2) * (-Ifges(4,5) * t244 + Ifges(4,6) * t185) / 0.2e1 + (Ifges(4,2) * t185 + t131 - t179) * t244 / 0.2e1 + (-t176 * t24 + t177 * t28 + t377 * t62 - t61 * t66 - t65 * t68) * m(6) + t377 * t119 + (t10 * t176 - t14 * t16 - t15 * t17 + t378 * t39) * m(7) + t378 * t69 - t137 * t240 + Ifges(4,6) * t235 - Ifges(4,5) * t236 + t221 + t376 * t117 + (-t137 * t141 - t138 * t142 - t200 * t212 + (t111 * t304 + t113 * t303) * pkin(2)) * m(4) - t150 * t90 - t142 * t155 - t141 * t156 - t113 * mrSges(4,2) + t111 * mrSges(4,1) - t68 * t91 - t17 * t71 - t16 * t72 + (-mrSges(3,1) * t268 + mrSges(3,2) * t274) * pkin(7) + t190 * t295 / 0.2e1 - t255 * t290 / 0.2e1 + (-Ifges(5,4) * t342 - Ifges(5,2) * t344 + Ifges(6,6) * t341 + Ifges(6,3) * t345 + t338 * t380 + t360) * t225 + t231 * t281; t278 + t279 + m(6) * t31 + m(5) * t144 + m(7) * (-qJD(6) * t260 + t1 * t218 - t2 * t215) + t218 * t33 - t215 * t32 - t185 * t156 + t226 + t155 * t244 + t370 * t72 + (-t292 - t300) * t71 + (-t137 * t185 + t138 * t244 + t208) * m(4) + t363 * t225 + (-m(7) * t39 + t366 - t69) * t134; (-(qJD(6) * t71 + t32) * t350 + (-t2 - t301) * mrSges(7,3)) * t218 + (-t1 * mrSges(7,3) - (-qJD(6) * t72 + t33) * t350) * t215 + t373 * t63 + ((Ifges(5,6) / 0.2e1 - Ifges(6,5) / 0.2e1) * t214 - (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t134 + t360 + (Ifges(5,4) / 0.2e1 + Ifges(6,6) / 0.2e1) * t225 + t259 * mrSges(7,3)) * t225 + (pkin(4) * t87 - qJ(5) * t88) * mrSges(6,1) + t372 * t64 - t305 * qJD(5) + t221 - t89 * t91 + t256 * t69 - t23 * t71 - t22 * t72 - ((Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * t214 + t395 + t396) * t134 - t350 * t229 + qJ(5) * t21 + (t10 * qJ(5) - t14 * t22 - t15 * t23 + t389 * t39) * m(7) + (-pkin(4) * t28 - qJ(5) * t24 + t374 * t62 - t61 * t64 - t65 * t89) * m(6); -t324 + t305 * t214 + t258 * qJD(6) + (t258 + t91) * t225 - m(7) * (t214 * t39 + t225 * t259) + (t214 * t62 + t225 * t65 + t28) * m(6) + t391; -t39 * (mrSges(7,1) * t114 + mrSges(7,2) * t112) + (Ifges(7,1) * t112 - t317) * t348 + t51 * t347 + (Ifges(7,5) * t112 - Ifges(7,6) * t114) * t346 - t14 * t71 + t15 * t72 + (t112 * t14 + t114 * t15) * mrSges(7,3) + t288 + (-Ifges(7,2) * t114 + t110 + t52) * t349 + t367;];
tauc  = t5(:);
