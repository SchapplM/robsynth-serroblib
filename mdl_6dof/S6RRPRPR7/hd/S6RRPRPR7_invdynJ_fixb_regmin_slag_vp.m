% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:47:48
% EndTime: 2019-03-09 10:48:00
% DurationCPUTime: 5.65s
% Computational Cost: add. (5792->458), mult. (12836->581), div. (0->0), fcn. (9185->12), ass. (0->254)
t207 = qJDD(2) - qJDD(4);
t208 = qJD(2) - qJD(4);
t220 = sin(qJ(6));
t224 = cos(qJ(6));
t225 = cos(qJ(4));
t226 = cos(qJ(2));
t221 = sin(qJ(4));
t222 = sin(qJ(2));
t331 = t222 * t221;
t132 = t226 * t225 + t331;
t116 = t132 * qJD(1);
t310 = qJD(1) * t226;
t311 = qJD(1) * t222;
t118 = -t221 * t310 + t225 * t311;
t216 = sin(pkin(10));
t217 = cos(pkin(10));
t261 = -t216 * t116 + t217 * t118;
t247 = t132 * qJD(4);
t303 = qJD(1) * qJD(2);
t292 = t222 * t303;
t302 = t226 * qJDD(1);
t369 = t292 - t302;
t195 = t222 * qJDD(1);
t291 = t226 * t303;
t370 = t291 + t195;
t50 = -qJD(1) * t247 + t221 * t369 + t225 * t370;
t306 = qJD(4) * t225;
t307 = qJD(4) * t221;
t308 = qJD(2) * t226;
t382 = t221 * t308 + t222 * t306 - t226 * t307;
t51 = qJD(1) * t382 + qJDD(1) * t132 - t225 * t292;
t30 = -t216 * t51 + t217 * t50;
t304 = qJD(6) * t224;
t305 = qJD(6) * t220;
t16 = -t220 * t207 - t208 * t304 + t224 * t30 - t261 * t305;
t262 = t220 * t208 - t224 * t261;
t17 = -qJD(6) * t262 + t224 * t207 + t220 * t30;
t66 = t217 * t116 + t216 * t118;
t384 = qJD(6) + t66;
t374 = t262 * t384;
t52 = t224 * t208 + t220 * t261;
t376 = t384 * t52;
t393 = (t16 - t376) * t224 + (-t17 + t374) * t220;
t350 = t262 * t261;
t29 = -t216 * t50 - t217 * t51;
t28 = qJDD(6) - t29;
t344 = t220 * t28;
t386 = t224 * t384;
t368 = -t384 * t386 - t344;
t390 = t350 - t368;
t189 = pkin(7) * t311;
t389 = -pkin(8) * t311 + qJD(3) + t189;
t346 = t16 * t220;
t388 = t262 * t386 - t346;
t190 = pkin(7) * t310;
t141 = -pkin(8) * t310 + t190;
t228 = -pkin(2) - pkin(3);
t278 = -t221 * qJ(3) + t225 * t228;
t385 = qJD(4) * t278 - t221 * t141 + t389 * t225;
t143 = t225 * qJ(3) + t221 * t228;
t383 = qJD(4) * t143 + t225 * t141 + t389 * t221;
t212 = qJD(2) * qJ(3);
t119 = t141 + t212;
t295 = t228 * qJD(2);
t96 = t295 + t389;
t263 = t225 * t119 + t221 * t96;
t184 = pkin(7) * t195;
t290 = pkin(7) * t291 + qJDD(3) + t184;
t73 = -pkin(8) * t370 + t228 * qJDD(2) + t290;
t185 = pkin(7) * t302;
t210 = qJDD(2) * qJ(3);
t211 = qJD(2) * qJD(3);
t95 = -pkin(7) * t292 + t185 + t210 + t211;
t75 = pkin(8) * t369 + t95;
t381 = -t263 * qJD(4) - t221 * t75 + t225 * t73;
t353 = t118 * pkin(4);
t380 = pkin(5) * t261 + t66 * pkin(9) + t353;
t283 = -t221 * t119 + t225 * t96;
t337 = t118 * qJ(5);
t46 = t283 - t337;
t43 = -t208 * pkin(4) + t46;
t339 = t116 * qJ(5);
t47 = t263 - t339;
t44 = t217 * t47;
t21 = t216 * t43 + t44;
t19 = -t208 * pkin(9) + t21;
t120 = -qJD(1) * pkin(1) - pkin(2) * t310 - qJ(3) * t311;
t94 = pkin(3) * t310 - t120;
t62 = t116 * pkin(4) + qJD(5) + t94;
t24 = pkin(5) * t66 - pkin(9) * t261 + t62;
t7 = -t220 * t19 + t224 * t24;
t379 = t7 * t261;
t8 = t224 * t19 + t220 * t24;
t378 = t8 * t261;
t377 = -t337 + t385;
t351 = t52 * t261;
t375 = -t116 * t208 + t50;
t373 = t384 * t261;
t372 = t339 - t383;
t223 = sin(qJ(1));
t324 = t226 * t221;
t258 = -t222 * t225 + t324;
t103 = t258 * t223;
t227 = cos(qJ(1));
t323 = t226 * t227;
t298 = t221 * t323;
t329 = t222 * t227;
t105 = -t225 * t329 + t298;
t354 = g(3) * t132;
t371 = -g(1) * t105 - g(2) * t103 + t94 * t118 - t354 - t381;
t314 = t226 * pkin(2) + t222 * qJ(3);
t366 = -pkin(1) - t314;
t364 = t118 * t208 + t51;
t356 = pkin(7) - pkin(8);
t149 = t356 * t222;
t150 = t356 * t226;
t317 = t221 * t149 + t225 * t150;
t206 = g(1) * t227;
t363 = g(2) * t223 + t206;
t176 = t216 * pkin(4) + pkin(9);
t12 = -t207 * pkin(4) - t50 * qJ(5) - t118 * qJD(5) + t381;
t250 = t119 * t307 - t221 * t73 - t225 * t75 - t96 * t306;
t15 = -t51 * qJ(5) - t116 * qJD(5) - t250;
t3 = t217 * t12 - t216 * t15;
t1 = t207 * pkin(5) - t3;
t209 = qJ(4) + pkin(10);
t193 = sin(t209);
t194 = cos(t209);
t259 = t222 * t193 + t226 * t194;
t325 = t226 * t193;
t330 = t222 * t223;
t244 = g(1) * (t193 * t323 - t194 * t329) + g(2) * (-t194 * t330 + t223 * t325) + g(3) * t259 - t1;
t360 = (qJD(6) * t176 + t380) * t384 - t244;
t181 = qJ(3) * t310;
t296 = qJD(1) * t228;
t102 = t222 * t296 + t181;
t138 = -pkin(4) + t278;
t83 = t216 * t138 + t217 * t143;
t81 = -pkin(9) + t83;
t359 = (qJD(6) * t81 + t102 - t380) * t384 + t244;
t342 = pkin(7) * qJDD(2);
t358 = (qJD(1) * t366 + t120) * qJD(2) - t342;
t309 = qJD(2) * t222;
t140 = t356 * t309;
t142 = qJD(2) * t150;
t237 = -qJD(4) * t317 + t221 * t140 + t225 * t142;
t85 = qJD(2) * t132 - t247;
t231 = -t85 * qJ(5) + qJD(5) * t258 + t237;
t248 = -t225 * t140 + t221 * t142 + t149 * t306 - t150 * t307;
t84 = -t225 * t309 + t382;
t32 = -t84 * qJ(5) - t132 * qJD(5) + t248;
t10 = t216 * t231 + t217 * t32;
t345 = t216 * t47;
t20 = t217 * t43 - t345;
t18 = t208 * pkin(5) - t20;
t4 = t216 * t12 + t217 * t15;
t288 = -t207 * pkin(9) + qJD(6) * t24 + t4;
t128 = t226 * pkin(3) - t366;
t254 = t132 * pkin(4) + t128;
t76 = t217 * t132 - t216 * t258;
t77 = -t216 * t132 - t217 * t258;
t36 = t76 * pkin(5) - t77 * pkin(9) + t254;
t275 = t225 * t149 - t221 * t150;
t251 = qJ(5) * t258 + t275;
t61 = -t132 * qJ(5) + t317;
t38 = t216 * t251 + t217 * t61;
t41 = -t216 * t84 + t217 * t85;
t357 = t1 * t77 + t18 * t41 - t38 * t28 - (qJD(6) * t36 + t10) * t384 - t288 * t76 + t206;
t205 = g(1) * t223;
t355 = g(2) * t227;
t352 = t18 * t77;
t348 = t216 * t377 - t217 * t372;
t347 = t216 * t372 + t217 * t377;
t343 = t220 * t384;
t26 = t224 * t28;
t341 = qJD(6) * t19;
t215 = qJDD(1) * pkin(1);
t340 = qJDD(2) * pkin(2);
t336 = t118 * t116;
t230 = qJD(1) ^ 2;
t328 = t222 * t230;
t327 = t223 * t226;
t183 = t225 * pkin(4) + pkin(3);
t326 = t226 * t183;
t131 = t216 * t225 + t217 * t221;
t320 = t208 * t131;
t130 = t216 * t221 - t217 * t225;
t319 = t208 * t130;
t196 = t222 * qJD(3);
t315 = qJ(3) * t308 + t196;
t213 = t222 ^ 2;
t214 = t226 ^ 2;
t312 = t213 - t214;
t301 = pkin(4) * t331;
t300 = t77 * t305;
t299 = t226 * t328;
t297 = -g(1) * t329 - g(2) * t330 + g(3) * t226;
t293 = t116 ^ 2 - t118 ^ 2;
t289 = t205 - t355;
t279 = -qJD(2) * pkin(2) + qJD(3);
t273 = t208 ^ 2;
t272 = t227 * pkin(1) + pkin(2) * t323 + t223 * pkin(7) + qJ(3) * t329;
t271 = -t184 - t297;
t270 = t222 * t295;
t90 = t259 * t223;
t269 = g(1) * t90 + t36 * t28;
t229 = qJD(2) ^ 2;
t268 = pkin(7) * t229 + t355;
t266 = t28 * t77 + t384 * t41;
t265 = -t341 - t355;
t264 = pkin(2) * t302 + qJ(3) * t370 + qJD(1) * t196 + t215;
t82 = t217 * t138 - t216 * t143;
t146 = t189 + t279;
t148 = t190 + t212;
t260 = t146 * t226 - t148 * t222;
t257 = -t305 * t384 - t343 * t66 + t26;
t256 = qJD(6) * t131 + t311;
t101 = t290 - t340;
t253 = -0.2e1 * pkin(1) * t303 - t342;
t87 = t270 + t315;
t245 = -t268 + 0.2e1 * t215;
t242 = g(2) * t90 + g(3) * (t222 * t194 - t325) - t288;
t23 = t217 * t46 - t345;
t241 = -t176 * t28 + (t18 + t23) * t384;
t240 = t84 * pkin(4) + t87;
t238 = -t81 * t28 + (-t18 - t347) * t384;
t109 = pkin(2) * t309 - t315;
t72 = pkin(2) * t292 - t264;
t236 = -qJD(1) * t109 - qJDD(1) * t366 - t268 - t72;
t55 = pkin(3) * t302 + qJD(1) * t270 + t264;
t235 = qJD(2) * t260 + t101 * t222 + t95 * t226;
t233 = t51 * pkin(4) + qJDD(5) + t55;
t104 = t132 * t223;
t106 = t132 * t227;
t232 = g(1) * t106 + g(2) * t104 - g(3) * t258 + t94 * t116 + t250;
t219 = -qJ(5) - pkin(8);
t201 = t227 * pkin(7);
t177 = -t217 * pkin(4) - pkin(5);
t172 = g(1) * t327;
t166 = qJ(3) * t323;
t164 = qJ(3) * t327;
t137 = pkin(2) * t311 - t181;
t92 = t259 * t227;
t80 = pkin(5) - t82;
t79 = -t223 * t220 + t92 * t224;
t78 = -t92 * t220 - t223 * t224;
t40 = t216 * t85 + t217 * t84;
t37 = t216 * t61 - t217 * t251;
t22 = t216 * t46 + t44;
t13 = t40 * pkin(5) - t41 * pkin(9) + t240;
t9 = t216 * t32 - t217 * t231;
t6 = -t29 * pkin(5) - t30 * pkin(9) + t233;
t5 = t224 * t6;
t2 = [qJDD(1), t289, t363, t213 * qJDD(1) + 0.2e1 * t222 * t291, 0.2e1 * t222 * t302 - 0.2e1 * t303 * t312, qJDD(2) * t222 + t229 * t226, qJDD(2) * t226 - t229 * t222, 0, t222 * t253 + t226 * t245 + t172, t253 * t226 + (-t245 - t205) * t222, t222 * t358 + t236 * t226 + t172 (t213 + t214) * qJDD(1) * pkin(7) + t235 - t363, -t358 * t226 + (t236 + t205) * t222, pkin(7) * t235 - g(1) * t201 - g(2) * t272 + t120 * t109 + (-t205 + t72) * t366, t118 * t85 - t258 * t50, -t85 * t116 - t118 * t84 - t50 * t132 + t258 * t51, t207 * t258 - t85 * t208, t132 * t207 + t84 * t208, 0, g(1) * t104 - g(2) * t106 + t87 * t116 + t128 * t51 + t55 * t132 - t207 * t275 - t208 * t237 + t94 * t84, -g(1) * t103 + g(2) * t105 + t87 * t118 + t128 * t50 + t207 * t317 + t208 * t248 - t258 * t55 + t94 * t85, -t10 * t66 - t20 * t41 - t21 * t40 + t261 * t9 + t38 * t29 - t3 * t77 + t37 * t30 - t4 * t76 + t363, t4 * t38 + t21 * t10 - t3 * t37 - t20 * t9 + t233 * t254 + t62 * t240 - g(1) * (t227 * t219 + t201) - g(2) * (t183 * t323 + t227 * t301 + t272) + (-g(1) * (t366 - t301 - t326) - g(2) * t219) * t223, t262 * t300 + (t16 * t77 - t262 * t41) * t224 (t220 * t262 - t224 * t52) * t41 + (-t346 - t17 * t224 + (t220 * t52 + t224 * t262) * qJD(6)) * t77, t16 * t76 + t224 * t266 - t262 * t40 - t300 * t384, -t304 * t384 * t77 - t17 * t76 - t220 * t266 - t52 * t40, t28 * t76 + t384 * t40, -g(2) * t79 + t37 * t17 + t7 * t40 + t5 * t76 + t9 * t52 + (t13 * t384 + (-t19 * t76 - t38 * t384 + t352) * qJD(6) + t269) * t224 + t357 * t220, -g(2) * t78 + t37 * t16 - t8 * t40 - t9 * t262 + (-(-qJD(6) * t38 + t13) * t384 - (t6 - t341) * t76 - qJD(6) * t352 - t269) * t220 + t357 * t224; 0, 0, 0, -t299, t312 * t230, t195, t302, qJDD(2), pkin(1) * t328 + t271, g(3) * t222 - t185 + (pkin(1) * t230 + t363) * t226, 0.2e1 * t340 - qJDD(3) + (-t120 * t222 + t137 * t226) * qJD(1) + t271 (-pkin(2) * t222 + qJ(3) * t226) * qJDD(1) + ((t148 - t212) * t222 + (-t146 + t279) * t226) * qJD(1), t185 + 0.2e1 * t210 + 0.2e1 * t211 + (qJD(1) * t137 - g(3)) * t222 + (qJD(1) * t120 - t363) * t226, t95 * qJ(3) + t148 * qJD(3) - t101 * pkin(2) - t120 * t137 - g(1) * (-pkin(2) * t329 + t166) - g(2) * (-pkin(2) * t330 + t164) - g(3) * t314 - t260 * qJD(1) * pkin(7), -t336, t293, -t375, t364, t207, -t102 * t116 - t278 * t207 + t208 * t383 + t371, -t102 * t118 + t143 * t207 + t208 * t385 - t232, t29 * t83 - t30 * t82 + (-t21 + t348) * t261 + (t20 - t347) * t66, t4 * t83 + t3 * t82 - t62 * (t181 - t353) - g(1) * (pkin(4) * t298 + t166) - g(2) * (pkin(4) * t223 * t324 + t164) - g(3) * (t314 + t326) + t347 * t21 - t348 * t20 + (-g(3) * pkin(4) * t221 - t62 * t296 + t363 * (pkin(2) + t183)) * t222, t388, -t393, -t390, t343 * t384 - t26 - t351, t373, t80 * t17 + t238 * t220 - t224 * t359 + t348 * t52 + t379, t80 * t16 + t220 * t359 + t238 * t224 - t262 * t348 - t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t299, t195, -t213 * t230 - t229, -t148 * qJD(2) + t120 * t311 + t101 + t297, 0, 0, 0, 0, 0, -t116 * t311 - t225 * t207 - t221 * t273, -t118 * t311 + t221 * t207 - t225 * t273, t130 * t30 + t131 * t29 - t261 * t320 - t319 * t66, -t3 * t130 + t4 * t131 + t20 * t320 + t21 * t319 - t311 * t62 + t297, 0, 0, 0, 0, 0, -t131 * t344 + t130 * t17 - t320 * t52 + (-t220 * t319 - t224 * t256) * t384, -t131 * t26 + t130 * t16 + t320 * t262 + (t220 * t256 - t224 * t319) * t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t336, -t293, t375, -t364, -t207, -t208 * t263 - t371, -t208 * t283 + t232 (t216 * t29 - t217 * t30) * pkin(4) - (t20 - t23) * t66 + (t21 - t22) * t261, t20 * t22 - t21 * t23 + (-t62 * t118 + t4 * t216 + t3 * t217 + t258 * t363 + t354) * pkin(4), -t388, t393, t390, t257 + t351, -t373, t177 * t17 - t22 * t52 + t241 * t220 - t224 * t360 - t379, t177 * t16 + t22 * t262 + t220 * t360 + t241 * t224 + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261 ^ 2 - t66 ^ 2, t20 * t261 + t21 * t66 + t233 + t289, 0, 0, 0, 0, 0, t257 - t351, t350 + t368; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262 * t52, t262 ^ 2 - t52 ^ 2, t16 + t376, -t17 - t374, t28, -g(1) * t78 + t18 * t262 + t220 * t242 + t224 * t265 + t384 * t8 + t5, g(1) * t79 + t18 * t52 + t7 * t384 + (-t265 - t6) * t220 + t242 * t224;];
tau_reg  = t2;
