% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:01:52
% EndTime: 2019-03-09 10:01:59
% DurationCPUTime: 4.86s
% Computational Cost: add. (6215->530), mult. (12768->648), div. (0->0), fcn. (7911->10), ass. (0->284)
t219 = sin(qJ(2));
t306 = qJD(1) * qJD(2);
t282 = t219 * t306;
t222 = cos(qJ(2));
t304 = t222 * qJDD(1);
t392 = t282 - t304;
t374 = pkin(3) + pkin(7);
t322 = qJD(1) * t219;
t174 = qJD(4) + t322;
t221 = cos(qJ(4));
t218 = sin(qJ(4));
t318 = qJD(2) * t218;
t321 = qJD(1) * t222;
t131 = t221 * t321 + t318;
t283 = t218 * t321;
t316 = qJD(2) * t221;
t133 = -t283 + t316;
t215 = sin(pkin(9));
t216 = cos(pkin(9));
t78 = t216 * t131 + t133 * t215;
t391 = t174 * t78;
t313 = qJD(4) * t221;
t284 = t216 * t313;
t289 = t221 * t322;
t314 = qJD(4) * t218;
t345 = t215 * t218;
t355 = t215 * t314 - t216 * t289 + t322 * t345 - t284;
t257 = t215 * t221 + t216 * t218;
t354 = -t215 * t313 - t216 * t314 - t257 * t322;
t220 = sin(qJ(1));
t223 = cos(qJ(1));
t264 = g(1) * t223 + g(2) * t220;
t250 = t264 * t222;
t365 = g(3) * t219;
t390 = t250 + t365;
t191 = pkin(7) * t322;
t389 = qJD(3) + t191;
t259 = -t131 * t215 + t216 * t133;
t388 = t259 ^ 2;
t224 = -pkin(2) - pkin(8);
t328 = qJ(5) - t224;
t276 = t328 * t221;
t102 = -qJD(4) * t276 - t218 * qJD(5);
t309 = t221 * qJD(5);
t239 = t314 * t328 - t309;
t196 = pkin(2) * t322;
t351 = qJ(3) * t222;
t261 = pkin(8) * t219 - t351;
t106 = qJD(1) * t261 + t196;
t192 = pkin(7) * t321;
t139 = pkin(3) * t321 + t192;
t122 = t221 * t139;
t343 = t218 * t219;
t54 = -t218 * t106 + t122 + (pkin(4) * t222 - qJ(5) * t343) * qJD(1);
t353 = t221 * t106 + t218 * t139;
t61 = qJ(5) * t289 + t353;
t361 = (t239 - t54) * t216 + (-t102 + t61) * t215;
t329 = t218 * pkin(4) + qJ(3);
t200 = t219 * qJ(3);
t205 = t222 * pkin(2);
t324 = t205 + t200;
t144 = -pkin(1) - t324;
t126 = -pkin(8) * t222 + t144;
t155 = t374 * t219;
t326 = t221 * t126 + t218 * t155;
t281 = t222 * t306;
t305 = t219 * qJDD(1);
t244 = t281 + t305;
t128 = qJDD(4) + t244;
t109 = t221 * t128;
t387 = -t174 * t314 + t109;
t217 = -qJ(5) - pkin(8);
t342 = t218 * t220;
t302 = pkin(4) * t342;
t339 = t219 * t220;
t386 = t217 * t339 + t222 * t302;
t340 = t218 * t223;
t177 = pkin(4) * t340;
t337 = t219 * t223;
t385 = t222 * t177 + t217 * t337;
t312 = qJD(4) * t222;
t285 = t218 * t312;
t384 = t219 * t316 + t285;
t73 = -qJD(4) * t283 + t218 * qJDD(2) + (qJD(2) * qJD(4) - t392) * t221;
t383 = pkin(4) * t73 + qJDD(5);
t72 = -qJD(4) * t131 + t221 * qJDD(2) + t218 * t392;
t39 = t215 * t72 + t216 * t73;
t40 = -t215 * t73 + t216 * t72;
t63 = t216 * t102 + t215 * t239;
t141 = t328 * t218;
t84 = -t141 * t215 + t216 * t276;
t85 = -t216 * t141 - t215 * t276;
t382 = -t85 * t39 + t40 * t84 - t63 * t78;
t381 = pkin(5) * t39 - qJ(6) * t40 - qJD(6) * t259;
t209 = qJ(4) + pkin(9);
t198 = cos(t209);
t208 = g(3) * t222;
t173 = pkin(2) * t282;
t310 = t219 * qJD(3);
t234 = qJD(2) * t261 - t310;
t279 = -pkin(1) - t200;
t243 = t224 * t222 + t279;
t57 = qJD(1) * t234 + qJDD(1) * t243 + t173;
t172 = pkin(7) * t281;
t188 = pkin(7) * t305;
t280 = qJDD(3) + t172 + t188;
t76 = t244 * pkin(3) + t224 * qJDD(2) + t280;
t277 = -t218 * t57 + t221 * t76;
t308 = pkin(3) * t322 + t389;
t101 = t224 * qJD(2) + t308;
t94 = t243 * qJD(1);
t59 = t101 * t218 + t221 * t94;
t231 = -qJD(4) * t59 + t277;
t12 = pkin(4) * t128 - qJ(5) * t72 - qJD(5) * t133 + t231;
t303 = -t101 * t313 - t218 * t76 - t221 * t57;
t15 = -qJ(5) * t73 - qJD(5) * t131 - t314 * t94 - t303;
t3 = t216 * t12 - t215 * t15;
t295 = -qJDD(6) + t3;
t212 = qJD(2) * qJ(3);
t114 = t212 + t139;
t83 = pkin(4) * t131 + qJD(5) + t114;
t31 = pkin(5) * t78 - qJ(6) * t259 + t83;
t197 = sin(t209);
t97 = t197 * t220 - t198 * t337;
t99 = t197 * t223 + t198 * t339;
t380 = g(1) * t97 - g(2) * t99 + t198 * t208 - t31 * t259 + t295;
t379 = t114 * t174 + t224 * t128;
t187 = pkin(4) * t221 + pkin(3);
t378 = pkin(4) * t313 + t187 * t322 + t389;
t4 = t215 * t12 + t216 * t15;
t300 = t128 * qJ(6) + t4;
t1 = qJD(6) * t174 + t300;
t49 = -qJ(5) * t131 + t59;
t357 = t215 * t49;
t58 = t221 * t101 - t218 * t94;
t48 = -qJ(5) * t133 + t58;
t44 = pkin(4) * t174 + t48;
t18 = t216 * t44 - t357;
t16 = -pkin(5) * t174 + qJD(6) - t18;
t46 = t216 * t49;
t19 = t215 * t44 + t46;
t17 = qJ(6) * t174 + t19;
t370 = pkin(5) * t128;
t2 = -t295 - t370;
t256 = -t216 * t221 + t345;
t291 = -g(1) * t337 - g(2) * t339 + t208;
t377 = t1 * t257 - t354 * t16 - t355 * t17 + t2 * t256 + t291;
t376 = t354 * t18 - t355 * t19 - t256 * t3 + t257 * t4 + t291;
t375 = t174 ^ 2;
t371 = pkin(2) * t219;
t369 = g(1) * t220;
t366 = g(2) * t223;
t21 = t215 * t48 + t46;
t364 = t21 * t259;
t315 = qJD(2) * t222;
t140 = t374 * t315;
t123 = t221 * t140;
t274 = qJ(5) * t222 - t126;
t317 = qJD(2) * t219;
t195 = pkin(2) * t317;
t88 = t195 + t234;
t27 = pkin(4) * t315 + t123 + t274 * t313 + (-qJ(5) * t317 - qJD(4) * t155 + qJD(5) * t222 - t88) * t218;
t245 = -t126 * t314 + t218 * t140 + t155 * t313 + t221 * t88;
t32 = qJ(5) * t384 - t222 * t309 + t245;
t9 = t215 * t27 + t216 * t32;
t362 = -t355 * pkin(5) - t354 * qJ(6) + qJD(6) * t256 + t378;
t29 = t215 * t54 + t216 * t61;
t360 = pkin(5) * t321 - t361;
t24 = qJ(6) * t321 + t29;
t359 = t63 - t24;
t135 = t221 * t155;
t64 = t219 * pkin(4) + t218 * t274 + t135;
t333 = t221 * t222;
t69 = -qJ(5) * t333 + t326;
t37 = t215 * t64 + t216 * t69;
t356 = t72 * t221;
t352 = pkin(7) * qJDD(2);
t350 = qJD(4) * t94;
t348 = qJDD(2) * pkin(2);
t347 = t131 * t174;
t346 = t133 * t174;
t344 = t218 * t128;
t341 = t218 * t222;
t338 = t219 * t221;
t226 = qJD(1) ^ 2;
t336 = t219 * t226;
t335 = t220 * t221;
t334 = t220 * t222;
t332 = t221 * t223;
t331 = t222 * t223;
t22 = t216 * t48 - t357;
t327 = qJD(6) - t22;
t298 = t219 * t335;
t325 = pkin(4) * t298 + t177;
t156 = t374 * t222;
t213 = t219 ^ 2;
t214 = t222 ^ 2;
t323 = t213 - t214;
t320 = qJD(2) * t131;
t319 = qJD(2) * t133;
t311 = qJD(4) * t224;
t301 = g(3) * t333;
t299 = t218 * t337;
t297 = t219 * t332;
t296 = t222 * t336;
t206 = t223 * pkin(7);
t294 = t223 * t187 + t217 * t334 + t206;
t189 = pkin(7) * t304;
t210 = qJDD(2) * qJ(3);
t211 = qJD(2) * qJD(3);
t293 = t189 + t210 + t211;
t176 = pkin(4) * t333;
t292 = t176 + t156;
t290 = t374 * qJD(2);
t288 = t218 * t317;
t275 = -qJD(2) * pkin(2) + qJD(3);
t273 = pkin(4) * t285;
t272 = t223 * pkin(1) + pkin(2) * t331 + t220 * pkin(7) + qJ(3) * t337;
t271 = -t188 - t291;
t270 = pkin(3) * t304 + t293;
t269 = qJD(1) * t290;
t169 = qJ(3) * t334;
t268 = -pkin(2) * t339 + t169;
t171 = qJ(3) * t331;
t267 = -pkin(2) * t337 + t171;
t180 = g(1) * t334;
t266 = -g(2) * t331 + t180;
t225 = qJD(2) ^ 2;
t265 = pkin(7) * t225 + t366;
t263 = -t78 ^ 2 - t388;
t262 = pkin(5) * t197 - qJ(6) * t198;
t8 = -t215 * t32 + t216 * t27;
t36 = -t215 * t69 + t216 * t64;
t143 = t191 + t275;
t149 = -t192 - t212;
t258 = t143 * t222 + t149 * t219;
t255 = qJD(2) * (-pkin(7) - t187);
t254 = t174 * t218;
t253 = t279 - t205;
t252 = pkin(4) * t343 - t217 * t222 + t324;
t115 = t253 * qJD(1);
t249 = t115 * t322 + qJDD(3) - t271;
t248 = -0.2e1 * pkin(1) * t306 - t352;
t247 = -t174 * t313 - t344;
t246 = -qJ(3) * t315 - t310;
t104 = t257 * t222;
t242 = 0.2e1 * qJDD(1) * pkin(1) - t265;
t238 = t352 + (-qJD(1) * t144 - t115) * qJD(2);
t235 = pkin(4) * t299 + t220 * t187 - t217 * t331 + t272;
t82 = -t219 * t269 + t270;
t233 = t82 - t390;
t232 = t256 * t40 - t257 * t39 - t259 * t354 + t355 * t78;
t108 = t195 + t246;
t74 = qJD(1) * t246 + qJDD(1) * t253 + t173;
t229 = qJD(1) * t108 + qJDD(1) * t144 + t265 + t74;
t105 = t280 - t348;
t95 = pkin(7) * t282 - t293;
t228 = qJD(2) * t258 + t105 * t219 - t95 * t222;
t43 = t82 + t383;
t227 = (-g(3) - t269) * t219 - t250 + t270 + t383;
t186 = -pkin(4) * t216 - pkin(5);
t182 = pkin(4) * t215 + qJ(6);
t153 = pkin(4) * t297;
t138 = t219 * t290;
t136 = -qJ(3) * t321 + t196;
t119 = -t218 * t339 + t332;
t118 = t298 + t340;
t117 = t299 + t335;
t116 = t297 - t342;
t103 = t256 * t222;
t100 = -t197 * t339 + t198 * t223;
t98 = t197 * t337 + t198 * t220;
t75 = pkin(5) * t257 + qJ(6) * t256 + t329;
t66 = -t215 * t384 - t216 * t288 + t222 * t284;
t65 = qJD(4) * t104 - t256 * t317;
t51 = -pkin(5) * t103 + qJ(6) * t104 + t292;
t38 = pkin(4) * t133 + pkin(5) * t259 + qJ(6) * t78;
t35 = -pkin(5) * t219 - t36;
t34 = qJ(6) * t219 + t37;
t20 = -pkin(5) * t65 + qJ(6) * t66 + qJD(6) * t104 + t219 * t255 - t273;
t7 = -pkin(5) * t315 - t8;
t6 = qJ(6) * t315 + qJD(6) * t219 + t9;
t5 = t43 + t381;
t10 = [qJDD(1), -t366 + t369, t264, qJDD(1) * t213 + 0.2e1 * t219 * t281, 0.2e1 * t219 * t304 - 0.2e1 * t306 * t323, qJDD(2) * t219 + t222 * t225, qJDD(2) * t222 - t219 * t225, 0, t219 * t248 + t222 * t242 + t180, t248 * t222 + (-t242 - t369) * t219 (t213 + t214) * qJDD(1) * pkin(7) + t228 - t264, t219 * t238 + t222 * t229 - t180, t238 * t222 + (-t229 + t369) * t219, pkin(7) * t228 - g(1) * t206 - g(2) * t272 + t115 * t108 + t74 * t144 - t253 * t369, -t72 * t341 + (-t221 * t312 + t288) * t133 (-t131 * t218 + t133 * t221) * t317 + (t218 * t73 - t356 + (t131 * t221 + t133 * t218) * qJD(4)) * t222 (t174 * t318 + t72) * t219 + (t247 + t319) * t222 (t174 * t316 - t73) * t219 + (-t320 - t387) * t222, t128 * t219 + t174 * t315 (-t218 * t88 + t123) * t174 + (-t126 * t218 + t135) * t128 + t277 * t219 - t138 * t131 + t156 * t73 + t82 * t333 - g(1) * t119 - g(2) * t117 + (-t114 * t338 + t222 * t58) * qJD(2) + (-t114 * t341 - t174 * t326 - t59 * t219) * qJD(4), -t245 * t174 - t326 * t128 - t138 * t133 + t156 * t72 + g(1) * t118 - g(2) * t116 + ((qJD(2) * t114 + t350) * t218 + t303) * t219 + (-qJD(2) * t59 - t114 * t313 - t82 * t218) * t222, t103 * t4 + t104 * t3 + t18 * t66 + t19 * t65 - t259 * t8 - t36 * t40 - t37 * t39 - t78 * t9 + t266, t4 * t37 + t19 * t9 + t3 * t36 + t18 * t8 + t43 * t292 - t83 * t273 - g(1) * (-pkin(1) * t220 - pkin(2) * t334 + t294) - g(2) * t235 + (t255 * t83 + t329 * t369) * t219, -g(1) * t100 - g(2) * t98 - t103 * t5 - t128 * t35 - t16 * t315 - t174 * t7 - t2 * t219 + t20 * t78 - t31 * t65 + t39 * t51, t1 * t103 - t104 * t2 - t16 * t66 + t17 * t65 + t259 * t7 - t34 * t39 + t35 * t40 - t6 * t78 + t266, -g(1) * t99 - g(2) * t97 + t1 * t219 + t104 * t5 + t128 * t34 + t17 * t315 + t174 * t6 - t20 * t259 + t31 * t66 - t40 * t51, t1 * t34 + t17 * t6 + t5 * t51 + t31 * t20 + t2 * t35 + t16 * t7 - g(1) * (pkin(5) * t100 + qJ(6) * t99 + t294) - g(2) * (pkin(5) * t98 + qJ(6) * t97 + t235) - (-t219 * t329 - pkin(1) - t205) * t369; 0, 0, 0, -t296, t323 * t226, t305, t304, qJDD(2), pkin(1) * t336 + t271, t365 - t189 + (pkin(1) * t226 + t264) * t222 (t351 - t371) * qJDD(1) + ((-t149 - t212) * t219 + (-t143 + t275) * t222) * qJD(1), -t136 * t321 + t249 - 0.2e1 * t348, t189 + 0.2e1 * t210 + 0.2e1 * t211 + (qJD(1) * t136 - g(3)) * t219 + (qJD(1) * t115 - t264) * t222, -t258 * qJD(1) * pkin(7) - t105 * pkin(2) - g(1) * t267 - g(2) * t268 - g(3) * t324 - t95 * qJ(3) - t149 * qJD(3) - t115 * t136, -t133 * t254 + t356 (-t73 - t346) * t221 + (-t72 + t347) * t218 (-t133 * t222 - t174 * t343) * qJD(1) + t387 (t131 * t222 - t174 * t338) * qJD(1) + t247, -t174 * t321, -t58 * t321 + qJ(3) * t73 - t122 * t174 + t308 * t131 + t379 * t221 + ((t106 - t311) * t174 + t233) * t218, qJ(3) * t72 + t353 * t174 + t59 * t321 + t308 * t133 - t379 * t218 + (-t174 * t311 + t233) * t221, -t259 * t361 + t29 * t78 - t376 + t382, t4 * t85 - t3 * t84 + t43 * t329 - g(1) * (t267 + t385) - g(2) * (t268 + t386) - g(3) * t252 + t378 * t83 + (t63 - t29) * t19 + t361 * t18, -t128 * t84 + t16 * t321 - t360 * t174 - t197 * t390 + t257 * t5 - t355 * t31 + t362 * t78 + t39 * t75, t24 * t78 + t259 * t360 - t377 + t382, t128 * t85 - t17 * t321 + t359 * t174 + t198 * t390 + t256 * t5 - t362 * t259 - t354 * t31 - t40 * t75, t1 * t85 + t5 * t75 + t2 * t84 - g(1) * (t171 + t385) - g(2) * (t169 + t386) - g(3) * (t219 * t262 + t252) + t362 * t31 + t359 * t17 + t360 * t16 + t264 * (-t222 * t262 + t371); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t305, qJDD(2) + t296, -t213 * t226 - t225, qJD(2) * t149 + t172 + t249 - t348, 0, 0, 0, 0, 0, -t174 * t254 + t109 - t320, -t375 * t221 - t319 - t344, t232, -qJD(2) * t83 + t376, -qJD(2) * t78 - t128 * t256 + t174 * t354, t232, qJD(2) * t259 + t128 * t257 - t174 * t355, -qJD(2) * t31 + t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 * t131, -t131 ^ 2 + t133 ^ 2, t72 + t347, t346 - t73, t128, -g(1) * t116 - g(2) * t118 - t114 * t133 + t174 * t59 + t231 + t301, g(1) * t117 - g(2) * t119 + t114 * t131 + t174 * t58 + (t350 - t208) * t218 + t303, t19 * t259 - t364 + (-t215 * t39 - t216 * t40) * pkin(4) + (-t18 + t22) * t78, -t19 * t22 + t18 * t21 - g(1) * t153 - g(2) * t325 + (g(1) * t342 - t133 * t83 + t215 * t4 + t216 * t3 + t301) * pkin(4), t174 * t21 - t38 * t78 + (pkin(5) - t186) * t128 + t380, t17 * t259 - t182 * t39 + t186 * t40 - t364 + (t16 - t327) * t78, t197 * t208 - g(1) * t98 + g(2) * t100 + t128 * t182 - t31 * t78 + t38 * t259 + (0.2e1 * qJD(6) - t22) * t174 + t300, t1 * t182 + t2 * t186 - t31 * t38 - t16 * t21 - g(1) * (-pkin(5) * t97 + qJ(6) * t98 + t153 - t302) - g(2) * (pkin(5) * t99 - qJ(6) * t100 + t325) - g(3) * (-t176 + (-pkin(5) * t198 - qJ(6) * t197) * t222) + t327 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, t18 * t259 + t19 * t78 + t227, t174 * t259 + t39, t263, -t40 + t391, -t16 * t259 + t17 * t78 + t227 + t381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259 * t78 - t128, t40 + t391, -t375 - t388, -t17 * t174 - t370 - t380;];
tau_reg  = t10;
