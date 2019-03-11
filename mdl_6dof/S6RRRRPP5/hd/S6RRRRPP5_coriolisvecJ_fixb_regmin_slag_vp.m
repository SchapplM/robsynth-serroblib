% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:09:13
% EndTime: 2019-03-09 21:09:25
% DurationCPUTime: 5.16s
% Computational Cost: add. (7652->505), mult. (18448->613), div. (0->0), fcn. (12498->6), ass. (0->239)
t215 = cos(qJ(3));
t214 = sin(qJ(2));
t216 = cos(qJ(2));
t309 = t215 * t216;
t241 = pkin(3) * t214 - pkin(9) * t309;
t340 = -pkin(9) - pkin(8);
t273 = qJD(3) * t340;
t246 = pkin(2) * t214 - pkin(8) * t216;
t155 = t246 * qJD(1);
t213 = sin(qJ(3));
t293 = qJD(1) * t214;
t265 = t213 * t293;
t296 = pkin(7) * t265 + t215 * t155;
t365 = -qJD(1) * t241 + t215 * t273 - t296;
t135 = t213 * t155;
t310 = t214 * t215;
t364 = -t135 - (-pkin(9) * t213 * t216 - pkin(7) * t310) * qJD(1) + t213 * t273;
t292 = qJD(1) * t216;
t193 = -qJD(3) + t292;
t177 = -qJD(4) + t193;
t287 = qJD(2) * t215;
t148 = -t265 + t287;
t289 = qJD(2) * t213;
t149 = t215 * t293 + t289;
t212 = sin(qJ(4));
t339 = cos(qJ(4));
t99 = -t339 * t148 + t149 * t212;
t324 = t177 * t99;
t284 = qJD(3) * t214;
t269 = t213 * t284;
t286 = qJD(2) * t216;
t229 = t215 * t286 - t269;
t279 = qJD(2) * qJD(3);
t222 = qJD(1) * t229 + t215 * t279;
t263 = t339 * qJD(4);
t261 = qJD(1) * t284;
t280 = qJD(1) * qJD(2);
t262 = t216 * t280;
t274 = t215 * t261 + (t262 + t279) * t213;
t282 = qJD(4) * t212;
t46 = -t148 * t263 + t149 * t282 + t212 * t274 - t339 * t222;
t225 = t46 + t324;
t342 = t99 ^ 2;
t235 = t212 * t148 + t339 * t149;
t97 = t235 ^ 2;
t363 = t97 - t342;
t271 = t213 * t292;
t272 = t339 * t215;
t313 = t212 * t213;
t343 = qJD(3) + qJD(4);
t344 = t339 * qJD(3) + t263;
t302 = -t212 * t271 - t215 * t344 + t272 * t292 + t313 * t343;
t327 = qJD(2) * pkin(2);
t169 = pkin(7) * t293 - t327;
t239 = t148 * pkin(3) - t169;
t227 = qJ(5) * t235 + t239;
t38 = pkin(4) * t99 - t227;
t362 = t38 * t99;
t361 = qJ(6) * t99;
t325 = t235 * t99;
t360 = t239 * t99;
t320 = t99 * qJ(5);
t171 = t340 * t213;
t172 = t340 * t215;
t359 = t171 * t263 + t172 * t282 + t212 * t365 + t364 * t339;
t151 = t212 * t215 + t339 * t213;
t112 = t343 * t151;
t301 = -t151 * t292 + t112;
t206 = pkin(7) * t292;
t285 = qJD(3) * t213;
t247 = -t206 + (-t271 + t285) * pkin(3);
t198 = t214 * t280;
t358 = qJ(5) * t198 - t177 * qJD(5);
t316 = t235 * t177;
t47 = qJD(4) * t235 + t212 * t222 + t339 * t274;
t357 = -t47 - t316;
t270 = t213 * t286;
t283 = qJD(3) * t215;
t356 = t214 * t283 + t270;
t341 = pkin(4) + pkin(5);
t20 = -t341 * t99 + qJD(6) + t227;
t352 = t47 * qJ(6) + t99 * qJD(6);
t355 = t20 * t99 + t352;
t354 = -0.2e1 * t280;
t353 = pkin(4) * t235;
t331 = qJ(5) * t293 - t359;
t326 = t235 * t38;
t351 = -t177 ^ 2 - t97;
t119 = t212 * t171 - t339 * t172;
t350 = qJD(4) * t119 + t364 * t212 - t339 * t365;
t349 = qJ(6) * t235;
t348 = t239 * t235;
t347 = t235 * t341;
t170 = qJD(2) * pkin(8) + t206;
t163 = -pkin(2) * t216 - pkin(8) * t214 - pkin(1);
t140 = t163 * qJD(1);
t312 = t213 * t140;
t108 = t215 * t170 + t312;
t82 = pkin(9) * t148 + t108;
t322 = t212 * t82;
t107 = t215 * t140 - t170 * t213;
t81 = -pkin(9) * t149 + t107;
t72 = -pkin(3) * t193 + t81;
t26 = t339 * t72 - t322;
t303 = qJD(5) - t26;
t346 = 0.2e1 * t358;
t345 = qJ(5) * t302 - t151 * qJD(5) + t247;
t210 = t214 ^ 2;
t281 = t210 * qJD(1);
t191 = pkin(4) * t198;
t231 = t241 * qJD(2);
t158 = t246 * qJD(2);
t141 = qJD(1) * t158;
t252 = pkin(7) * t198;
t299 = t215 * t141 + t213 * t252;
t37 = qJD(1) * t231 - qJD(3) * t82 + t299;
t300 = t140 * t283 + t213 * t141;
t223 = -t170 * t285 - t215 * t252 + t300;
t48 = -pkin(9) * t274 + t223;
t259 = t212 * t48 + t82 * t263 + t72 * t282 - t339 * t37;
t5 = -t191 + t259;
t2 = -pkin(5) * t198 + t46 * qJ(6) - qJD(6) * t235 + t5;
t224 = -t20 * t235 + t2;
t338 = pkin(7) * t213;
t337 = pkin(8) * t193;
t335 = -t301 * t341 - t345;
t276 = t214 * t341;
t334 = -t302 * qJ(6) - qJD(1) * t276 + t151 * qJD(6) - t350;
t150 = -t272 + t313;
t333 = -t301 * qJ(6) - qJD(6) * t150 + t331;
t332 = t301 * pkin(4) + t345;
t330 = pkin(4) * t293 + t350;
t33 = t339 * t81 - t322;
t79 = t339 * t82;
t27 = t212 * t72 + t79;
t328 = qJ(5) * t47;
t197 = pkin(3) * t212 + qJ(5);
t323 = t197 * t47;
t321 = t33 * t177;
t147 = t215 * t163;
t106 = -pkin(9) * t310 + t147 + (-pkin(3) - t338) * t216;
t195 = pkin(7) * t309;
t295 = t213 * t163 + t195;
t311 = t213 * t214;
t113 = -pkin(9) * t311 + t295;
t319 = t212 * t106 + t339 * t113;
t18 = t33 + t349;
t185 = pkin(3) * t263 + qJD(5);
t318 = -t18 + t185;
t317 = t185 - t33;
t315 = t149 * t193;
t314 = t169 * t213;
t308 = t216 * t193;
t219 = qJD(1) ^ 2;
t307 = t216 * t219;
t218 = qJD(2) ^ 2;
t306 = t218 * t214;
t305 = t218 * t216;
t15 = t26 + t349;
t304 = qJD(5) - t15;
t298 = t213 * t158 + t163 * t283;
t288 = qJD(2) * t214;
t297 = t215 * t158 + t288 * t338;
t96 = t274 * pkin(3) + pkin(7) * t262;
t159 = pkin(3) * t311 + t214 * pkin(7);
t294 = -t216 ^ 2 + t210;
t118 = -t339 * t171 - t212 * t172;
t291 = qJD(2) * t118;
t290 = qJD(2) * t119;
t278 = pkin(3) * t282;
t275 = t212 * t311;
t121 = t356 * pkin(3) + pkin(7) * t286;
t204 = -t215 * pkin(3) - pkin(2);
t267 = t193 * t283;
t266 = t177 * t282;
t260 = -t212 * t37 - t72 * t263 + t82 * t282 - t339 * t48;
t258 = -t198 + t325;
t257 = pkin(1) * t354;
t256 = t193 + t292;
t255 = -t148 + t287;
t254 = -t149 + t289;
t253 = qJD(3) + t292;
t203 = -t339 * pkin(3) - pkin(4);
t16 = t27 + t361;
t251 = t339 * t286;
t32 = t212 * t81 + t79;
t17 = t32 + t361;
t250 = -t17 + t278;
t249 = -t32 + t278;
t55 = -qJ(5) * t216 + t319;
t8 = t47 * pkin(4) + t46 * qJ(5) - qJD(5) * t235 + t96;
t129 = t214 * t272 - t275;
t244 = t129 * qJ(5) - t159;
t243 = t339 * t106 - t212 * t113;
t242 = -t149 * pkin(3) - t320;
t4 = -t260 + t358;
t240 = t151 * qJ(5) - t204;
t56 = t216 * pkin(4) - t243;
t238 = -t177 * t26 + t260;
t237 = -t177 * t27 - t259;
t236 = -t32 * t177 - t259;
t59 = t231 + (-t195 + (pkin(9) * t214 - t163) * t213) * qJD(3) + t297;
t62 = -t356 * pkin(9) + (-t214 * t287 - t216 * t285) * pkin(7) + t298;
t234 = t106 * t282 + t113 * t263 + t212 * t62 - t339 * t59;
t233 = t106 * t263 - t113 * t282 + t212 * t59 + t339 * t62;
t232 = t253 * t289;
t3 = -pkin(5) * t47 - t8;
t230 = -t185 * t177 + t197 * t198 + t4;
t63 = t112 * t214 + t212 * t270 - t215 * t251;
t228 = -t63 * qJ(5) + t129 * qJD(5) - t121;
t10 = qJ(5) * t288 - qJD(5) * t216 + t233;
t196 = -pkin(5) + t203;
t162 = t177 * qJ(5);
t144 = pkin(3) * t266;
t128 = t151 * t214;
t95 = pkin(4) * t150 - t240;
t85 = qJ(6) * t150 + t119;
t84 = -t151 * qJ(6) + t118;
t80 = -t341 * t150 + t240;
t73 = pkin(4) * t128 - t244;
t64 = t213 * t251 - t212 * t269 - qJD(4) * t275 + (t212 * t286 + t214 * t344) * t215;
t58 = t320 + t353;
t57 = -t341 * t128 + t244;
t50 = -t242 + t353;
t43 = qJ(6) * t128 + t55;
t36 = t216 * pkin(5) - t129 * qJ(6) + t56;
t31 = -t320 - t347;
t23 = -t162 + t27;
t22 = pkin(4) * t177 + t303;
t21 = t242 - t347;
t14 = t16 - t162;
t13 = t341 * t177 + t304;
t12 = pkin(4) * t64 - t228;
t11 = -pkin(4) * t288 + t234;
t9 = -t341 * t64 + t228;
t7 = qJ(6) * t64 + qJD(6) * t128 + t10;
t6 = t63 * qJ(6) - qJD(2) * t276 - t129 * qJD(6) + t234;
t1 = t4 + t352;
t19 = [0, 0, 0, 0.2e1 * t216 * t198, t294 * t354, t305, -t306, 0, -pkin(7) * t305 + t214 * t257, pkin(7) * t306 + t216 * t257, t149 * t229 + t222 * t310 (t148 * t215 - t149 * t213) * t286 + ((-t148 + t265) * t285 + (-t149 * qJD(3) - t232 - t274) * t215) * t214, t256 * t269 + (t149 * t214 + (t281 + (-t193 - t253) * t216) * t215) * qJD(2), t214 * t267 + t274 * t216 + (t148 * t214 + (-t281 + t308) * t213) * qJD(2), -t256 * t288 -(-t163 * t285 + t297) * t193 + (pkin(7) * t274 + t169 * t283 + (qJD(1) * t147 + t107) * qJD(2)) * t214 + ((-pkin(7) * t148 + t314) * qJD(2) + (t312 + (pkin(7) * t193 + t170) * t215) * qJD(3) - t299) * t216, t298 * t193 + t300 * t216 + (-t169 * t214 - t170 * t216 + (-t308 - t281) * pkin(7)) * t285 + ((pkin(7) * t149 + t169 * t215) * t216 + (-t295 * qJD(1) - t108 + (-t193 + t253) * pkin(7) * t215) * t214) * qJD(2), -t129 * t46 - t235 * t63, t128 * t46 - t129 * t47 - t235 * t64 + t63 * t99, t63 * t177 + t46 * t216 + (qJD(1) * t129 + t235) * t288, t64 * t177 + t47 * t216 + (-qJD(1) * t128 - t99) * t288 (-t177 - t292) * t288, t234 * t177 + t259 * t216 + t121 * t99 + t159 * t47 + t96 * t128 - t239 * t64 + (qJD(1) * t243 + t26) * t288, t233 * t177 - t260 * t216 + t121 * t235 - t159 * t46 + t96 * t129 + t239 * t63 + (-t319 * qJD(1) - t27) * t288, t11 * t177 + t12 * t99 + t8 * t128 + t5 * t216 + t38 * t64 + t73 * t47 + (-qJD(1) * t56 - t22) * t288, -t10 * t99 + t11 * t235 - t128 * t4 + t129 * t5 - t22 * t63 - t23 * t64 - t46 * t56 - t47 * t55, -t10 * t177 - t12 * t235 - t8 * t129 - t4 * t216 + t38 * t63 + t73 * t46 + (qJD(1) * t55 + t23) * t288, t10 * t23 + t11 * t22 + t12 * t38 + t4 * t55 + t5 * t56 + t73 * t8, -t3 * t128 + t6 * t177 + t2 * t216 - t20 * t64 - t57 * t47 - t9 * t99 + (-qJD(1) * t36 - t13) * t288, -t1 * t216 + t9 * t235 + t3 * t129 - t7 * t177 - t20 * t63 - t57 * t46 + (qJD(1) * t43 + t14) * t288, t1 * t128 - t129 * t2 + t13 * t63 + t14 * t64 - t235 * t6 + t36 * t46 + t43 * t47 + t7 * t99, t1 * t43 + t13 * t6 + t14 * t7 + t2 * t36 + t20 * t9 + t3 * t57; 0, 0, 0, -t214 * t307, t294 * t219, 0, 0, 0, t219 * pkin(1) * t214, pkin(1) * t307, -t213 ^ 2 * t261 + (t232 - t315) * t215 (-t274 + t315) * t213 + ((t148 + t287) * qJD(3) + (t216 * t255 - t269) * qJD(1)) * t215, -t267 + (t214 * t254 + t215 * t308) * qJD(1), t193 * t285 + (-t213 * t308 + t214 * t255) * qJD(1), t193 * t293, -pkin(2) * t274 + t296 * t193 + (t215 * t337 + t314) * qJD(3) + ((-pkin(8) * t289 - t107) * t214 + (-pkin(7) * t255 - t314) * t216) * qJD(1), -t135 * t193 + (-t213 * t337 + (t169 - t327) * t215) * qJD(3) + ((-t169 - t327) * t309 + (pkin(2) * t285 - pkin(8) * t287 + t108) * t214 + (t193 * t310 + t216 * t254) * pkin(7)) * qJD(1), -t46 * t151 - t235 * t302, t46 * t150 - t151 * t47 - t235 * t301 + t302 * t99, t302 * t177 + (qJD(2) * t151 - t235) * t293, t301 * t177 + (-qJD(2) * t150 + t99) * t293, t177 * t293, t96 * t150 + t204 * t47 + t247 * t99 + t350 * t177 - t301 * t239 + (-t26 - t291) * t293, t96 * t151 - t204 * t46 + t359 * t177 + t302 * t239 + t247 * t235 + (t27 - t290) * t293, t8 * t150 + t95 * t47 + t332 * t99 + t301 * t38 + t330 * t177 + (t22 - t291) * t293, -t118 * t46 - t119 * t47 - t4 * t150 + t5 * t151 - t302 * t22 - t301 * t23 + t235 * t330 + t331 * t99, -t8 * t151 + t95 * t46 + t302 * t38 + t331 * t177 - t332 * t235 + (-t23 + t290) * t293, t5 * t118 + t4 * t119 + t330 * t22 - t331 * t23 + t332 * t38 + t8 * t95, -t3 * t150 - t80 * t47 - t335 * t99 - t301 * t20 - t334 * t177 + (-qJD(2) * t84 + t13) * t293, t3 * t151 - t80 * t46 - t302 * t20 + t333 * t177 + t335 * t235 + (qJD(2) * t85 - t14) * t293, t1 * t150 + t302 * t13 + t301 * t14 - t2 * t151 + t235 * t334 - t333 * t99 + t84 * t46 + t85 * t47, t1 * t85 - t334 * t13 - t333 * t14 + t2 * t84 + t335 * t20 + t3 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149 * t148, -t148 ^ 2 + t149 ^ 2, t148 * t193 + t222, -t274 - t315, t198, -t169 * t149 + t299 + (-qJD(3) - t193) * t108, -t107 * t193 - t148 * t169 - t223, t325, t363, -t225, t357, t198, t348 + (-t149 * t99 + t339 * t198 + t266) * pkin(3) + t236, -t360 - t321 + (-t149 * t235 + t177 * t263 - t198 * t212) * pkin(3) + t260, -t198 * t203 - t50 * t99 + t144 + t191 + t236 - t326, -t203 * t46 - t323 + (t23 + t249) * t235 + (t22 - t317) * t99, t235 * t50 + t230 + t321 - t362, t4 * t197 + t5 * t203 + t249 * t22 + t317 * t23 - t38 * t50, -t17 * t177 - t196 * t198 + t21 * t99 + t144 - t224, t177 * t18 - t21 * t235 + t230 + t355, t196 * t46 + t323 + (-t14 - t250) * t235 + (-t13 + t318) * t99, t1 * t197 + t13 * t250 + t14 * t318 + t2 * t196 - t20 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t325, t363, -t225, t357, t198, t237 + t348, t238 - t360, -t58 * t99 + 0.2e1 * t191 + t237 - t326, pkin(4) * t46 - t328 + (t23 - t27) * t235 + (t22 - t303) * t99, t235 * t58 - t238 + t346 - t362, -t5 * pkin(4) + t4 * qJ(5) - t22 * t27 + t303 * t23 - t38 * t58, -t16 * t177 + t198 * t341 + t31 * t99 - t224, t15 * t177 - t235 * t31 - t260 + t346 + t355, t328 - t341 * t46 + (-t14 + t16) * t235 + (-t13 + t304) * t99, t1 * qJ(5) - t13 * t16 + t14 * t304 - t2 * t341 - t20 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t258, -t225, t351, t177 * t23 + t326 + t5, t258, t351, t225, t14 * t177 + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 + t316, -t46 + t324, -t97 - t342, t13 * t235 - t14 * t99 + t3;];
tauc_reg  = t19;
