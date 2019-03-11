% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPR8
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
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:53:20
% EndTime: 2019-03-09 10:53:32
% DurationCPUTime: 4.37s
% Computational Cost: add. (5526->443), mult. (13818->592), div. (0->0), fcn. (10115->8), ass. (0->227)
t214 = sin(qJ(6));
t217 = cos(qJ(6));
t212 = sin(pkin(10));
t216 = sin(qJ(2));
t299 = qJD(1) * t216;
t276 = t212 * t299;
t213 = cos(pkin(10));
t294 = qJD(2) * t213;
t164 = -t276 + t294;
t279 = t213 * t299;
t295 = qJD(2) * t212;
t165 = t279 + t295;
t215 = sin(qJ(4));
t328 = cos(qJ(4));
t110 = -t328 * t164 + t165 * t215;
t238 = -t215 * t164 - t328 * t165;
t46 = t110 * t214 - t217 * t238;
t280 = t328 * t213;
t218 = cos(qJ(2));
t292 = qJD(2) * t218;
t245 = t280 * t292;
t287 = qJD(1) * qJD(2);
t274 = t218 * t287;
t262 = t212 * t274;
t275 = qJD(4) * t328;
t62 = t215 * (qJD(4) * t165 + t262) - qJD(1) * t245 - t164 * t275;
t170 = t328 * t212 + t215 * t213;
t228 = t218 * t170;
t226 = qJD(2) * t228;
t332 = qJD(4) * t238;
t63 = qJD(1) * t226 - t332;
t14 = qJD(6) * t46 - t214 * t62 - t217 * t63;
t298 = qJD(1) * t218;
t198 = -qJD(4) + t298;
t288 = -qJD(6) - t198;
t345 = t288 * t46;
t348 = t14 + t345;
t315 = t110 * t198;
t347 = t62 - t315;
t237 = -t215 * t212 + t280;
t291 = qJD(4) * t215;
t333 = -t212 * t291 + t213 * t275;
t302 = -t237 * t298 + t333;
t346 = t110 ^ 2;
t330 = t238 ^ 2;
t253 = pkin(2) * t216 - qJ(3) * t218;
t172 = t253 * qJD(1);
t155 = t212 * t172;
t311 = t213 * t216;
t312 = t212 * t218;
t234 = -pkin(7) * t311 - pkin(8) * t312;
t117 = qJD(1) * t234 + t155;
t327 = pkin(8) + qJ(3);
t180 = t327 * t212;
t181 = t327 * t213;
t130 = pkin(7) * t276 + t213 * t172;
t310 = t213 * t218;
t243 = pkin(3) * t216 - pkin(8) * t310;
t98 = qJD(1) * t243 + t130;
t344 = qJD(3) * t280 - t328 * t117 - t180 * t275 + (-qJD(3) * t212 - qJD(4) * t181 - t98) * t215;
t343 = t198 * t238;
t248 = -t217 * t110 - t214 * t238;
t342 = t248 * t288;
t154 = t170 * qJD(4);
t301 = -qJD(1) * t228 + t154;
t341 = -t248 ^ 2 + t46 ^ 2;
t205 = pkin(7) * t299;
t320 = qJD(2) * pkin(2);
t269 = -qJD(3) + t320;
t257 = -t205 + t269;
t233 = pkin(3) * t164 + t257;
t225 = -qJ(5) * t238 + t233;
t329 = pkin(4) + pkin(5);
t26 = -t329 * t110 + t225;
t184 = t198 * qJD(5);
t202 = t216 * t287;
t194 = qJ(5) * t202;
t230 = t243 * qJD(2);
t150 = qJD(2) * t253 - qJD(3) * t216;
t140 = t150 * qJD(1);
t175 = (qJD(3) - t205) * qJD(2);
t96 = t213 * t140 - t175 * t212;
t66 = qJD(1) * t230 + t96;
t177 = -pkin(2) * t218 - qJ(3) * t216 - pkin(1);
t157 = t177 * qJD(1);
t206 = pkin(7) * t298;
t183 = qJD(2) * qJ(3) + t206;
t120 = t213 * t157 - t183 * t212;
t284 = pkin(3) * t298;
t71 = -pkin(8) * t165 + t120 - t284;
t97 = t212 * t140 + t213 * t175;
t78 = -pkin(8) * t262 + t97;
t121 = t212 * t157 + t213 * t183;
t81 = pkin(8) * t164 + t121;
t239 = -t215 * t66 - t71 * t275 + t81 * t291 - t328 * t78;
t8 = -t184 + t194 - t239;
t6 = pkin(9) * t63 + t8;
t282 = t216 * t329;
t264 = qJD(2) * t282;
t272 = t215 * t78 + t81 * t275 + t71 * t291 - t328 * t66;
t7 = pkin(9) * t62 - qJD(1) * t264 + t272;
t281 = -t214 * t6 + t217 * t7;
t340 = -t26 * t46 + t281;
t338 = -0.2e1 * t287;
t326 = qJ(5) * t299 - t344;
t337 = t26 * t248;
t336 = t46 * t248;
t128 = -t215 * t180 + t328 * t181;
t335 = t170 * qJD(3) + t128 * qJD(4) - t215 * t117 + t328 * t98;
t32 = -t215 * t81 + t328 * t71;
t305 = qJD(5) - t32;
t158 = t212 * t284 + t206;
t334 = qJ(5) * t302 + qJD(5) * t170 + t158;
t289 = qJD(6) * t217;
t290 = qJD(6) * t214;
t270 = -t110 * t289 - t214 * t63 + t217 * t62 - t238 * t290;
t331 = t270 + t342;
t325 = pkin(4) * t299 + t335;
t246 = -t170 * t214 - t217 * t237;
t324 = qJD(6) * t246 + t214 * t301 + t217 * t302;
t116 = t170 * t217 - t214 * t237;
t323 = qJD(6) * t116 + t214 * t302 - t217 * t301;
t322 = -t301 * t329 + t334;
t321 = -pkin(4) * t301 + t334;
t33 = t215 * t71 + t328 * t81;
t37 = pkin(4) * t110 - t225;
t319 = t238 * t37;
t185 = t198 * qJ(5);
t23 = pkin(9) * t110 + t33;
t21 = -t185 + t23;
t318 = t214 * t21;
t316 = qJ(5) * t110;
t314 = t238 * t110;
t313 = t212 * t216;
t221 = qJD(1) ^ 2;
t309 = t218 * t221;
t220 = qJD(2) ^ 2;
t308 = t220 * t216;
t307 = t220 * t218;
t306 = -pkin(9) * t238 - t305;
t163 = t213 * t177;
t119 = -pkin(8) * t311 + t163 + (-pkin(7) * t212 - pkin(3)) * t218;
t196 = pkin(7) * t310;
t136 = t212 * t177 + t196;
t126 = -pkin(8) * t313 + t136;
t304 = t215 * t119 + t328 * t126;
t293 = qJD(2) * t216;
t283 = pkin(7) * t293;
t124 = t213 * t150 + t212 * t283;
t197 = pkin(7) * t274;
t149 = pkin(3) * t262 + t197;
t207 = pkin(7) * t292;
t278 = t212 * t292;
t159 = pkin(3) * t278 + t207;
t173 = pkin(3) * t313 + t216 * pkin(7);
t300 = t216 ^ 2 - t218 ^ 2;
t127 = t328 * t180 + t215 * t181;
t297 = qJD(2) * t127;
t296 = qJD(2) * t128;
t20 = t329 * t198 - t306;
t286 = t20 * t289 + t214 * t7 + t217 * t6;
t285 = pkin(7) * t312;
t201 = -pkin(3) * t213 - pkin(2);
t271 = pkin(1) * t338;
t268 = t288 ^ 2;
t267 = -t164 + t294;
t266 = -t165 + t295;
t265 = pkin(4) * t202;
t91 = -pkin(9) * t237 + t128;
t261 = pkin(9) * t302 - qJD(1) * t282 + qJD(6) * t91 - t335;
t90 = -t170 * pkin(9) + t127;
t260 = -pkin(9) * t301 - qJD(6) * t90 + t326;
t258 = t257 - t269;
t49 = -qJ(5) * t218 + t304;
t147 = t237 * t216;
t255 = qJ(5) * t147 - t173;
t254 = t328 * t119 - t215 * t126;
t2 = t214 * t20 + t217 * t21;
t50 = t218 * pkin(4) - t254;
t36 = t218 * pkin(5) - t147 * pkin(9) + t50;
t146 = t170 * t216;
t38 = pkin(9) * t146 + t49;
t252 = -t214 * t38 + t217 * t36;
t251 = t214 * t36 + t217 * t38;
t250 = qJ(5) * t217 - t214 * t329;
t249 = qJ(5) * t214 + t217 * t329;
t247 = t217 * t146 - t147 * t214;
t87 = t146 * t214 + t147 * t217;
t242 = qJ(5) * t170 - t201;
t241 = -t21 * t290 + t286;
t240 = -t198 * t33 - t272;
t141 = t212 * t150;
t101 = qJD(2) * t234 + t141;
t88 = t230 + t124;
t236 = t215 * t101 + t119 * t291 + t126 * t275 - t328 * t88;
t235 = t328 * t101 + t119 * t275 - t126 * t291 + t215 * t88;
t232 = -qJ(5) * t62 - qJD(5) * t238 - t149;
t92 = t154 * t216 + t215 * t278 - t245;
t231 = -qJ(5) * t92 + qJD(5) * t147 - t159;
t10 = -t265 + t272;
t227 = -t198 * t32 + t239;
t15 = pkin(4) * t63 - t232;
t16 = qJ(5) * t293 - qJD(5) * t218 + t235;
t222 = t63 + t343;
t135 = t163 - t285;
t131 = -pkin(7) * t279 + t155;
t125 = -t213 * t283 + t141;
t108 = -pkin(4) * t237 - t242;
t93 = t333 * t216 + t226;
t77 = t237 * t329 + t242;
t72 = pkin(4) * t146 - t255;
t52 = -pkin(4) * t238 + t316;
t51 = -t329 * t146 + t255;
t35 = -t62 - t315;
t34 = t238 * t329 - t316;
t29 = -t185 + t33;
t28 = pkin(4) * t198 + t305;
t27 = pkin(4) * t93 - t231;
t25 = qJD(6) * t87 - t214 * t92 - t217 * t93;
t24 = qJD(6) * t247 + t214 * t93 - t217 * t92;
t19 = -t329 * t93 + t231;
t18 = -pkin(4) * t293 + t236;
t12 = pkin(9) * t93 + t16;
t11 = t92 * pkin(9) + t236 - t264;
t9 = -t329 * t63 + t232;
t1 = t20 * t217 - t318;
t3 = [0, 0, 0, 0.2e1 * t218 * t202, t300 * t338, t307, -t308, 0, -pkin(7) * t307 + t216 * t271, pkin(7) * t308 + t218 * t271 (-qJD(1) * t124 - t96) * t218 + ((-pkin(7) * t164 - t212 * t257) * t218 + (t120 + (t135 + 0.2e1 * t285) * qJD(1)) * t216) * qJD(2) (qJD(1) * t125 + t97) * t218 + ((pkin(7) * t165 - t213 * t257) * t218 + (-t121 + (-t136 + 0.2e1 * t196) * qJD(1)) * t216) * qJD(2), -t124 * t165 + t125 * t164 + (-t212 * t97 - t213 * t96) * t216 + (-t120 * t213 - t121 * t212 + (-t135 * t213 - t136 * t212) * qJD(1)) * t292, t120 * t124 + t121 * t125 + t135 * t96 + t136 * t97 + (-t257 + t205) * t207, -t147 * t62 + t238 * t92, t110 * t92 + t146 * t62 - t147 * t63 + t238 * t93, t198 * t92 + t218 * t62 + (qJD(1) * t147 - t238) * t293, t198 * t93 + t218 * t63 + (-qJD(1) * t146 - t110) * t293 (-t198 - t298) * t293, t236 * t198 + t272 * t218 + t159 * t110 + t173 * t63 + t149 * t146 - t233 * t93 + (qJD(1) * t254 + t32) * t293, t235 * t198 - t239 * t218 - t159 * t238 - t173 * t62 + t149 * t147 + t233 * t92 + (-t304 * qJD(1) - t33) * t293, t10 * t218 + t110 * t27 + t146 * t15 + t18 * t198 + t37 * t93 + t63 * t72 + (-qJD(1) * t50 - t28) * t293, t10 * t147 - t110 * t16 - t146 * t8 - t18 * t238 - t28 * t92 - t29 * t93 - t49 * t63 - t50 * t62, t238 * t27 - t147 * t15 - t16 * t198 - t218 * t8 + t37 * t92 + t62 * t72 + (qJD(1) * t49 + t29) * t293, t10 * t50 + t15 * t72 + t16 * t29 + t18 * t28 + t27 * t37 + t49 * t8, t24 * t46 - t270 * t87, -t14 * t87 - t24 * t248 - t247 * t270 - t25 * t46, -t270 * t218 - t288 * t24 + (-qJD(1) * t87 - t46) * t293, -t14 * t218 + t288 * t25 + (-qJD(1) * t247 + t248) * t293 (t288 - t298) * t293 -(t11 * t217 - t12 * t214) * t288 + t281 * t218 + t19 * t248 + t51 * t14 - t9 * t247 + t26 * t25 + (-t2 * t218 + t251 * t288) * qJD(6) + (-qJD(1) * t252 - t1) * t293 (qJD(6) * t252 + t11 * t214 + t12 * t217) * t288 - t241 * t218 + t19 * t46 - t51 * t270 + t9 * t87 + t26 * t24 + (qJD(1) * t251 + t2) * t293; 0, 0, 0, -t216 * t309, t300 * t221, 0, 0, 0, t221 * pkin(1) * t216, pkin(1) * t309 ((-qJ(3) * t295 - t120) * t216 + (-pkin(7) * t267 + t212 * t258 + t130) * t218) * qJD(1) ((-qJ(3) * t294 + t121) * t216 + (pkin(7) * t266 + t213 * t258 - t131) * t218) * qJD(1), t130 * t165 - t131 * t164 + (qJD(3) * t164 + t120 * t298 + t97) * t213 + (qJD(3) * t165 + t121 * t298 - t96) * t212, -t120 * t130 - t121 * t131 + (-t120 * t212 + t121 * t213) * qJD(3) + (-t212 * t96 + t213 * t97) * qJ(3) + (t257 - t320) * t206, -t170 * t62 - t238 * t302, -t302 * t110 - t170 * t63 - t237 * t62 + t238 * t301, -t302 * t198 + (qJD(2) * t170 + t238) * t299, t301 * t198 + (qJD(2) * t237 + t110) * t299, t198 * t299, -t158 * t110 - t149 * t237 + t201 * t63 + t335 * t198 - t301 * t233 + (-t32 - t297) * t299, t158 * t238 + t149 * t170 - t201 * t62 + t344 * t198 - t302 * t233 + (t33 - t296) * t299, t108 * t63 - t15 * t237 + t301 * t37 + t325 * t198 - t321 * t110 + (t28 - t297) * t299, t10 * t170 + t326 * t110 - t127 * t62 - t128 * t63 + t237 * t8 - t238 * t325 + t302 * t28 - t301 * t29, t108 * t62 - t15 * t170 - t302 * t37 + t326 * t198 - t321 * t238 + (-t29 + t296) * t299, t10 * t127 + t108 * t15 + t128 * t8 + t325 * t28 - t326 * t29 - t321 * t37, -t116 * t270 + t324 * t46, -t116 * t14 - t246 * t270 - t248 * t324 - t323 * t46, -t324 * t288 + (-qJD(2) * t116 + t46) * t299, t323 * t288 + (-qJD(2) * t246 - t248) * t299, -t288 * t299, -t9 * t246 + t77 * t14 + t322 * t248 + t323 * t26 - (t214 * t260 - t217 * t261) * t288 + (-(-t214 * t91 + t217 * t90) * qJD(2) + t1) * t299, t9 * t116 - t77 * t270 + t322 * t46 + t324 * t26 - (t214 * t261 + t217 * t260) * t288 + ((t214 * t90 + t217 * t91) * qJD(2) - t2) * t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266 * t298, t267 * t298, -t164 ^ 2 - t165 ^ 2, t120 * t165 - t121 * t164 + t197, 0, 0, 0, 0, 0, t222, -t347, t222, -t330 - t346, t347, t110 * t29 + t238 * t28 + t15, 0, 0, 0, 0, 0, -t14 + t345, t270 - t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t314, t330 - t346, t35, -t170 * t274 + t332 + t343, t202, -t233 * t238 + t240, -t110 * t233 + t227, -t110 * t52 + t240 + 0.2e1 * t265 + t319, pkin(4) * t62 - qJ(5) * t63 - (t29 - t33) * t238 + (t28 - t305) * t110, -t110 * t37 - t238 * t52 - 0.2e1 * t184 + 0.2e1 * t194 - t227, -pkin(4) * t10 + qJ(5) * t8 - t28 * t33 + t305 * t29 - t37 * t52, -t336, -t341, t331, t348, t202, t249 * t202 - t34 * t248 - (t214 * t306 - t217 * t23) * t288 + (t250 * t288 + t2) * qJD(6) - t340, t250 * t202 - t34 * t46 - t337 - (t214 * t23 + t217 * t306) * t288 + (-t249 * t288 - t318) * qJD(6) + t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202 - t314, t35, -t198 ^ 2 - t330, t198 * t29 + t10 - t319, 0, 0, 0, 0, 0, -t202 * t217 - t214 * t268 + t238 * t248, t202 * t214 - t217 * t268 + t238 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t336, t341, -t331, -t348, -t202 (-qJD(6) - t288) * t2 + t340, -t1 * t288 - t241 + t337;];
tauc_reg  = t3;
