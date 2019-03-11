% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRP11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:49:03
% EndTime: 2019-03-09 17:49:22
% DurationCPUTime: 6.83s
% Computational Cost: add. (6825->487), mult. (17641->646), div. (0->0), fcn. (13231->8), ass. (0->232)
t200 = sin(qJ(3));
t198 = sin(pkin(6));
t204 = cos(qJ(2));
t286 = qJD(1) * t204;
t266 = t198 * t286;
t243 = t200 * t266;
t283 = qJD(3) * t200;
t336 = t243 - t283;
t178 = -qJD(3) + t266;
t201 = sin(qJ(2));
t304 = cos(pkin(6));
t255 = t304 * qJD(1);
t240 = pkin(1) * t255;
t186 = t201 * t240;
t147 = pkin(8) * t266 + t186;
t348 = -pkin(3) * t336 - t200 * qJD(4) - t147;
t203 = cos(qJ(3));
t347 = -t348 + t178 * (pkin(10) * t200 - qJ(4) * t203);
t287 = qJD(1) * t198;
t267 = t201 * t287;
t144 = -pkin(8) * t267 + t204 * t240;
t226 = (pkin(2) * t201 - pkin(9) * t204) * t198;
t145 = qJD(1) * t226;
t252 = -t200 * t144 + t203 * t145;
t281 = qJD(3) * t203;
t295 = t203 * t204;
t327 = pkin(4) + pkin(9);
t328 = pkin(3) + pkin(10);
t346 = t327 * t281 - (pkin(4) * t295 - t201 * t328) * t287 + t252;
t337 = t203 * t266 - t281;
t179 = t327 * t200;
t345 = -qJD(5) * t179 + t347;
t230 = t255 + qJD(2);
t106 = pkin(9) * t230 + t147;
t140 = (-pkin(2) * t204 - pkin(9) * t201 - pkin(1)) * t198;
t119 = qJD(1) * t140;
t69 = t200 * t106 - t203 * t119;
t293 = -qJD(4) - t69;
t199 = sin(qJ(5));
t202 = cos(qJ(5));
t113 = t199 * t267 - t202 * t243;
t282 = qJD(3) * t202;
t344 = t200 * t282 + t113;
t218 = t203 * t230;
t125 = t200 * t267 - t218;
t274 = qJD(1) * qJD(2);
t259 = t198 * t274;
t181 = t201 * t259;
t279 = qJD(5) * t202;
t280 = qJD(5) * t199;
t239 = t204 * t259;
t127 = t200 * t230 + t203 * t267;
t276 = t127 * qJD(3);
t338 = t200 * t239 + t276;
t44 = -t125 * t279 - t178 * t280 - t202 * t181 - t199 * t338;
t334 = t127 * pkin(4) - t293;
t39 = t178 * t328 + t334;
t105 = -pkin(2) * t230 - t144;
t209 = -t127 * qJ(4) + t105;
t43 = t125 * t328 + t209;
t12 = t199 * t39 + t202 * t43;
t285 = qJD(2) * t201;
t265 = t198 * t285;
t232 = t328 * t265;
t146 = qJD(2) * t226;
t135 = qJD(1) * t146;
t256 = t204 * t304;
t299 = t198 * t201;
t335 = pkin(1) * t256 - pkin(8) * t299;
t148 = t335 * qJD(2);
t136 = qJD(1) * t148;
t246 = t106 * t281 + t119 * t283 - t203 * t135 + t200 * t136;
t270 = t200 * t299;
t241 = qJD(3) * t270;
t92 = qJD(1) * t241 - qJD(3) * t218 - t203 * t239;
t21 = -t92 * pkin(4) - qJD(1) * t232 + t246;
t137 = pkin(8) * t239 + qJD(2) * t186;
t33 = pkin(3) * t338 + t92 * qJ(4) - t127 * qJD(4) + t137;
t23 = pkin(10) * t338 + t33;
t6 = -qJD(5) * t12 - t199 * t23 + t202 * t21;
t89 = t125 * t199 - t178 * t202;
t1 = -t92 * pkin(5) + t44 * qJ(6) - t89 * qJD(6) + t6;
t87 = -t202 * t125 - t178 * t199;
t10 = -qJ(6) * t87 + t12;
t120 = qJD(5) + t127;
t45 = qJD(5) * t89 + t181 * t199 - t202 * t338;
t5 = t199 * t21 + t202 * t23 + t39 * t279 - t280 * t43;
t2 = -qJ(6) * t45 - qJD(6) * t87 + t5;
t11 = -t199 * t43 + t202 * t39;
t9 = -qJ(6) * t89 + t11;
t7 = pkin(5) * t120 + t9;
t343 = -(t120 * t7 - t2) * t199 + (t10 * t120 + t1) * t202;
t306 = qJ(4) * t337 + t348;
t342 = t346 * t202;
t318 = t120 * t87;
t341 = t44 - t318;
t309 = t89 * t120;
t340 = -t45 + t309;
t258 = -t200 * qJ(4) - pkin(2);
t156 = -t203 * t328 + t258;
t339 = t156 * t280 - t179 * t279 - t346 * t199 + t202 * t347;
t167 = t178 * qJ(4);
t70 = t203 * t106 + t200 * t119;
t53 = -pkin(4) * t125 + t70;
t46 = -t167 + t53;
t331 = t120 * t46 + t328 * t92;
t330 = t89 ^ 2;
t329 = t127 ^ 2;
t206 = qJD(1) ^ 2;
t326 = t7 - t9;
t325 = pkin(5) * t202;
t324 = pkin(9) * t200;
t193 = t203 * pkin(9);
t297 = t200 * t204;
t114 = (t199 * t297 + t201 * t202) * t287;
t254 = qJ(6) * t203 - t156;
t323 = qJ(6) * t114 + t254 * t279 + t342 + (-qJ(6) * t283 + qJD(6) * t203 + t345) * t199 - t337 * pkin(5);
t278 = qJD(5) * t203;
t262 = t199 * t278;
t275 = t202 * qJD(6);
t322 = -t203 * t275 - t339 + (t262 + t344) * qJ(6);
t153 = t200 * t304 + t203 * t299;
t298 = t198 * t204;
t245 = t304 * t201 * pkin(1);
t139 = pkin(8) * t298 + pkin(9) * t304 + t245;
t253 = -t200 * t139 + t140 * t203;
t74 = pkin(3) * t298 - t253;
t54 = pkin(4) * t153 + pkin(10) * t298 + t74;
t152 = -t203 * t304 + t270;
t138 = -pkin(2) * t304 - t335;
t210 = -t153 * qJ(4) + t138;
t62 = t152 * t328 + t210;
t321 = t199 * t54 + t202 * t62;
t302 = t125 * qJ(4);
t66 = t127 * t328 + t302;
t320 = t199 * t53 + t202 * t66;
t61 = t125 * pkin(3) + t209;
t317 = t127 * t61;
t316 = t178 * t87;
t315 = t178 * t89;
t173 = qJ(4) * t181;
t247 = t106 * t283 - t119 * t281 - t200 * t135 - t203 * t136;
t236 = qJD(4) * t178 + t247;
t31 = -t173 + t236;
t19 = -pkin(4) * t338 - t31;
t314 = t19 * t199;
t313 = t19 * t202;
t312 = t199 * t92;
t311 = t202 * t44;
t83 = t202 * t92;
t85 = t92 * t200;
t294 = qJ(6) + t328;
t303 = qJ(6) * t127;
t48 = t202 * t53;
t308 = t280 * t294 - t275 + t125 * pkin(5) - t48 - (-t66 - t303) * t199;
t170 = t294 * t202;
t307 = -qJD(5) * t170 - t199 * qJD(6) - t202 * t303 - t320;
t290 = t203 * t144 + t200 * t145;
t76 = -qJ(4) * t267 - t290;
t214 = pkin(4) * t243 + t76;
t305 = -t327 * t283 + t214;
t301 = t127 * t125;
t195 = t198 ^ 2;
t300 = t195 * t206;
t296 = t202 * t203;
t291 = t203 * t139 + t200 * t140;
t289 = t202 * t156 + t199 * t179;
t261 = qJD(2) * t298;
t149 = pkin(8) * t261 + qJD(2) * t245;
t180 = t203 * pkin(4) + t193;
t288 = t201 ^ 2 - t204 ^ 2;
t284 = qJD(2) * t203;
t277 = qJD(5) * t328;
t273 = t178 * t324;
t272 = pkin(9) * t284;
t271 = t201 * t300;
t269 = t202 * t298;
t263 = t178 * t281;
t260 = t195 * t274;
t257 = -t199 * t62 + t202 * t54;
t251 = t120 * t199;
t248 = 0.2e1 * t260;
t238 = -0.2e1 * pkin(1) * t260;
t237 = t199 * t283 - t114;
t233 = pkin(3) * t181;
t34 = -t233 + t246;
t234 = t34 * t200 - t31 * t203;
t229 = 0.2e1 * t255 + qJD(2);
t73 = qJ(4) * t298 - t291;
t228 = -t139 * t281 - t140 * t283 + t203 * t146 - t200 * t148;
t225 = -t120 * t251 - t83;
t99 = t152 * t202 + t199 * t298;
t98 = -t241 + (qJD(3) * t304 + t261) * t203;
t27 = t98 * pkin(4) - t228 - t232;
t216 = -t98 * qJ(4) - t153 * qJD(4) + t149;
t97 = qJD(3) * t153 + t200 * t261;
t32 = t328 * t97 + t216;
t224 = t199 * t27 + t202 * t32 + t54 * t279 - t280 * t62;
t223 = -t178 * t70 - t246;
t222 = t178 * t203;
t220 = -t139 * t283 + t140 * t281 + t200 * t146 + t203 * t148;
t215 = -t120 ^ 2 * t202 + t312;
t63 = -pkin(4) * t152 - t73;
t213 = -qJD(5) * t321 - t199 * t32 + t202 * t27;
t212 = -t125 * t178 - t92;
t211 = pkin(9) * t263 - t181 * t324;
t36 = -qJ(4) * t265 + qJD(4) * t298 - t220;
t28 = -pkin(4) * t97 - t36;
t8 = t45 * pkin(5) + t19;
t172 = -pkin(3) * t203 + t258;
t169 = t294 * t199;
t160 = t202 * t179;
t100 = -t152 * t199 + t269;
t91 = -qJ(6) * t296 + t289;
t86 = t87 ^ 2;
t80 = t200 * pkin(5) + t199 * t254 + t160;
t79 = pkin(3) * t127 + t302;
t78 = -pkin(3) * t267 - t252;
t75 = t92 * t153;
t72 = t152 * pkin(3) + t210;
t67 = t167 - t70;
t65 = pkin(3) * t178 - t293;
t59 = qJD(5) * t99 + t97 * t199 + t202 * t265;
t58 = -qJD(5) * t269 - t97 * t202 + (qJD(5) * t152 + t265) * t199;
t42 = pkin(3) * t97 + t216;
t40 = -pkin(3) * t265 - t228;
t29 = pkin(5) * t87 + qJD(6) + t46;
t15 = qJ(6) * t99 + t321;
t14 = pkin(5) * t153 + qJ(6) * t100 + t257;
t4 = -qJ(6) * t58 + qJD(6) * t99 + t224;
t3 = t98 * pkin(5) - t59 * qJ(6) + t100 * qJD(6) + t213;
t13 = [0, 0, 0, t201 * t204 * t248, -t288 * t248, t229 * t261, -t229 * t265, 0, -t137 * t304 - t149 * t230 + t201 * t238, -t136 * t304 - t148 * t230 + t204 * t238, t127 * t98 - t75, -t98 * t125 - t127 * t97 + t92 * t152 - t153 * t338, -t98 * t178 + (t204 * t92 + (qJD(1) * t153 + t127) * t285) * t198, -t125 * t265 - t152 * t181 + t97 * t178 + t298 * t338 (-t178 * t198 - t195 * t286) * t285, t105 * t97 + t149 * t125 + t137 * t152 + t138 * t338 - t178 * t228 + t181 * t253 + t246 * t298 - t265 * t69, t220 * t178 + t149 * t127 - t138 * t92 + t137 * t153 + t105 * t98 + (-t247 * t204 + (-qJD(1) * t291 - t70) * t285) * t198, t36 * t125 + t40 * t127 + t31 * t152 + t34 * t153 + t338 * t73 + t65 * t98 + t67 * t97 - t74 * t92, -t42 * t125 - t33 * t152 - t40 * t178 + t181 * t74 + t265 * t65 - t298 * t34 - t338 * t72 - t61 * t97, -t42 * t127 - t33 * t153 + t36 * t178 - t61 * t98 + t72 * t92 + (t204 * t31 + (-qJD(1) * t73 - t67) * t285) * t198, t31 * t73 + t33 * t72 + t34 * t74 + t36 * t67 + t40 * t65 + t42 * t61, t100 * t44 + t59 * t89, t100 * t45 - t44 * t99 - t58 * t89 - t59 * t87, t100 * t92 + t120 * t59 - t153 * t44 + t89 * t98, -t120 * t58 - t153 * t45 - t87 * t98 - t92 * t99, t120 * t98 - t75, t11 * t98 + t120 * t213 + t6 * t153 - t19 * t99 - t257 * t92 + t28 * t87 + t63 * t45 + t46 * t58, -t19 * t100 - t12 * t98 - t120 * t224 - t5 * t153 + t28 * t89 + t321 * t92 - t63 * t44 + t46 * t59, t1 * t100 - t10 * t58 + t14 * t44 - t15 * t45 + t2 * t99 - t3 * t89 - t4 * t87 - t59 * t7, t2 * t15 + t10 * t4 + t1 * t14 + t7 * t3 + t8 * (-pkin(5) * t99 + t63) + t29 * (pkin(5) * t58 + t28); 0, 0, 0, -t204 * t271, t288 * t300, -t198 * t206 * t256, t230 * t267 - t181, 0, pkin(1) * t271 + t147 * t230 - t137, pkin(8) * t181 + t144 * t230 + (-qJD(2) * t255 + t300) * t204 * pkin(1), -t127 * t222 - t85, t127 * t243 - t92 * t203 + (-t338 - t276) * t200 + t337 * t125, -t263 + (t178 * t295 + (t200 * qJD(2) - t127) * t201) * t287, t178 * t283 + (-t178 * t297 + (t125 + t284) * t201) * t287, t178 * t267, -pkin(2) * t338 - t105 * t336 - t147 * t125 - t137 * t203 + t178 * t252 + t267 * t69 + t211, pkin(2) * t92 + t137 * t200 - t290 * t178 - t147 * t127 + (t105 * t203 - t273) * qJD(3) + (-t105 * t295 + (t70 - t272) * t201) * t287, -t76 * t125 - t78 * t127 + t234 - t336 * t67 - t337 * t65 + (-t338 + t276) * t193 + (t125 * t283 - t85) * pkin(9), -t306 * t125 - t172 * t338 + t78 * t178 + t33 * t203 - t267 * t65 + t336 * t61 - t211, t172 * t92 - t76 * t178 - t33 * t200 - t306 * t127 + (-t203 * t61 + t273) * qJD(3) + (t61 * t295 + (t67 + t272) * t201) * t287, t33 * t172 - t65 * t78 - t67 * t76 + t306 * t61 + ((t200 * t67 + t203 * t65) * qJD(3) + t234) * pkin(9), t44 * t199 * t203 + (-t202 * t278 + t237) * t89, t89 * t113 + t114 * t87 + (-t199 * t87 + t202 * t89) * t283 + (t199 * t45 + t311 + (t199 * t89 + t202 * t87) * qJD(5)) * t203, -t44 * t200 + t237 * t120 + (-t120 * t279 + t312 - t315) * t203, -t45 * t200 + t344 * t120 + (t120 * t280 + t316 + t83) * t203, -t120 * t222 - t85 -(-t156 * t199 + t160) * t92 + t180 * t45 - t46 * t113 + t305 * t87 + (-t282 * t46 + t6) * t200 + (-t156 * t279 + t199 * t345 + t342) * t120 + (-t11 * t178 - t280 * t46 + t313) * t203, t289 * t92 - t180 * t44 - t46 * t114 + t305 * t89 + (qJD(3) * t199 * t46 - t5) * t200 + t339 * t120 + (t12 * t178 - t279 * t46 - t314) * t203, t10 * t113 + t7 * t114 + t80 * t44 - t91 * t45 - t323 * t89 - t322 * t87 + (t10 * t202 - t199 * t7) * t283 + (t1 * t199 - t2 * t202 + (t10 * t199 + t202 * t7) * qJD(5)) * t203, t2 * t91 + t1 * t80 + t8 * (pkin(5) * t296 + t180) + t323 * t7 + t322 * t10 + ((-t113 - t262) * pkin(5) + t214 + (-t325 - t327) * t283) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301, -t125 ^ 2 + t329, t212, -t127 * t178 - t338, t181, -t105 * t127 + t223, t105 * t125 + t178 * t69 + t247, pkin(3) * t92 - qJ(4) * t338 + (-t67 - t70) * t127 + (t65 + t293) * t125, t125 * t79 - t223 - 0.2e1 * t233 + t317, -t61 * t125 + t79 * t127 + t178 * t293 + 0.2e1 * t173 - t236, -pkin(3) * t34 - qJ(4) * t31 + t293 * t67 - t61 * t79 - t65 * t70, -t251 * t89 - t311 (-t45 - t309) * t202 + (t44 + t318) * t199, t89 * t125 + t225, -t87 * t125 + t215, t120 * t125, qJ(4) * t45 + t11 * t125 + t314 + t334 * t87 + (-t48 + (t66 + t277) * t199) * t120 + t331 * t202, -qJ(4) * t44 - t12 * t125 + t313 + t334 * t89 + (t202 * t277 + t320) * t120 - t331 * t199, t169 * t45 - t170 * t44 - t307 * t87 - t308 * t89 - t343, -t2 * t169 - t1 * t170 + t8 * (pkin(5) * t199 + qJ(4)) + t308 * t7 + (t120 * t325 + t334) * t29 + t307 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t181 - t301, -t178 ^ 2 - t329, -t178 * t67 + t317 + t34, 0, 0, 0, 0, 0, t225 + t316, t215 + t315, t340 * t199 + t202 * t341, t29 * t178 + t343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t87, -t86 + t330, -t341, t340, -t92, t12 * t120 - t46 * t89 + t6, t11 * t120 + t46 * t87 - t5, pkin(5) * t44 - t326 * t87, t326 * t10 + (-t29 * t89 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86 - t330, t10 * t87 + t7 * t89 + t8;];
tauc_reg  = t13;
