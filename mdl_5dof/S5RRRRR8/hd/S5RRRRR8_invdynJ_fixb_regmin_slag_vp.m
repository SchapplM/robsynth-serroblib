% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:51
% EndTime: 2019-12-31 22:26:04
% DurationCPUTime: 5.33s
% Computational Cost: add. (5364->411), mult. (11920->562), div. (0->0), fcn. (9010->14), ass. (0->236)
t201 = qJD(2) + qJD(3);
t208 = sin(qJ(3));
t213 = cos(qJ(2));
t323 = cos(qJ(3));
t267 = qJD(1) * t323;
t209 = sin(qJ(2));
t288 = qJD(1) * t209;
t341 = -t208 * t288 + t213 * t267;
t345 = t341 * t201;
t205 = qJ(2) + qJ(3);
t198 = cos(t205);
t320 = g(3) * t198;
t196 = sin(t205);
t210 = sin(qJ(1));
t214 = cos(qJ(1));
t249 = g(1) * t214 + g(2) * t210;
t343 = t249 * t196;
t344 = -t343 + t320;
t326 = qJD(4) + qJD(5);
t342 = t341 - t326;
t325 = pkin(6) + pkin(7);
t169 = t325 * t213;
t154 = qJD(1) * t169;
t141 = t208 * t154;
t167 = t325 * t209;
t152 = qJD(1) * t167;
t100 = -t152 * t323 - t141;
t266 = qJD(3) * t323;
t335 = pkin(2) * t266 - t100;
t207 = sin(qJ(4));
t212 = cos(qJ(4));
t292 = t208 * t213;
t140 = -qJD(1) * t292 - t209 * t267;
t93 = -pkin(3) * t140 - pkin(8) * t341;
t82 = pkin(2) * t288 + t93;
t340 = -t207 * t335 - t212 * t82;
t112 = -t140 * t207 - t201 * t212;
t206 = sin(qJ(5));
t211 = cos(qJ(5));
t239 = t140 * t212 - t201 * t207;
t240 = t112 * t206 + t211 * t239;
t58 = t112 * t211 - t206 * t239;
t339 = t240 * t58;
t148 = t206 * t207 - t211 * t212;
t311 = t342 * t148;
t295 = t206 * t212;
t150 = t207 * t211 + t295;
t338 = t342 * t150;
t200 = qJDD(2) + qJDD(3);
t282 = qJD(1) * qJD(2);
t264 = t213 * t282;
t281 = t209 * qJDD(1);
t109 = qJDD(2) * pkin(2) + t325 * (-t264 - t281);
t265 = t209 * t282;
t280 = t213 * qJDD(1);
t111 = t325 * (-t265 + t280);
t315 = qJD(2) * pkin(2);
t143 = -t152 + t315;
t287 = qJD(3) * t208;
t256 = t109 * t323 - t111 * t208 - t143 * t287 - t154 * t266;
t36 = -pkin(3) * t200 - t256;
t337 = t36 + t320;
t333 = t240 ^ 2 - t58 ^ 2;
t134 = qJD(4) - t341;
t127 = qJD(5) + t134;
t283 = qJD(5) * t211;
t284 = qJD(5) * t206;
t285 = qJD(4) * t212;
t286 = qJD(4) * t207;
t260 = qJDD(1) * t323;
t68 = t208 * t280 + t209 * t260 + t345;
t45 = t140 * t286 + t200 * t207 + t201 * t285 + t212 * t68;
t46 = -qJD(4) * t239 - t200 * t212 + t207 * t68;
t13 = -t112 * t283 - t206 * t46 + t211 * t45 + t239 * t284;
t332 = t127 * t58 + t13;
t204 = qJ(4) + qJ(5);
t195 = sin(t204);
t197 = cos(t204);
t297 = t198 * t210;
t120 = t195 * t214 - t197 * t297;
t296 = t198 * t214;
t122 = t195 * t210 + t197 * t296;
t188 = g(3) * t196;
t194 = -pkin(2) * t213 - pkin(1);
t165 = t194 * qJD(1);
t80 = -pkin(3) * t341 + pkin(8) * t140 + t165;
t142 = t323 * t154;
t97 = t143 * t208 + t142;
t84 = pkin(8) * t201 + t97;
t48 = t207 * t80 + t212 * t84;
t29 = -pkin(9) * t112 + t48;
t26 = t29 * t284;
t96 = t143 * t323 - t141;
t83 = -pkin(3) * t201 - t96;
t55 = pkin(4) * t112 + t83;
t331 = g(1) * t122 - g(2) * t120 + t188 * t197 + t55 * t58 + t26;
t119 = t195 * t297 + t197 * t214;
t121 = -t195 * t296 + t197 * t210;
t135 = pkin(2) * t265 + qJDD(1) * t194;
t151 = t209 * t323 + t292;
t105 = t201 * t151;
t243 = t208 * t281 - t213 * t260;
t69 = qJD(1) * t105 + t243;
t27 = pkin(3) * t69 - pkin(8) * t68 + t135;
t25 = t212 * t27;
t220 = t208 * t109 + t111 * t323 + t143 * t266 - t154 * t287;
t35 = pkin(8) * t200 + t220;
t67 = qJDD(4) + t69;
t3 = pkin(4) * t67 - pkin(9) * t45 - qJD(4) * t48 - t207 * t35 + t25;
t234 = t207 * t27 + t212 * t35 + t285 * t80 - t286 * t84;
t5 = -pkin(9) * t46 + t234;
t274 = -t206 * t5 + t211 * t3;
t47 = -t207 * t84 + t212 * t80;
t28 = pkin(9) * t239 + t47;
t23 = pkin(4) * t134 + t28;
t313 = t211 * t29;
t9 = t206 * t23 + t313;
t330 = -g(1) * t121 + g(2) * t119 - qJD(5) * t9 + t188 * t195 + t240 * t55 + t274;
t223 = qJD(5) * t240 - t206 * t45 - t211 * t46;
t329 = -t127 * t240 + t223;
t99 = -t152 * t208 + t142;
t254 = pkin(2) * t287 - t99;
t88 = t150 * t151;
t301 = t341 * t207;
t328 = (t286 - t301) * pkin(4);
t327 = t207 * t82 - t212 * t335;
t324 = -pkin(8) - pkin(9);
t319 = t212 * pkin(4);
t199 = t212 * pkin(9);
t191 = pkin(2) * t208 + pkin(8);
t318 = -pkin(9) - t191;
t317 = t207 * t93 + t212 * t96;
t314 = t207 * t45;
t312 = t83 * t341;
t116 = -t167 * t208 + t169 * t323;
t110 = t212 * t116;
t232 = -t208 * t209 + t213 * t323;
t95 = -pkin(3) * t232 - pkin(8) * t151 + t194;
t309 = t207 * t95 + t110;
t308 = t328 + t254;
t104 = t201 * t232;
t307 = t104 * t207;
t306 = t104 * t212;
t305 = t112 * t134;
t304 = t239 * t134;
t303 = t127 * t140;
t302 = t134 * t140;
t300 = t140 * t341;
t299 = t151 * t207;
t298 = t151 * t212;
t294 = t207 * t210;
t293 = t207 * t214;
t291 = t210 * t212;
t290 = t212 * t214;
t202 = t209 ^ 2;
t289 = -t213 ^ 2 + t202;
t279 = pkin(9) * t301;
t278 = t209 * t315;
t276 = qJD(4) * pkin(8) * t134;
t74 = t83 * t286;
t273 = qJD(2) * t325;
t272 = qJD(4) * t324;
t270 = t151 * t286;
t269 = t151 * t285;
t263 = qJD(5) * t23 + t5;
t261 = qJD(4) * t318;
t259 = -qJD(4) * t80 - t35;
t257 = t134 * t212;
t192 = -pkin(2) * t323 - pkin(3);
t253 = -t97 + t328;
t252 = -pkin(4) * t140 - t199 * t341;
t251 = -t285 * t84 + t25;
t250 = -pkin(8) * t67 - t312;
t248 = g(1) * t210 - g(2) * t214;
t144 = t318 * t207;
t247 = -qJD(5) * t144 - t207 * t261 - t279 + t327;
t145 = t191 * t212 + t199;
t246 = qJD(5) * t145 - t212 * t261 + t252 - t340;
t166 = t324 * t207;
t245 = -qJD(5) * t166 - t207 * t272 - t279 + t317;
t168 = pkin(8) * t212 + t199;
t87 = t212 * t93;
t244 = qJD(5) * t168 - t207 * t96 - t212 * t272 + t252 + t87;
t242 = -t191 * t67 - t312;
t238 = -t48 * t140 + t207 * t337 + t83 * t285;
t237 = t47 * t140 + t212 * t343 + t74;
t235 = -0.2e1 * pkin(1) * t282 - pkin(6) * qJDD(2);
t233 = -t167 * t323 - t169 * t208;
t231 = t269 + t307;
t230 = -t270 + t306;
t54 = pkin(3) * t105 - pkin(8) * t104 + t278;
t153 = t209 * t273;
t155 = t213 * t273;
t62 = qJD(3) * t233 - t153 * t323 - t208 * t155;
t229 = -t116 * t286 + t207 * t54 + t212 * t62 + t285 * t95;
t215 = qJD(2) ^ 2;
t225 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t215 + t248;
t216 = qJD(1) ^ 2;
t224 = pkin(1) * t216 - pkin(6) * qJDD(1) + t249;
t222 = t165 * t140 + t256 - t344;
t17 = pkin(4) * t46 + t36;
t8 = -t206 * t29 + t211 * t23;
t221 = t8 * t140 + t17 * t148 - t197 * t344 - t338 * t55;
t63 = qJD(3) * t116 - t208 * t153 + t155 * t323;
t219 = -t9 * t140 + t17 * t150 + t195 * t344 + t311 * t55;
t217 = g(1) * t296 + g(2) * t297 - t165 * t341 + t188 - t220;
t193 = -pkin(3) - t319;
t164 = t192 - t319;
t132 = t198 * t290 + t294;
t131 = -t198 * t293 + t291;
t130 = -t198 * t291 + t293;
t129 = t198 * t294 + t290;
t91 = t212 * t95;
t89 = t148 * t151;
t81 = pkin(4) * t299 - t233;
t71 = t140 ^ 2 - t341 ^ 2;
t65 = qJDD(5) + t67;
t53 = -t243 + (-qJD(1) * t151 - t140) * t201;
t52 = t68 - t345;
t51 = t212 * t54;
t49 = -pkin(9) * t299 + t309;
t38 = -pkin(4) * t232 - pkin(9) * t298 - t116 * t207 + t91;
t32 = pkin(4) * t231 + t63;
t22 = t104 * t295 - t206 * t270 - t284 * t299 + (t298 * t326 + t307) * t211;
t21 = -t148 * t104 - t326 * t88;
t20 = t134 * t257 - t140 * t239 + t207 * t67;
t19 = -t134 ^ 2 * t207 - t112 * t140 + t212 * t67;
t18 = -t239 * t257 + t314;
t12 = t127 * t338 - t140 * t58 - t148 * t65;
t11 = t127 * t311 - t140 * t240 + t150 * t65;
t10 = -pkin(9) * t231 + t229;
t7 = -pkin(9) * t306 + pkin(4) * t105 - t207 * t62 + t51 + (-t110 + (pkin(9) * t151 - t95) * t207) * qJD(4);
t6 = (t45 - t305) * t212 + (-t46 + t304) * t207;
t4 = t13 * t150 - t240 * t311;
t1 = -t13 * t148 + t150 * t223 - t240 * t338 - t311 * t58;
t2 = [qJDD(1), t248, t249, qJDD(1) * t202 + 0.2e1 * t209 * t264, 0.2e1 * t209 * t280 - 0.2e1 * t282 * t289, qJDD(2) * t209 + t213 * t215, qJDD(2) * t213 - t209 * t215, 0, t209 * t235 + t213 * t225, -t209 * t225 + t213 * t235, -t104 * t140 + t151 * t68, t104 * t341 + t105 * t140 - t151 * t69 + t232 * t68, t104 * t201 + t151 * t200, -t105 * t201 + t200 * t232, 0, t105 * t165 - t135 * t232 + t194 * t69 + t198 * t248 + t200 * t233 - t201 * t63 - t278 * t341, t104 * t165 - t116 * t200 + t135 * t151 - t140 * t278 + t194 * t68 - t196 * t248 - t201 * t62, -t230 * t239 + t298 * t45, (-t112 * t212 + t207 * t239) * t104 + (-t314 - t212 * t46 + (t112 * t207 + t212 * t239) * qJD(4)) * t151, -t105 * t239 + t134 * t230 - t232 * t45 + t298 * t67, -t105 * t112 - t134 * t231 + t232 * t46 - t299 * t67, t105 * t134 - t232 * t67, t63 * t112 - t233 * t46 + t83 * t269 + (-t116 * t285 + t51) * t134 + t91 * t67 - t251 * t232 + t47 * t105 - g(1) * t130 - g(2) * t132 + (t36 * t151 + t83 * t104 + (-qJD(4) * t95 - t62) * t134 - t116 * t67 - t259 * t232) * t207, -t63 * t239 - t233 * t45 + t83 * t306 - t229 * t134 - t309 * t67 + t234 * t232 - t48 * t105 - g(1) * t129 - g(2) * t131 + (t36 * t212 - t74) * t151, -t13 * t89 - t21 * t240, -t13 * t88 - t21 * t58 + t22 * t240 - t223 * t89, -t105 * t240 + t127 * t21 - t13 * t232 - t65 * t89, -t105 * t58 - t127 * t22 - t223 * t232 - t65 * t88, t105 * t127 - t232 * t65, (-t10 * t206 + t211 * t7) * t127 + (-t206 * t49 + t211 * t38) * t65 - t274 * t232 + t8 * t105 + t32 * t58 - t81 * t223 + t17 * t88 + t55 * t22 - g(1) * t120 - g(2) * t122 + ((-t206 * t38 - t211 * t49) * t127 + t9 * t232) * qJD(5), -g(1) * t119 - g(2) * t121 - t9 * t105 + t81 * t13 - t26 * t232 - t17 * t89 + t55 * t21 - t32 * t240 + (-(-qJD(5) * t49 + t7) * t127 - t38 * t65 + t3 * t232) * t206 + (-(qJD(5) * t38 + t10) * t127 - t49 * t65 + t263 * t232) * t211; 0, 0, 0, -t209 * t216 * t213, t289 * t216, t281, t280, qJDD(2), -g(3) * t213 + t209 * t224, g(3) * t209 + t213 * t224, t300, t71, t52, t53, t200, t99 * t201 + (t200 * t323 - t201 * t287 + t288 * t341) * pkin(2) + t222, t100 * t201 + (t140 * t288 - t200 * t208 - t201 * t266) * pkin(2) + t217, t18, t6, t20, t19, t302, t192 * t46 - t337 * t212 + t242 * t207 + t254 * t112 + (-t191 * t285 + t340) * t134 + t237, t192 * t45 + t242 * t212 - t207 * t343 - t254 * t239 + (t191 * t286 + t327) * t134 + t238, t4, t1, t11, t12, t303, (t144 * t211 - t145 * t206) * t65 - t164 * t223 + t308 * t58 + (t206 * t247 - t211 * t246) * t127 + t221, -(t144 * t206 + t145 * t211) * t65 + t164 * t13 - t308 * t240 + (t206 * t246 + t211 * t247) * t127 + t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t300, t71, t52, t53, t200, t201 * t97 + t222, t201 * t96 + t217, t18, t6, t20, t19, t302, -pkin(3) * t46 - t97 * t112 - t87 * t134 + (t134 * t96 + t250) * t207 + (-t337 - t276) * t212 + t237, -pkin(3) * t45 + t97 * t239 + t317 * t134 + t250 * t212 + (-t343 + t276) * t207 + t238, t4, t1, t11, t12, t303, (t166 * t211 - t168 * t206) * t65 - t193 * t223 + t253 * t58 + (t206 * t245 - t211 * t244) * t127 + t221, -(t166 * t206 + t168 * t211) * t65 + t193 * t13 - t253 * t240 + (t206 * t244 + t211 * t245) * t127 + t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t239 * t112, -t112 ^ 2 + t239 ^ 2, t45 + t305, -t46 - t304, t67, -g(1) * t131 + g(2) * t129 + t239 * t83 + t134 * t48 + (t259 + t188) * t207 + t251, g(1) * t132 - g(2) * t130 + t112 * t83 + t134 * t47 + t188 * t212 - t234, -t339, t333, t332, t329, t65, -(-t206 * t28 - t313) * t127 + (-t127 * t284 + t211 * t65 + t239 * t58) * pkin(4) + t330, (-t127 * t29 - t3) * t206 + (t127 * t28 - t263) * t211 + (-t127 * t283 - t206 * t65 - t239 * t240) * pkin(4) + t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t339, t333, t332, t329, t65, t127 * t9 + t330, t127 * t8 - t206 * t3 - t263 * t211 + t331;];
tau_reg = t2;
