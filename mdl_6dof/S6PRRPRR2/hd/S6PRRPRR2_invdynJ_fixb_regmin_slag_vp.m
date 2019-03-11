% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:19
% EndTime: 2019-03-08 22:00:33
% DurationCPUTime: 5.50s
% Computational Cost: add. (5101->476), mult. (12107->687), div. (0->0), fcn. (9948->18), ass. (0->249)
t194 = sin(pkin(12));
t203 = sin(qJ(3));
t278 = qJD(2) * t203;
t197 = cos(pkin(12));
t207 = cos(qJ(3));
t288 = t197 * t207;
t149 = qJD(2) * t288 - t194 * t278;
t332 = qJD(5) + qJD(6);
t343 = t149 - t332;
t159 = t194 * t207 + t197 * t203;
t150 = t159 * qJD(3);
t158 = t194 * t203 - t288;
t153 = t158 * qJD(3);
t196 = sin(pkin(6));
t204 = sin(qJ(2));
t279 = qJD(1) * t204;
t262 = t196 * t279;
t320 = qJD(3) * pkin(3);
t266 = t203 * t320;
t342 = -pkin(4) * t150 - pkin(9) * t153 + t262 - t266;
t208 = cos(qJ(2));
t290 = t196 * t208;
t261 = qJD(1) * t290;
t124 = t158 * t261;
t200 = -qJ(4) - pkin(8);
t250 = qJD(3) * t200;
t143 = qJD(4) * t207 + t203 * t250;
t144 = -qJD(4) * t203 + t207 * t250;
t80 = t143 * t197 + t144 * t194;
t308 = t80 + t124;
t151 = t159 * qJD(2);
t202 = sin(qJ(5));
t206 = cos(qJ(5));
t272 = t206 * qJD(3);
t127 = t151 * t202 - t272;
t129 = qJD(3) * t202 + t151 * t206;
t201 = sin(qJ(6));
t205 = cos(qJ(6));
t233 = t127 * t201 - t205 * t129;
t61 = t205 * t127 + t129 * t201;
t341 = t233 * t61;
t285 = t201 * t206;
t162 = t202 * t205 + t285;
t311 = t343 * t162;
t276 = qJD(5) * t202;
t302 = t149 * t202;
t340 = t276 - t302;
t339 = t233 ^ 2 - t61 ^ 2;
t140 = qJD(5) - t149;
t138 = qJD(6) + t140;
t273 = qJD(6) * t205;
t274 = qJD(6) * t201;
t270 = qJD(2) * qJD(3);
t255 = t207 * t270;
t256 = t203 * t270;
t103 = qJDD(2) * t159 - t194 * t256 + t197 * t255;
t46 = qJD(5) * t272 + t202 * qJDD(3) + t206 * t103 - t151 * t276;
t47 = qJD(5) * t129 - t206 * qJDD(3) + t103 * t202;
t8 = -t127 * t273 - t129 * t274 - t201 * t47 + t205 * t46;
t338 = t138 * t61 + t8;
t195 = sin(pkin(11));
t198 = cos(pkin(11));
t199 = cos(pkin(6));
t287 = t199 * t204;
t146 = t195 * t208 + t198 * t287;
t190 = qJ(3) + pkin(12);
t183 = sin(t190);
t184 = cos(t190);
t293 = t196 * t198;
t105 = t146 * t184 - t183 * t293;
t148 = -t195 * t287 + t198 * t208;
t294 = t195 * t196;
t107 = t148 * t184 + t183 * t294;
t292 = t196 * t204;
t136 = t183 * t199 + t184 * t292;
t286 = t199 * t208;
t145 = t195 * t204 - t198 * t286;
t147 = t195 * t286 + t198 * t204;
t245 = -qJD(2) * t200 + t262;
t280 = qJD(1) * t199;
t119 = -t203 * t245 + t207 * t280;
t112 = t119 + t320;
t120 = t203 * t280 + t207 * t245;
t289 = t197 * t120;
t50 = t194 * t112 + t289;
t45 = qJD(3) * pkin(9) + t50;
t182 = pkin(3) * t207 + pkin(2);
t139 = -t182 * qJD(2) + qJD(4) - t261;
t66 = -pkin(4) * t149 - pkin(9) * t151 + t139;
t25 = t202 * t66 + t206 * t45;
t19 = -pkin(10) * t127 + t25;
t17 = t19 * t274;
t193 = qJ(5) + qJ(6);
t188 = sin(t193);
t189 = cos(t193);
t108 = t194 * t120;
t49 = t112 * t197 - t108;
t44 = -qJD(3) * pkin(4) - t49;
t32 = pkin(5) * t127 + t44;
t337 = t32 * t61 - g(1) * (-t107 * t189 - t147 * t188) - g(2) * (-t105 * t189 - t145 * t188) - g(3) * (-t136 * t189 + t188 * t290) + t17;
t269 = t199 * qJDD(1);
t173 = t207 * t269;
t271 = qJD(1) * qJD(2);
t131 = qJDD(2) * pkin(8) + (qJDD(1) * t204 + t208 * t271) * t196;
t215 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t280 + t131;
t232 = t245 * qJD(3);
t39 = qJDD(3) * pkin(3) - t203 * t215 - t207 * t232 + t173;
t40 = (-t232 + t269) * t203 + t215 * t207;
t16 = t194 * t39 + t197 * t40;
t14 = qJDD(3) * pkin(9) + t16;
t257 = t204 * t271;
t171 = t196 * t257;
t218 = pkin(3) * t256 - t182 * qJDD(2) + qJDD(4) + t171;
t253 = qJDD(1) * t290;
t101 = t218 - t253;
t267 = t207 * qJDD(2);
t268 = t203 * qJDD(2);
t237 = -t194 * t268 + t197 * t267;
t102 = -qJD(3) * t151 + t237;
t30 = -pkin(4) * t102 - pkin(9) * t103 + t101;
t29 = t206 * t30;
t217 = -qJD(5) * t25 - t202 * t14 + t29;
t98 = qJD(2) * t150 + qJDD(5) - t237;
t2 = pkin(5) * t98 - pkin(10) * t46 + t217;
t275 = qJD(5) * t206;
t227 = -t206 * t14 - t202 * t30 - t66 * t275 + t45 * t276;
t3 = -pkin(10) * t47 - t227;
t263 = t205 * t2 - t201 * t3;
t24 = -t202 * t45 + t206 * t66;
t18 = -pkin(10) * t129 + t24;
t11 = pkin(5) * t140 + t18;
t316 = t19 * t205;
t5 = t11 * t201 + t316;
t336 = t32 * t233 - g(1) * (-t107 * t188 + t147 * t189) - g(2) * (-t105 * t188 + t145 * t189) - g(3) * (-t136 * t188 - t189 * t290) - qJD(6) * t5 + t263;
t214 = qJD(6) * t233 - t201 * t46 - t205 * t47;
t335 = -t138 * t233 + t214;
t334 = -t124 * t202 - t342 * t206;
t100 = pkin(4) * t158 - pkin(9) * t159 - t182;
t166 = t200 * t203;
t167 = t200 * t207;
t126 = t166 * t194 - t167 * t197;
t333 = -t100 * t275 + t126 * t276 + t342 * t202 - t308 * t206;
t309 = t143 * t194 - t197 * t144 - t159 * t261;
t89 = t162 * t159;
t161 = t201 * t202 - t205 * t206;
t312 = t343 * t161;
t94 = qJDD(6) + t98;
t331 = -t312 * t138 - t162 * t94;
t209 = qJD(3) ^ 2;
t242 = g(1) * t147 + g(2) * t145;
t330 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t209 + t196 * (-g(3) * t208 + t257) - t171 + t242 + t253;
t252 = qJD(6) * t11 + t3;
t329 = t201 * t2 + t205 * t252;
t113 = t206 * t126;
t301 = t153 * t206;
t328 = pkin(10) * t301 + pkin(5) * t150 - t202 * t80 + (-t113 + (pkin(10) * t159 - t100) * t202) * qJD(5) + t334;
t154 = t199 * t207 - t203 * t292;
t327 = g(3) * t154;
t326 = g(3) * t196;
t178 = pkin(3) * t194 + pkin(9);
t324 = pkin(10) + t178;
t284 = t202 * t153;
t226 = t159 * t275 - t284;
t323 = pkin(10) * t226 + t333;
t54 = t119 * t197 - t108;
t81 = pkin(3) * t278 + pkin(4) * t151 - pkin(9) * t149;
t322 = t202 * t81 + t206 * t54;
t321 = qJD(2) * pkin(2);
t319 = t151 * t61;
t318 = t151 * t233;
t315 = t202 * t46;
t314 = t202 * t98;
t313 = t202 * t100 + t113;
t310 = pkin(5) * t226 + t309;
t306 = t127 * t140;
t305 = t127 * t151;
t304 = t129 * t140;
t303 = t129 * t151;
t300 = t159 * t202;
t299 = t159 * t206;
t298 = t184 * t188;
t297 = t184 * t189;
t296 = t184 * t202;
t295 = t184 * t208;
t291 = t196 * t207;
t283 = t202 * t208;
t282 = qJDD(1) - g(3);
t191 = t203 ^ 2;
t281 = -t207 ^ 2 + t191;
t277 = qJD(2) * t204;
t265 = t196 * t283;
t264 = t206 * t290;
t179 = -pkin(3) * t197 - pkin(4);
t260 = t196 * t277;
t259 = qJD(2) * t290;
t258 = t159 * t276;
t254 = t208 * t270;
t15 = -t194 * t40 + t197 * t39;
t249 = qJD(5) * t324;
t248 = t196 * t282;
t52 = t119 * t194 + t289;
t125 = -t197 * t166 - t167 * t194;
t246 = t140 * t206;
t244 = t138 * t311 - t161 * t94;
t243 = pkin(5) * t340 - t52;
t241 = g(1) * t148 + g(2) * t146;
t240 = g(1) * t195 - g(2) * t198;
t156 = t324 * t202;
t239 = -pkin(10) * t302 + qJD(6) * t156 + t202 * t249 + t322;
t157 = t324 * t206;
t74 = t206 * t81;
t238 = pkin(5) * t151 + qJD(6) * t157 - t202 * t54 + t74 + (-pkin(10) * t149 + t249) * t206;
t92 = t206 * t100;
t31 = pkin(5) * t158 - pkin(10) * t299 - t202 * t126 + t92;
t33 = -pkin(10) * t300 + t313;
t236 = t201 * t31 + t205 * t33;
t155 = t199 * t203 + t204 * t291;
t87 = t154 * t194 + t155 * t197;
t229 = -t206 * t87 + t265;
t67 = -t202 * t87 - t264;
t235 = t201 * t229 + t205 * t67;
t234 = t201 * t67 - t205 * t229;
t210 = qJD(2) ^ 2;
t231 = qJDD(2) * t208 - t204 * t210;
t230 = -t140 * t340 + t206 * t98;
t13 = -qJDD(3) * pkin(4) - t15;
t225 = -t258 - t301;
t223 = t140 * t44 - t178 * t98;
t222 = g(1) * (-t148 * t183 + t184 * t294) + g(2) * (-t146 * t183 - t184 * t293) + g(3) * (-t183 * t292 + t184 * t199);
t220 = -g(3) * t292 - t241;
t164 = -t261 - t321;
t219 = -qJD(2) * t164 - t131 + t241;
t213 = qJD(5) * t140 * t178 + t13 + t222;
t212 = -pkin(8) * qJDD(3) + (t164 + t261 - t321) * qJD(3);
t165 = -pkin(5) * t206 + t179;
t118 = -qJD(3) * t155 - t203 * t259;
t117 = qJD(3) * t154 + t207 * t259;
t90 = t161 * t159;
t86 = -t197 * t154 + t155 * t194;
t78 = pkin(5) * t300 + t125;
t53 = t117 * t197 + t118 * t194;
t51 = t117 * t194 - t197 * t118;
t27 = -t153 * t285 - t201 * t258 - t274 * t300 + (t299 * t332 - t284) * t205;
t26 = t161 * t153 - t332 * t89;
t22 = qJD(5) * t229 - t202 * t53 + t206 * t260;
t21 = qJD(5) * t67 + t202 * t260 + t206 * t53;
t6 = pkin(5) * t47 + t13;
t4 = t11 * t205 - t19 * t201;
t1 = [t282, 0, t231 * t196 (-qJDD(2) * t204 - t208 * t210) * t196, 0, 0, 0, 0, 0, qJD(3) * t118 + qJDD(3) * t154 + (-t203 * t254 + t207 * t231) * t196, -qJD(3) * t117 - qJDD(3) * t155 + (-t203 * t231 - t207 * t254) * t196, t102 * t87 + t103 * t86 + t149 * t53 + t151 * t51, -t15 * t86 + t16 * t87 - t49 * t51 + t50 * t53 - g(3) + (-t101 * t208 + t139 * t277) * t196, 0, 0, 0, 0, 0, t127 * t51 + t140 * t22 + t47 * t86 + t67 * t98, t129 * t51 - t140 * t21 + t229 * t98 + t46 * t86, 0, 0, 0, 0, 0 (-qJD(6) * t234 - t201 * t21 + t205 * t22) * t138 + t235 * t94 + t51 * t61 - t86 * t214 -(qJD(6) * t235 + t201 * t22 + t205 * t21) * t138 - t234 * t94 - t51 * t233 + t86 * t8; 0, qJDD(2), t282 * t290 + t242, -t204 * t248 + t241, qJDD(2) * t191 + 0.2e1 * t203 * t255, 0.2e1 * t203 * t267 - 0.2e1 * t281 * t270, qJDD(3) * t203 + t207 * t209, qJDD(3) * t207 - t203 * t209, 0, t212 * t203 + t207 * t330, -t203 * t330 + t212 * t207, t102 * t126 + t103 * t125 + t308 * t149 - t15 * t159 - t150 * t50 + t309 * t151 + t153 * t49 - t158 * t16 + t220, t16 * t126 - t15 * t125 - t101 * t182 + t139 * t266 - g(1) * (-t147 * t182 - t148 * t200) - g(2) * (-t145 * t182 - t146 * t200) + t308 * t50 - t309 * t49 + (-t139 * t279 - g(3) * (t182 * t208 - t200 * t204)) * t196, t129 * t225 + t46 * t299 -(-t127 * t206 - t129 * t202) * t153 + (-t315 - t206 * t47 + (t127 * t202 - t129 * t206) * qJD(5)) * t159, t129 * t150 + t140 * t225 + t158 * t46 + t98 * t299, -t127 * t150 - t140 * t226 - t158 * t47 - t98 * t300, t140 * t150 + t158 * t98, t125 * t47 + t24 * t150 + t29 * t158 + t92 * t98 + t334 * t140 + t309 * t127 + ((-g(3) * t290 + t242) * t184 + (-t126 * t140 - t158 * t45 + t159 * t44) * qJD(5)) * t206 + ((-qJD(5) * t100 - t80) * t140 - t126 * t98 + (-qJD(5) * t66 - t14) * t158 + t13 * t159 - t44 * t153 + t220) * t202, -t313 * t98 + t227 * t158 - t25 * t150 + t125 * t46 - t44 * t301 - g(1) * (t147 * t296 + t148 * t206) - g(2) * (t145 * t296 + t146 * t206) - (-t184 * t283 + t204 * t206) * t326 + (t13 * t206 - t276 * t44) * t159 + t333 * t140 + t309 * t129, -t233 * t26 - t8 * t90, -t214 * t90 + t233 * t27 - t26 * t61 - t8 * t89, t138 * t26 - t150 * t233 + t158 * t8 - t90 * t94, -t138 * t27 - t150 * t61 + t158 * t214 - t89 * t94, t138 * t150 + t158 * t94 (-t201 * t33 + t205 * t31) * t94 + t263 * t158 + t4 * t150 - t78 * t214 + t6 * t89 + t32 * t27 - g(1) * (-t147 * t297 + t148 * t188) - g(2) * (-t145 * t297 + t146 * t188) + t310 * t61 - (t188 * t204 + t189 * t295) * t326 + (t323 * t201 + t205 * t328) * t138 + (-t138 * t236 - t158 * t5) * qJD(6), -t236 * t94 - (-t17 + t329) * t158 - t5 * t150 + t78 * t8 - t6 * t90 + t32 * t26 - g(1) * (t147 * t298 + t148 * t189) - g(2) * (t145 * t298 + t146 * t189) - t310 * t233 - (-t188 * t295 + t189 * t204) * t326 + ((-qJD(6) * t31 + t323) * t205 + (qJD(6) * t33 - t328) * t201) * t138; 0, 0, 0, 0, -t203 * t210 * t207, t281 * t210, t268, t267, qJDD(3), t203 * t219 - t240 * t291 + t173 - t327, g(3) * t155 + (t196 * t240 - t269) * t203 + t219 * t207 (t50 - t52) * t151 + (t49 - t54) * t149 + (t102 * t194 - t103 * t197) * pkin(3), t49 * t52 - t50 * t54 + (t16 * t194 + t15 * t197 - t139 * t278 - g(1) * (-t148 * t203 + t195 * t291) - g(2) * (-t146 * t203 - t198 * t291) - t327) * pkin(3), t129 * t246 + t315 (t46 - t306) * t206 + (-t47 - t304) * t202, t140 * t246 - t303 + t314, t230 + t305, -t140 * t151, -t52 * t127 - t74 * t140 - t24 * t151 + t179 * t47 + (t54 * t140 + t223) * t202 - t213 * t206, -t52 * t129 + t322 * t140 + t25 * t151 + t179 * t46 + t213 * t202 + t223 * t206, t162 * t8 - t233 * t312, -t161 * t8 + t162 * t214 - t233 * t311 - t312 * t61, t318 - t331, t244 + t319, -t138 * t151 (-t156 * t205 - t157 * t201) * t94 - t165 * t214 + t6 * t161 - t4 * t151 + t243 * t61 - t311 * t32 + (t201 * t239 - t205 * t238) * t138 - t222 * t189 -(-t156 * t201 + t157 * t205) * t94 + t165 * t8 + t6 * t162 + t5 * t151 - t243 * t233 + t312 * t32 + (t201 * t238 + t205 * t239) * t138 + t222 * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149 ^ 2 - t151 ^ 2, -t149 * t50 + t151 * t49 - t208 * t248 + t218 - t242, 0, 0, 0, 0, 0, t230 - t305, -t140 ^ 2 * t206 - t303 - t314, 0, 0, 0, 0, 0, t244 - t319, t318 + t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129 * t127, -t127 ^ 2 + t129 ^ 2, t46 + t306, t304 - t47, t98, t25 * t140 - t44 * t129 - g(1) * (-t107 * t202 + t147 * t206) - g(2) * (-t105 * t202 + t145 * t206) - g(3) * (-t136 * t202 - t264) + t217, t24 * t140 + t44 * t127 - g(1) * (-t107 * t206 - t147 * t202) - g(2) * (-t105 * t206 - t145 * t202) - g(3) * (-t136 * t206 + t265) + t227, -t341, t339, t338, t335, t94 -(-t18 * t201 - t316) * t138 + (-t129 * t61 - t138 * t274 + t205 * t94) * pkin(5) + t336 (-t138 * t19 - t2) * t201 + (t138 * t18 - t252) * t205 + (t129 * t233 - t138 * t273 - t201 * t94) * pkin(5) + t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t341, t339, t338, t335, t94, t138 * t5 + t336, t138 * t4 - t329 + t337;];
tau_reg  = t1;
