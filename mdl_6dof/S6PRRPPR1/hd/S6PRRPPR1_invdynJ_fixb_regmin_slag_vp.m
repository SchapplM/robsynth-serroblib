% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:04
% EndTime: 2021-01-16 02:06:24
% DurationCPUTime: 5.93s
% Computational Cost: add. (5081->489), mult. (12208->690), div. (0->0), fcn. (9902->18), ass. (0->244)
t202 = sin(qJ(3));
t203 = sin(qJ(2));
t197 = sin(pkin(6));
t287 = qJD(1) * t197;
t268 = t203 * t287;
t321 = qJD(3) * pkin(3);
t343 = t202 * t321 - t268;
t195 = sin(pkin(11));
t205 = cos(qJ(3));
t311 = cos(pkin(11));
t254 = t311 * t202;
t155 = t195 * t205 + t254;
t143 = t155 * qJD(3);
t253 = t311 * t205;
t300 = t195 * t202;
t225 = t253 - t300;
t146 = t225 * qJD(3);
t342 = pkin(4) * t143 - qJ(5) * t146 - qJD(5) * t155 + t343;
t200 = qJ(4) + pkin(8);
t259 = qJD(3) * t200;
t135 = qJD(4) * t205 - t202 * t259;
t221 = -qJD(4) * t202 - t205 * t259;
t206 = cos(qJ(2));
t267 = t206 * t287;
t313 = t135 * t311 + t195 * t221 - t225 * t267;
t173 = qJD(2) * t253;
t285 = qJD(2) * t202;
t141 = t195 * t285 - t173;
t134 = qJD(6) + t141;
t341 = t134 - qJD(6);
t144 = t155 * qJD(2);
t194 = sin(pkin(12));
t198 = cos(pkin(12));
t121 = qJD(3) * t194 + t144 * t198;
t122 = qJD(3) * t198 - t144 * t194;
t201 = sin(qJ(6));
t204 = cos(qJ(6));
t60 = t121 * t201 - t122 * t204;
t340 = t134 * t60;
t156 = t194 * t204 + t198 * t201;
t148 = t156 * qJD(6);
t316 = t156 * t141 + t148;
t234 = -t121 * t204 - t122 * t201;
t339 = t134 * t234;
t199 = cos(pkin(6));
t296 = t197 * t203;
t149 = t199 * t205 - t202 * t296;
t324 = -t194 * t313 + t198 * t342;
t323 = t194 * t342 + t198 * t313;
t277 = t199 * qJDD(1);
t172 = t205 * t277;
t279 = qJD(1) * qJD(2);
t127 = qJDD(2) * pkin(8) + (qJDD(1) * t203 + t206 * t279) * t197;
t286 = qJD(1) * t199;
t211 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t286 + t127;
t251 = qJD(2) * t200 + t268;
t232 = t251 * qJD(3);
t43 = qJDD(3) * pkin(3) - t202 * t211 - t205 * t232 + t172;
t44 = (-t232 + t277) * t202 + t211 * t205;
t13 = -t195 * t44 + t311 * t43;
t12 = -qJDD(3) * pkin(4) + qJDD(5) - t13;
t312 = cos(pkin(10));
t257 = t312 * t203;
t196 = sin(pkin(10));
t297 = t196 * t206;
t139 = t199 * t257 + t297;
t191 = qJ(3) + pkin(11);
t187 = sin(t191);
t189 = cos(t191);
t258 = t197 * t312;
t100 = t139 * t187 + t189 * t258;
t256 = t312 * t206;
t298 = t196 * t203;
t137 = t199 * t298 - t256;
t299 = t196 * t197;
t102 = t137 * t187 + t189 * t299;
t129 = t187 * t296 - t189 * t199;
t223 = g(1) * t102 - g(2) * t100 - g(3) * t129;
t220 = t12 + t223;
t113 = t202 * t286 + t205 * t251;
t104 = t195 * t113;
t112 = -t202 * t251 + t205 * t286;
t56 = t112 * t311 - t104;
t273 = pkin(3) * t285;
t83 = pkin(4) * t144 + qJ(5) * t141 + t273;
t28 = t194 * t83 + t198 * t56;
t338 = qJD(5) * t198 - t28;
t27 = -t194 * t56 + t198 * t83;
t337 = -qJD(5) * t194 - t27;
t314 = t135 * t195 - t155 * t267 - t221 * t311;
t138 = -t199 * t256 + t298;
t140 = t199 * t297 + t257;
t246 = g(1) * t140 + g(2) * t138;
t294 = t197 * t206;
t335 = g(3) * t294 - t246;
t336 = t335 * t187;
t154 = t194 * t201 - t198 * t204;
t317 = t134 * t154;
t276 = t202 * qJDD(2);
t239 = -qJDD(2) * t253 + t195 * t276;
t96 = qJD(2) * t143 + t239;
t92 = qJDD(6) + t96;
t334 = t134 * t317 - t156 * t92;
t207 = qJD(3) ^ 2;
t261 = qJDD(1) * t294;
t264 = t203 * t279;
t240 = t197 * t264 - t261;
t333 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t207 + t197 * (-g(3) * t206 + t264) - t240 + t246;
t278 = qJD(2) * qJD(3);
t263 = t202 * t278;
t213 = qJDD(2) * t155 - t195 * t263;
t97 = qJD(3) * t173 + t213;
t79 = -qJDD(3) * t198 + t194 * t97;
t80 = qJDD(3) * t194 + t198 * t97;
t17 = -qJD(6) * t234 + t201 * t80 + t204 * t79;
t136 = t141 ^ 2;
t332 = pkin(3) * t195;
t331 = pkin(3) * t202;
t330 = pkin(9) * t198;
t328 = g(3) * t197;
t178 = qJ(5) + t332;
t327 = pkin(9) + t178;
t14 = t195 * t43 + t311 * t44;
t11 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t14;
t183 = pkin(3) * t205 + pkin(2);
t95 = pkin(3) * t263 - qJDD(2) * t183 + qJDD(4) + t240;
t26 = pkin(4) * t96 - qJ(5) * t97 - qJD(5) * t144 + t95;
t7 = t11 * t198 + t194 * t26;
t326 = pkin(5) * t143 - t146 * t330 + t324;
t308 = t146 * t194;
t325 = pkin(9) * t308 - t323;
t108 = t112 + t321;
t255 = t311 * t113;
t52 = t108 * t195 + t255;
t47 = qJD(3) * qJ(5) + t52;
t133 = -qJD(2) * t183 + qJD(4) - t267;
t66 = pkin(4) * t141 - qJ(5) * t144 + t133;
t22 = t194 * t66 + t198 * t47;
t322 = qJD(2) * pkin(2);
t320 = t144 * t60;
t319 = t144 * t234;
t165 = t200 * t205;
t118 = t165 * t311 - t200 * t300;
t93 = -pkin(4) * t225 - qJ(5) * t155 - t183;
t49 = t118 * t198 + t194 * t93;
t315 = pkin(5) * t308 + t314;
t309 = t141 * t194;
t307 = t155 * t194;
t306 = t155 * t198;
t305 = t178 * t194;
t304 = t178 * t198;
t190 = pkin(12) + qJ(6);
t186 = sin(t190);
t303 = t186 * t189;
t188 = cos(t190);
t302 = t188 * t189;
t301 = t189 * t206;
t295 = t197 * t205;
t51 = t108 * t311 - t104;
t46 = -qJD(3) * pkin(4) + qJD(5) - t51;
t292 = -qJD(5) + t46;
t291 = qJDD(1) - g(3);
t290 = -t137 * t200 - t140 * t183;
t289 = t183 * t294 + t200 * t296;
t192 = t202 ^ 2;
t288 = -t205 ^ 2 + t192;
t284 = qJD(2) * t203;
t281 = qJD(6) * t201;
t280 = qJD(6) * t204;
t275 = t205 * qJDD(2);
t274 = g(3) * t296;
t6 = -t11 * t194 + t198 * t26;
t2 = pkin(5) * t96 - pkin(9) * t80 + t6;
t3 = -pkin(9) * t79 + t7;
t270 = t2 * t204 - t201 * t3;
t269 = t311 * pkin(3);
t266 = t197 * t284;
t265 = qJD(2) * t294;
t262 = t205 * t278;
t21 = -t194 * t47 + t198 * t66;
t48 = -t118 * t194 + t198 * t93;
t54 = t112 * t195 + t255;
t252 = -t138 * t183 + t139 * t200;
t99 = t137 * t189 - t187 * t299;
t117 = t165 * t195 + t200 * t254;
t250 = t202 * t265;
t249 = -t134 * t316 - t154 * t92;
t248 = pkin(3) * t196 * t295 + t137 * t331;
t182 = -t269 - pkin(4);
t247 = g(1) * t137 - g(2) * t139;
t244 = -t194 * t7 - t198 * t6;
t243 = t201 * t2 + t204 * t3;
t15 = pkin(9) * t122 + t22;
t9 = pkin(5) * t141 - pkin(9) * t121 + t21;
t5 = t15 * t204 + t201 * t9;
t242 = t15 * t201 - t204 * t9;
t241 = pkin(4) * t189 + qJ(5) * t187;
t31 = -pkin(5) * t225 - pkin(9) * t306 + t48;
t33 = -pkin(9) * t307 + t49;
t238 = -t201 * t33 + t204 * t31;
t237 = t201 * t31 + t204 * t33;
t150 = t199 * t202 + t203 * t295;
t89 = t149 * t195 + t150 * t311;
t67 = -t194 * t89 - t198 * t294;
t68 = -t194 * t294 + t198 * t89;
t236 = -t201 * t68 + t204 * t67;
t235 = t201 * t67 + t204 * t68;
t233 = t149 * pkin(3);
t208 = qJD(2) ^ 2;
t231 = qJDD(2) * t206 - t203 * t208;
t229 = -g(1) * t196 + g(2) * t312;
t151 = t327 * t194;
t227 = pkin(9) * t309 + qJD(6) * t151 - t338;
t152 = t327 * t198;
t226 = pkin(5) * t144 + qJD(6) * t152 + t141 * t330 - t337;
t16 = -t121 * t281 + t122 * t280 - t201 * t79 + t204 * t80;
t101 = t139 * t189 - t187 * t258;
t130 = t187 * t199 + t189 * t296;
t224 = g(1) * t99 - g(2) * t101 - g(3) * t130;
t222 = t231 * t197;
t219 = t150 * qJD(3);
t162 = -t267 - t322;
t217 = -qJD(2) * t162 - t127 - t247;
t216 = (-t139 * t202 - t205 * t258) * pkin(3);
t215 = t335 * t189;
t214 = t12 * t155 + t146 * t46 + t247;
t210 = -pkin(8) * qJDD(3) + (t162 + t267 - t322) * qJD(3);
t209 = -t219 - t250;
t163 = -pkin(5) * t198 + t182;
t111 = qJD(3) * t149 + t205 * t265;
t88 = -t149 * t311 + t150 * t195;
t87 = t154 * t155;
t86 = t156 * t155;
t78 = pkin(5) * t307 + t117;
t55 = t111 * t311 + t195 * t209;
t53 = t111 * t195 - t209 * t311;
t40 = t194 * t266 + t198 * t55;
t39 = -t194 * t55 + t198 * t266;
t36 = -pkin(5) * t309 + t54;
t35 = t146 * t156 + t280 * t306 - t281 * t307;
t34 = -t146 * t154 - t148 * t155;
t32 = -pkin(5) * t122 + t46;
t8 = pkin(5) * t79 + t12;
t1 = [t291, 0, t222, (-qJDD(2) * t203 - t206 * t208) * t197, 0, 0, 0, 0, 0, t149 * qJDD(3) + t205 * t222 + (-t219 - 0.2e1 * t250) * qJD(3), -qJD(3) * t111 - qJDD(3) * t150 + (-t202 * t231 - t206 * t262) * t197, -qJD(3) * t53 - qJDD(3) * t88 + (t141 * t284 - t206 * t96) * t197, -qJD(3) * t55 - qJDD(3) * t89 + (t144 * t284 - t206 * t97) * t197, -t141 * t55 + t144 * t53 + t88 * t97 - t89 * t96, -t13 * t88 + t14 * t89 - t51 * t53 + t52 * t55 - g(3) + (t133 * t284 - t206 * t95) * t197, -t122 * t53 + t141 * t39 + t67 * t96 + t79 * t88, t121 * t53 - t141 * t40 - t68 * t96 + t80 * t88, -t121 * t39 + t122 * t40 - t67 * t80 - t68 * t79, t12 * t88 + t21 * t39 + t22 * t40 + t46 * t53 + t6 * t67 + t68 * t7 - g(3), 0, 0, 0, 0, 0, (-qJD(6) * t235 - t201 * t40 + t204 * t39) * t134 + t236 * t92 + t53 * t60 + t88 * t17, -(qJD(6) * t236 + t201 * t39 + t204 * t40) * t134 - t235 * t92 - t53 * t234 + t88 * t16; 0, qJDD(2), -t335 + t261, -t291 * t296 - t247, qJDD(2) * t192 + 0.2e1 * t202 * t262, 0.2e1 * t202 * t275 - 0.2e1 * t278 * t288, qJDD(3) * t202 + t205 * t207, qJDD(3) * t205 - t202 * t207, 0, t202 * t210 + t205 * t333, -t202 * t333 + t205 * t210, -t141 * t268 - qJDD(3) * t117 + t133 * t143 - t225 * t95 - t183 * t96 - t215 + (t141 * t331 - t314) * qJD(3), -t144 * t268 - qJDD(3) * t118 + t133 * t146 + t155 * t95 - t183 * t97 + t336 + (t144 * t331 - t313) * qJD(3), t117 * t97 - t118 * t96 - t13 * t155 + t14 * t225 - t141 * t313 - t143 * t52 + t144 * t314 - t146 * t51 + t247 - t274, -g(1) * t290 - g(2) * t252 - g(3) * t289 - t13 * t117 + t14 * t118 + t133 * t343 - t95 * t183 + t313 * t52 - t314 * t51, t117 * t79 + t21 * t143 - t6 * t225 + t48 * t96 - t198 * t215 + (t214 - t274) * t194 + t324 * t141 - t314 * t122, t117 * t80 - t22 * t143 + t7 * t225 - t49 * t96 - t246 * t194 * t189 + t214 * t198 - (-t194 * t301 + t198 * t203) * t328 - t323 * t141 + t314 * t121, -t48 * t80 - t49 * t79 + t244 * t155 + (-t194 * t22 - t198 * t21) * t146 - t324 * t121 + t323 * t122 - t336, t7 * t49 + t6 * t48 + t12 * t117 - g(1) * (-t140 * t241 + t290) - g(2) * (-t138 * t241 + t252) - g(3) * (t241 * t294 + t289) + t314 * t46 + t323 * t22 + t324 * t21, -t16 * t87 - t234 * t34, -t16 * t86 + t17 * t87 + t234 * t35 - t34 * t60, t134 * t34 - t143 * t234 - t16 * t225 - t87 * t92, -t134 * t35 - t143 * t60 + t17 * t225 - t86 * t92, t134 * t143 - t225 * t92, t238 * t92 - t270 * t225 - t242 * t143 + t78 * t17 + t8 * t86 + t32 * t35 - g(1) * (-t137 * t186 - t140 * t302) - g(2) * (-t138 * t302 + t139 * t186) + t315 * t60 - (t186 * t203 + t188 * t301) * t328 + (t201 * t325 + t204 * t326) * t134 + (-t134 * t237 + t225 * t5) * qJD(6), -t237 * t92 + t243 * t225 - t5 * t143 + t78 * t16 - t8 * t87 + t32 * t34 - g(1) * (-t137 * t188 + t140 * t303) - g(2) * (t138 * t303 + t139 * t188) - t315 * t234 - (-t186 * t301 + t188 * t203) * t328 + (-t201 * t326 + t204 * t325) * t134 + (-t134 * t238 - t225 * t242) * qJD(6); 0, 0, 0, 0, -t202 * t208 * t205, t288 * t208, t276, t275, qJDD(3), -g(3) * t149 + t202 * t217 + t229 * t295 + t172, g(3) * t150 + (-t197 * t229 - t277) * t202 + t217 * t205, t54 * qJD(3) - t133 * t144 + (qJDD(3) * t311 - t141 * t285) * pkin(3) - t223 + t13, qJD(3) * t56 + t133 * t141 + (-qJDD(3) * t195 - t144 * t285) * pkin(3) - t224 - t14, (t52 - t54) * t144 + (-t51 + t56) * t141 + (-t195 * t96 - t311 * t97) * pkin(3), -g(1) * t248 - g(2) * t216 - g(3) * t233 + t13 * t269 - t133 * t273 + t14 * t332 + t51 * t54 - t52 * t56, -t96 * t305 + t122 * t54 - t144 * t21 + t182 * t79 + (t194 * t292 - t27) * t141 - t220 * t198, -t96 * t304 - t121 * t54 + t144 * t22 + t182 * t80 + (t198 * t292 + t28) * t141 + t220 * t194, -t122 * t28 + t121 * t27 + (qJD(5) * t122 - t141 * t21 - t178 * t79 + t7) * t198 + (qJD(5) * t121 - t141 * t22 + t178 * t80 - t6) * t194 + t224, t7 * t304 - t6 * t305 + t12 * t182 - t46 * t54 - g(1) * (pkin(4) * t102 - qJ(5) * t99 + t248) - g(2) * (-pkin(4) * t100 + qJ(5) * t101 + t216) - g(3) * (-pkin(4) * t129 + qJ(5) * t130 + t233) + t338 * t22 + t337 * t21, t156 * t16 + t234 * t317, -t154 * t16 - t156 * t17 + t234 * t316 + t317 * t60, t319 - t334, t249 + t320, -t134 * t144, (-t151 * t204 - t152 * t201) * t92 + t163 * t17 + t8 * t154 + t242 * t144 - t36 * t60 + t316 * t32 + (t201 * t227 - t204 * t226) * t134 - t223 * t188, -(-t151 * t201 + t152 * t204) * t92 + t163 * t16 + t8 * t156 + t5 * t144 + t36 * t234 - t317 * t32 + (t201 * t226 + t204 * t227) * t134 + t223 * t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t144 * qJD(3) + t239, (t173 - t141) * qJD(3) + t213, -t144 ^ 2 - t136, t141 * t52 + t144 * t51 + t335 + t95, t122 * t144 - t136 * t194 + t198 * t96, -t121 * t144 - t136 * t198 - t194 * t96, -t194 * t79 - t198 * t80 + (t121 * t194 + t122 * t198) * t141, -t144 * t46 + (-t194 * t21 + t198 * t22) * t141 + t335 - t244, 0, 0, 0, 0, 0, t249 - t320, t319 + t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t141 + t79, t122 * t141 + t80, -t121 ^ 2 - t122 ^ 2, t121 * t21 - t122 * t22 + t220, 0, 0, 0, 0, 0, t17 - t339, t16 - t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234 * t60, t234 ^ 2 - t60 ^ 2, t16 + t340, -t17 - t339, t92, t32 * t234 - g(1) * (t140 * t188 + t186 * t99) - g(2) * (-t101 * t186 + t138 * t188) - g(3) * (-t130 * t186 - t188 * t294) + t270 + t341 * t5, t32 * t60 - g(1) * (-t140 * t186 + t188 * t99) - g(2) * (-t101 * t188 - t138 * t186) - g(3) * (-t130 * t188 + t186 * t294) - t243 - t341 * t242;];
tau_reg = t1;
