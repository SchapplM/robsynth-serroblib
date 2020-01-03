% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR14_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:40
% EndTime: 2019-12-31 20:38:56
% DurationCPUTime: 5.84s
% Computational Cost: add. (6123->476), mult. (15566->679), div. (0->0), fcn. (12844->14), ass. (0->244)
t204 = sin(qJ(5));
t208 = cos(qJ(5));
t202 = cos(pkin(10));
t201 = sin(pkin(5));
t210 = cos(qJ(2));
t285 = qJD(1) * qJD(2);
t268 = t210 * t285;
t206 = sin(qJ(2));
t283 = qJDD(1) * t206;
t231 = t268 + t283;
t219 = t231 * t201;
t203 = cos(pkin(5));
t284 = qJDD(1) * t203;
t186 = qJDD(2) + t284;
t200 = sin(pkin(10));
t315 = t186 * t200;
t100 = t202 * t219 + t315;
t205 = sin(qJ(4));
t209 = cos(qJ(4));
t297 = qJD(1) * t203;
t187 = qJD(2) + t297;
t298 = qJD(1) * t201;
t273 = t206 * t298;
t130 = t187 * t202 - t200 * t273;
t131 = t187 * t200 + t202 * t273;
t238 = t130 * t205 + t131 * t209;
t162 = t202 * t186;
t99 = t200 * t219 - t162;
t34 = t238 * qJD(4) + t100 * t205 + t209 * t99;
t32 = qJDD(5) + t34;
t78 = -t209 * t130 + t131 * t205;
t75 = qJD(5) + t78;
t344 = t75 ^ 2;
t345 = -t344 * t204 + t208 * t32;
t296 = qJD(1) * t210;
t272 = t201 * t296;
t172 = -qJD(4) + t272;
t343 = t172 * t78;
t156 = t200 * t205 - t209 * t202;
t310 = t201 * t210;
t223 = t156 * t310;
t117 = qJD(1) * t223;
t149 = t156 * qJD(4);
t302 = t117 - t149;
t157 = t200 * t209 + t202 * t205;
t224 = t157 * t310;
t301 = -qJD(1) * t224 + t157 * qJD(4);
t207 = sin(qJ(1));
t305 = t207 * t210;
t211 = cos(qJ(1));
t306 = t206 * t211;
t152 = t203 * t306 + t305;
t197 = pkin(10) + qJ(4);
t194 = sin(t197);
t195 = cos(t197);
t309 = t201 * t211;
t102 = t152 * t195 - t194 * t309;
t303 = t210 * t211;
t307 = t206 * t207;
t151 = -t203 * t303 + t307;
t342 = t102 * t204 - t151 * t208;
t341 = t102 * t208 + t151 * t204;
t196 = t201 ^ 2;
t281 = 0.2e1 * t196;
t53 = t208 * t172 + t204 * t238;
t340 = t238 * t53;
t239 = t172 * t204 - t208 * t238;
t339 = t238 * t239;
t324 = pkin(8) + qJ(3);
t167 = t324 * t200;
t168 = t324 * t202;
t115 = -t167 * t205 + t168 * t209;
t228 = (-pkin(8) * t202 * t210 + pkin(3) * t206) * t201;
t247 = pkin(2) * t206 - qJ(3) * t210;
t143 = t247 * t298;
t280 = pkin(1) * t297;
t144 = -pkin(7) * t273 + t210 * t280;
t87 = t202 * t143 - t144 * t200;
t59 = qJD(1) * t228 + t87;
t257 = t200 * t272;
t88 = t200 * t143 + t202 * t144;
t71 = -pkin(8) * t257 + t88;
t338 = t157 * qJD(3) + t115 * qJD(4) - t205 * t71 + t209 * t59;
t237 = -t167 * t209 - t168 * t205;
t337 = -t156 * qJD(3) + t237 * qJD(4) - t205 * t59 - t209 * t71;
t336 = t172 * t238;
t335 = -pkin(2) * t186 + qJDD(3);
t154 = -t203 * t307 + t303;
t311 = t201 * t207;
t104 = -t154 * t194 + t195 * t311;
t312 = t201 * t206;
t317 = t152 * t194;
t229 = g(1) * t104 + g(2) * (-t195 * t309 - t317) + g(3) * (-t194 * t312 + t195 * t203);
t282 = qJDD(1) * t210;
t185 = t201 * t282;
t269 = t206 * t285;
t255 = t201 * t269;
t142 = qJDD(4) - t185 + t255;
t180 = pkin(7) * t272;
t145 = t206 * t280 + t180;
t118 = qJ(3) * t187 + t145;
t234 = -pkin(2) * t210 - qJ(3) * t206 - pkin(1);
t139 = t234 * t201;
t123 = qJD(1) * t139;
t61 = -t118 * t200 + t202 * t123;
t40 = -pkin(3) * t272 - pkin(8) * t131 + t61;
t62 = t202 * t118 + t200 * t123;
t43 = pkin(8) * t130 + t62;
t16 = t205 * t40 + t209 * t43;
t230 = t269 - t282;
t279 = pkin(1) * qJD(2) * t203;
t259 = qJD(1) * t279;
t278 = pkin(1) * t284;
t274 = -pkin(7) * t185 - t206 * t278 - t210 * t259;
t217 = -pkin(7) * t255 - t274;
t66 = qJ(3) * t186 + qJD(3) * t187 + t217;
t216 = t247 * qJD(2) - qJD(3) * t206;
t72 = (qJD(1) * t216 + qJDD(1) * t234) * t201;
t36 = -t200 * t66 + t202 * t72;
t19 = pkin(3) * t201 * t230 - pkin(8) * t100 + t36;
t37 = t200 * t72 + t202 * t66;
t23 = -pkin(8) * t99 + t37;
t215 = -t16 * qJD(4) + t209 * t19 - t205 * t23;
t4 = -pkin(4) * t142 - t215;
t334 = (pkin(4) * t238 + t75 * pkin(9)) * t75 + t229 + t4;
t148 = t200 * t203 + t202 * t312;
t331 = pkin(1) * t206;
t300 = pkin(7) * t310 + t203 * t331;
t138 = qJ(3) * t203 + t300;
t83 = -t138 * t200 + t202 * t139;
t48 = -pkin(3) * t310 - pkin(8) * t148 + t83;
t147 = t200 * t312 - t203 * t202;
t84 = t202 * t138 + t200 * t139;
t58 = -pkin(8) * t147 + t84;
t241 = t205 * t48 + t209 * t58;
t121 = t216 * t201;
t293 = qJD(2) * t206;
t271 = t201 * t293;
t236 = -pkin(7) * t271 + t210 * t279;
t127 = qJD(3) * t203 + t236;
t67 = t202 * t121 - t127 * t200;
t46 = qJD(2) * t228 + t67;
t292 = qJD(2) * t210;
t270 = t201 * t292;
t256 = t200 * t270;
t68 = t200 * t121 + t202 * t127;
t56 = -pkin(8) * t256 + t68;
t333 = -t241 * qJD(4) - t205 * t56 + t209 * t46;
t290 = qJD(4) * t209;
t291 = qJD(4) * t205;
t33 = t209 * t100 + t130 * t290 - t131 * t291 - t205 * t99;
t12 = -t239 * qJD(5) - t208 * t142 + t204 * t33;
t14 = -pkin(9) * t172 + t16;
t111 = -pkin(2) * t187 + qJD(3) - t144;
t76 = -pkin(3) * t130 + t111;
t22 = pkin(4) * t78 - pkin(9) * t238 + t76;
t246 = t14 * t204 - t208 * t22;
t233 = t205 * t19 + t209 * t23 + t40 * t290 - t43 * t291;
t3 = pkin(9) * t142 + t233;
t188 = pkin(7) * t312;
t258 = qJD(2) * t180 + qJDD(1) * t188 + t206 * t259 - t210 * t278;
t82 = t258 + t335;
t42 = pkin(3) * t99 + t82;
t8 = pkin(4) * t34 - pkin(9) * t33 + t42;
t1 = -t246 * qJD(5) + t204 * t8 + t208 * t3;
t109 = pkin(3) * t257 + t145;
t153 = t203 * t305 + t306;
t253 = g(1) * t153 + g(2) * t151;
t193 = -pkin(3) * t202 - pkin(2);
t98 = pkin(4) * t156 - pkin(9) * t157 + t193;
t332 = t195 * t253 - (-t301 * pkin(4) + t302 * pkin(9) + qJD(5) * t115 + t109) * t75 + t98 * t32;
t329 = g(1) * t211;
t328 = g(3) * t201;
t327 = t53 * t75;
t326 = t239 * t75;
t322 = pkin(4) * t273 + t338;
t288 = qJD(5) * t208;
t289 = qJD(5) * t204;
t11 = t204 * t142 - t172 * t288 + t208 * t33 - t238 * t289;
t321 = t11 * t204;
t320 = t204 * t32;
t316 = t157 * t208;
t314 = t186 * t203;
t313 = t196 * qJD(1) ^ 2;
t308 = t204 * t210;
t304 = t208 * t210;
t146 = pkin(7) * t270 + t206 * t279;
t198 = t206 ^ 2;
t299 = -t210 ^ 2 + t198;
t295 = qJD(2) * t200;
t294 = qJD(2) * t202;
t287 = qJD(2) - t187;
t286 = qJ(3) * qJDD(1);
t277 = t210 * t313;
t276 = t201 * t308;
t275 = t201 * t304;
t110 = pkin(3) * t256 + t146;
t262 = t208 * t75;
t261 = t187 + t297;
t260 = t186 + t284;
t252 = -g(1) * t151 + g(2) * t153;
t251 = g(1) * t154 + g(2) * t152;
t250 = g(1) * t152 - g(2) * t154;
t249 = g(2) * t207 + t329;
t6 = t14 * t208 + t204 * t22;
t25 = -pkin(9) * t310 + t241;
t91 = t209 * t147 + t148 * t205;
t92 = -t147 * t205 + t148 * t209;
t141 = t188 + (-pkin(1) * t210 - pkin(2)) * t203;
t96 = pkin(3) * t147 + t141;
t35 = pkin(4) * t91 - pkin(9) * t92 + t96;
t245 = t204 * t35 + t208 * t25;
t244 = -t204 * t25 + t208 * t35;
t15 = -t205 * t43 + t209 * t40;
t242 = -t205 * t58 + t209 * t48;
t235 = t253 - t82;
t73 = t204 * t92 + t275;
t232 = t205 * t46 + t209 * t56 + t48 * t290 - t58 * t291;
t89 = -t117 * t204 - t208 * t273;
t227 = -t204 * t149 + t157 * t288 - t89;
t90 = -t117 * t208 + t204 * t273;
t226 = -t208 * t149 - t157 * t289 - t90;
t222 = g(3) * t310 - t253;
t221 = -g(3) * t312 - t251;
t220 = -qJ(3) * t293 + (qJD(3) - t111) * t210;
t13 = pkin(4) * t172 - t15;
t218 = -pkin(9) * t32 + (t13 + t15) * t75;
t2 = -t6 * qJD(5) - t204 * t3 + t208 * t8;
t214 = -t222 - t258;
t213 = -t115 * t32 + t4 * t157 + (pkin(9) * t273 - qJD(5) * t98 - t337) * t75 - t251;
t137 = t194 * t203 + t195 * t312;
t105 = t154 * t195 + t194 * t311;
t74 = t208 * t92 - t276;
t70 = t105 * t208 + t153 * t204;
t69 = -t105 * t204 + t153 * t208;
t51 = qJD(2) * t224 + t92 * qJD(4);
t50 = -qJD(2) * t223 - t91 * qJD(4);
t27 = -qJD(5) * t276 + t204 * t50 - t208 * t271 + t92 * t288;
t26 = -qJD(5) * t73 + t204 * t271 + t208 * t50;
t24 = pkin(4) * t310 - t242;
t20 = pkin(4) * t51 - pkin(9) * t50 + t110;
t10 = -pkin(4) * t271 - t333;
t9 = pkin(9) * t271 + t232;
t5 = [qJDD(1), g(1) * t207 - g(2) * t211, t249, (qJDD(1) * t198 + 0.2e1 * t206 * t268) * t196, (t206 * t282 - t299 * t285) * t281, (t260 * t206 + t261 * t292) * t201, (t260 * t210 - t261 * t293) * t201, t314, -t146 * t187 - t188 * t186 - t258 * t203 + (t210 * t314 - t230 * t281) * pkin(1) + t250, -t231 * pkin(1) * t281 - t300 * t186 - t236 * t187 - t217 * t203 + t252, -t146 * t130 + t141 * t99 + t82 * t147 + t250 * t202 + (-t249 * t200 + (qJD(1) * t83 + t61) * t293 + (-qJD(1) * t67 - qJDD(1) * t83 + t111 * t295 - t36) * t210) * t201, t141 * t100 + t146 * t131 + t82 * t148 - t250 * t200 + (-t249 * t202 + (-qJD(1) * t84 - t62) * t293 + (qJD(1) * t68 + qJDD(1) * t84 + t111 * t294 + t37) * t210) * t201, -t100 * t83 + t130 * t68 - t131 * t67 - t147 * t37 - t148 * t36 - t84 * t99 + (-t200 * t62 - t202 * t61) * t270 - t252, t37 * t84 + t62 * t68 + t36 * t83 + t61 * t67 + t82 * t141 + t111 * t146 - g(1) * (-pkin(1) * t207 - pkin(2) * t152 + pkin(7) * t309 - qJ(3) * t151) - g(2) * (pkin(1) * t211 + pkin(2) * t154 + pkin(7) * t311 + qJ(3) * t153), t238 * t50 + t33 * t92, -t238 * t51 - t33 * t91 - t34 * t92 - t50 * t78, t142 * t92 - t172 * t50 + (-t210 * t33 + t238 * t293) * t201, -t142 * t91 + t172 * t51 + (t210 * t34 - t78 * t293) * t201, (-t142 * t210 - t172 * t293) * t201, -t333 * t172 + t242 * t142 + t110 * t78 + t96 * t34 + t42 * t91 + t76 * t51 + g(1) * t102 - g(2) * t105 + (t15 * t293 - t210 * t215) * t201, t232 * t172 - t241 * t142 + t110 * t238 + t96 * t33 + t42 * t92 + t76 * t50 - g(1) * t317 - g(2) * t104 + (-t16 * t293 - t195 * t329 + t210 * t233) * t201, t11 * t74 - t239 * t26, -t11 * t73 - t12 * t74 + t239 * t27 - t26 * t53, t11 * t91 - t239 * t51 + t26 * t75 + t32 * t74, -t12 * t91 - t27 * t75 - t32 * t73 - t51 * t53, t32 * t91 + t51 * t75, (-qJD(5) * t245 + t20 * t208 - t204 * t9) * t75 + t244 * t32 + t2 * t91 - t246 * t51 + t10 * t53 + t24 * t12 + t4 * t73 + t13 * t27 + g(1) * t341 - g(2) * t70, -(qJD(5) * t244 + t20 * t204 + t208 * t9) * t75 - t245 * t32 - t1 * t91 - t6 * t51 - t10 * t239 + t24 * t11 + t4 * t74 + t13 * t26 - g(1) * t342 - g(2) * t69; 0, 0, 0, -t206 * t277, t299 * t313, (t287 * t296 + t283) * t201, -t287 * t273 + t185, t186, t145 * t187 + t313 * t331 + t214, pkin(1) * t277 + t144 * t187 + (pkin(7) * t285 + g(3)) * t312 + t251 + t274, -pkin(2) * t99 + t130 * t145 + t235 * t202 + ((-g(3) * t202 + t200 * t286) * t210 + (t200 * t220 - t206 * t61 + t210 * t87) * qJD(1)) * t201, -pkin(2) * t100 - t131 * t145 - t235 * t200 + ((g(3) * t200 + t202 * t286) * t210 + (t202 * t220 + t206 * t62 - t210 * t88) * qJD(1)) * t201, -t88 * t130 + t87 * t131 + (-qJ(3) * t99 + qJD(3) * t130 + t61 * t272 + t37) * t202 + (qJ(3) * t100 + qJD(3) * t131 + t62 * t272 - t36) * t200 + t221, -t111 * t145 - t61 * t87 - t62 * t88 + (-t200 * t61 + t202 * t62) * qJD(3) + (-t222 - t82) * pkin(2) + (-t200 * t36 + t202 * t37 + t221) * qJ(3), t157 * t33 + t238 * t302, -t156 * t33 - t157 * t34 - t238 * t301 - t302 * t78, t142 * t157 - t302 * t172 - t238 * t273, -t142 * t156 + t301 * t172 + t78 * t273, t172 * t273, -t109 * t78 + t142 * t237 - t15 * t273 + t42 * t156 + t338 * t172 + t193 * t34 - t222 * t195 + t301 * t76, -t109 * t238 - t115 * t142 + t42 * t157 + t16 * t273 + t337 * t172 + t193 * t33 + t222 * t194 + t302 * t76, t11 * t316 - t226 * t239, t53 * t90 - t239 * t89 - (t204 * t239 - t208 * t53) * t149 + (-t321 - t12 * t208 + (t204 * t53 + t208 * t239) * qJD(5)) * t157, t11 * t156 + t226 * t75 - t239 * t301 + t32 * t316, -t12 * t156 - t157 * t320 - t227 * t75 - t301 * t53, t156 * t32 + t301 * t75, -t237 * t12 + t2 * t156 + t322 * t53 - t301 * t246 + t332 * t208 + t213 * t204 - (t195 * t304 + t204 * t206) * t328 + t227 * t13, -t1 * t156 - t237 * t11 - t301 * t6 - t322 * t239 - t332 * t204 + t213 * t208 - (-t195 * t308 + t206 * t208) * t328 + t226 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162 + (t200 * t283 + (-t131 + t295) * t296) * t201, t315 + (t202 * t283 + (-t130 + t294) * t296) * t201, -t130 ^ 2 - t131 ^ 2, -t130 * t62 + t131 * t61 - t214 + t335, 0, 0, 0, 0, 0, t34 - t336, t33 + t343, 0, 0, 0, 0, 0, -t340 + t345, -t208 * t344 - t320 + t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238 * t78, t238 ^ 2 - t78 ^ 2, t33 - t343, -t34 - t336, t142, -t16 * t172 - t238 * t76 + t215 - t229, g(1) * t105 + g(2) * t102 + g(3) * t137 - t15 * t172 + t76 * t78 - t233, -t239 * t262 + t321, (t11 - t327) * t208 + (-t12 + t326) * t204, t262 * t75 + t320 + t339, t340 + t345, -t75 * t238, -pkin(4) * t12 - t16 * t53 + t218 * t204 - t208 * t334 + t238 * t246, -pkin(4) * t11 + t16 * t239 + t204 * t334 + t218 * t208 + t6 * t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t239 * t53, t239 ^ 2 - t53 ^ 2, t11 + t327, -t12 - t326, t32, t6 * t75 + t13 * t239 - g(1) * t69 + g(2) * t342 - g(3) * (-t137 * t204 - t275) + t2, -t246 * t75 + t13 * t53 + g(1) * t70 + g(2) * t341 - g(3) * (-t137 * t208 + t276) - t1;];
tau_reg = t5;
