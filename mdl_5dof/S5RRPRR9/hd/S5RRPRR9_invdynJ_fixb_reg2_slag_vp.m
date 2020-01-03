% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:22:04
% EndTime: 2019-12-31 20:22:14
% DurationCPUTime: 5.47s
% Computational Cost: add. (8853->519), mult. (20693->679), div. (0->0), fcn. (15173->14), ass. (0->245)
t203 = cos(qJ(2));
t297 = cos(pkin(9));
t245 = t297 * t203;
t173 = qJD(1) * t245;
t196 = sin(pkin(9));
t200 = sin(qJ(2));
t273 = qJD(1) * t200;
t145 = t196 * t273 - t173;
t136 = qJD(4) + t145;
t178 = pkin(2) * t196 + pkin(7);
t192 = qJ(2) + pkin(9);
t185 = sin(t192);
t186 = cos(t192);
t201 = sin(qJ(1));
t204 = cos(qJ(1));
t236 = g(1) * t204 + g(2) * t201;
t217 = -g(3) * t186 + t185 * t236;
t197 = -qJ(3) - pkin(6);
t249 = qJD(2) * t197;
t143 = qJD(3) * t203 + t200 * t249;
t168 = t197 * t203;
t109 = qJD(1) * t143 - qJDD(1) * t168;
t218 = -qJD(3) * t200 + t203 * t249;
t257 = t197 * t200;
t99 = qJDD(2) * pkin(2) + qJD(1) * t218 + qJDD(1) * t257;
t310 = t196 * t109 - t297 * t99;
t46 = -qJDD(2) * pkin(3) + t310;
t348 = qJD(4) * t178 * t136 - t217 + t46;
t246 = t297 * t200;
t158 = t196 * t203 + t246;
t148 = t158 * qJD(1);
t199 = sin(qJ(4));
t202 = cos(qJ(4));
t121 = -t202 * qJD(2) + t148 * t199;
t123 = qJD(2) * t199 + t148 * t202;
t198 = sin(qJ(5));
t327 = cos(qJ(5));
t227 = -t198 * t121 + t327 * t123;
t61 = t327 * t121 + t123 * t198;
t314 = t61 * t227;
t260 = t327 * t202;
t284 = t198 * t199;
t226 = t260 - t284;
t333 = qJD(4) + qJD(5);
t255 = t327 * qJD(5);
t335 = t327 * qJD(4) + t255;
t300 = -t226 * t145 - t335 * t202 + t333 * t284;
t261 = t327 * t199;
t161 = t198 * t202 + t261;
t114 = t333 * t161;
t299 = t161 * t145 + t114;
t313 = pkin(8) + t178;
t248 = qJD(4) * t313;
t282 = t199 * t145;
t163 = qJD(1) * t168;
t151 = t196 * t163;
t162 = qJD(1) * t257;
t111 = t162 * t297 + t151;
t81 = pkin(2) * t273 + pkin(3) * t148 + pkin(7) * t145;
t44 = t202 * t111 + t199 * t81;
t347 = pkin(8) * t282 + t199 * t248 + t44;
t43 = -t111 * t199 + t202 * t81;
t346 = -pkin(4) * t148 - t43 + (-pkin(8) * t145 - t248) * t202;
t269 = qJD(1) * qJD(2);
t254 = t200 * t269;
t212 = qJDD(1) * t158 - t196 * t254;
t108 = qJD(2) * t173 + t212;
t345 = qJD(2) * qJD(4) + t108;
t344 = -g(1) * t201 + g(2) * t204;
t272 = qJD(4) * t199;
t343 = t272 + t282;
t342 = t227 ^ 2 - t61 ^ 2;
t133 = qJD(5) + t136;
t271 = qJD(4) * t202;
t264 = t148 * t271 + t199 * t345;
t229 = qJDD(2) * t202 - t264;
t270 = qJD(5) * t198;
t56 = -t199 * qJDD(2) + t148 * t272 - t202 * t345;
t18 = t121 * t255 + t123 * t270 - t198 * t229 + t327 * t56;
t341 = t133 * t61 - t18;
t325 = pkin(2) * t203;
t184 = pkin(1) + t325;
t166 = -qJD(1) * t184 + qJD(3);
t70 = pkin(3) * t145 - pkin(7) * t148 + t166;
t309 = qJD(2) * pkin(2);
t154 = t162 + t309;
t247 = t297 * t163;
t106 = t196 * t154 - t247;
t92 = qJD(2) * pkin(7) + t106;
t37 = -t199 * t92 + t202 * t70;
t31 = -pkin(8) * t123 + t37;
t26 = pkin(4) * t136 + t31;
t38 = t199 * t70 + t202 * t92;
t32 = -pkin(8) * t121 + t38;
t147 = t158 * qJD(2);
t267 = t200 * qJDD(1);
t233 = -qJDD(1) * t245 + t196 * t267;
t107 = qJD(1) * t147 + t233;
t100 = qJDD(4) + t107;
t135 = pkin(2) * t254 - qJDD(1) * t184 + qJDD(3);
t39 = pkin(3) * t107 - pkin(7) * t108 + t135;
t36 = t202 * t39;
t51 = t297 * t109 + t196 * t99;
t47 = qJDD(2) * pkin(7) + t51;
t9 = -qJD(4) * t38 - t199 * t47 + t36;
t6 = pkin(4) * t100 + pkin(8) * t56 + t9;
t8 = t199 * t39 + t202 * t47 + t70 * t271 - t272 * t92;
t7 = pkin(8) * t229 + t8;
t1 = t198 * t6 + t26 * t255 - t32 * t270 + t327 * t7;
t195 = qJ(4) + qJ(5);
t190 = cos(t195);
t286 = t190 * t201;
t189 = sin(t195);
t287 = t189 * t204;
t126 = -t186 * t286 + t287;
t285 = t190 * t204;
t288 = t189 * t201;
t128 = t186 * t285 + t288;
t318 = g(3) * t185;
t105 = t154 * t297 + t151;
t91 = -qJD(2) * pkin(3) - t105;
t59 = t121 * pkin(4) + t91;
t340 = g(1) * t128 - g(2) * t126 + t190 * t318 + t59 * t61 - t1;
t234 = -t136 * t37 + t8;
t244 = t136 * t199;
t337 = t123 * t244;
t276 = t202 * t204;
t280 = t199 * t201;
t137 = t186 * t280 + t276;
t278 = t201 * t202;
t279 = t199 * t204;
t139 = -t186 * t279 + t278;
t334 = -g(1) * t139 + g(2) * t137;
t125 = t186 * t288 + t285;
t127 = -t186 * t287 + t286;
t263 = t327 * t32;
t11 = t198 * t26 + t263;
t2 = -qJD(5) * t11 - t198 * t7 + t327 * t6;
t332 = -g(1) * t127 + g(2) * t125 + t189 * t318 - t59 * t227 + t2;
t19 = qJD(5) * t227 - t198 * t56 - t327 * t229;
t331 = t133 * t227 - t19;
t95 = qJDD(5) + t100;
t330 = -t133 * t300 + t161 * t95;
t215 = -t236 * t186 - t318;
t329 = -t18 * t226 - t227 * t299;
t328 = t148 ^ 2;
t326 = pkin(2) * t200;
t171 = t204 * t184;
t320 = g(2) * t171;
t316 = g(3) * t203;
t315 = t202 * pkin(4);
t155 = t313 * t199;
t156 = t313 * t202;
t101 = -t327 * t155 - t198 * t156;
t312 = t101 * qJD(5) + t198 * t346 - t327 * t347;
t102 = -t198 * t155 + t327 * t156;
t311 = -t102 * qJD(5) + t198 * t347 + t327 * t346;
t307 = t136 * t38;
t306 = t148 * t61;
t305 = t148 * t227;
t302 = t198 * t32;
t301 = t56 * t199;
t223 = -t196 * t200 + t245;
t104 = -pkin(3) * t223 - pkin(7) * t158 - t184;
t120 = -t168 * t297 + t196 * t257;
t112 = t202 * t120;
t58 = t199 * t104 + t112;
t53 = t199 * t229;
t298 = -t121 * t271 + t53;
t296 = pkin(6) * qJDD(1);
t295 = t121 * t145;
t294 = t121 * t148;
t293 = t123 * t121;
t292 = t123 * t148;
t291 = t148 * t145;
t290 = t158 * t199;
t289 = t158 * t202;
t283 = t199 * t100;
t150 = t223 * qJD(2);
t281 = t199 * t150;
t86 = t202 * t100;
t277 = t202 * t150;
t193 = t200 ^ 2;
t194 = t203 ^ 2;
t275 = t193 - t194;
t274 = t193 + t194;
t266 = t203 * qJDD(1);
t265 = t200 * t309;
t207 = qJD(1) ^ 2;
t262 = t200 * t207 * t203;
t258 = t158 * t272;
t253 = pkin(4) * t199 - t197;
t80 = t297 * t143 + t196 * t218;
t82 = pkin(3) * t147 - pkin(7) * t150 + t265;
t250 = -t199 * t80 + t202 * t82;
t57 = t202 * t104 - t120 * t199;
t79 = t143 * t196 - t297 * t218;
t110 = t162 * t196 - t247;
t119 = -t168 * t196 - t197 * t246;
t243 = t136 * t202;
t242 = -t161 * t19 + t300 * t61;
t241 = -t133 * t299 + t226 * t95;
t240 = t203 * t254;
t239 = pkin(4) * t343 - t110;
t238 = t344 * t185;
t179 = -pkin(2) * t297 - pkin(3);
t237 = pkin(3) * t186 + pkin(7) * t185;
t232 = -t199 * t38 - t202 * t37;
t183 = pkin(3) + t315;
t205 = -pkin(8) - pkin(7);
t231 = t183 * t186 - t185 * t205;
t230 = -t136 * t343 + t86;
t33 = -pkin(4) * t223 - pkin(8) * t289 + t57;
t41 = -pkin(8) * t290 + t58;
t21 = -t198 * t41 + t327 * t33;
t22 = t198 * t33 + t327 * t41;
t228 = -0.2e1 * pkin(1) * t269 - pkin(6) * qJDD(2);
t225 = t158 * t271 + t281;
t224 = -t258 + t277;
t23 = t104 * t271 - t120 * t272 + t199 * t82 + t202 * t80;
t222 = -t178 * t100 + t136 * t91;
t220 = t229 * t202;
t206 = qJD(2) ^ 2;
t214 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t206 - t344;
t213 = pkin(1) * t207 + t236 - t296;
t167 = t179 - t315;
t144 = t145 ^ 2;
t140 = t186 * t276 + t280;
t138 = -t186 * t278 + t279;
t88 = t226 * t158;
t87 = t161 * t158;
t78 = pkin(4) * t290 + t119;
t45 = pkin(4) * t225 + t79;
t28 = t150 * t261 - t198 * t258 - t270 * t290 + (t150 * t198 + t335 * t158) * t202;
t27 = t114 * t158 - t150 * t260 + t198 * t281;
t25 = -pkin(4) * t229 + t46;
t24 = -qJD(4) * t58 + t250;
t20 = -pkin(8) * t225 + t23;
t17 = -pkin(8) * t277 + pkin(4) * t147 + (-t112 + (pkin(8) * t158 - t104) * t199) * qJD(4) + t250;
t13 = t327 * t31 - t302;
t12 = -t198 * t31 - t263;
t10 = t327 * t26 - t302;
t4 = -t22 * qJD(5) + t327 * t17 - t198 * t20;
t3 = t21 * qJD(5) + t198 * t17 + t327 * t20;
t5 = [0, 0, 0, 0, 0, qJDD(1), -t344, t236, 0, 0, qJDD(1) * t193 + 0.2e1 * t240, 0.2e1 * t200 * t266 - 0.2e1 * t269 * t275, qJDD(2) * t200 + t203 * t206, qJDD(1) * t194 - 0.2e1 * t240, qJDD(2) * t203 - t200 * t206, 0, t200 * t228 + t203 * t214, -t200 * t214 + t203 * t228, 0.2e1 * t274 * t296 - t236, -g(1) * (-pkin(1) * t201 + pkin(6) * t204) - g(2) * (pkin(1) * t204 + pkin(6) * t201) + (pkin(6) ^ 2 * t274 + pkin(1) ^ 2) * qJDD(1), t108 * t158 + t148 * t150, -t107 * t158 + t108 * t223 - t145 * t150 - t147 * t148, qJD(2) * t150 + qJDD(2) * t158, -t107 * t223 + t145 * t147, -qJD(2) * t147 + qJDD(2) * t223, 0, -qJDD(2) * t119 - t107 * t184 - t135 * t223 + t147 * t166 - t344 * t186 + (t145 * t326 - t79) * qJD(2), -qJDD(2) * t120 - t108 * t184 + t135 * t158 + t150 * t166 + (t148 * t326 - t80) * qJD(2) + t238, -t105 * t150 - t106 * t147 - t107 * t120 + t108 * t119 - t145 * t80 + t148 * t79 + t158 * t310 + t223 * t51 - t236, t51 * t120 + t106 * t80 + t310 * t119 - t105 * t79 - t135 * t184 + t166 * t265 - g(1) * (-t184 * t201 - t197 * t204) - g(2) * (-t197 * t201 + t171), t123 * t224 - t289 * t56, (-t121 * t202 - t123 * t199) * t150 + (t220 + t301 + (t121 * t199 - t123 * t202) * qJD(4)) * t158, t123 * t147 + t136 * t224 + t158 * t86 + t223 * t56, t121 * t225 - t158 * t53, -t121 * t147 - t136 * t225 - t158 * t283 - t223 * t229, -t100 * t223 + t136 * t147, t24 * t136 + t57 * t100 - t9 * t223 + t37 * t147 + t79 * t121 - t119 * t229 + t91 * t281 - g(1) * t138 - g(2) * t140 + (t199 * t46 + t271 * t91) * t158, t91 * t277 - g(1) * t137 - g(2) * t139 - t100 * t58 - t119 * t56 + t123 * t79 - t136 * t23 - t147 * t38 + t223 * t8 + (t46 * t202 - t272 * t91) * t158, -t23 * t121 + t58 * t229 - t24 * t123 + t57 * t56 + t232 * t150 + (-t8 * t199 - t9 * t202 + (t199 * t37 - t202 * t38) * qJD(4)) * t158 - t238, -t320 + t46 * t119 + t38 * t23 + t37 * t24 + t9 * t57 + t8 * t58 + t91 * t79 + (g(1) * t197 - g(2) * t237) * t204 + (-g(1) * (-t184 - t237) + g(2) * t197) * t201, -t18 * t88 - t227 * t27, t18 * t87 - t19 * t88 - t227 * t28 + t27 * t61, -t133 * t27 + t147 * t227 + t18 * t223 + t88 * t95, t19 * t87 + t28 * t61, -t133 * t28 - t147 * t61 + t19 * t223 - t87 * t95, t133 * t147 - t223 * t95, -g(1) * t126 - g(2) * t128 + t10 * t147 + t133 * t4 + t19 * t78 - t2 * t223 + t21 * t95 + t25 * t87 + t28 * t59 + t45 * t61, -g(1) * t125 - g(2) * t127 + t1 * t223 - t11 * t147 - t133 * t3 - t18 * t78 - t22 * t95 + t227 * t45 + t25 * t88 - t27 * t59, -t1 * t87 + t10 * t27 - t11 * t28 + t18 * t21 - t19 * t22 - t2 * t88 - t227 * t4 - t3 * t61 - t238, -t320 + t1 * t22 + t10 * t4 + t11 * t3 + t2 * t21 + t25 * t78 + t59 * t45 + (-g(1) * t253 - g(2) * t231) * t204 + (-g(1) * (-t184 - t231) - g(2) * t253) * t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, t275 * t207, t267, t262, t266, qJDD(2), t200 * t213 - t316, g(3) * t200 + t203 * t213, 0, 0, t291, -t144 + t328, (t173 + t145) * qJD(2) + t212, -t291, -t233, qJDD(2), t110 * qJD(2) - t166 * t148 + (qJDD(2) * t297 - t145 * t273) * pkin(2) + t217 - t310, qJD(2) * t111 + t145 * t166 + (-qJDD(2) * t196 - t148 * t273) * pkin(2) - t51 - t215, (t106 - t110) * t148 + (-t105 + t111) * t145 + (-t107 * t196 - t108 * t297) * pkin(2), t105 * t110 - t106 * t111 + (-t297 * t310 - t316 + t196 * t51 + (-qJD(1) * t166 + t236) * t200) * pkin(2), t123 * t243 - t301, (-t56 - t295) * t202 - t337 + t298, t136 * t243 + t283 - t292, t121 * t244 + t220, t230 + t294, -t136 * t148, t179 * t264 - t43 * t136 - t37 * t148 - t110 * t121 + t222 * t199 + (-t179 * qJDD(2) - t348) * t202, -t110 * t123 + t136 * t44 + t148 * t38 - t179 * t56 + t348 * t199 + t222 * t202, t44 * t121 + t43 * t123 + ((qJD(4) * t123 + t229) * t178 + t234) * t202 + (-t38 * t145 - t178 * t56 - t9 + (t121 * t178 - t38) * qJD(4)) * t199 + t215, t46 * t179 - t38 * t44 - t37 * t43 - t91 * t110 - g(3) * (t237 + t325) + (qJD(4) * t232 - t9 * t199 + t8 * t202) * t178 + t236 * (pkin(3) * t185 - pkin(7) * t186 + t326), -t161 * t18 - t227 * t300, t242 + t329, -t305 + t330, -t19 * t226 + t299 * t61, t241 + t306, -t133 * t148, -t10 * t148 + t101 * t95 + t311 * t133 + t167 * t19 + t217 * t190 - t226 * t25 + t239 * t61 + t299 * t59, -t102 * t95 + t11 * t148 - t133 * t312 + t161 * t25 - t167 * t18 - t189 * t217 + t227 * t239 - t300 * t59, t1 * t226 + t10 * t300 + t101 * t18 - t102 * t19 - t11 * t299 - t161 * t2 - t227 * t311 - t312 * t61 + t215, t1 * t102 + t2 * t101 + t25 * t167 - g(3) * (t231 + t325) + t239 * t59 + t312 * t11 + t311 * t10 + t236 * (t183 * t185 + t186 * t205 + t326); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t148 * qJD(2) + t233, (t173 - t145) * qJD(2) + t212, -t144 - t328, t105 * t148 + t106 * t145 + t135 + t344, 0, 0, 0, 0, 0, 0, t230 - t294, -t136 ^ 2 * t202 - t283 - t292, (t56 - t295) * t202 + t337 + t298, -t148 * t91 + (t9 + t307) * t202 + t234 * t199 + t344, 0, 0, 0, 0, 0, 0, t241 - t306, -t305 - t330, t242 - t329, t1 * t161 - t10 * t299 - t11 * t300 - t148 * t59 + t2 * t226 + t344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, -t121 ^ 2 + t123 ^ 2, t121 * t136 - t56, -t293, t123 * t136 + t229, t100, -t92 * t271 - t123 * t91 + t307 + t36 + (-qJD(4) * t70 + t318 - t47) * t199 + t334, g(1) * t140 - g(2) * t138 + t121 * t91 + t202 * t318 - t234, 0, 0, t314, t342, t341, -t314, t331, t95, -t12 * t133 + (-t123 * t61 - t133 * t270 + t327 * t95) * pkin(4) + t332, t13 * t133 + (-t123 * t227 - t133 * t255 - t198 * t95) * pkin(4) + t340, -t10 * t61 + t11 * t227 + t12 * t227 + t13 * t61 + (t327 * t18 - t19 * t198 + (t198 * t227 - t327 * t61) * qJD(5)) * pkin(4), -t10 * t12 - t11 * t13 + (t1 * t198 + t2 * t327 - t59 * t123 + t199 * t318 + (-t10 * t198 + t11 * t327) * qJD(5) + t334) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, t342, t341, -t314, t331, t95, t11 * t133 + t332, t10 * t133 + t340, 0, 0;];
tau_reg = t5;
