% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 01:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:00:24
% EndTime: 2019-05-05 01:00:36
% DurationCPUTime: 6.15s
% Computational Cost: add. (32905->468), mult. (74048->711), div. (0->0), fcn. (58032->14), ass. (0->275)
t243 = sin(qJ(6));
t238 = sin(pkin(12));
t240 = cos(pkin(12));
t245 = sin(qJ(4));
t249 = cos(qJ(4));
t267 = t238 * t249 + t240 * t245;
t216 = t267 * qJD(2);
t244 = sin(qJ(5));
t248 = cos(qJ(5));
t198 = -t248 * qJD(4) + t216 * t244;
t200 = qJD(4) * t244 + t216 * t248;
t247 = cos(qJ(6));
t168 = t247 * t198 + t200 * t243;
t170 = -t198 * t243 + t200 * t247;
t132 = t170 * t168;
t278 = t240 * qJDD(2);
t279 = t238 * qJDD(2);
t268 = t245 * t279 - t249 * t278;
t284 = t216 * qJD(4);
t191 = -t268 - t284;
t184 = qJDD(5) - t191;
t183 = qJDD(6) + t184;
t321 = -t132 + t183;
t330 = t243 * t321;
t174 = t200 * t198;
t319 = -t174 + t184;
t329 = t244 * t319;
t286 = qJD(2) * t240;
t289 = t238 * t245;
t214 = qJD(2) * t289 - t249 * t286;
t194 = t216 * t214;
t317 = qJDD(4) - t194;
t328 = t245 * t317;
t327 = t247 * t321;
t326 = t248 * t319;
t325 = t249 * t317;
t233 = t238 ^ 2;
t234 = t240 ^ 2;
t315 = qJD(2) ^ 2;
t221 = (t233 + t234) * t315;
t239 = sin(pkin(6));
t241 = cos(pkin(6));
t302 = sin(pkin(11));
t303 = cos(pkin(11));
t261 = t302 * g(1) - t303 * g(2);
t260 = t241 * t261;
t287 = -g(3) + qJDD(1);
t324 = t239 * t287 + t260;
t246 = sin(qJ(2));
t250 = cos(qJ(2));
t262 = t303 * g(1) + t302 * g(2);
t182 = t246 * t324 - t250 * t262;
t255 = qJDD(2) * qJ(3) + t182;
t258 = -t239 * t261 + t241 * t287;
t313 = 2 * qJD(3);
t148 = t240 * (-t315 * pkin(2) + t255) + t238 * t258 + t286 * t313;
t283 = t234 * t315;
t141 = -pkin(3) * t283 + pkin(8) * t278 + t148;
t256 = t240 * t258;
t323 = qJ(3) + pkin(8);
t253 = t256 + (-t323 * qJDD(2) + (-(2 * qJD(3)) + (t240 * pkin(3) + pkin(2)) * qJD(2)) * qJD(2) - t182) * t238;
t108 = t249 * t141 + t245 * t253;
t213 = t267 * qJDD(2);
t285 = t214 * qJD(4);
t193 = t213 - t285;
t266 = -t244 * qJDD(4) - t248 * t193;
t160 = -qJD(5) * t198 - t266;
t269 = -t248 * qJDD(4) + t244 * t193;
t264 = qJD(5) * t200 + t269;
t105 = -t168 * qJD(6) + t247 * t160 - t243 * t264;
t210 = qJD(5) + t214;
t206 = qJD(6) + t210;
t156 = t206 * t168;
t320 = -t156 + t105;
t180 = t210 * t198;
t137 = t160 + t180;
t270 = t243 * t160 + t247 * t264;
t83 = (qJD(6) - t206) * t170 + t270;
t133 = (qJD(5) - t210) * t200 + t269;
t166 = t168 ^ 2;
t167 = t170 ^ 2;
t316 = t198 ^ 2;
t197 = t200 ^ 2;
t205 = t206 ^ 2;
t209 = t210 ^ 2;
t211 = t214 ^ 2;
t212 = t216 ^ 2;
t314 = qJD(4) ^ 2;
t181 = t246 * t262 + t250 * t324;
t237 = qJDD(2) * pkin(2);
t176 = -t315 * qJ(3) + qJDD(3) - t181 - t237;
t162 = -pkin(3) * t278 + t176 + (-t233 * t315 - t283) * pkin(8);
t111 = (-t193 + t285) * pkin(9) + (-t191 + t284) * pkin(4) + t162;
t185 = pkin(4) * t214 - pkin(9) * t216;
t96 = -t314 * pkin(4) + qJDD(4) * pkin(9) - t214 * t185 + t108;
t66 = -t248 * t111 + t244 * t96;
t47 = pkin(5) * t319 - pkin(10) * t137 - t66;
t177 = pkin(5) * t210 - pkin(10) * t200;
t67 = t244 * t111 + t248 * t96;
t49 = -t316 * pkin(5) - pkin(10) * t264 - t210 * t177 + t67;
t26 = t243 * t49 - t247 * t47;
t27 = t243 * t47 + t247 * t49;
t13 = t243 * t27 - t247 * t26;
t312 = pkin(5) * t13;
t86 = t156 + t105;
t56 = -t243 * t83 - t247 * t86;
t311 = pkin(5) * t56;
t310 = t13 * t244;
t309 = t13 * t248;
t107 = t245 * t141 - t249 * t253;
t71 = -t107 * t249 + t108 * t245;
t308 = t238 * t71;
t95 = -qJDD(4) * pkin(4) - t314 * pkin(9) + t216 * t185 + t107;
t68 = pkin(5) * t264 - t316 * pkin(10) + t200 * t177 + t95;
t307 = t243 * t68;
t306 = t244 * t95;
t305 = t247 * t68;
t304 = t248 * t95;
t123 = t132 + t183;
t301 = t123 * t243;
t300 = t123 * t247;
t144 = t174 + t184;
t299 = t144 * t244;
t298 = t144 * t248;
t297 = t162 * t245;
t296 = t162 * t249;
t188 = qJDD(4) + t194;
t295 = t188 * t245;
t294 = t188 * t249;
t293 = t206 * t243;
t292 = t206 * t247;
t291 = t210 * t244;
t290 = t210 * t248;
t281 = qJD(5) + t210;
t277 = t245 * t132;
t276 = t249 * t132;
t275 = t245 * t174;
t274 = t249 * t174;
t273 = -pkin(4) * t249 - pkin(3);
t14 = t243 * t26 + t247 * t27;
t39 = t244 * t66 + t248 * t67;
t271 = -t176 + t237;
t72 = t107 * t245 + t249 * t108;
t147 = -t256 + ((-pkin(2) * qJD(2) + t313) * qJD(2) + t255) * t238;
t113 = t147 * t238 + t240 * t148;
t31 = t245 * t39 - t249 * t95;
t32 = t245 * t95 + t249 * t39;
t16 = -t238 * t31 + t240 * t32;
t38 = t244 * t67 - t248 * t66;
t127 = -t205 - t166;
t76 = t127 * t243 + t327;
t265 = pkin(5) * t76 - t26;
t140 = -t167 - t205;
t93 = t140 * t247 - t301;
t263 = pkin(5) * t93 - t27;
t230 = t234 * qJDD(2);
t229 = t233 * qJDD(2);
t220 = t230 + t229;
t219 = t240 * t221;
t218 = t238 * t221;
t204 = -t212 - t314;
t203 = -t212 + t314;
t202 = t211 - t314;
t192 = t213 - 0.2e1 * t285;
t190 = t268 + 0.2e1 * t284;
t186 = -t314 - t211;
t179 = -t197 + t209;
t178 = -t209 + t316;
t173 = -t211 - t212;
t172 = t197 - t316;
t165 = -t197 - t209;
t164 = -t204 * t245 - t294;
t163 = t204 * t249 - t295;
t161 = -t209 - t316;
t155 = t197 + t316;
t154 = t213 * t245 - t249 * t268;
t153 = -t213 * t249 - t245 * t268;
t152 = -t167 + t205;
t151 = t166 - t205;
t150 = t186 * t249 - t328;
t149 = t186 * t245 + t325;
t142 = (-t198 * t248 + t200 * t244) * t210;
t138 = t281 * t198 + t266;
t136 = t160 - t180;
t135 = -t281 * t200 - t269;
t131 = t167 - t166;
t130 = t160 * t248 - t200 * t291;
t129 = t198 * t290 + t244 * t264;
t128 = -t163 * t238 + t164 * t240;
t126 = t178 * t248 - t299;
t125 = -t179 * t244 + t326;
t121 = -t165 * t244 - t298;
t120 = t165 * t248 - t299;
t119 = (-t168 * t247 + t170 * t243) * t206;
t118 = (-t168 * t243 - t170 * t247) * t206;
t117 = -t153 * t238 + t154 * t240;
t116 = t161 * t248 - t329;
t115 = t161 * t244 + t326;
t114 = -t149 * t238 + t150 * t240;
t112 = -t166 - t167;
t104 = -qJD(6) * t170 - t270;
t103 = t135 * t248 - t136 * t244;
t102 = -t133 * t248 + t137 * t244;
t101 = -t133 * t244 - t137 * t248;
t100 = t151 * t247 - t301;
t99 = -t152 * t243 + t327;
t98 = t151 * t243 + t300;
t97 = t152 * t247 + t330;
t94 = -t140 * t243 - t300;
t91 = t121 * t249 - t138 * t245;
t90 = t121 * t245 + t138 * t249;
t89 = t116 * t249 - t135 * t245;
t88 = t116 * t245 + t135 * t249;
t82 = (qJD(6) + t206) * t170 + t270;
t81 = t105 * t247 - t170 * t293;
t80 = t105 * t243 + t170 * t292;
t79 = -t104 * t243 + t168 * t292;
t78 = t104 * t247 + t168 * t293;
t77 = t127 * t247 - t330;
t75 = t102 * t249 - t155 * t245;
t74 = t102 * t245 + t155 * t249;
t73 = -t118 * t244 + t119 * t248;
t70 = -pkin(9) * t120 + t304;
t69 = -pkin(9) * t115 + t306;
t64 = t100 * t248 - t244 * t98;
t63 = -t244 * t97 + t248 * t99;
t62 = -t244 * t93 + t248 * t94;
t61 = t244 * t94 + t248 * t93;
t60 = -t238 * t90 + t240 * t91;
t59 = -t243 * t320 - t247 * t82;
t58 = t243 * t86 - t247 * t83;
t57 = -t243 * t82 + t247 * t320;
t55 = -t238 * t88 + t240 * t89;
t54 = -t244 * t80 + t248 * t81;
t53 = -t244 * t78 + t248 * t79;
t52 = -pkin(4) * t120 + t67;
t51 = -t244 * t76 + t248 * t77;
t50 = t244 * t77 + t248 * t76;
t48 = -pkin(4) * t115 + t66;
t46 = -t238 * t74 + t240 * t75;
t44 = -pkin(10) * t93 + t305;
t43 = -pkin(10) * t76 + t307;
t42 = t240 * t72 - t308;
t41 = t245 * t320 + t249 * t62;
t40 = t245 * t62 - t249 * t320;
t37 = t245 * t82 + t249 * t51;
t36 = t245 * t51 - t249 * t82;
t35 = -pkin(5) * t320 + pkin(10) * t94 + t307;
t34 = -pkin(5) * t82 + pkin(10) * t77 - t305;
t33 = -pkin(9) * t101 - t38;
t30 = -t244 * t57 + t248 * t59;
t29 = -t244 * t56 + t248 * t58;
t28 = t244 * t58 + t248 * t56;
t24 = t112 * t245 + t249 * t29;
t23 = -t112 * t249 + t245 * t29;
t22 = -t238 * t40 + t240 * t41;
t21 = -pkin(4) * t28 - t311;
t20 = -t238 * t36 + t240 * t37;
t19 = -pkin(4) * t61 - t263;
t18 = -pkin(4) * t50 - t265;
t17 = -pkin(9) * t61 - t244 * t35 + t248 * t44;
t15 = -pkin(9) * t50 - t244 * t34 + t248 * t43;
t12 = -t23 * t238 + t24 * t240;
t11 = -pkin(5) * t68 + pkin(10) * t14;
t10 = -pkin(10) * t56 - t13;
t9 = -pkin(5) * t112 + pkin(10) * t58 + t14;
t8 = t14 * t248 - t310;
t7 = t14 * t244 + t309;
t6 = t245 * t68 + t249 * t8;
t5 = t245 * t8 - t249 * t68;
t4 = -pkin(9) * t28 + t10 * t248 - t244 * t9;
t3 = -pkin(4) * t7 - t312;
t2 = -pkin(9) * t7 - pkin(10) * t309 - t11 * t244;
t1 = -t238 * t5 + t240 * t6;
t25 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t287, 0, 0, 0, 0, 0, 0, (qJDD(2) * t250 - t315 * t246) * t239, (-qJDD(2) * t246 - t315 * t250) * t239, 0, t241 ^ 2 * t287 + (t250 * t181 + t246 * t182 - t260) * t239, 0, 0, 0, 0, 0, 0, (-t219 * t246 + t250 * t278) * t239, (t218 * t246 - t250 * t279) * t239, (t220 * t246 + t221 * t250) * t239, t241 * (-t147 * t240 + t148 * t238) + (t113 * t246 - t176 * t250) * t239, 0, 0, 0, 0, 0, 0, t241 * (t149 * t240 + t150 * t238) + (t114 * t246 - t190 * t250) * t239, t241 * (t163 * t240 + t164 * t238) + (t128 * t246 - t192 * t250) * t239, t241 * (t153 * t240 + t154 * t238) + (t117 * t246 - t173 * t250) * t239, t241 * (t238 * t72 + t240 * t71) + (-t162 * t250 + t246 * t42) * t239, 0, 0, 0, 0, 0, 0, t241 * (t238 * t89 + t240 * t88) + (-t115 * t250 + t246 * t55) * t239, t241 * (t238 * t91 + t240 * t90) + (-t120 * t250 + t246 * t60) * t239, t241 * (t238 * t75 + t240 * t74) + (-t101 * t250 + t246 * t46) * t239, t241 * (t238 * t32 + t240 * t31) + (t16 * t246 - t250 * t38) * t239, 0, 0, 0, 0, 0, 0, t241 * (t238 * t37 + t240 * t36) + (t20 * t246 - t250 * t50) * t239, t241 * (t238 * t41 + t240 * t40) + (t22 * t246 - t250 * t61) * t239, t241 * (t23 * t240 + t238 * t24) + (t12 * t246 - t250 * t28) * t239, t241 * (t238 * t6 + t240 * t5) + (t1 * t246 - t250 * t7) * t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t181, -t182, 0, 0, t229, 0.2e1 * t238 * t278, 0, t230, 0, 0, -qJ(3) * t219 + t240 * t271, qJ(3) * t218 - t238 * t271, pkin(2) * t221 + qJ(3) * t220 + t113, -pkin(2) * t176 + qJ(3) * t113, t238 * (t193 * t249 - t245 * t284) + t240 * (t193 * t245 + t249 * t284), t238 * (-t190 * t249 - t192 * t245) + t240 * (-t190 * t245 + t192 * t249), t238 * (-t203 * t245 + t325) + t240 * (t203 * t249 + t328), t238 * (-t191 * t245 + t249 * t285) + t240 * (t191 * t249 + t245 * t285), t238 * (t202 * t249 - t295) + t240 * (t202 * t245 + t294), (t238 * (-t214 * t249 + t216 * t245) + t240 * (-t214 * t245 - t216 * t249)) * qJD(4), t238 * (-pkin(8) * t149 + t297) + t240 * (-pkin(3) * t190 + pkin(8) * t150 - t296) - pkin(2) * t190 + qJ(3) * t114, t238 * (-pkin(8) * t163 + t296) + t240 * (-pkin(3) * t192 + pkin(8) * t164 + t297) - pkin(2) * t192 + qJ(3) * t128, t238 * (-pkin(8) * t153 - t71) + t240 * (-pkin(3) * t173 + pkin(8) * t154 + t72) - pkin(2) * t173 + qJ(3) * t117, -pkin(8) * t308 + t240 * (-pkin(3) * t162 + pkin(8) * t72) - pkin(2) * t162 + qJ(3) * t42, t238 * (t130 * t249 + t275) + t240 * (t130 * t245 - t274), t238 * (t103 * t249 + t172 * t245) + t240 * (t103 * t245 - t172 * t249), t238 * (t125 * t249 + t137 * t245) + t240 * (t125 * t245 - t137 * t249), t238 * (t129 * t249 - t275) + t240 * (t129 * t245 + t274), t238 * (t126 * t249 - t133 * t245) + t240 * (t126 * t245 + t133 * t249), t238 * (t142 * t249 + t184 * t245) + t240 * (t142 * t245 - t184 * t249), t238 * (-pkin(8) * t88 - t245 * t48 + t249 * t69) + t240 * (-pkin(3) * t115 + pkin(8) * t89 + t245 * t69 + t249 * t48) - pkin(2) * t115 + qJ(3) * t55, t238 * (-pkin(8) * t90 - t245 * t52 + t249 * t70) + t240 * (-pkin(3) * t120 + pkin(8) * t91 + t245 * t70 + t249 * t52) - pkin(2) * t120 + qJ(3) * t60, t238 * (-pkin(8) * t74 + t249 * t33) + t240 * (pkin(8) * t75 + t245 * t33) + qJ(3) * t46 + (pkin(4) * t289 + t240 * t273 - pkin(2)) * t101, (t238 * (pkin(4) * t245 - pkin(9) * t249) + t240 * (-pkin(9) * t245 + t273) - pkin(2)) * t38 + t323 * t16, t238 * (t249 * t54 + t277) + t240 * (t245 * t54 - t276), t238 * (t131 * t245 + t249 * t30) + t240 * (-t131 * t249 + t245 * t30), t238 * (t245 * t86 + t249 * t63) + t240 * (t245 * t63 - t249 * t86), t238 * (t249 * t53 - t277) + t240 * (t245 * t53 + t276), t238 * (-t245 * t83 + t249 * t64) + t240 * (t245 * t64 + t249 * t83), t238 * (t183 * t245 + t249 * t73) + t240 * (-t183 * t249 + t245 * t73), t238 * (-pkin(8) * t36 + t15 * t249 - t18 * t245) + t240 * (-pkin(3) * t50 + pkin(8) * t37 + t15 * t245 + t18 * t249) - pkin(2) * t50 + qJ(3) * t20, t238 * (-pkin(8) * t40 + t17 * t249 - t19 * t245) + t240 * (-pkin(3) * t61 + pkin(8) * t41 + t17 * t245 + t19 * t249) - pkin(2) * t61 + qJ(3) * t22, t238 * (-pkin(8) * t23 - t21 * t245 + t249 * t4) + t240 * (-pkin(3) * t28 + pkin(8) * t24 + t21 * t249 + t245 * t4) - pkin(2) * t28 + qJ(3) * t12, t238 * (-pkin(8) * t5 + t2 * t249 - t245 * t3) + t240 * (-pkin(3) * t7 + pkin(8) * t6 + t2 * t245 + t249 * t3) - pkin(2) * t7 + qJ(3) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, t279, -t221, t176, 0, 0, 0, 0, 0, 0, t190, t192, t173, t162, 0, 0, 0, 0, 0, 0, t115, t120, t101, t38, 0, 0, 0, 0, 0, 0, t50, t61, t28, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, t212 - t211, t213, -t194, -t268, qJDD(4), -t107, -t108, 0, 0, t160 * t244 + t200 * t290, t135 * t244 + t136 * t248, t179 * t248 + t329, t198 * t291 - t248 * t264, t178 * t244 + t298, (-t198 * t244 - t200 * t248) * t210, pkin(4) * t135 + pkin(9) * t116 - t304, pkin(4) * t138 + pkin(9) * t121 + t306, pkin(4) * t155 + pkin(9) * t102 + t39, -pkin(4) * t95 + pkin(9) * t39, t244 * t81 + t248 * t80, t244 * t59 + t248 * t57, t244 * t99 + t248 * t97, t244 * t79 + t248 * t78, t100 * t244 + t248 * t98, t118 * t248 + t119 * t244, -pkin(4) * t82 + pkin(9) * t51 + t244 * t43 + t248 * t34, -pkin(4) * t320 + pkin(9) * t62 + t244 * t44 + t248 * t35, -pkin(4) * t112 + pkin(9) * t29 + t10 * t244 + t248 * t9, -pkin(4) * t68 + pkin(9) * t8 - pkin(10) * t310 + t11 * t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t172, t137, -t174, -t133, t184, -t66, -t67, 0, 0, t132, t131, t86, -t132, -t83, t183, t265, t263, t311, t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t131, t86, -t132, -t83, t183, -t26, -t27, 0, 0;];
tauJ_reg  = t25;
