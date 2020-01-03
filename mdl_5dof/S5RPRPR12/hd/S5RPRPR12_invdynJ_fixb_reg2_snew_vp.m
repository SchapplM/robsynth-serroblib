% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR12_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:17
% EndTime: 2019-12-31 18:30:29
% DurationCPUTime: 4.99s
% Computational Cost: add. (19524->405), mult. (48534->596), div. (0->0), fcn. (36226->10), ass. (0->244)
t216 = sin(pkin(8));
t211 = t216 ^ 2;
t218 = cos(pkin(8));
t212 = t218 ^ 2;
t255 = t211 + t212;
t283 = sin(qJ(1));
t284 = cos(qJ(1));
t229 = t284 * g(1) + t283 * g(2);
t304 = (2 * qJD(2) * qJD(1)) - t229;
t215 = sin(pkin(9));
t222 = cos(qJ(3));
t220 = sin(qJ(3));
t265 = t218 * t220;
t231 = t216 * t222 + t265;
t201 = t231 * qJD(1);
t217 = cos(pkin(9));
t184 = -t217 * qJD(3) + t215 * t201;
t186 = t215 * qJD(3) + t217 * t201;
t151 = t186 * t184;
t249 = t218 * qJDD(1);
t250 = t216 * qJDD(1);
t233 = t220 * t250 - t222 * t249;
t253 = t201 * qJD(3);
t172 = t233 + t253;
t293 = -t151 + t172;
t303 = t215 * t293;
t302 = t217 * t293;
t219 = sin(qJ(5));
t221 = cos(qJ(5));
t146 = t221 * t184 + t219 * t186;
t148 = -t219 * t184 + t221 * t186;
t114 = t148 * t146;
t165 = qJDD(5) + t172;
t295 = -t114 + t165;
t301 = t219 * t295;
t267 = t216 * t220;
t199 = (-t218 * t222 + t267) * qJD(1);
t177 = t201 * t199;
t291 = qJDD(3) - t177;
t300 = t220 * t291;
t299 = t221 * t295;
t298 = t222 * t291;
t223 = qJD(1) ^ 2;
t297 = -(t223 * pkin(1)) + qJDD(1) * qJ(2) + t304;
t234 = t283 * g(1) - t284 * g(2);
t230 = -qJDD(2) + t234;
t241 = pkin(2) * t218 + pkin(1);
t167 = t241 * qJDD(1) + (t255 * pkin(6) + qJ(2)) * t223 + t230;
t296 = pkin(6) + qJ(2);
t198 = t231 * qJDD(1);
t254 = t199 * qJD(3);
t174 = t198 - t254;
t157 = t215 * qJDD(3) + t217 * t174;
t237 = -t217 * qJDD(3) + t215 * t174;
t100 = -t146 * qJD(5) + t221 * t157 - t219 * t237;
t193 = qJD(5) + t199;
t138 = t193 * t146;
t294 = -t138 + t100;
t161 = t199 * t184;
t127 = -t157 - t161;
t292 = t157 - t161;
t256 = t223 * qJ(2);
t274 = qJDD(1) * pkin(1);
t194 = t230 + t256 + t274;
t290 = t255 * t256 - t194 - t274;
t238 = t219 * t157 + t221 * t237;
t82 = (qJD(5) - t193) * t148 + t238;
t144 = t146 ^ 2;
t145 = t148 ^ 2;
t287 = t184 ^ 2;
t183 = t186 ^ 2;
t192 = t193 ^ 2;
t286 = t199 ^ 2;
t197 = t201 ^ 2;
t285 = qJD(3) ^ 2;
t282 = t218 * g(3);
t164 = t199 * pkin(3) - t201 * qJ(4);
t227 = (-t296 * qJDD(1) + t241 * t223 - t304) * t216;
t226 = t227 - t282;
t235 = -t216 * g(3) + t297 * t218;
t160 = -t212 * t223 * pkin(2) + pkin(6) * t249 + t235;
t259 = t222 * t160;
t101 = -t285 * pkin(3) + qJDD(3) * qJ(4) - t199 * t164 + t220 * t226 + t259;
t109 = (-t174 + t254) * qJ(4) + (t172 + t253) * pkin(3) - t167;
t57 = 0.2e1 * qJD(4) * t186 + t215 * t101 - t217 * t109;
t44 = pkin(4) * t293 + t127 * pkin(7) - t57;
t156 = t199 * pkin(4) - t186 * pkin(7);
t58 = -0.2e1 * qJD(4) * t184 + t217 * t101 + t215 * t109;
t51 = -t287 * pkin(4) - t237 * pkin(7) - t199 * t156 + t58;
t21 = t219 * t51 - t221 * t44;
t22 = t219 * t44 + t221 * t51;
t12 = -t221 * t21 + t219 * t22;
t281 = t215 * t12;
t116 = t220 * t160 - t222 * t226;
t98 = -qJDD(3) * pkin(3) - t285 * qJ(4) + t201 * t164 + qJDD(4) + t116;
t280 = t215 * t98;
t117 = -g(3) * t265 + t220 * t227 + t259;
t87 = -t222 * t116 + t220 * t117;
t279 = t216 * t87;
t278 = t217 * t12;
t277 = t217 * t98;
t62 = t237 * pkin(4) - t287 * pkin(7) + t186 * t156 + t98;
t276 = t219 * t62;
t275 = t221 * t62;
t273 = t193 * t219;
t272 = t193 * t221;
t271 = t199 * t186;
t270 = t199 * t215;
t269 = t199 * t217;
t129 = t151 + t172;
t268 = t215 * t129;
t266 = t217 * t129;
t105 = t114 + t165;
t264 = t219 * t105;
t263 = t220 * t167;
t169 = qJDD(3) + t177;
t262 = t220 * t169;
t261 = t220 * t172;
t260 = t221 * t105;
t258 = t222 * t167;
t257 = t222 * t169;
t245 = t220 * t114;
t244 = t220 * t151;
t243 = t222 * t114;
t242 = t222 * t151;
t240 = -pkin(3) * t222 - pkin(2);
t13 = t219 * t21 + t221 * t22;
t36 = t215 * t57 + t217 * t58;
t88 = t220 * t116 + t222 * t117;
t236 = t216 * (t297 * t216 + t282) + t218 * t235;
t35 = t215 * t58 - t217 * t57;
t123 = t237 - t271;
t208 = t212 * qJDD(1);
t206 = t211 * qJDD(1);
t202 = t255 * t223;
t189 = -t197 - t285;
t188 = -t197 + t285;
t187 = t286 - t285;
t173 = t198 - 0.2e1 * t254;
t171 = t233 + 0.2e1 * t253;
t166 = -t286 - t285;
t163 = t222 * t172;
t159 = -t183 + t286;
t158 = -t286 + t287;
t152 = -t286 - t197;
t149 = -t183 + t287;
t143 = -t183 - t286;
t142 = -t220 * t189 - t257;
t141 = t222 * t189 - t262;
t140 = -t286 - t287;
t137 = t220 * t198 - t222 * t233;
t136 = -t222 * t198 - t220 * t233;
t135 = -t145 + t192;
t134 = t144 - t192;
t133 = t222 * t166 - t300;
t132 = t220 * t166 + t298;
t131 = -t183 - t287;
t122 = t237 + t271;
t121 = (-t184 * t217 + t186 * t215) * t199;
t120 = t217 * t157 - t186 * t270;
t119 = t184 * t269 + t215 * t237;
t118 = -t145 - t192;
t113 = t145 - t144;
t112 = -t192 - t144;
t111 = t217 * t158 - t268;
t110 = -t215 * t159 + t302;
t103 = -t215 * t143 - t266;
t102 = t217 * t143 - t268;
t99 = -t148 * qJD(5) - t238;
t97 = (-t146 * t221 + t148 * t219) * t193;
t96 = (-t146 * t219 - t148 * t221) * t193;
t94 = t217 * t140 - t303;
t93 = t215 * t140 + t302;
t92 = -t144 - t145;
t91 = -t123 * t217 - t215 * t127;
t90 = -t217 * t122 - t215 * t292;
t89 = -t123 * t215 + t217 * t127;
t85 = t138 + t100;
t81 = (qJD(5) + t193) * t148 + t238;
t80 = t221 * t134 - t264;
t79 = -t219 * t135 + t299;
t78 = t219 * t134 + t260;
t77 = t221 * t135 + t301;
t76 = t221 * t100 - t148 * t273;
t75 = t219 * t100 + t148 * t272;
t74 = t146 * t272 - t219 * t99;
t73 = t146 * t273 + t221 * t99;
t72 = t222 * t103 + t220 * t292;
t71 = t220 * t103 - t222 * t292;
t70 = -t219 * t118 - t260;
t69 = t221 * t118 - t264;
t68 = t220 * t122 + t222 * t94;
t67 = -t222 * t122 + t220 * t94;
t66 = t220 * t131 + t222 * t91;
t65 = -t222 * t131 + t220 * t91;
t64 = t221 * t112 - t301;
t63 = t219 * t112 + t299;
t61 = -qJ(4) * t102 + t277;
t60 = -t215 * t96 + t217 * t97;
t59 = -qJ(4) * t93 + t280;
t55 = t219 * t85 - t221 * t82;
t54 = -t219 * t294 - t221 * t81;
t53 = -t219 * t82 - t221 * t85;
t52 = -t219 * t81 + t221 * t294;
t50 = -t215 * t78 + t217 * t80;
t49 = -t215 * t77 + t217 * t79;
t48 = -t215 * t75 + t217 * t76;
t47 = -t215 * t73 + t217 * t74;
t46 = -pkin(3) * t102 + t58;
t45 = -pkin(3) * t93 + t57;
t42 = -t215 * t69 + t217 * t70;
t41 = t215 * t70 + t217 * t69;
t40 = -pkin(7) * t69 + t275;
t39 = -t215 * t63 + t217 * t64;
t38 = t215 * t64 + t217 * t63;
t37 = -pkin(7) * t63 + t276;
t34 = t220 * t294 + t222 * t42;
t33 = t220 * t42 - t222 * t294;
t32 = -pkin(4) * t294 + pkin(7) * t70 + t276;
t31 = -pkin(4) * t81 + pkin(7) * t64 - t275;
t30 = t220 * t81 + t222 * t39;
t29 = t220 * t39 - t222 * t81;
t26 = -qJ(4) * t89 - t35;
t25 = -t215 * t53 + t217 * t55;
t24 = -t215 * t52 + t217 * t54;
t23 = t215 * t55 + t217 * t53;
t19 = t220 * t92 + t222 * t25;
t18 = t220 * t25 - t222 * t92;
t17 = -pkin(3) * t23 - pkin(4) * t53;
t16 = -pkin(3) * t41 - pkin(4) * t69 + t22;
t15 = -qJ(4) * t41 - t215 * t32 + t217 * t40;
t14 = -pkin(3) * t38 - pkin(4) * t63 + t21;
t11 = -qJ(4) * t38 - t215 * t31 + t217 * t37;
t10 = -pkin(4) * t62 + pkin(7) * t13;
t9 = -pkin(7) * t53 - t12;
t8 = -pkin(4) * t92 + pkin(7) * t55 + t13;
t7 = t217 * t13 - t281;
t6 = t215 * t13 + t278;
t5 = t220 * t62 + t222 * t7;
t4 = t220 * t7 - t222 * t62;
t3 = -pkin(3) * t6 - pkin(4) * t12;
t2 = -qJ(4) * t23 - t215 * t8 + t217 * t9;
t1 = -pkin(7) * t278 - qJ(4) * t6 - t215 * t10;
t20 = [0, 0, 0, 0, 0, qJDD(1), t234, t229, 0, 0, t206, 0.2e1 * t216 * t249, 0, t208, 0, 0, -t290 * t218, t290 * t216, pkin(1) * t202 + qJ(2) * (t208 + t206) + t236, pkin(1) * t194 + qJ(2) * t236, t216 * (t222 * t174 - t220 * t253) + t218 * (t220 * t174 + t222 * t253), t216 * (-t222 * t171 - t220 * t173) + t218 * (-t220 * t171 + t222 * t173), t216 * (-t220 * t188 + t298) + t218 * (t222 * t188 + t300), t216 * (t222 * t254 + t261) + t218 * (t220 * t254 - t163), t216 * (t222 * t187 - t262) + t218 * (t220 * t187 + t257), (t216 * (-t199 * t222 + t201 * t220) + t218 * (-t199 * t220 - t201 * t222)) * qJD(3), t216 * (-pkin(6) * t132 - t263) + t218 * (-pkin(2) * t171 + pkin(6) * t133 + t258) - pkin(1) * t171 + qJ(2) * (-t216 * t132 + t218 * t133), t216 * (-pkin(6) * t141 - t258) + t218 * (-pkin(2) * t173 + pkin(6) * t142 - t263) - pkin(1) * t173 + qJ(2) * (-t216 * t141 + t218 * t142), t216 * (-pkin(6) * t136 - t87) + t218 * (-pkin(2) * t152 + pkin(6) * t137 + t88) - pkin(1) * t152 + qJ(2) * (-t216 * t136 + t218 * t137), -pkin(6) * t279 + t218 * (pkin(2) * t167 + pkin(6) * t88) + pkin(1) * t167 + qJ(2) * (t218 * t88 - t279), t216 * (t222 * t120 + t244) + t218 * (t220 * t120 - t242), t216 * (-t220 * t149 + t222 * t90) + t218 * (t222 * t149 + t220 * t90), t216 * (t222 * t110 - t220 * t127) + t218 * (t220 * t110 + t222 * t127), t216 * (t222 * t119 - t244) + t218 * (t220 * t119 + t242), t216 * (t222 * t111 - t220 * t123) + t218 * (t220 * t111 + t222 * t123), t216 * (t222 * t121 + t261) + t218 * (t220 * t121 - t163), t216 * (-pkin(6) * t67 - t220 * t45 + t222 * t59) + t218 * (-pkin(2) * t93 + pkin(6) * t68 + t220 * t59 + t222 * t45) - pkin(1) * t93 + qJ(2) * (-t216 * t67 + t218 * t68), t216 * (-pkin(6) * t71 - t220 * t46 + t222 * t61) + t218 * (-pkin(2) * t102 + pkin(6) * t72 + t220 * t61 + t222 * t46) - pkin(1) * t102 + qJ(2) * (-t216 * t71 + t218 * t72), t216 * (-pkin(6) * t65 + t222 * t26) + t218 * (pkin(6) * t66 + t220 * t26) + qJ(2) * (-t216 * t65 + t218 * t66) + (pkin(3) * t267 + t218 * t240 - pkin(1)) * t89, (t216 * (pkin(3) * t220 - qJ(4) * t222) + t218 * (-qJ(4) * t220 + t240) - pkin(1)) * t35 + t296 * (-t216 * (t220 * t36 - t222 * t98) + t218 * (t220 * t98 + t222 * t36)), t216 * (t222 * t48 + t245) + t218 * (t220 * t48 - t243), t216 * (t220 * t113 + t222 * t24) + t218 * (-t222 * t113 + t220 * t24), t216 * (t220 * t85 + t222 * t49) + t218 * (t220 * t49 - t222 * t85), t216 * (t222 * t47 - t245) + t218 * (t220 * t47 + t243), t216 * (-t220 * t82 + t222 * t50) + t218 * (t220 * t50 + t222 * t82), t216 * (t220 * t165 + t222 * t60) + t218 * (-t222 * t165 + t220 * t60), t216 * (-pkin(6) * t29 + t222 * t11 - t220 * t14) + t218 * (-pkin(2) * t38 + pkin(6) * t30 + t220 * t11 + t222 * t14) - pkin(1) * t38 + qJ(2) * (-t216 * t29 + t218 * t30), t216 * (-pkin(6) * t33 + t222 * t15 - t220 * t16) + t218 * (-pkin(2) * t41 + pkin(6) * t34 + t220 * t15 + t222 * t16) - pkin(1) * t41 + qJ(2) * (-t216 * t33 + t218 * t34), t216 * (-pkin(6) * t18 - t220 * t17 + t222 * t2) + t218 * (-pkin(2) * t23 + pkin(6) * t19 + t222 * t17 + t220 * t2) - pkin(1) * t23 + qJ(2) * (-t216 * t18 + t218 * t19), t216 * (-pkin(6) * t4 + t222 * t1 - t220 * t3) + t218 * (-pkin(2) * t6 + pkin(6) * t5 + t220 * t1 + t222 * t3) - pkin(1) * t6 + qJ(2) * (-t216 * t4 + t218 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t249, t250, -t202, -t194, 0, 0, 0, 0, 0, 0, t171, t173, t152, -t167, 0, 0, 0, 0, 0, 0, t93, t102, t89, t35, 0, 0, 0, 0, 0, 0, t38, t41, t23, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t197 - t286, t198, -t177, -t233, qJDD(3), -t116, -t117, 0, 0, t215 * t157 + t186 * t269, -t215 * t122 + t217 * t292, t217 * t159 + t303, t184 * t270 - t217 * t237, t215 * t158 + t266, (-t184 * t215 - t186 * t217) * t199, -pkin(3) * t122 + qJ(4) * t94 - t277, -pkin(3) * t292 + qJ(4) * t103 + t280, -pkin(3) * t131 + qJ(4) * t91 + t36, -pkin(3) * t98 + qJ(4) * t36, t215 * t76 + t217 * t75, t215 * t54 + t217 * t52, t215 * t79 + t217 * t77, t215 * t74 + t217 * t73, t215 * t80 + t217 * t78, t215 * t97 + t217 * t96, -pkin(3) * t81 + qJ(4) * t39 + t215 * t37 + t217 * t31, -pkin(3) * t294 + qJ(4) * t42 + t215 * t40 + t217 * t32, -pkin(3) * t92 + qJ(4) * t25 + t215 * t9 + t217 * t8, -pkin(3) * t62 - pkin(7) * t281 + qJ(4) * t7 + t217 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t292, t131, t98, 0, 0, 0, 0, 0, 0, t81, t294, t92, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t113, t85, -t114, -t82, t165, -t21, -t22, 0, 0;];
tauJ_reg = t20;
