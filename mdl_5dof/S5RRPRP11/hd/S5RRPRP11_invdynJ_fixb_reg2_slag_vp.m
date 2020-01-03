% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:10
% EndTime: 2019-12-31 20:14:16
% DurationCPUTime: 3.30s
% Computational Cost: add. (3017->425), mult. (6301->497), div. (0->0), fcn. (3628->6), ass. (0->230)
t144 = sin(qJ(2));
t227 = qJD(1) * qJD(2);
t213 = t144 * t227;
t147 = cos(qJ(2));
t225 = t147 * qJDD(1);
t305 = -t213 + t225;
t280 = pkin(3) + pkin(6);
t146 = cos(qJ(4));
t233 = qJD(4) * t146;
t143 = sin(qJ(4));
t237 = qJD(1) * t147;
t214 = t143 * t237;
t33 = -qJD(4) * t214 + t143 * qJDD(2) + (qJD(2) * qJD(4) + t305) * t146;
t27 = t33 * t143;
t231 = t143 * qJD(2);
t80 = t146 * t237 + t231;
t304 = t80 * t233 + t27;
t145 = sin(qJ(1));
t148 = cos(qJ(1));
t201 = g(1) * t148 + g(2) * t145;
t274 = g(3) * t144;
t303 = -t201 * t147 - t274;
t236 = qJD(2) * t144;
t228 = t146 * qJD(2);
t82 = -t214 + t228;
t268 = t146 * t82;
t32 = qJD(4) * t80 - t146 * qJDD(2) + t143 * t305;
t269 = t146 * t32;
t270 = t143 * t80;
t302 = (qJD(4) * (t143 * t82 + t146 * t80) + t27 + t269) * t147 - (-t268 + t270) * t236;
t230 = t144 * qJD(1);
t111 = qJD(4) + t230;
t184 = t82 * t111;
t185 = t80 * t111;
t301 = (t32 + t185) * t143 - (t33 + t184) * t146;
t212 = t147 * t227;
t226 = t144 * qJDD(1);
t172 = t212 + t226;
t79 = qJDD(4) + t172;
t279 = t79 * pkin(4);
t110 = pkin(2) * t213;
t259 = qJ(3) * t147;
t196 = pkin(7) * t144 - t259;
t229 = t144 * qJD(3);
t160 = qJD(2) * t196 - t229;
t127 = t144 * qJ(3);
t210 = -pkin(1) - t127;
t281 = pkin(2) + pkin(7);
t171 = -t147 * t281 + t210;
t20 = qJD(1) * t160 + qJDD(1) * t171 + t110;
t234 = qJD(4) * t143;
t109 = pkin(6) * t212;
t118 = pkin(6) * t226;
t211 = qJDD(3) + t109 + t118;
t37 = pkin(3) * t172 - qJDD(2) * t281 + t211;
t52 = t171 * qJD(1);
t121 = pkin(6) * t230;
t289 = qJD(3) + t121;
t244 = pkin(3) * t230 + t289;
t57 = -qJD(2) * t281 + t244;
t4 = -t143 * t20 + t146 * t37 - t233 * t52 - t234 * t57;
t2 = qJDD(5) - t4 - t279;
t22 = t143 * t57 + t146 * t52;
t16 = qJ(5) * t111 + t22;
t264 = t16 * t111;
t300 = -t2 + t264;
t262 = t22 * t111;
t299 = t4 + t262;
t154 = qJD(4) * t111 * t281 + t303;
t119 = pkin(6) * t225;
t138 = qJDD(2) * qJ(3);
t139 = qJD(2) * qJD(3);
t220 = t119 + t138 + t139;
t87 = t280 * t236;
t38 = pkin(3) * t225 - qJD(1) * t87 + t220;
t5 = t33 * pkin(4) + t32 * qJ(5) - t82 * qJD(5) + t38;
t298 = t5 + t154;
t64 = t143 * t79;
t297 = -t111 * t233 - t64;
t295 = t33 - t184;
t254 = t111 * t143;
t258 = qJD(2) * t80;
t65 = t146 * t79;
t293 = -t111 * t254 - t258 + t65;
t292 = -t111 * t234 + t65;
t266 = t281 * t79;
t140 = qJD(2) * qJ(3);
t122 = pkin(6) * t237;
t88 = pkin(3) * t237 + t122;
t66 = t140 + t88;
t288 = t111 * t66 - t266;
t137 = g(3) * t147;
t224 = t143 * t37 + t146 * t20 + t233 * t57;
t246 = t148 * t143;
t249 = t145 * t146;
t69 = t144 * t246 + t249;
t245 = t148 * t146;
t250 = t145 * t143;
t71 = -t144 * t250 + t245;
t287 = -g(1) * t69 + g(2) * t71 + (-qJD(4) * t52 + t137) * t143 + t224;
t248 = t145 * t147;
t106 = qJ(3) * t248;
t247 = t147 * t148;
t108 = qJ(3) * t247;
t133 = t147 * pkin(2);
t242 = t133 + t127;
t219 = pkin(7) * t147 + t242;
t286 = t144 * t201 * t281 - g(1) * t108 - g(2) * t106 - g(3) * t219;
t74 = -pkin(1) - t219;
t98 = t280 * t144;
t271 = t143 * t98 + t146 * t74;
t124 = pkin(2) * t236;
t50 = t124 + t160;
t235 = qJD(2) * t147;
t89 = t280 * t235;
t12 = -qJD(4) * t271 - t143 * t50 + t146 * t89;
t284 = t144 * (t111 * t228 - t33) + t147 * (-t258 - t292);
t283 = t82 ^ 2;
t282 = t111 ^ 2;
t278 = g(1) * t145;
t275 = g(2) * t148;
t273 = t82 * t80;
t125 = pkin(2) * t230;
t60 = qJD(1) * t196 + t125;
t31 = t143 * t88 + t146 * t60;
t198 = pkin(4) * t146 + qJ(5) * t143;
t181 = -pkin(3) - t198;
t272 = qJD(4) * t198 - t146 * qJD(5) - t181 * t230 + t289;
t267 = t281 * t32;
t265 = t281 * t82;
t21 = -t143 * t52 + t146 * t57;
t263 = t21 * t111;
t261 = t79 * qJ(5);
t260 = pkin(6) * qJDD(2);
t257 = qJD(2) * t82;
t255 = qJDD(2) * pkin(2);
t253 = t144 * t145;
t252 = t144 * t148;
t151 = qJD(1) ^ 2;
t251 = t144 * t151;
t243 = qJD(5) - t21;
t99 = t280 * t147;
t134 = t148 * pkin(6);
t241 = pkin(3) * t148 + t134;
t240 = pkin(1) * t148 + pkin(6) * t145;
t141 = t144 ^ 2;
t142 = t147 ^ 2;
t239 = t141 - t142;
t238 = t141 + t142;
t232 = qJD(4) * t147;
t217 = t146 * t230;
t223 = t217 * t80 + t304;
t222 = t111 * t217 - t297;
t221 = t80 ^ 2 - t283;
t218 = -g(1) * t252 - g(2) * t253 + t137;
t215 = t111 * t237;
t209 = -qJD(2) * pkin(2) + qJD(3);
t208 = pkin(2) * t247 + qJ(3) * t252 + t240;
t207 = -t118 - t218;
t206 = t144 * t212;
t68 = -t144 * t245 + t250;
t70 = t144 * t249 + t246;
t205 = -g(1) * t70 - g(2) * t68;
t204 = -g(1) * t71 - g(2) * t69;
t203 = t238 * qJDD(1) * pkin(6);
t150 = qJD(2) ^ 2;
t202 = pkin(6) * t150 + t275;
t199 = t281 * t304 - t218;
t197 = pkin(4) * t143 - qJ(5) * t146;
t14 = -pkin(4) * t111 + t243;
t194 = t14 * t143 + t146 * t16;
t193 = -t143 * t21 + t146 * t22;
t30 = -t143 * t60 + t146 * t88;
t39 = -t143 * t74 + t146 * t98;
t91 = t121 + t209;
t97 = -t122 - t140;
t187 = t144 * t97 + t147 * t91;
t182 = t210 - t133;
t178 = pkin(3) * t145 + pkin(7) * t247 + t208;
t67 = t182 * qJD(1);
t176 = t230 * t67 + qJDD(3) - t207;
t175 = -0.2e1 * pkin(1) * t227 - t260;
t174 = t237 * t80 - t222;
t3 = -t234 * t52 + t224;
t11 = t143 * t89 + t146 * t50 + t233 * t98 - t234 * t74;
t173 = -qJ(3) * t235 - t229;
t169 = 0.2e1 * qJDD(1) * pkin(1) - t202;
t93 = -pkin(1) - t242;
t168 = t260 + (-qJD(1) * t93 - t67) * qJD(2);
t23 = pkin(4) * t80 - qJ(5) * t82 + t66;
t167 = t111 * t23 - t266;
t164 = g(1) * t68 - g(2) * t70 + t137 * t146 + t4;
t163 = t171 * t278;
t34 = qJD(1) * t173 + qJDD(1) * t182 + t110;
t63 = t124 + t173;
t159 = qJD(1) * t63 + qJDD(1) * t93 + t202 + t34;
t158 = t143 * t184 - t223 + t269;
t53 = pkin(6) * t213 - t220;
t59 = t211 - t255;
t157 = qJD(2) * t187 + t59 * t144 - t53 * t147;
t156 = t23 * t82 + qJDD(5) - t164;
t155 = -t232 * t270 + (t147 * t33 - t236 * t80) * t146;
t153 = t38 + t154;
t116 = g(1) * t248;
t105 = t147 * t251;
t96 = t239 * t151;
t95 = qJDD(2) * t147 - t144 * t150;
t94 = qJDD(2) * t144 + t147 * t150;
t92 = qJ(3) + t197;
t85 = -qJ(3) * t237 + t125;
t76 = qJDD(1) * t142 - 0.2e1 * t206;
t75 = qJDD(1) * t141 + 0.2e1 * t206;
t49 = 0.2e1 * t144 * t225 - 0.2e1 * t227 * t239;
t48 = t147 * t198 + t99;
t42 = t111 * t235 + t144 * t79;
t41 = pkin(4) * t82 + qJ(5) * t80;
t36 = -pkin(4) * t144 - t39;
t35 = qJ(5) * t144 + t271;
t25 = -pkin(4) * t237 - t30;
t24 = qJ(5) * t237 + t31;
t17 = (-qJD(4) * t197 + qJD(5) * t143) * t147 + (-pkin(6) + t181) * t236;
t15 = t185 - t32;
t13 = (-t144 * t254 - t147 * t82) * qJD(1) + t292;
t10 = -t254 * t82 - t269;
t9 = -pkin(4) * t235 - t12;
t8 = -t232 * t268 + (t147 * t32 + t236 * t82) * t143;
t7 = qJ(5) * t235 + qJD(5) * t144 + t11;
t6 = (t111 * t231 - t32) * t144 + (t257 + t297) * t147;
t1 = t111 * qJD(5) + t261 + t3;
t18 = [0, 0, 0, 0, 0, qJDD(1), -t275 + t278, t201, 0, 0, t75, t49, t94, t76, t95, 0, t144 * t175 + t147 * t169 + t116, t175 * t147 + (-t169 - t278) * t144, -t201 + 0.2e1 * t203, -g(1) * (-t145 * pkin(1) + t134) - g(2) * t240 + (pkin(6) ^ 2 * t238 + pkin(1) ^ 2) * qJDD(1), 0, -t94, -t95, t75, t49, t76, t203 + t157 - t201, t144 * t168 + t147 * t159 - t116, t168 * t147 + (-t159 + t278) * t144, pkin(6) * t157 - g(1) * t134 - g(2) * t208 - t182 * t278 + t34 * t93 + t67 * t63, t8, t302, t6, t155, t284, t42, t12 * t111 + t99 * t33 + t39 * t79 - t87 * t80 + (-t228 * t66 + t4) * t144 + (qJD(2) * t21 + t38 * t146 - t234 * t66) * t147 + t204, -t11 * t111 - t99 * t32 - t271 * t79 - t87 * t82 + (t231 * t66 - t3) * t144 + (-qJD(2) * t22 - t38 * t143 - t233 * t66) * t147 - t205, -t11 * t80 - t12 * t82 + t39 * t32 - t271 * t33 + t116 + t193 * t236 + (-t275 + t143 * t4 - t146 * t3 + (t143 * t22 + t146 * t21) * qJD(4)) * t147, -g(1) * t241 - g(2) * t178 + t22 * t11 + t21 * t12 + t271 * t3 + t38 * t99 + t4 * t39 - t66 * t87 - t163, t8, t6, -t302, t42, -t284, t155, -t9 * t111 + t17 * t80 + t48 * t33 - t36 * t79 + (-t228 * t23 - t2) * t144 + (-qJD(2) * t14 + t5 * t146 - t23 * t234) * t147 + t204, -t36 * t32 - t35 * t33 - t7 * t80 + t9 * t82 + t116 + t194 * t236 + (-t275 - t1 * t146 - t143 * t2 + (-t14 * t146 + t143 * t16) * qJD(4)) * t147, t7 * t111 - t17 * t82 + t48 * t32 + t35 * t79 + (-t23 * t231 + t1) * t144 + (qJD(2) * t16 + t5 * t143 + t23 * t233) * t147 + t205, t1 * t35 + t16 * t7 + t5 * t48 + t23 * t17 + t2 * t36 + t14 * t9 - g(1) * (t71 * pkin(4) + t70 * qJ(5) + t241) - g(2) * (pkin(4) * t69 + qJ(5) * t68 + t178) - t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t96, t226, t105, t225, qJDD(2), pkin(1) * t251 + t207, t274 - t119 + (pkin(1) * t151 + t201) * t147, 0, 0, qJDD(2), -t226, -t225, -t105, t96, t105, (-pkin(2) * t144 + t259) * qJDD(1) + ((-t97 - t140) * t144 + (t209 - t91) * t147) * qJD(1), -t237 * t85 + t176 - 0.2e1 * t255, t119 + 0.2e1 * t138 + 0.2e1 * t139 + (qJD(1) * t85 - g(3)) * t144 + (qJD(1) * t67 - t201) * t147, -t53 * qJ(3) - t97 * qJD(3) - t59 * pkin(2) - t67 * t85 - g(1) * (-pkin(2) * t252 + t108) - g(2) * (-pkin(2) * t253 + t106) - g(3) * t242 - t187 * qJD(1) * pkin(6), t10, t301, t13, t223, t174, -t215, qJ(3) * t33 - t30 * t111 + t153 * t143 + t146 * t288 - t21 * t237 + t244 * t80, -qJ(3) * t32 + t31 * t111 - t143 * t288 + t153 * t146 + t22 * t237 + t244 * t82, t30 * t82 + t31 * t80 + (-t267 - t299) * t146 + (t21 * t230 - t3 + (t21 - t265) * qJD(4)) * t143 + t199, t38 * qJ(3) - t22 * t31 - t21 * t30 + t244 * t66 - (qJD(4) * t193 + t3 * t143 + t4 * t146) * t281 + t286, t10, t13, -t301, -t215, -t174, t146 * t185 + t27, t25 * t111 + t14 * t237 + t143 * t298 + t146 * t167 + t272 * t80 + t92 * t33, t24 * t80 - t25 * t82 + (-t267 - t300) * t146 + (-t14 * t230 - t1 + (-t14 - t265) * qJD(4)) * t143 + t199, -t24 * t111 + t143 * t167 - t146 * t298 - t16 * t237 - t272 * t82 + t92 * t32, t5 * t92 - t16 * t24 - t14 * t25 + t272 * t23 - (qJD(4) * t194 + t1 * t143 - t2 * t146) * t281 + t303 * t197 + t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, qJDD(2) + t105, -t141 * t151 - t150, qJD(2) * t97 + t109 + t176 - t255, 0, 0, 0, 0, 0, 0, t293, -t146 * t282 - t257 - t64, t158, -t66 * qJD(2) + t299 * t146 + (t3 - t263) * t143 + t218, 0, 0, 0, 0, 0, 0, t293, t158, t222 + t257, -t23 * qJD(2) + t300 * t146 + (t111 * t14 + t1) * t143 + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, -t221, t15, -t273, -t295, t79, -t66 * t82 + t164 + t262, t66 * t80 + t263 - t287, 0, 0, t273, t15, t221, t79, t295, -t273, -t41 * t80 - t156 + t262 + 0.2e1 * t279, pkin(4) * t32 - t33 * qJ(5) + (t16 - t22) * t82 + (t14 - t243) * t80, 0.2e1 * t261 - t23 * t80 + t41 * t82 + (0.2e1 * qJD(5) - t21) * t111 + t287, t1 * qJ(5) - t2 * pkin(4) - t23 * t41 - t14 * t22 - g(1) * (-pkin(4) * t68 + qJ(5) * t69) - g(2) * (pkin(4) * t70 - qJ(5) * t71) + t243 * t16 + t198 * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 + t273, t15, -t282 - t283, t156 - t264 - t279;];
tau_reg = t18;
