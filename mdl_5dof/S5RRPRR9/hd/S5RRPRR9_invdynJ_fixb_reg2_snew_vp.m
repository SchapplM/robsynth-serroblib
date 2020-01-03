% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:22:01
% EndTime: 2019-12-31 20:22:13
% DurationCPUTime: 5.14s
% Computational Cost: add. (25012->445), mult. (57388->633), div. (0->0), fcn. (40876->10), ass. (0->261)
t225 = sin(pkin(9));
t226 = cos(pkin(9));
t234 = cos(qJ(2));
t271 = qJD(1) * t234;
t230 = sin(qJ(2));
t272 = qJD(1) * t230;
t199 = t225 * t272 - t226 * t271;
t201 = t225 * t271 + t226 * t272;
t175 = t201 * t199;
t303 = qJDD(2) - t175;
t312 = t225 * t303;
t311 = t226 * t303;
t228 = sin(qJ(5));
t229 = sin(qJ(4));
t233 = cos(qJ(4));
t182 = -t233 * qJD(2) + t229 * t201;
t184 = t229 * qJD(2) + t233 * t201;
t232 = cos(qJ(5));
t152 = t232 * t182 + t228 * t184;
t154 = -t228 * t182 + t232 * t184;
t119 = t154 * t152;
t215 = t230 * qJDD(1);
t263 = qJD(1) * qJD(2);
t256 = t234 * t263;
t206 = t215 + t256;
t217 = t234 * qJDD(1);
t257 = t230 * t263;
t207 = t217 - t257;
t249 = t225 * t206 - t226 * t207;
t174 = qJDD(4) + t249;
t173 = qJDD(5) + t174;
t305 = -t119 + t173;
t310 = t228 * t305;
t159 = t184 * t182;
t304 = -t159 + t174;
t309 = t229 * t304;
t308 = t232 * t305;
t307 = t233 * t304;
t270 = qJD(2) * t201;
t160 = t249 + t270;
t196 = qJD(4) + t199;
t191 = qJD(5) + t196;
t138 = t191 * t152;
t177 = t226 * t206 + t225 * t207;
t247 = -t229 * qJDD(2) - t233 * t177;
t143 = -t182 * qJD(4) - t247;
t251 = -t233 * qJDD(2) + t229 * t177;
t242 = t184 * qJD(4) + t251;
t93 = -t152 * qJD(5) + t232 * t143 - t228 * t242;
t306 = -t138 + t93;
t167 = t196 * t182;
t124 = t143 + t167;
t194 = qJD(2) * t199;
t162 = t177 - t194;
t236 = qJD(1) ^ 2;
t231 = sin(qJ(1));
t298 = cos(qJ(1));
t245 = t298 * g(1) + t231 * g(2);
t289 = qJDD(1) * pkin(6);
t239 = -t236 * pkin(1) - t245 + t289;
t186 = -t230 * g(3) + t234 * t239;
t223 = t234 ^ 2;
t220 = t223 * t236;
t243 = qJD(2) * pkin(2) - qJ(3) * t272;
t155 = -pkin(2) * t220 + t207 * qJ(3) - qJD(2) * t243 + t186;
t238 = t230 * t239;
t277 = t230 * t236;
t237 = -t238 - t206 * qJ(3) + qJDD(2) * pkin(2) + (pkin(2) * t277 + qJ(3) * t263 - g(3)) * t234;
t112 = -0.2e1 * qJD(3) * t199 + t226 * t155 + t225 * t237;
t253 = t228 * t143 + t232 * t242;
t76 = (qJD(5) - t191) * t154 + t253;
t120 = (qJD(4) - t196) * t184 + t251;
t254 = t231 * g(1) - t298 * g(2);
t244 = qJDD(1) * pkin(1) + t254;
t157 = t207 * pkin(2) + (qJ(3) * t223 + pkin(6)) * t236 - t243 * t272 - qJDD(3) + t244;
t150 = t152 ^ 2;
t151 = t154 ^ 2;
t302 = t182 ^ 2;
t181 = t184 ^ 2;
t190 = t191 ^ 2;
t195 = t196 ^ 2;
t197 = t199 ^ 2;
t198 = t201 ^ 2;
t301 = 0.2e1 * qJD(3);
t104 = pkin(3) * t160 - t162 * pkin(7) - t157;
t168 = t199 * pkin(3) - t201 * pkin(7);
t235 = qJD(2) ^ 2;
t96 = -t235 * pkin(3) + qJDD(2) * pkin(7) - t199 * t168 + t112;
t59 = -t233 * t104 + t229 * t96;
t48 = pkin(4) * t304 - t124 * pkin(8) - t59;
t164 = t196 * pkin(4) - t184 * pkin(8);
t60 = t229 * t104 + t233 * t96;
t50 = -t302 * pkin(4) - pkin(8) * t242 - t196 * t164 + t60;
t24 = t228 * t50 - t232 * t48;
t25 = t228 * t48 + t232 * t50;
t12 = t228 * t25 - t232 * t24;
t300 = pkin(4) * t12;
t79 = t138 + t93;
t44 = -t228 * t76 - t232 * t79;
t299 = pkin(4) * t44;
t297 = pkin(3) * t225;
t252 = t225 * t155 - t226 * t237;
t95 = -qJDD(2) * pkin(3) - t235 * pkin(7) + (t301 + t168) * t201 + t252;
t56 = pkin(4) * t242 - t302 * pkin(8) + t184 * t164 + t95;
t296 = t228 * t56;
t295 = t229 * t12;
t294 = t229 * t95;
t111 = t201 * t301 + t252;
t71 = -t226 * t111 + t225 * t112;
t293 = t230 * t71;
t292 = t232 * t56;
t291 = t233 * t12;
t290 = t233 * t95;
t288 = t191 * t228;
t287 = t191 * t232;
t286 = t196 * t229;
t285 = t196 * t233;
t284 = t225 * t157;
t171 = qJDD(2) + t175;
t283 = t225 * t171;
t282 = t226 * t157;
t281 = t226 * t171;
t108 = t119 + t173;
t280 = t228 * t108;
t131 = t159 + t174;
t279 = t229 * t131;
t214 = t234 * t277;
t278 = t230 * (qJDD(2) + t214);
t276 = t232 * t108;
t275 = t233 * t131;
t274 = t234 * (qJDD(2) - t214);
t269 = qJD(2) * t225;
t268 = qJD(2) * t226;
t265 = qJD(4) + t196;
t262 = t225 * t119;
t261 = t225 * t159;
t260 = t226 * t119;
t259 = t226 * t159;
t258 = -pkin(3) * t226 - pkin(2);
t13 = t228 * t24 + t232 * t25;
t36 = t229 * t59 + t233 * t60;
t72 = t225 * t111 + t226 * t112;
t185 = t234 * g(3) + t238;
t250 = t230 * t185 + t234 * t186;
t35 = t229 * t60 - t233 * t59;
t115 = -t190 - t150;
t65 = t228 * t115 + t308;
t246 = pkin(4) * t65 - t24;
t161 = -t249 + t270;
t128 = -t151 - t190;
t83 = t232 * t128 - t280;
t241 = pkin(4) * t83 - t25;
t222 = t230 ^ 2;
t218 = t222 * t236;
t208 = t217 - 0.2e1 * t257;
t205 = t215 + 0.2e1 * t256;
t203 = t236 * pkin(6) + t244;
t189 = -t198 - t235;
t188 = -t198 + t235;
t187 = t197 - t235;
t169 = -t235 - t197;
t166 = -t181 + t195;
t165 = -t195 + t302;
t163 = t177 + t194;
t158 = -t197 - t198;
t156 = t181 - t302;
t146 = -t181 - t195;
t145 = -t225 * t189 - t281;
t144 = t226 * t189 - t283;
t141 = -t195 - t302;
t137 = t181 + t302;
t136 = -t151 + t190;
t135 = t150 - t190;
t134 = t226 * t169 - t312;
t133 = t225 * t169 + t311;
t129 = (-t182 * t233 + t184 * t229) * t196;
t127 = t226 * t161 + t225 * t163;
t126 = t225 * t161 - t226 * t163;
t125 = t265 * t182 + t247;
t123 = t143 - t167;
t121 = -t265 * t184 - t251;
t118 = t151 - t150;
t117 = t233 * t143 - t184 * t286;
t116 = t182 * t285 + t229 * t242;
t114 = t233 * t165 - t279;
t113 = -t229 * t166 + t307;
t106 = -t229 * t146 - t275;
t105 = t233 * t146 - t279;
t101 = (-t152 * t232 + t154 * t228) * t191;
t100 = (-t152 * t228 - t154 * t232) * t191;
t99 = t233 * t141 - t309;
t98 = t229 * t141 + t307;
t97 = -t150 - t151;
t92 = -t154 * qJD(5) - t253;
t91 = -t120 * t233 + t229 * t124;
t90 = t233 * t121 - t229 * t123;
t89 = -t120 * t229 - t233 * t124;
t88 = t232 * t135 - t280;
t87 = -t228 * t136 + t308;
t86 = t228 * t135 + t276;
t85 = t232 * t136 + t310;
t84 = -t228 * t128 - t276;
t82 = t226 * t106 - t225 * t125;
t81 = t225 * t106 + t226 * t125;
t75 = (qJD(5) + t191) * t154 + t253;
t74 = -t225 * t121 + t226 * t99;
t73 = t226 * t121 + t225 * t99;
t70 = -t154 * t288 + t232 * t93;
t69 = t154 * t287 + t228 * t93;
t68 = t152 * t287 - t228 * t92;
t67 = t152 * t288 + t232 * t92;
t66 = t232 * t115 - t310;
t64 = -t225 * t137 + t226 * t91;
t63 = t226 * t137 + t225 * t91;
t62 = -t229 * t100 + t233 * t101;
t61 = -pkin(7) * t105 + t290;
t57 = -pkin(7) * t98 + t294;
t55 = -t229 * t86 + t233 * t88;
t54 = -t229 * t85 + t233 * t87;
t53 = -pkin(3) * t105 + t60;
t52 = -t229 * t83 + t233 * t84;
t51 = t229 * t84 + t233 * t83;
t49 = -pkin(3) * t98 + t59;
t46 = t228 * t79 - t232 * t76;
t45 = -t228 * t306 - t232 * t75;
t43 = -t228 * t75 + t232 * t306;
t42 = -t229 * t69 + t233 * t70;
t41 = -t229 * t67 + t233 * t68;
t40 = -t229 * t65 + t233 * t66;
t39 = t229 * t66 + t233 * t65;
t38 = -pkin(8) * t83 + t292;
t37 = -pkin(8) * t65 + t296;
t34 = t225 * t306 + t226 * t52;
t33 = t225 * t52 - t226 * t306;
t32 = t225 * t75 + t226 * t40;
t31 = t225 * t40 - t226 * t75;
t30 = -pkin(4) * t306 + pkin(8) * t84 + t296;
t28 = t225 * t36 - t226 * t95;
t27 = -pkin(4) * t75 + pkin(8) * t66 - t292;
t26 = -pkin(7) * t89 - t35;
t22 = -t229 * t44 + t233 * t46;
t21 = -t229 * t43 + t233 * t45;
t20 = t229 * t46 + t233 * t44;
t19 = t226 * t22 + t225 * t97;
t18 = t225 * t22 - t226 * t97;
t17 = -pkin(3) * t20 - t299;
t16 = -pkin(3) * t51 - t241;
t15 = -pkin(3) * t39 - t246;
t14 = -pkin(7) * t51 - t229 * t30 + t233 * t38;
t11 = -pkin(7) * t39 - t229 * t27 + t233 * t37;
t10 = -pkin(4) * t56 + pkin(8) * t13;
t9 = -pkin(8) * t44 - t12;
t8 = -pkin(4) * t97 + pkin(8) * t46 + t13;
t7 = t233 * t13 - t295;
t6 = t229 * t13 + t291;
t5 = t225 * t56 + t226 * t7;
t4 = t225 * t7 - t226 * t56;
t3 = -pkin(3) * t6 - t300;
t2 = -pkin(7) * t20 - t229 * t8 + t233 * t9;
t1 = -pkin(7) * t6 - pkin(8) * t291 - t229 * t10;
t23 = [0, 0, 0, 0, 0, qJDD(1), t254, t245, 0, 0, (t206 + t256) * t230, t234 * t205 + t230 * t208, t278 + t234 * (-t218 + t235), (t207 - t257) * t234, t230 * (t220 - t235) + t274, 0, t234 * t203 + pkin(1) * t208 + pkin(6) * (t234 * (-t220 - t235) - t278), -t230 * t203 - pkin(1) * t205 + pkin(6) * (-t274 - t230 * (-t218 - t235)), pkin(1) * (t218 + t220) + (t222 + t223) * t289 + t250, pkin(1) * t203 + pkin(6) * t250, t230 * (t226 * t177 - t201 * t269) + t234 * (t225 * t177 + t201 * t268), t230 * (-t226 * t160 - t225 * t162) + t234 * (-t225 * t160 + t226 * t162), t230 * (-t225 * t188 + t311) + t234 * (t226 * t188 + t312), t230 * (t199 * t268 + t225 * t249) + t234 * (t199 * t269 - t226 * t249), t230 * (t226 * t187 - t283) + t234 * (t225 * t187 + t281), (t230 * (-t199 * t226 + t201 * t225) + t234 * (-t199 * t225 - t201 * t226)) * qJD(2), t230 * (-qJ(3) * t133 - t284) + t234 * (-pkin(2) * t160 + qJ(3) * t134 + t282) - pkin(1) * t160 + pkin(6) * (-t230 * t133 + t234 * t134), t230 * (-qJ(3) * t144 - t282) + t234 * (-pkin(2) * t162 + qJ(3) * t145 - t284) - pkin(1) * t162 + pkin(6) * (-t230 * t144 + t234 * t145), t230 * (-qJ(3) * t126 - t71) + t234 * (-pkin(2) * t158 + qJ(3) * t127 + t72) - pkin(1) * t158 + pkin(6) * (-t230 * t126 + t234 * t127), -qJ(3) * t293 + t234 * (pkin(2) * t157 + qJ(3) * t72) + pkin(1) * t157 + pkin(6) * (t234 * t72 - t293), t230 * (t226 * t117 + t261) + t234 * (t225 * t117 - t259), t230 * (t225 * t156 + t226 * t90) + t234 * (-t226 * t156 + t225 * t90), t230 * (t226 * t113 + t225 * t124) + t234 * (t225 * t113 - t226 * t124), t230 * (t226 * t116 - t261) + t234 * (t225 * t116 + t259), t230 * (t226 * t114 - t225 * t120) + t234 * (t225 * t114 + t226 * t120), t230 * (t226 * t129 + t225 * t174) + t234 * (t225 * t129 - t226 * t174), t230 * (-qJ(3) * t73 - t225 * t49 + t226 * t57) + t234 * (-pkin(2) * t98 + qJ(3) * t74 + t225 * t57 + t226 * t49) - pkin(1) * t98 + pkin(6) * (-t230 * t73 + t234 * t74), t230 * (-qJ(3) * t81 - t225 * t53 + t226 * t61) + t234 * (-pkin(2) * t105 + qJ(3) * t82 + t225 * t61 + t226 * t53) - pkin(1) * t105 + pkin(6) * (-t230 * t81 + t234 * t82), t230 * (-qJ(3) * t63 + t226 * t26) + t234 * (qJ(3) * t64 + t225 * t26) + pkin(6) * (-t230 * t63 + t234 * t64) + (t230 * t297 + t234 * t258 - pkin(1)) * t89, (t230 * (-pkin(7) * t226 + t297) + t234 * (-pkin(7) * t225 + t258) - pkin(1)) * t35 + (pkin(6) + qJ(3)) * (-t230 * t28 + t234 * (t225 * t95 + t226 * t36)), t230 * (t226 * t42 + t262) + t234 * (t225 * t42 - t260), t230 * (t225 * t118 + t226 * t21) + t234 * (-t226 * t118 + t225 * t21), t230 * (t225 * t79 + t226 * t54) + t234 * (t225 * t54 - t226 * t79), t230 * (t226 * t41 - t262) + t234 * (t225 * t41 + t260), t230 * (-t225 * t76 + t226 * t55) + t234 * (t225 * t55 + t226 * t76), t230 * (t225 * t173 + t226 * t62) + t234 * (-t226 * t173 + t225 * t62), t230 * (-qJ(3) * t31 + t226 * t11 - t225 * t15) + t234 * (-pkin(2) * t39 + qJ(3) * t32 + t225 * t11 + t226 * t15) - pkin(1) * t39 + pkin(6) * (-t230 * t31 + t234 * t32), t230 * (-qJ(3) * t33 + t226 * t14 - t225 * t16) + t234 * (-pkin(2) * t51 + qJ(3) * t34 + t225 * t14 + t226 * t16) - pkin(1) * t51 + pkin(6) * (-t230 * t33 + t234 * t34), t230 * (-qJ(3) * t18 - t225 * t17 + t226 * t2) + t234 * (-pkin(2) * t20 + qJ(3) * t19 + t226 * t17 + t225 * t2) - pkin(1) * t20 + pkin(6) * (-t230 * t18 + t234 * t19), t230 * (-qJ(3) * t4 + t226 * t1 - t225 * t3) + t234 * (-pkin(2) * t6 + qJ(3) * t5 + t225 * t1 + t226 * t3) - pkin(1) * t6 + pkin(6) * (-t230 * t4 + t234 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214, t218 - t220, t215, t214, t217, qJDD(2), -t185, -t186, 0, 0, t175, t198 - t197, t163, -t175, t161, qJDD(2), pkin(2) * t133 - t111, pkin(2) * t144 - t112, pkin(2) * t126, pkin(2) * t71, t229 * t143 + t184 * t285, t229 * t121 + t233 * t123, t233 * t166 + t309, t182 * t286 - t233 * t242, t229 * t165 + t275, (-t182 * t229 - t184 * t233) * t196, pkin(2) * t73 + pkin(3) * t121 + pkin(7) * t99 - t290, pkin(2) * t81 + pkin(3) * t125 + pkin(7) * t106 + t294, pkin(2) * t63 + pkin(3) * t137 + pkin(7) * t91 + t36, pkin(2) * t28 - pkin(3) * t95 + pkin(7) * t36, t229 * t70 + t233 * t69, t229 * t45 + t233 * t43, t229 * t87 + t233 * t85, t229 * t68 + t233 * t67, t229 * t88 + t233 * t86, t233 * t100 + t229 * t101, pkin(2) * t31 - pkin(3) * t75 + pkin(7) * t40 + t229 * t37 + t233 * t27, pkin(2) * t33 - pkin(3) * t306 + pkin(7) * t52 + t229 * t38 + t233 * t30, pkin(2) * t18 - pkin(3) * t97 + pkin(7) * t22 + t229 * t9 + t233 * t8, pkin(2) * t4 - pkin(3) * t56 + pkin(7) * t7 - pkin(8) * t295 + t233 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t162, t158, -t157, 0, 0, 0, 0, 0, 0, t98, t105, t89, t35, 0, 0, 0, 0, 0, 0, t39, t51, t20, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t156, t124, -t159, -t120, t174, -t59, -t60, 0, 0, t119, t118, t79, -t119, -t76, t173, t246, t241, t299, t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t118, t79, -t119, -t76, t173, -t24, -t25, 0, 0;];
tauJ_reg = t23;
