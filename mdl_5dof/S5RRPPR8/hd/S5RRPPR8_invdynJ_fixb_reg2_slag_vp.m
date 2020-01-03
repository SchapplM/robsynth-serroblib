% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:20
% EndTime: 2019-12-31 19:39:25
% DurationCPUTime: 2.97s
% Computational Cost: add. (3292->409), mult. (7252->512), div. (0->0), fcn. (4922->10), ass. (0->214)
t195 = cos(pkin(8));
t202 = cos(qJ(2));
t194 = sin(pkin(8));
t199 = sin(qJ(2));
t273 = t199 * t194;
t121 = t202 * t195 + t273;
t102 = t121 * qJD(1);
t255 = qJD(1) * t202;
t241 = t194 * t255;
t256 = qJD(1) * t199;
t105 = t195 * t256 - t241;
t198 = sin(qJ(5));
t201 = cos(qJ(5));
t52 = -t198 * t102 + t201 * t105;
t306 = t52 ^ 2;
t295 = t201 * t102 + t198 * t105;
t305 = t295 ^ 2;
t187 = qJD(2) - qJD(5);
t304 = t52 * t187;
t282 = t52 * t295;
t303 = t187 * t295;
t173 = t199 * qJDD(1);
t249 = qJD(1) * qJD(2);
t238 = t202 * t249;
t302 = -t238 - t173;
t301 = t305 - t306;
t251 = qJD(5) * t201;
t252 = qJD(5) * t198;
t239 = t199 * t249;
t212 = qJDD(1) * t121 - t195 * t239;
t58 = t194 * t238 + t212;
t111 = t121 * qJD(2);
t248 = t202 * qJDD(1);
t222 = t195 * t173 - t194 * t248;
t59 = qJD(1) * t111 + t222;
t12 = t102 * t251 + t105 * t252 + t198 * t58 - t201 * t59;
t300 = t12 + t303;
t163 = pkin(6) * t173;
t237 = pkin(6) * t238 + qJDD(3) + t163;
t250 = t199 * qJD(4);
t290 = pkin(2) + pkin(3);
t42 = t302 * qJ(4) - qJD(1) * t250 - t290 * qJDD(2) + t237;
t164 = pkin(6) * t248;
t188 = qJDD(2) * qJ(3);
t189 = qJD(2) * qJD(3);
t245 = t164 + t188 + t189;
t254 = qJD(2) * t199;
t280 = pkin(6) - qJ(4);
t98 = -t202 * qJD(4) - t280 * t254;
t44 = -qJ(4) * t248 + qJD(1) * t98 + t245;
t235 = t194 * t44 - t195 * t42;
t10 = -qJDD(2) * pkin(4) - t59 * pkin(7) - t235;
t19 = t194 * t42 + t195 * t44;
t11 = -t58 * pkin(7) + t19;
t283 = t105 * pkin(7);
t168 = pkin(6) * t255;
t130 = -qJ(4) * t255 + t168;
t190 = qJD(2) * qJ(3);
t115 = t130 + t190;
t167 = pkin(6) * t256;
t128 = qJ(4) * t256 - t167;
t242 = t290 * qJD(2);
t83 = qJD(3) - t242 - t128;
t36 = -t194 * t115 + t195 * t83;
t24 = -qJD(2) * pkin(4) - t283 + t36;
t284 = t102 * pkin(7);
t37 = t195 * t115 + t194 * t83;
t25 = t37 - t284;
t1 = t198 * t10 + t201 * t11 + t24 * t251 - t25 * t252;
t186 = pkin(8) + qJ(5);
t171 = sin(t186);
t172 = cos(t186);
t218 = t202 * t171 - t199 * t172;
t116 = -qJD(1) * pkin(1) - pkin(2) * t255 - qJ(3) * t256;
t76 = pkin(3) * t255 + qJD(4) - t116;
t45 = t102 * pkin(4) + t76;
t200 = sin(qJ(1));
t217 = t199 * t171 + t202 * t172;
t73 = t217 * t200;
t203 = cos(qJ(1));
t75 = t217 * t203;
t299 = g(1) * t75 + g(2) * t73 - g(3) * t218 + t295 * t45 - t1;
t6 = t198 * t24 + t201 * t25;
t2 = -t6 * qJD(5) + t201 * t10 - t198 * t11;
t72 = t218 * t200;
t265 = t202 * t203;
t271 = t199 * t203;
t74 = t171 * t265 - t172 * t271;
t298 = g(1) * t74 + g(2) * t72 + g(3) * t217 - t45 * t52 + t2;
t261 = t202 * pkin(2) + t199 * qJ(3);
t297 = -pkin(1) - t261;
t296 = g(1) * t203 + g(2) * t200;
t277 = pkin(6) * qJDD(2);
t293 = (qJD(1) * t297 + t116) * qJD(2) - t277;
t292 = t102 ^ 2;
t291 = t105 ^ 2;
t286 = t58 * pkin(4);
t183 = g(1) * t200;
t285 = g(2) * t203;
t177 = t202 * pkin(3);
t120 = -t198 * t194 + t201 * t195;
t62 = -t194 * t128 + t195 * t130;
t29 = t62 - t284;
t63 = t195 * t128 + t194 * t130;
t30 = t63 + t283;
t131 = -t194 * qJ(3) - t195 * t290;
t127 = -pkin(4) + t131;
t132 = t195 * qJ(3) - t194 * t290;
t64 = t201 * t127 - t198 * t132;
t279 = t120 * qJD(3) + t64 * qJD(5) - t198 * t29 - t201 * t30;
t123 = t201 * t194 + t198 * t195;
t65 = t198 * t127 + t201 * t132;
t278 = -t123 * qJD(3) - t65 * qJD(5) + t198 * t30 - t201 * t29;
t141 = t280 * t202;
t99 = qJD(2) * t141 - t250;
t39 = t194 * t99 + t195 * t98;
t193 = qJDD(1) * pkin(1);
t276 = qJDD(2) * pkin(2);
t275 = t105 * t102;
t272 = t199 * t200;
t206 = qJD(1) ^ 2;
t270 = t199 * t206;
t269 = t200 * t202;
t156 = t195 * pkin(4) + pkin(3);
t267 = t202 * t156;
t266 = t202 * t194;
t264 = t187 * t123;
t263 = t187 * t120;
t140 = t280 * t199;
t68 = t194 * t140 + t195 * t141;
t174 = t199 * qJD(3);
t253 = qJD(2) * t202;
t262 = qJ(3) * t253 + t174;
t260 = t203 * pkin(1) + t200 * pkin(6);
t191 = t199 ^ 2;
t192 = t202 ^ 2;
t258 = -t191 + t192;
t257 = t191 + t192;
t247 = pkin(4) * t266;
t246 = pkin(4) * t273;
t244 = t177 + t261;
t243 = -g(1) * t271 - g(2) * t272 + g(3) * t202;
t236 = t183 - t285;
t38 = -t194 * t98 + t195 * t99;
t234 = t198 * t59 + t201 * t58;
t233 = -qJD(2) * pkin(2) + qJD(3);
t232 = t194 * qJD(3) + t62;
t231 = t195 * qJD(3) - t63;
t67 = t195 * t140 - t194 * t141;
t117 = pkin(1) + t244;
t229 = pkin(2) * t265 + qJ(3) * t271 + t260;
t228 = -t163 - t243;
t227 = t199 * t242;
t226 = t199 * t238;
t225 = t257 * qJDD(1) * pkin(6);
t205 = qJD(2) ^ 2;
t224 = pkin(6) * t205 + t285;
t221 = t36 * t194 - t37 * t195;
t216 = -t199 * t195 + t266;
t40 = pkin(7) * t216 + t67;
t41 = -t121 * pkin(7) + t68;
t16 = -t198 * t41 + t201 * t40;
t17 = t198 * t40 + t201 * t41;
t220 = pkin(2) * t248 - t302 * qJ(3) + qJD(1) * t174 + t193;
t61 = -t198 * t121 - t201 * t216;
t134 = t167 + t233;
t139 = t168 + t190;
t219 = t134 * t202 - t139 * t199;
t161 = qJ(3) * t255;
t95 = -t290 * t256 + t161;
t93 = t237 - t276;
t214 = -0.2e1 * pkin(1) * t249 - t277;
t71 = -t227 + t262;
t213 = -t224 + 0.2e1 * t193;
t13 = t52 * qJD(5) + t234;
t100 = pkin(2) * t254 - t262;
t56 = pkin(2) * t239 - t220;
t211 = -qJD(1) * t100 - qJDD(1) * t297 - t224 - t56;
t81 = -pkin(6) * t239 + t245;
t210 = t219 * qJD(2) + t93 * t199 + t81 * t202;
t28 = pkin(3) * t248 - qJD(1) * t227 + qJDD(4) + t220;
t208 = t236 + t28;
t197 = -pkin(7) - qJ(4);
t185 = qJDD(2) - qJDD(5);
t179 = t203 * pkin(6);
t154 = g(1) * t269;
t150 = qJ(3) * t265;
t148 = qJ(3) * t269;
t147 = t202 * t270;
t138 = t258 * t206;
t137 = qJDD(2) * t202 - t205 * t199;
t136 = qJDD(2) * t199 + t205 * t202;
t129 = pkin(2) * t256 - t161;
t119 = t192 * qJDD(1) - 0.2e1 * t226;
t118 = t191 * qJDD(1) + 0.2e1 * t226;
t110 = t194 * t253 - t195 * t254;
t92 = t121 * t203;
t91 = t216 * t203;
t90 = t121 * t200;
t89 = t216 * t200;
t82 = t199 * t248 + t258 * t249;
t66 = t121 * pkin(4) + t117;
t60 = t201 * t121 - t198 * t216;
t57 = -t105 * pkin(4) + t95;
t43 = t110 * pkin(4) + t71;
t27 = -t110 * pkin(7) + t39;
t26 = -t111 * pkin(7) + t38;
t22 = t61 * qJD(5) + t201 * t110 + t198 * t111;
t21 = t198 * t110 - t201 * t111 + t121 * t251 - t216 * t252;
t20 = t28 + t286;
t5 = -t198 * t25 + t201 * t24;
t4 = -t17 * qJD(5) - t198 * t27 + t201 * t26;
t3 = t16 * qJD(5) + t198 * t26 + t201 * t27;
t7 = [0, 0, 0, 0, 0, qJDD(1), t236, t296, 0, 0, t118, 0.2e1 * t82, t136, t119, t137, 0, t199 * t214 + t202 * t213 + t154, t214 * t202 + (-t213 - t183) * t199, 0.2e1 * t225 - t296, -g(1) * (-t200 * pkin(1) + t179) - g(2) * t260 + (t257 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t118, t136, -0.2e1 * t82, 0, -t137, t119, t293 * t199 + t211 * t202 + t154, t225 + t210 - t296, -t293 * t202 + (t211 + t183) * t199, t210 * pkin(6) - g(1) * t179 - g(2) * t229 + t116 * t100 + (-t183 + t56) * t297, t105 * t111 - t216 * t59, -t111 * t102 - t105 * t110 - t59 * t121 + t216 * t58, -t111 * qJD(2) + qJDD(2) * t216, t102 * t110 + t58 * t121, t110 * qJD(2) + t121 * qJDD(2), 0, g(1) * t90 - g(2) * t92 - t38 * qJD(2) - t67 * qJDD(2) + t71 * t102 + t76 * t110 + t117 * t58 + t28 * t121, -g(1) * t89 + g(2) * t91 + t39 * qJD(2) + t68 * qJDD(2) + t71 * t105 + t76 * t111 + t117 * t59 - t216 * t28, -t39 * t102 - t38 * t105 - t37 * t110 - t36 * t111 - t19 * t121 - t216 * t235 - t68 * t58 - t67 * t59 + t296, t19 * t68 + t37 * t39 - t235 * t67 + t36 * t38 + t28 * t117 + t76 * t71 - g(1) * (-t203 * qJ(4) + t179) - g(2) * (pkin(3) * t265 + t229) + (-g(1) * (t297 - t177) + g(2) * qJ(4)) * t200, -t12 * t61 - t52 * t21, t12 * t60 - t61 * t13 + t21 * t295 - t52 * t22, -t61 * t185 + t21 * t187, t13 * t60 + t22 * t295, t60 * t185 + t22 * t187, 0, g(1) * t73 - g(2) * t75 + t66 * t13 - t16 * t185 - t4 * t187 + t20 * t60 + t45 * t22 + t295 * t43, -g(1) * t72 + g(2) * t74 - t66 * t12 + t17 * t185 + t3 * t187 + t20 * t61 - t45 * t21 + t43 * t52, -t1 * t60 + t16 * t12 - t17 * t13 - t2 * t61 + t5 * t21 - t6 * t22 - t295 * t3 - t4 * t52 + t296, t1 * t17 + t6 * t3 + t2 * t16 + t5 * t4 + t20 * t66 + t45 * t43 - g(1) * (t203 * t197 + t179) - g(2) * (t156 * t265 + t203 * t246 + t229) + (-g(1) * (t297 - t246 - t267) - g(2) * t197) * t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, -t138, t173, t147, t248, qJDD(2), pkin(1) * t270 + t228, g(3) * t199 - t164 + (pkin(1) * t206 + t296) * t202, 0, 0, -t147, t173, t138, qJDD(2), -t248, t147, 0.2e1 * t276 - qJDD(3) + (-t116 * t199 + t129 * t202) * qJD(1) + t228, (-pkin(2) * t199 + qJ(3) * t202) * qJDD(1) + ((t139 - t190) * t199 + (-t134 + t233) * t202) * qJD(1), t164 + 0.2e1 * t188 + 0.2e1 * t189 + (qJD(1) * t129 - g(3)) * t199 + (qJD(1) * t116 - t296) * t202, t81 * qJ(3) + t139 * qJD(3) - t93 * pkin(2) - t116 * t129 - g(1) * (-pkin(2) * t271 + t150) - g(2) * (-pkin(2) * t272 + t148) - g(3) * t261 - t219 * qJD(1) * pkin(6), -t275, -t291 + t292, -t222, t275, (t105 + t241) * qJD(2) + t212, qJDD(2), -g(1) * t91 - g(2) * t89 - g(3) * t121 + t232 * qJD(2) - t131 * qJDD(2) - t95 * t102 + t76 * t105 + t235, -g(1) * t92 - g(2) * t90 + g(3) * t216 + t231 * qJD(2) + t132 * qJDD(2) - t76 * t102 - t95 * t105 + t19, -t131 * t59 - t132 * t58 + (t232 - t37) * t105 + (-t231 + t36) * t102, t296 * t199 * t290 - g(1) * t150 - g(2) * t148 - g(3) * t244 - t221 * qJD(3) - t131 * t235 + t19 * t132 - t36 * t62 - t37 * t63 - t76 * t95, -t282, t301, t300, t282, t13 + t304, t185, -t64 * t185 - t278 * t187 - t295 * t57 - t298, t65 * t185 + t279 * t187 - t57 * t52 - t299, t64 * t12 - t65 * t13 + (-t278 - t6) * t52 + (-t279 + t5) * t295, t1 * t65 + t2 * t64 - t45 * t57 - g(1) * (t203 * t247 + t150) - g(2) * (t200 * t247 + t148) - g(3) * (t261 + t267) + t279 * t6 + t278 * t5 + (-g(3) * pkin(4) * t194 + t296 * (pkin(2) + t156)) * t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t147, t173, -t191 * t206 - t205, -t139 * qJD(2) + t116 * t256 + t243 + t93, 0, 0, 0, 0, 0, 0, -t195 * qJDD(2) - t102 * t256 - t194 * t205, t194 * qJDD(2) - t105 * t256 - t195 * t205, -t194 * t58 - t195 * t59 + (t102 * t195 - t105 * t194) * qJD(2), t221 * qJD(2) + t19 * t194 - t195 * t235 - t76 * t256 + t243, 0, 0, 0, 0, 0, 0, -t120 * t185 - t264 * t187 - t256 * t295, t123 * t185 - t263 * t187 - t52 * t256, t120 * t12 - t123 * t13 + t263 * t295 - t264 * t52, t1 * t123 + t2 * t120 - t256 * t45 - t263 * t6 + t264 * t5 + t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t105 + t241) * qJD(2) + t212, 0.2e1 * t102 * qJD(2) + t222, -t291 - t292, t37 * t102 + t36 * t105 + t208, 0, 0, 0, 0, 0, 0, t13 - t304, -t12 + t303, -t305 - t306, t295 * t6 + t5 * t52 + t208 + t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, -t301, -t300, -t282, -t234 + (-qJD(5) - t187) * t52, -t185, -t6 * t187 + t298, -t5 * t187 + t299, 0, 0;];
tau_reg = t7;
