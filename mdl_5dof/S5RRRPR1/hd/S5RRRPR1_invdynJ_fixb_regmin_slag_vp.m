% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:48:58
% EndTime: 2021-01-15 22:49:12
% DurationCPUTime: 3.42s
% Computational Cost: add. (4734->338), mult. (11561->445), div. (0->0), fcn. (8650->16), ass. (0->191)
t200 = sin(qJ(5));
t204 = cos(qJ(5));
t205 = cos(qJ(3));
t206 = cos(qJ(2));
t257 = qJD(1) * t206;
t247 = t205 * t257;
t201 = sin(qJ(3));
t202 = sin(qJ(2));
t258 = qJD(1) * t202;
t248 = t201 * t258;
t124 = -t247 + t248;
t126 = -t201 * t257 - t205 * t258;
t198 = sin(pkin(9));
t199 = cos(pkin(9));
t229 = -t198 * t124 - t199 * t126;
t97 = -t199 * t124 + t198 * t126;
t51 = t200 * t97 + t204 * t229;
t273 = t204 * t97;
t54 = -t200 * t229 + t273;
t278 = t51 * t54;
t194 = qJD(2) + qJD(3);
t187 = qJD(5) + t194;
t271 = t54 * t187;
t254 = qJD(5) * t200;
t251 = t206 * qJDD(1);
t253 = qJD(1) * qJD(2);
t245 = t206 * t253;
t252 = t202 * qJDD(1);
t290 = t245 + t252;
t81 = qJD(3) * t247 - t194 * t248 + t201 * t251 + t290 * t205;
t136 = t201 * t206 + t205 * t202;
t106 = t194 * t136;
t234 = t201 * t252 - t205 * t251;
t82 = qJD(1) * t106 + t234;
t39 = t198 * t81 + t199 * t82;
t40 = -t198 * t82 + t199 * t81;
t9 = qJD(5) * t273 - t200 * t39 + t204 * t40 - t229 * t254;
t4 = t9 - t271;
t214 = -qJD(5) * t51 - t200 * t40 - t204 * t39;
t272 = t51 * t187;
t5 = t214 + t272;
t8 = t51 ^ 2 - t54 ^ 2;
t197 = qJ(2) + qJ(3);
t188 = pkin(9) + t197;
t180 = qJ(5) + t188;
t169 = sin(t180);
t170 = cos(t180);
t203 = sin(qJ(1));
t207 = cos(qJ(1));
t237 = g(1) * t207 + g(2) * t203;
t284 = t97 * pkin(8);
t285 = pkin(6) + pkin(7);
t157 = t285 * t206;
t143 = qJD(1) * t157;
t131 = t205 * t143;
t156 = t285 * t202;
t141 = qJD(1) * t156;
t275 = qJD(2) * pkin(2);
t133 = -t141 + t275;
t228 = -t201 * t133 - t131;
t266 = t124 * qJ(4);
t80 = -t228 - t266;
t274 = t199 * t80;
t120 = t126 * qJ(4);
t127 = t201 * t143;
t242 = t205 * t133 - t127;
t79 = t120 + t242;
t69 = t194 * pkin(3) + t79;
t33 = t198 * t69 + t274;
t25 = t33 + t284;
t191 = t206 * pkin(2);
t277 = pkin(1) + t191;
t155 = t277 * qJD(1);
t108 = t124 * pkin(3) + qJD(4) - t155;
t64 = -pkin(4) * t97 + t108;
t224 = g(3) * t169 + t237 * t170 + t25 * t254 - t64 * t54;
t192 = qJDD(2) + qJDD(3);
t107 = qJDD(2) * pkin(2) - t285 * t290;
t246 = t202 * t253;
t109 = t285 * (-t246 + t251);
t213 = t228 * qJD(3) + t205 * t107 - t201 * t109;
t20 = t192 * pkin(3) - t81 * qJ(4) + t126 * qJD(4) + t213;
t256 = qJD(3) * t201;
t286 = (qJD(3) * t133 + t109) * t205 + t201 * t107 - t143 * t256;
t22 = -t82 * qJ(4) - t124 * qJD(4) + t286;
t6 = -t198 * t22 + t199 * t20;
t2 = t192 * pkin(4) - t40 * pkin(8) + t6;
t7 = t198 * t20 + t199 * t22;
t3 = -t39 * pkin(8) + t7;
t216 = -g(3) * t170 + t237 * t169 + t204 * t2 - t200 * t3 - t64 * t51;
t189 = sin(t197);
t190 = cos(t197);
t292 = -g(3) * t190 + t237 * t189;
t92 = t229 * pkin(8);
t263 = t199 * t201;
t276 = pkin(2) * qJD(3);
t241 = t201 * t141 - t131;
t83 = t241 + t266;
t262 = -t205 * t141 - t127;
t84 = t120 + t262;
t269 = t198 * t84 - t199 * t83 + (-t198 * t205 - t263) * t276;
t264 = t198 * t201;
t267 = -t198 * t83 - t199 * t84 + (t199 * t205 - t264) * t276;
t261 = -t201 * t156 + t205 * t157;
t287 = qJD(5) - t187;
t283 = pkin(3) * t198;
t279 = t126 * pkin(3);
t135 = t201 * t202 - t205 * t206;
t249 = qJD(2) * t285;
t142 = t202 * t249;
t144 = t206 * t249;
t255 = qJD(3) * t205;
t221 = -t205 * t142 - t201 * t144 - t156 * t255 - t157 * t256;
t45 = -t106 * qJ(4) - t135 * qJD(4) + t221;
t105 = t194 * t135;
t212 = -t261 * qJD(3) + t201 * t142 - t205 * t144;
t46 = t105 * qJ(4) - t136 * qJD(4) + t212;
t16 = t198 * t46 + t199 * t45;
t70 = t198 * t80;
t38 = t199 * t79 - t70;
t240 = -t205 * t156 - t201 * t157;
t93 = -t136 * qJ(4) + t240;
t94 = -t135 * qJ(4) + t261;
t50 = t198 * t93 + t199 * t94;
t270 = t284 + t269;
t268 = -t92 - t267;
t265 = t126 * t124;
t260 = pkin(3) * t190 + t191;
t195 = t202 ^ 2;
t259 = -t206 ^ 2 + t195;
t185 = t202 * t275;
t99 = t106 * pkin(3) + t185;
t15 = -t198 * t45 + t199 * t46;
t32 = t199 * t69 - t70;
t37 = -t198 * t79 - t274;
t49 = -t198 * t94 + t199 * t93;
t111 = t135 * pkin(3) - t277;
t182 = t205 * pkin(2) + pkin(3);
t118 = -pkin(2) * t264 + t199 * t182;
t68 = pkin(4) * t229 - t279;
t236 = g(1) * t203 - g(2) * t207;
t235 = t229 * t33 + t32 * t97;
t23 = t194 * pkin(4) + t32 - t92;
t233 = -t200 * t23 - t204 * t25;
t103 = -t198 * t135 + t199 * t136;
t30 = -t103 * pkin(8) + t49;
t102 = t199 * t135 + t198 * t136;
t31 = -t102 * pkin(8) + t50;
t232 = -t200 * t31 + t204 * t30;
t231 = t200 * t30 + t204 * t31;
t56 = t204 * t102 + t200 * t103;
t57 = -t200 * t102 + t204 * t103;
t177 = t199 * pkin(3) + pkin(4);
t227 = t200 * t177 + t204 * t283;
t226 = t204 * t177 - t200 * t283;
t225 = -0.2e1 * pkin(1) * t253 - pkin(6) * qJDD(2);
t121 = pkin(2) * t246 - qJDD(1) * t277;
t175 = sin(t188);
t176 = cos(t188);
t220 = g(3) * t175 - t108 * t97 + t237 * t176 - t7;
t208 = qJD(2) ^ 2;
t218 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t208 + t236;
t209 = qJD(1) ^ 2;
t217 = pkin(1) * t209 - pkin(6) * qJDD(1) + t237;
t59 = t82 * pkin(3) + qJDD(4) + t121;
t215 = -g(3) * t176 - t108 * t229 + t237 * t175 + t6;
t211 = g(3) * t189 - t155 * t124 + t237 * t190 - t286;
t210 = -t155 * t126 + t213 + t292;
t193 = -qJ(4) - t285;
t186 = qJDD(5) + t192;
t184 = pkin(2) * t258;
t140 = pkin(1) + t260;
t119 = pkin(2) * t263 + t198 * t182;
t114 = pkin(4) + t118;
t110 = t184 - t279;
t85 = -t124 ^ 2 + t126 ^ 2;
t74 = t102 * pkin(4) + t111;
t65 = t184 + t68;
t63 = -t199 * t105 - t198 * t106;
t62 = -t198 * t105 + t199 * t106;
t61 = -t234 + (-qJD(1) * t136 - t126) * t194;
t60 = t124 * t194 + t81;
t44 = t62 * pkin(4) + t99;
t27 = -t92 + t38;
t26 = t37 - t284;
t17 = t39 * pkin(4) + t59;
t14 = t57 * qJD(5) + t200 * t63 + t204 * t62;
t13 = -t56 * qJD(5) - t200 * t62 + t204 * t63;
t12 = -t62 * pkin(8) + t16;
t11 = -t63 * pkin(8) + t15;
t1 = [qJDD(1), t236, t237, t195 * qJDD(1) + 0.2e1 * t202 * t245, 0.2e1 * t202 * t251 - 0.2e1 * t259 * t253, qJDD(2) * t202 + t208 * t206, qJDD(2) * t206 - t208 * t202, 0, t225 * t202 + t218 * t206, -t218 * t202 + t225 * t206, t126 * t105 + t81 * t136, t105 * t124 + t126 * t106 - t81 * t135 - t136 * t82, -t105 * t194 + t136 * t192, -t106 * t194 - t135 * t192, 0, -t155 * t106 + t121 * t135 + t124 * t185 + t236 * t190 + t240 * t192 + t212 * t194 - t277 * t82, t155 * t105 + t121 * t136 - t126 * t185 - t236 * t189 - t261 * t192 - t221 * t194 - t277 * t81, t59 * t102 + t108 * t62 + t111 * t39 + t15 * t194 + t236 * t176 + t49 * t192 - t97 * t99, t59 * t103 + t108 * t63 + t111 * t40 - t16 * t194 - t236 * t175 - t50 * t192 + t229 * t99, -t7 * t102 - t6 * t103 - t15 * t229 + t16 * t97 - t32 * t63 - t33 * t62 - t50 * t39 - t49 * t40 - t237, t7 * t50 + t33 * t16 + t6 * t49 + t32 * t15 + t59 * t111 + t108 * t99 - g(1) * (-t203 * t140 - t207 * t193) - g(2) * (t207 * t140 - t203 * t193), t13 * t51 + t9 * t57, t13 * t54 - t14 * t51 + t214 * t57 - t9 * t56, t13 * t187 + t57 * t186, -t14 * t187 - t56 * t186, 0, -t44 * t54 - t74 * t214 + t17 * t56 + t64 * t14 + (-qJD(5) * t231 + t204 * t11 - t200 * t12) * t187 + t232 * t186 + t236 * t170, t44 * t51 + t74 * t9 + t17 * t57 + t64 * t13 - (qJD(5) * t232 + t200 * t11 + t204 * t12) * t187 - t231 * t186 - t236 * t169; 0, 0, 0, -t202 * t209 * t206, t259 * t209, t252, t251, qJDD(2), -g(3) * t206 + t202 * t217, g(3) * t202 + t206 * t217, -t265, t85, t60, t61, t192, -t241 * t194 + (-t124 * t258 + t205 * t192 - t194 * t256) * pkin(2) + t210, t262 * t194 + (t126 * t258 - t201 * t192 - t194 * t255) * pkin(2) + t211, t110 * t97 + t118 * t192 + t269 * t194 + t215, -t110 * t229 - t119 * t192 - t267 * t194 + t220, -t118 * t40 - t119 * t39 - t229 * t269 + t267 * t97 + t235, t7 * t119 + t6 * t118 - t108 * t110 - g(3) * t260 + t267 * t33 + t269 * t32 - t237 * (-t202 * pkin(2) - pkin(3) * t189), -t278, t8, t4, t5, t186, (t204 * t114 - t200 * t119) * t186 + t65 * t54 + (t268 * t200 + t270 * t204) * t187 + ((-t200 * t114 - t204 * t119) * t187 + t233) * qJD(5) + t216, -t65 * t51 + (-t114 * t186 - t2 + (qJD(5) * t119 - t270) * t187) * t200 + (-qJD(5) * t23 - t119 * t186 - t3 + (-qJD(5) * t114 + t268) * t187) * t204 + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t265, t85, t60, t61, t192, -t228 * t194 + t210, t242 * t194 + t211, -t37 * t194 + (-t126 * t97 + t192 * t199) * pkin(3) + t215, t38 * t194 + (t126 * t229 - t192 * t198) * pkin(3) + t220, t37 * t229 - t38 * t97 + (-t198 * t39 - t199 * t40) * pkin(3) + t235, -t32 * t37 - t33 * t38 + (t108 * t126 + t198 * t7 + t199 * t6 + t292) * pkin(3), -t278, t8, t4, t5, t186, t226 * t186 + t68 * t54 - (-t200 * t27 + t204 * t26) * t187 + (-t187 * t227 + t233) * qJD(5) + t216, -t227 * t186 - t204 * t3 - t200 * t2 - t68 * t51 + (t200 * t26 + t204 * t27) * t187 + (-t187 * t226 - t204 * t23) * qJD(5) + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194 * t229 + t39, t97 * t194 + t40, -t229 ^ 2 - t97 ^ 2, t229 * t32 - t33 * t97 - t236 + t59, 0, 0, 0, 0, 0, -t214 + t272, t9 + t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, t8, t4, t5, t186, t287 * t233 + t216, (-t25 * t187 - t2) * t200 + (-t287 * t23 - t3) * t204 + t224;];
tau_reg = t1;
