% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:29:58
% EndTime: 2021-01-15 20:30:11
% DurationCPUTime: 2.68s
% Computational Cost: add. (3620->388), mult. (8442->485), div. (0->0), fcn. (5985->10), ass. (0->196)
t136 = sin(qJ(4));
t139 = cos(qJ(4));
t133 = sin(pkin(8));
t140 = cos(qJ(2));
t137 = sin(qJ(2));
t210 = cos(pkin(8));
t176 = t210 * t137;
t100 = t133 * t140 + t176;
t90 = t100 * qJD(1);
t74 = qJD(2) * t136 + t139 * t90;
t255 = t136 * t74;
t175 = t210 * t140;
t115 = qJD(1) * t175;
t193 = qJD(1) * qJD(2);
t181 = t137 * t193;
t150 = qJDD(1) * t100 - t133 * t181;
t146 = qJD(2) * t115 + t150;
t189 = pkin(2) * t181 + qJDD(3);
t190 = t140 * qJDD(1);
t208 = qJDD(1) * pkin(1);
t191 = t137 * qJDD(1);
t165 = -qJDD(1) * t175 + t133 * t191;
t89 = t100 * qJD(2);
t63 = qJD(1) * t89 + t165;
t18 = -pkin(2) * t190 + t63 * pkin(3) - pkin(7) * t146 + t189 - t208;
t15 = t139 * t18;
t229 = t140 * pkin(2);
t123 = pkin(1) + t229;
t106 = -qJD(1) * t123 + qJD(3);
t196 = qJD(1) * t137;
t87 = t133 * t196 - t115;
t33 = t87 * pkin(3) - t90 * pkin(7) + t106;
t135 = -qJ(3) - pkin(6);
t108 = t135 * t140;
t103 = qJD(1) * t108;
t177 = t210 * t103;
t182 = t135 * t137;
t102 = qJD(1) * t182;
t218 = qJD(2) * pkin(2);
t96 = t102 + t218;
t62 = t133 * t96 - t177;
t52 = qJD(2) * pkin(7) + t62;
t17 = t136 * t33 + t139 * t52;
t178 = qJD(2) * t135;
t154 = -t137 * qJD(3) + t140 * t178;
t57 = qJDD(2) * pkin(2) + qJD(1) * t154 + qJDD(1) * t182;
t86 = t140 * qJD(3) + t137 * t178;
t64 = qJD(1) * t86 - qJDD(1) * t108;
t25 = t133 * t57 + t210 * t64;
t23 = qJDD(2) * pkin(7) + t25;
t149 = -qJD(4) * t17 - t136 * t23 + t15;
t192 = qJD(4) * qJD(2);
t195 = qJD(4) * t136;
t157 = t136 * qJDD(2) - t90 * t195 + (t146 + t192) * t139;
t212 = t157 * qJ(5);
t58 = qJDD(4) + t63;
t245 = t58 * pkin(4);
t1 = -t74 * qJD(5) + t149 - t212 + t245;
t72 = -qJD(2) * t139 + t136 * t90;
t11 = -qJ(5) * t72 + t17;
t81 = qJD(4) + t87;
t232 = t11 * t81;
t254 = t1 + t232;
t253 = t81 * t255;
t252 = t139 * t58 - t81 * t195;
t130 = qJ(2) + pkin(8);
t124 = sin(t130);
t138 = sin(qJ(1));
t141 = cos(qJ(1));
t167 = g(1) * t141 + g(2) * t138;
t160 = t167 * t124;
t239 = g(1) * t138;
t251 = -g(2) * t141 + t239;
t234 = g(3) * t136;
t125 = cos(t130);
t200 = t141 * t139;
t204 = t138 * t136;
t82 = t125 * t204 + t200;
t201 = t141 * t136;
t203 = t138 * t139;
t84 = -t125 * t201 + t203;
t250 = -g(1) * t84 + g(2) * t82 + t124 * t234;
t236 = g(3) * t124;
t249 = t125 * t167 + t236;
t248 = t74 ^ 2;
t244 = t72 * pkin(4);
t16 = -t136 * t52 + t139 * t33;
t10 = -qJ(5) * t74 + t16;
t8 = pkin(4) * t81 + t10;
t243 = t10 - t8;
t231 = t133 * pkin(2);
t118 = pkin(7) + t231;
t199 = qJ(5) + t118;
t172 = qJD(4) * t199;
t202 = t139 * qJ(5);
t43 = pkin(2) * t196 + pkin(3) * t90 + pkin(7) * t87;
t37 = t139 * t43;
t93 = t133 * t103;
t66 = t102 * t210 + t93;
t242 = -t90 * pkin(4) - t139 * t172 - t87 * t202 - t37 + (-qJD(5) + t66) * t136;
t241 = pkin(2) * t137;
t240 = pkin(4) * t136;
t235 = g(3) * t125;
t233 = g(3) * t140;
t230 = t139 * pkin(4);
t228 = t72 * t81;
t227 = t72 * t87;
t226 = t74 * t81;
t225 = t74 * t90;
t224 = t90 * t72;
t194 = qJD(4) * t139;
t145 = -qJDD(2) * t139 + t136 * t146;
t209 = qJD(4) * t74;
t28 = t145 + t209;
t223 = -t136 * t28 - t194 * t72;
t24 = -t133 * t64 + t210 * t57;
t222 = t136 * t43 + t139 * t66;
t156 = -t133 * t137 + t175;
t60 = -pkin(3) * t156 - pkin(7) * t100 - t123;
t70 = -t108 * t210 + t133 * t182;
t67 = t139 * t70;
t221 = t136 * t60 + t67;
t205 = t136 * qJ(5);
t220 = t139 * qJD(5) - t136 * t172 - t205 * t87 - t222;
t219 = qJ(5) * t28;
t217 = t136 * t58;
t215 = t136 * t87;
t92 = t156 * qJD(2);
t214 = t136 * t92;
t213 = t139 * t92;
t211 = t157 * t136;
t207 = t124 * t141;
t122 = pkin(3) + t230;
t206 = t125 * t122;
t198 = (g(1) * t200 + g(2) * t203) * t124;
t131 = t137 ^ 2;
t197 = -t140 ^ 2 + t131;
t42 = t133 * t154 + t210 * t86;
t187 = t137 * t218;
t44 = pkin(3) * t89 - pkin(7) * t92 + t187;
t188 = t136 * t44 + t139 * t42 + t194 * t60;
t22 = -qJDD(2) * pkin(3) - t24;
t5 = t28 * pkin(4) + qJDD(5) + t22;
t185 = -t5 - t235;
t184 = t210 * pkin(2);
t183 = t100 * t194;
t180 = qJD(5) + t244;
t41 = t133 * t86 - t154 * t210;
t174 = t136 * t18 + t139 * t23 + t194 * t33 - t195 * t52;
t65 = t102 * t133 - t177;
t173 = t139 * t81;
t69 = -t108 * t133 - t135 * t176;
t171 = -g(1) * t82 - g(2) * t84;
t83 = -t125 * t203 + t201;
t85 = t125 * t200 + t204;
t170 = -g(1) * t83 - g(2) * t85;
t119 = -t184 - pkin(3);
t169 = qJD(4) * t118 * t81 + t22;
t2 = -qJD(5) * t72 + t174 - t219;
t168 = -t8 * t81 + t2;
t61 = t210 * t96 + t93;
t164 = -qJ(5) * t92 - qJD(5) * t100;
t134 = -qJ(5) - pkin(7);
t163 = t124 * t134 - t206;
t162 = -t215 * t81 + t252;
t159 = -0.2e1 * pkin(1) * t193 - pkin(6) * qJDD(2);
t51 = -qJD(2) * pkin(3) - t61;
t158 = -t118 * t58 + t51 * t81;
t80 = -qJDD(1) * t123 + t189;
t153 = g(1) * t85 - g(2) * t83 + t139 * t236 - t174;
t142 = qJD(2) ^ 2;
t152 = -pkin(6) * t142 + 0.2e1 * t208 + t251;
t143 = qJD(1) ^ 2;
t151 = pkin(1) * t143 - pkin(6) * qJDD(1) + t167;
t148 = t149 + t250;
t144 = t136 * t192 + t194 * t90 + t145;
t121 = t135 * t141;
t113 = t141 * t123;
t112 = t125 * t234;
t107 = t119 - t230;
t98 = t199 * t139;
t97 = t199 * t136;
t71 = t72 ^ 2;
t50 = t139 * t60;
t39 = t100 * t240 + t69;
t38 = t139 * t44;
t30 = -pkin(4) * t215 + t65;
t29 = t180 + t51;
t21 = (t183 + t214) * pkin(4) + t41;
t19 = -t100 * t205 + t221;
t12 = -pkin(4) * t156 - t100 * t202 - t136 * t70 + t50;
t7 = -t139 * t81 ^ 2 - t217 - t225;
t6 = t162 - t224;
t4 = -qJ(5) * t183 + (-qJD(4) * t70 + t164) * t136 + t188;
t3 = t89 * pkin(4) - t136 * t42 + t38 + t164 * t139 + (-t67 + (qJ(5) * t100 - t60) * t136) * qJD(4);
t9 = [qJDD(1), t251, t167, qJDD(1) * t131 + 0.2e1 * t140 * t181, 0.2e1 * t137 * t190 - 0.2e1 * t193 * t197, qJDD(2) * t137 + t140 * t142, qJDD(2) * t140 - t137 * t142, 0, t137 * t159 + t140 * t152, -t137 * t152 + t140 * t159, -t69 * qJDD(2) + t106 * t89 - t123 * t63 - t80 * t156 + t251 * t125 + (t241 * t87 - t41) * qJD(2), g(2) * t207 - t42 * qJD(2) - t70 * qJDD(2) + t80 * t100 + t106 * t92 - t123 * t146 - t124 * t239 + t187 * t90, -t24 * t100 + t146 * t69 + t156 * t25 + t41 * t90 - t42 * t87 - t61 * t92 - t62 * t89 - t70 * t63 - t167, t25 * t70 + t62 * t42 - t24 * t69 - t61 * t41 - t80 * t123 + t106 * t187 - g(1) * (-t123 * t138 - t121) - g(2) * (-t138 * t135 + t113), t74 * t213 + (t139 * t157 - t195 * t74) * t100, (-t139 * t72 - t255) * t92 + (-t211 - t139 * t28 + (t136 * t72 - t139 * t74) * qJD(4)) * t100, t100 * t252 - t156 * t157 + t213 * t81 + t74 * t89, -t81 * t214 + t28 * t156 - t72 * t89 + (-t194 * t81 - t217) * t100, -t156 * t58 + t81 * t89, (-t70 * t194 + t38) * t81 + t50 * t58 - (-t52 * t194 + t15) * t156 + t16 * t89 + t41 * t72 + t69 * t28 + t51 * t183 + ((-qJD(4) * t60 - t42) * t81 - t70 * t58 - (-qJD(4) * t33 - t23) * t156 + t22 * t100 + t51 * t92) * t136 + t170, -(-t195 * t70 + t188) * t81 - t221 * t58 + t174 * t156 - t17 * t89 + t41 * t74 + t69 * t157 + t51 * t213 + (t22 * t139 - t195 * t51) * t100 + t171, t29 * t214 - t1 * t156 + t12 * t58 + t21 * t72 + t39 * t28 + t3 * t81 + t8 * t89 + (t5 * t136 + t194 * t29) * t100 + t170, t29 * t213 - t11 * t89 - t19 * t58 + t2 * t156 + t21 * t74 + t39 * t157 - t4 * t81 + (t5 * t139 - t195 * t29) * t100 + t171, -t12 * t157 - t19 * t28 - t3 * t74 - t4 * t72 + (-t11 * t136 - t139 * t8) * t92 + t251 * t124 + (-t1 * t139 - t2 * t136 + (-t11 * t139 + t136 * t8) * qJD(4)) * t100, t2 * t19 + t11 * t4 + t1 * t12 + t8 * t3 + t5 * t39 + t29 * t21 - g(1) * (pkin(4) * t201 - t121) - g(2) * (-t134 * t207 + t141 * t206 + t113) + (-g(1) * (-t123 + t163) - g(2) * (-t135 + t240)) * t138; 0, 0, 0, -t137 * t143 * t140, t197 * t143, t191, t190, qJDD(2), t137 * t151 - t233, g(3) * t137 + t140 * t151, -t235 + t65 * qJD(2) - t106 * t90 + t160 + (qJDD(2) * t210 - t196 * t87) * pkin(2) + t24, t66 * qJD(2) + t106 * t87 + (-qJDD(2) * t133 - t196 * t90) * pkin(2) - t25 + t249, -t146 * t184 - t63 * t231 - (-t62 + t65) * t90 + (t66 - t61) * t87, t61 * t65 - t62 * t66 + (t210 * t24 - t233 + t133 * t25 + (-qJD(1) * t106 + t167) * t137) * pkin(2), t173 * t74 + t211, (t157 - t227) * t139 - t253 + t223, t173 * t81 + t217 - t225, t162 + t224, -t81 * t90, t119 * t28 - t16 * t90 - t37 * t81 - t65 * t72 + (-t169 - t235) * t139 + (t66 * t81 + t158) * t136 + t198, t119 * t157 + t222 * t81 + t17 * t90 - t65 * t74 + t112 + t158 * t139 + (-t160 + t169) * t136, t107 * t28 - t30 * t72 - t97 * t58 - t8 * t90 + t242 * t81 + t185 * t139 + (t29 * t87 + (t29 + t244) * qJD(4)) * t136 + t198, t107 * t157 + t11 * t90 - t30 * t74 - t98 * t58 + t112 - t220 * t81 + t29 * t173 + (pkin(4) * t209 - t160 + t5) * t136, -t136 * t254 + t139 * t168 + t157 * t97 - t220 * t72 - t242 * t74 - t98 * t28 - t249, t2 * t98 - t1 * t97 + t5 * t107 - g(3) * (-t163 + t229) + t242 * t8 + (pkin(4) * t195 - t30) * t29 + t220 * t11 + t167 * (t122 * t124 + t125 * t134 + t241); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * t90 + t165, (t115 - t87) * qJD(2) + t150, -t87 ^ 2 - t90 ^ 2, t61 * t90 + t62 * t87 - t251 + t80, 0, 0, 0, 0, 0, t6, t7, t6, t7, (-t157 - t227) * t139 + t253 + t223, t136 * t168 + t139 * t254 - t29 * t90 - t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t72, -t71 + t248, t157 + t228, -t144 + t226, t58, t17 * t81 - t51 * t74 + t148, t16 * t81 + t51 * t72 + t153, 0.2e1 * t245 - t212 + t232 + (-t180 - t29) * t74 + t148, -pkin(4) * t248 + t219 + t10 * t81 + (qJD(5) + t29) * t72 + t153, -pkin(4) * t157 + t243 * t72, -t243 * t11 + (-t29 * t74 + t1 + t250) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144 + t226, t157 - t228, -t71 - t248, t11 * t72 + t8 * t74 - t160 - t185;];
tau_reg = t9;
