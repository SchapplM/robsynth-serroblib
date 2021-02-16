% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP7
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
% Datum: 2021-01-15 20:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:41:23
% EndTime: 2021-01-15 20:41:36
% DurationCPUTime: 2.87s
% Computational Cost: add. (3670->398), mult. (8469->495), div. (0->0), fcn. (5973->10), ass. (0->180)
t127 = qJ(2) + pkin(8);
t121 = sin(t127);
t134 = sin(qJ(1));
t137 = cos(qJ(1));
t237 = -g(1) * t134 + g(2) * t137;
t240 = t237 * t121;
t136 = cos(qJ(2));
t207 = cos(pkin(8));
t175 = t207 * t136;
t111 = qJD(1) * t175;
t130 = sin(pkin(8));
t133 = sin(qJ(2));
t196 = qJD(1) * t133;
t86 = t130 * t196 - t111;
t79 = qJD(4) + t86;
t236 = -g(1) * t137 - g(2) * t134;
t239 = t236 * t121;
t176 = t207 * t133;
t98 = t130 * t136 + t176;
t89 = t98 * qJD(1);
t238 = t89 * qJD(2);
t192 = qJD(1) * qJD(2);
t181 = t133 * t192;
t145 = qJDD(1) * t98 - t130 * t181;
t141 = qJD(2) * t111 + t145;
t131 = -qJ(3) - pkin(6);
t107 = t131 * t136;
t101 = qJD(1) * t107;
t92 = t130 * t101;
t182 = t131 * t133;
t100 = qJD(1) * t182;
t212 = qJD(2) * pkin(2);
t95 = t100 + t212;
t59 = t207 * t95 + t92;
t49 = -qJD(2) * pkin(3) - t59;
t132 = sin(qJ(4));
t135 = cos(qJ(4));
t71 = -qJD(2) * t135 + t132 * t89;
t73 = qJD(2) * t132 + t135 * t89;
t16 = pkin(4) * t71 - qJ(5) * t73 + t49;
t223 = t130 * pkin(2);
t114 = pkin(7) + t223;
t190 = t133 * qJDD(1);
t162 = -qJDD(1) * t175 + t130 * t190;
t88 = t98 * qJD(2);
t61 = qJD(1) * t88 + t162;
t56 = qJDD(4) + t61;
t211 = t114 * t56;
t235 = t16 * t79 - t211;
t122 = cos(t127);
t226 = g(3) * t121;
t234 = t122 * t236 - t226;
t233 = t73 ^ 2;
t232 = t79 ^ 2;
t231 = pkin(4) * t56;
t222 = t136 * pkin(2);
t120 = pkin(1) + t222;
t105 = -qJD(1) * t120 + qJD(3);
t34 = t86 * pkin(3) - t89 * pkin(7) + t105;
t177 = t207 * t101;
t60 = t130 * t95 - t177;
t50 = qJD(2) * pkin(7) + t60;
t15 = t132 * t34 + t135 * t50;
t9 = qJ(5) * t79 + t15;
t230 = t79 * t9;
t225 = g(3) * t122;
t224 = g(3) * t136;
t221 = t15 * t79;
t220 = t71 * t86;
t219 = t73 * t71;
t178 = t73 * t79;
t218 = t73 * t89;
t217 = t89 * t71;
t193 = qJD(4) * t135;
t140 = -qJDD(2) * t135 + t132 * t141;
t29 = qJD(4) * t73 + t140;
t216 = -t132 * t29 - t193 * t71;
t179 = qJD(2) * t131;
t151 = -t133 * qJD(3) + t136 * t179;
t55 = qJDD(2) * pkin(2) + qJD(1) * t151 + qJDD(1) * t182;
t84 = t136 * qJD(3) + t133 * t179;
t62 = qJD(1) * t84 - qJDD(1) * t107;
t24 = -t130 * t62 + t207 * t55;
t25 = t130 * t55 + t207 * t62;
t41 = pkin(2) * t196 + pkin(3) * t89 + pkin(7) * t86;
t64 = t100 * t207 + t92;
t215 = t132 * t41 + t135 * t64;
t152 = -t130 * t133 + t175;
t58 = -pkin(3) * t152 - pkin(7) * t98 - t120;
t70 = -t107 * t207 + t130 * t182;
t214 = t132 * t58 + t135 * t70;
t163 = pkin(4) * t132 - qJ(5) * t135;
t63 = t130 * t100 - t177;
t213 = -t132 * qJD(5) + t163 * t79 - t63;
t46 = t132 * t56;
t210 = t132 * t79;
t47 = t135 * t56;
t191 = qJD(4) * qJD(2);
t194 = qJD(4) * t132;
t28 = -t132 * qJDD(2) + t89 * t194 + (-t141 - t191) * t135;
t209 = t28 * t132;
t208 = t56 * qJ(5);
t206 = qJD(4) * t98;
t205 = qJDD(1) * pkin(1);
t204 = t134 * t131;
t203 = t134 * t132;
t202 = t134 * t135;
t201 = t137 * t132;
t200 = t137 * t135;
t14 = -t132 * t50 + t135 * t34;
t199 = qJD(5) - t14;
t198 = (g(1) * t200 + g(2) * t202) * t121;
t128 = t133 ^ 2;
t197 = -t136 ^ 2 + t128;
t195 = qJD(4) * t114;
t189 = t136 * qJDD(1);
t188 = pkin(2) * t181 + qJDD(3);
t187 = t133 * t212;
t186 = t79 * t195;
t185 = t98 * t194;
t184 = t98 * t193;
t183 = t207 * pkin(2);
t39 = t130 * t84 - t151 * t207;
t17 = -pkin(2) * t189 + t61 * pkin(3) - pkin(7) * t141 + t188 - t205;
t23 = qJDD(2) * pkin(7) + t25;
t174 = t132 * t23 - t135 * t17 + t193 * t50 + t194 * t34;
t69 = -t130 * t107 - t131 * t176;
t80 = t122 * t203 + t200;
t82 = t122 * t201 - t202;
t172 = -g(1) * t80 + g(2) * t82;
t81 = t122 * t202 - t201;
t83 = t122 * t200 + t203;
t171 = g(1) * t81 - g(2) * t83;
t22 = -qJDD(2) * pkin(3) - t24;
t3 = pkin(4) * t29 + qJ(5) * t28 - qJD(5) * t73 + t22;
t91 = t152 * qJD(2);
t170 = t16 * t91 + t3 * t98;
t8 = -pkin(4) * t79 + t199;
t167 = -t132 * t9 + t135 * t8;
t166 = t22 * t98 + t49 * t91;
t165 = t56 * t98 + t79 * t91;
t164 = -pkin(4) * t135 - qJ(5) * t132;
t161 = t46 + (t135 * t86 + t193) * t79;
t160 = -t194 * t79 - t210 * t86 + t47;
t159 = pkin(3) - t164;
t157 = t186 + t225;
t156 = -0.2e1 * pkin(1) * t192 - pkin(6) * qJDD(2);
t155 = t132 * t17 + t135 * t23 + t193 * t34 - t194 * t50;
t40 = t130 * t151 + t207 * t84;
t42 = pkin(3) * t88 - pkin(7) * t91 + t187;
t154 = t132 * t42 + t135 * t40 + t193 * t58 - t194 * t70;
t153 = t49 * t79 - t211;
t78 = -qJDD(1) * t120 + t188;
t149 = g(1) * t82 + g(2) * t80 + t132 * t226 - t174;
t148 = -t225 - t239;
t138 = qJD(2) ^ 2;
t147 = -pkin(6) * t138 + 0.2e1 * t205 - t237;
t139 = qJD(1) ^ 2;
t146 = pkin(1) * t139 - pkin(6) * qJDD(1) - t236;
t144 = t16 * t73 + qJDD(5) - t149;
t143 = -g(1) * t83 - g(2) * t81 - t135 * t226 + t155;
t117 = t131 * t137;
t115 = -t183 - pkin(3);
t106 = -pkin(3) * t130 + pkin(7) * t207;
t104 = pkin(3) * t207 + pkin(7) * t130 + pkin(2);
t96 = -t183 - t159;
t67 = t104 * t136 + t106 * t133 + pkin(1);
t31 = pkin(4) * t73 + qJ(5) * t71;
t30 = t163 * t98 + t69;
t19 = pkin(4) * t152 + t132 * t70 - t135 * t58;
t18 = -qJ(5) * t152 + t214;
t11 = -pkin(4) * t89 + t132 * t64 - t135 * t41;
t10 = qJ(5) * t89 + t215;
t7 = t71 * t79 - t28;
t6 = (pkin(4) * t91 + qJ(5) * t206) * t132 + (-qJ(5) * t91 + (pkin(4) * qJD(4) - qJD(5)) * t98) * t135 + t39;
t5 = -t88 * pkin(4) + qJD(4) * t214 + t132 * t40 - t135 * t42;
t4 = qJ(5) * t88 - qJD(5) * t152 + t154;
t2 = qJDD(5) + t174 - t231;
t1 = qJD(5) * t79 + t155 + t208;
t12 = [qJDD(1), -t237, -t236, qJDD(1) * t128 + 0.2e1 * t136 * t181, 0.2e1 * t133 * t189 - 0.2e1 * t192 * t197, qJDD(2) * t133 + t136 * t138, qJDD(2) * t136 - t133 * t138, 0, t133 * t156 + t136 * t147, -t133 * t147 + t136 * t156, -t69 * qJDD(2) + t105 * t88 - t120 * t61 - t78 * t152 - t237 * t122 + (pkin(2) * t133 * t86 - t39) * qJD(2), -t40 * qJD(2) - t70 * qJDD(2) + t105 * t91 - t120 * t141 + t89 * t187 + t78 * t98 + t240, t141 * t69 + t152 * t25 - t24 * t98 + t39 * t89 - t40 * t86 - t59 * t91 - t60 * t88 - t70 * t61 + t236, t25 * t70 + t60 * t40 - t24 * t69 - t59 * t39 - t78 * t120 + t105 * t187 - g(1) * (-t120 * t134 - t117) - g(2) * (t120 * t137 - t204), -t73 * t185 + (-t28 * t98 + t73 * t91) * t135, (-t132 * t73 - t135 * t71) * t91 + (t209 - t135 * t29 + (t132 * t71 - t135 * t73) * qJD(4)) * t98, t135 * t165 + t152 * t28 - t185 * t79 + t73 * t88, -t132 * t165 + t152 * t29 - t184 * t79 - t71 * t88, -t152 * t56 + t79 * t88, t174 * t152 + t14 * t88 + t39 * t71 + t69 * t29 + ((-qJD(4) * t70 + t42) * t79 + t58 * t56 + t49 * t206) * t135 + ((-qJD(4) * t58 - t40) * t79 - t70 * t56 + t166) * t132 + t171, t135 * t166 - t15 * t88 + t152 * t155 - t154 * t79 - t185 * t49 - t214 * t56 - t69 * t28 + t39 * t73 + t172, t132 * t170 + t152 * t2 + t16 * t184 - t19 * t56 + t30 * t29 - t5 * t79 + t6 * t71 - t8 * t88 + t171, -t18 * t29 - t19 * t28 - t4 * t71 + t5 * t73 + t167 * t91 - t240 + (-t1 * t132 + t2 * t135 + (-t132 * t8 - t135 * t9) * qJD(4)) * t98, -t1 * t152 - t135 * t170 + t16 * t185 + t18 * t56 + t30 * t28 + t4 * t79 - t6 * t73 + t9 * t88 - t172, t1 * t18 + t9 * t4 + t3 * t30 + t16 * t6 + t2 * t19 + t8 * t5 - g(1) * (-pkin(4) * t81 - qJ(5) * t80 - t134 * t67 - t117) - g(2) * (pkin(4) * t83 + qJ(5) * t82 + t137 * t67 - t204); 0, 0, 0, -t133 * t139 * t136, t197 * t139, t190, t189, qJDD(2), t133 * t146 - t224, g(3) * t133 + t136 * t146, t63 * qJD(2) - t105 * t89 + (qJDD(2) * t207 - t196 * t86) * pkin(2) + t148 + t24, t64 * qJD(2) + t105 * t86 + (-qJDD(2) * t130 - t196 * t89) * pkin(2) - t25 - t234, -t141 * t183 - t61 * t223 - (-t60 + t63) * t89 + (t64 - t59) * t86, t59 * t63 - t60 * t64 + (t207 * t24 - t224 + t130 * t25 + (-qJD(1) * t105 - t236) * t133) * pkin(2), t135 * t178 - t209, (-t28 - t220) * t135 - t73 * t210 + t216, t161 - t218, t160 + t217, -t79 * t89, t115 * t29 - t14 * t89 - t63 * t71 + (-t225 - t22 + (-t41 - t195) * t79) * t135 + (t64 * t79 + t153) * t132 + t198, -t115 * t28 + t215 * t79 + t15 * t89 - t63 * t73 + t153 * t135 + (t157 + t22 + t239) * t132, t11 * t79 + t96 * t29 + t8 * t89 + t213 * t71 + (-t157 - t3) * t135 + t235 * t132 + t198, t10 * t71 - t11 * t73 + (-t114 * t29 + t8 * t86 + t1 + (t114 * t73 + t8) * qJD(4)) * t135 + (-t114 * t28 - t86 * t9 + t2 + (t114 * t71 - t9) * qJD(4)) * t132 + t234, -t10 * t79 + t96 * t28 - t9 * t89 - t213 * t73 - t235 * t135 + (t148 - t3 - t186) * t132, t3 * t96 - t9 * t10 - t8 * t11 - g(3) * (t121 * pkin(7) + t222) + t213 * t16 - t159 * t225 + (qJD(4) * t167 + t1 * t135 + t2 * t132) * t114 + t236 * (-t104 * t133 + t106 * t136 + t121 * t164); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162 + 0.2e1 * t238, (t111 - t86) * qJD(2) + t145, -t86 ^ 2 - t89 ^ 2, t59 * t89 + t60 * t86 + t237 + t78, 0, 0, 0, 0, 0, t160 - t217, -t135 * t232 - t218 - t46, -t210 * t79 - t217 + t47, (t28 - t220) * t135 + t132 * t178 + t216, t161 + t218, -t16 * t89 + (-t2 + t230) * t135 + (t79 * t8 + t1) * t132 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219, -t71 ^ 2 + t233, t7, -t132 * t191 - t193 * t89 - t140 + t178, t56, -t49 * t73 + t149 + t221, t14 * t79 + t49 * t71 - t143, -t31 * t71 - t144 + t221 + 0.2e1 * t231, pkin(4) * t28 - qJ(5) * t29 + (-t15 + t9) * t73 + (t8 - t199) * t71, 0.2e1 * t208 - t16 * t71 + t31 * t73 + (0.2e1 * qJD(5) - t14) * t79 + t143, t1 * qJ(5) - t2 * pkin(4) - t16 * t31 - t8 * t15 - g(1) * (-pkin(4) * t82 + qJ(5) * t83) - g(2) * (-pkin(4) * t80 + qJ(5) * t81) + t199 * t9 + t163 * t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t162 + t219 - t238, t7, -t232 - t233, t144 - t230 - t231;];
tau_reg = t12;
