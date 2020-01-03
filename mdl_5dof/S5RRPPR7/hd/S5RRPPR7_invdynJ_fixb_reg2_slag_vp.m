% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:27
% EndTime: 2019-12-31 19:36:33
% DurationCPUTime: 3.41s
% Computational Cost: add. (3882->405), mult. (8937->501), div. (0->0), fcn. (6289->10), ass. (0->206)
t129 = sin(pkin(8));
t130 = cos(pkin(8));
t133 = sin(qJ(2));
t136 = cos(qJ(2));
t91 = t129 * t136 + t130 * t133;
t252 = t91 * qJD(1);
t258 = qJD(5) + t252;
t135 = cos(qJ(5));
t132 = sin(qJ(5));
t200 = t132 * qJD(2);
t215 = t130 * t136;
t189 = qJD(1) * t215;
t203 = qJD(1) * t133;
t80 = t129 * t203 - t189;
t62 = -t135 * t80 + t200;
t184 = t258 * t62;
t201 = qJD(5) * t135;
t196 = t136 * qJDD(1);
t197 = t133 * qJDD(1);
t162 = t129 * t197 - t130 * t196;
t82 = t91 * qJD(2);
t55 = qJD(1) * t82 + t162;
t22 = qJD(5) * t200 - qJDD(2) * t135 - t132 * t55 - t201 * t80;
t259 = t22 - t184;
t233 = t136 * pkin(2);
t118 = pkin(1) + t233;
t95 = -qJD(1) * t118 + qJD(3);
t149 = -qJ(4) * t252 + t95;
t244 = pkin(3) + pkin(7);
t24 = t244 * t80 + t149;
t227 = qJ(3) + pkin(6);
t97 = t227 * t136;
t94 = qJD(1) * t97;
t86 = t129 * t94;
t96 = t227 * t133;
t93 = qJD(1) * t96;
t89 = qJD(2) * pkin(2) - t93;
t53 = t130 * t89 - t86;
t174 = qJD(4) - t53;
t241 = t252 * pkin(4);
t25 = -qJD(2) * t244 + t174 + t241;
t6 = -t132 * t24 + t135 * t25;
t261 = t6 * t258;
t198 = qJD(1) * qJD(2);
t188 = t133 * t198;
t150 = qJDD(1) * t91 - t129 * t188;
t187 = t136 * t198;
t56 = t130 * t187 + t150;
t221 = t56 * qJ(4);
t69 = pkin(2) * t188 - qJDD(1) * t118 + qJDD(3);
t143 = -qJD(4) * t252 - t221 + t69;
t5 = t244 * t55 + t143;
t7 = t132 * t25 + t135 * t24;
t185 = qJD(2) * t227;
t76 = -t133 * qJD(3) - t136 * t185;
t50 = qJDD(2) * pkin(2) + qJD(1) * t76 - qJDD(1) * t96;
t75 = t136 * qJD(3) - t133 * t185;
t57 = qJD(1) * t75 + qJDD(1) * t97;
t226 = t129 * t57 - t130 * t50;
t191 = qJDD(4) + t226;
t9 = t56 * pkin(4) - qJDD(2) * t244 + t191;
t2 = -qJD(5) * t7 - t132 * t5 + t135 * t9;
t249 = t258 * t7 + t2;
t253 = t135 * t258;
t64 = qJD(2) * t135 + t132 * t80;
t260 = t64 * t253;
t182 = t132 * t258;
t51 = qJDD(5) + t56;
t158 = t135 * t51 - t182 * t258;
t257 = 0.2e1 * qJD(2) * t252 + t162;
t245 = t80 ^ 2;
t79 = t252 ^ 2;
t255 = -t245 - t79;
t254 = -t245 + t79;
t59 = -t130 * t93 - t86;
t208 = -qJD(4) + t59;
t125 = qJ(2) + pkin(8);
t121 = sin(t125);
t122 = cos(t125);
t251 = -pkin(3) * t122 - t121 * qJ(4);
t134 = sin(qJ(1));
t137 = cos(qJ(1));
t250 = g(1) * t134 - g(2) * t137;
t237 = g(2) * t134;
t171 = g(1) * t137 + t237;
t117 = -pkin(2) * t130 - pkin(3);
t110 = -pkin(7) + t117;
t242 = t80 * pkin(4);
t224 = t130 * t94;
t54 = t129 * t89 + t224;
t45 = -qJD(2) * qJ(4) - t54;
t30 = -t45 - t242;
t248 = t110 * t51 + t258 * t30;
t61 = -t129 * t96 + t130 * t97;
t247 = -t61 * qJDD(2) - t121 * t250;
t148 = -g(3) * t121 - t122 * t171;
t246 = qJD(2) * (-t80 + t189) + t150;
t243 = t55 * pkin(3);
t240 = pkin(2) * t133;
t115 = g(3) * t122;
t235 = g(3) * t136;
t234 = t122 * pkin(7);
t232 = t64 * t62;
t231 = t64 * t80;
t230 = t80 * t62;
t229 = t80 * t252;
t228 = pkin(4) + t227;
t19 = t129 * t50 + t130 * t57;
t223 = t135 * t22;
t179 = t132 * qJDD(2) - t135 * t55;
t23 = qJD(5) * t64 + t179;
t222 = t23 * t132;
t220 = pkin(6) * qJDD(1);
t219 = qJDD(2) * pkin(3);
t218 = t121 * t137;
t217 = t122 * t134;
t216 = t122 * t137;
t214 = t134 * t132;
t213 = t134 * t135;
t212 = t137 * t227;
t211 = t137 * t132;
t210 = t137 * t135;
t58 = -t129 * t93 + t224;
t209 = t58 * qJD(2);
t207 = t241 - t208;
t127 = t133 ^ 2;
t128 = t136 ^ 2;
t205 = t127 - t128;
t204 = t127 + t128;
t202 = qJD(2) * t133;
t195 = qJDD(2) * qJ(4) + t19;
t120 = pkin(2) * t202;
t90 = t129 * t133 - t215;
t194 = qJD(5) * t132 * t90;
t193 = t90 * t201;
t139 = qJD(1) ^ 2;
t192 = t133 * t139 * t136;
t190 = -g(1) * t218 - t121 * t237 + t115;
t37 = t129 * t75 - t130 * t76;
t60 = t129 * t97 + t130 * t96;
t183 = pkin(2) * t203 + qJ(4) * t80;
t106 = t137 * t118;
t178 = g(2) * (pkin(3) * t216 + qJ(4) * t218 + t106);
t177 = t133 * t187;
t176 = g(1) * t217 - g(2) * t216;
t175 = t233 - t251;
t1 = qJD(5) * t6 + t132 * t9 + t135 * t5;
t173 = t1 - t261;
t172 = -pkin(3) * t121 - t240;
t16 = -qJD(2) * qJD(4) - t195;
t10 = -pkin(4) * t55 - t16;
t168 = t10 * t90 + t30 * t82;
t167 = -t132 * t6 + t135 * t7;
t166 = t258 * t82 + t51 * t90;
t165 = t55 * t90 + t80 * t82;
t85 = qJD(2) * t215 - t129 * t202;
t164 = t252 * t85 + t56 * t91;
t163 = -t190 - t226;
t38 = t129 * t76 + t130 * t75;
t159 = -t91 * qJ(4) - t118;
t32 = t244 * t90 + t159;
t39 = pkin(4) * t91 + t60;
t15 = t132 * t39 + t135 * t32;
t14 = -t132 * t32 + t135 * t39;
t161 = qJD(2) * t82 + qJDD(2) * t90;
t160 = qJD(2) * t85 + qJDD(2) * t91;
t157 = -t118 + t251;
t156 = -t85 * qJ(4) - t91 * qJD(4) + t120;
t155 = -0.2e1 * pkin(1) * t198 - pkin(6) * qJDD(2);
t153 = -t60 * qJDD(2) + t176;
t152 = -t132 * t51 - t253 * t258;
t35 = t80 * pkin(3) + t149;
t151 = t252 * t35 + qJDD(4) - t163;
t147 = -t252 * t82 - t55 * t91 - t56 * t90 - t80 * t85;
t138 = qJD(2) ^ 2;
t146 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t138 + t250;
t145 = pkin(1) * t139 + t171 - t220;
t144 = -t250 + t69;
t142 = -qJD(5) * t110 * t258 + t10 + t148;
t141 = t252 * t37 - t38 * t80 - t55 * t61 + t56 * t60 - t171;
t31 = (t80 + t189) * qJD(2) + t150;
t113 = pkin(2) * t129 + qJ(4);
t100 = qJ(4) * t216;
t98 = qJ(4) * t217;
t74 = -t121 * t214 + t210;
t73 = t121 * t213 + t211;
t72 = t121 * t211 + t213;
t71 = t121 * t210 - t214;
t52 = pkin(3) * t90 + t159;
t41 = -qJD(2) * pkin(3) + t174;
t40 = -pkin(4) * t90 + t61;
t36 = pkin(3) * t252 + t183;
t33 = t58 - t242;
t29 = pkin(3) * t82 + t156;
t28 = -pkin(4) * t82 + t38;
t27 = pkin(4) * t85 + t37;
t26 = t244 * t252 + t183;
t21 = t135 * t23;
t20 = t244 * t82 + t156;
t17 = t191 - t219;
t13 = t143 + t243;
t12 = t132 * t33 + t135 * t26;
t11 = -t132 * t26 + t135 * t33;
t4 = -qJD(5) * t15 - t132 * t20 + t135 * t27;
t3 = qJD(5) * t14 + t132 * t27 + t135 * t20;
t8 = [0, 0, 0, 0, 0, qJDD(1), t250, t171, 0, 0, qJDD(1) * t127 + 0.2e1 * t177, 0.2e1 * t133 * t196 - 0.2e1 * t198 * t205, qJDD(2) * t133 + t136 * t138, qJDD(1) * t128 - 0.2e1 * t177, qJDD(2) * t136 - t133 * t138, 0, t133 * t155 + t136 * t146, -t133 * t146 + t136 * t155, 0.2e1 * t204 * t220 - t171, -g(1) * (-pkin(1) * t134 + pkin(6) * t137) - g(2) * (pkin(1) * t137 + pkin(6) * t134) + (pkin(6) ^ 2 * t204 + pkin(1) ^ 2) * qJDD(1), t164, t147, t160, t165, -t161, 0, -t118 * t55 + t69 * t90 + t95 * t82 + (t240 * t80 - t37) * qJD(2) + t153, -t118 * t56 + t69 * t91 + t95 * t85 + (t240 * t252 - t38) * qJD(2) + t247, -t19 * t90 + t226 * t91 - t53 * t85 - t54 * t82 + t141, t19 * t61 + t54 * t38 + t226 * t60 - t53 * t37 - t69 * t118 + t95 * t120 - g(1) * (-t134 * t118 + t212) - g(2) * (t134 * t227 + t106), 0, -t160, t161, t164, t147, t165, t16 * t90 + t17 * t91 + t41 * t85 + t45 * t82 + t141, qJD(2) * t37 - t13 * t90 - t29 * t80 - t35 * t82 - t52 * t55 - t153, t38 * qJD(2) - t13 * t91 - t252 * t29 - t35 * t85 - t52 * t56 - t247, t13 * t52 + t35 * t29 - t16 * t61 - t45 * t38 + t17 * t60 + t41 * t37 - g(1) * t212 - t178 + (-g(1) * t157 - g(2) * t227) * t134, t64 * t193 + (-t22 * t90 + t64 * t82) * t132, (-t132 * t62 + t135 * t64) * t82 + (-t222 - t223 + (-t132 * t64 - t135 * t62) * qJD(5)) * t90, t132 * t166 + t193 * t258 - t22 * t91 + t64 * t85, t62 * t194 + (-t23 * t90 - t62 * t82) * t135, t135 * t166 - t194 * t258 - t23 * t91 - t62 * t85, t258 * t85 + t51 * t91, -g(1) * t74 - g(2) * t72 - t135 * t168 + t14 * t51 + t194 * t30 + t2 * t91 + t40 * t23 + t258 * t4 + t28 * t62 + t6 * t85, g(1) * t73 - g(2) * t71 - t1 * t91 + t132 * t168 - t15 * t51 + t193 * t30 - t40 * t22 - t258 * t3 + t28 * t64 - t7 * t85, t14 * t22 - t15 * t23 - t3 * t62 - t4 * t64 + t167 * t82 + (t1 * t135 - t2 * t132 + (-t132 * t7 - t135 * t6) * qJD(5)) * t90 + t176, t1 * t15 + t7 * t3 + t2 * t14 + t6 * t4 + t10 * t40 + t30 * t28 - t178 + (-g(1) * t228 - g(2) * t234) * t137 + (-g(1) * (t157 - t234) - g(2) * t228) * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, t205 * t139, t197, t192, t196, qJDD(2), t133 * t145 - t235, g(3) * t133 + t136 * t145, 0, 0, t229, t254, t31, -t229, -t162, qJDD(2), t209 - t95 * t252 + (qJDD(2) * t130 - t203 * t80) * pkin(2) + t163, t59 * qJD(2) + t95 * t80 + (-qJDD(2) * t129 - t203 * t252) * pkin(2) - t19 - t148, (t54 - t58) * t252 + (-t53 + t59) * t80 + (-t129 * t55 - t130 * t56) * pkin(2), t53 * t58 - t54 * t59 + (-t235 + t129 * t19 - t130 * t226 + (-qJD(1) * t95 + t171) * t133) * pkin(2), qJDD(2), -t31, t162, t229, t254, -t229, -t113 * t55 + t117 * t56 + (-t45 - t58) * t252 + (t41 + t208) * t80, -t209 + t36 * t80 + (-pkin(3) + t117) * qJDD(2) + t151, t113 * qJDD(2) - t35 * t80 + t36 * t252 + (0.2e1 * qJD(4) - t59) * qJD(2) + t148 + t195, -t16 * t113 + t17 * t117 - t35 * t36 - t41 * t58 - g(1) * (t137 * t172 + t100) - g(2) * (t134 * t172 + t98) - g(3) * t175 + t208 * t45, -t182 * t64 - t223, -t21 - t260 + (t22 + t184) * t132, t158 + t231, t253 * t62 + t222, t152 - t230, t258 * t80, -t11 * t258 + t113 * t23 + t132 * t142 + t135 * t248 + t207 * t62 + t6 * t80, -t113 * t22 + t12 * t258 - t132 * t248 + t135 * t142 + t207 * t64 - t7 * t80, t11 * t64 + t12 * t62 + (t110 * t22 - t7 * t252 - t2 + (-t110 * t62 - t7) * qJD(5)) * t135 + (-t110 * t23 + t6 * t252 - t1 + (t110 * t64 + t6) * qJD(5)) * t132 - t190, t10 * t113 - t7 * t12 - t6 * t11 - g(1) * (-t137 * t240 + t100) - g(2) * (-t134 * t240 + t98) - g(3) * (t175 + t234) + t207 * t30 + (qJD(5) * t167 + t1 * t132 + t2 * t135) * t110 + t171 * t121 * t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257, t246, t255, t252 * t53 + t54 * t80 + t144, 0, 0, 0, 0, 0, 0, t255, -t257, -t246, t243 - t221 - t45 * t80 + (-qJD(4) - t41) * t252 + t144, 0, 0, 0, 0, 0, 0, t152 + t230, t231 - t158, -t132 * t259 - t21 + t260, -t132 * t249 + t135 * t173 + t30 * t80 - t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, qJDD(2) - t229, -t79 - t138, qJD(2) * t45 + t151 - t219, 0, 0, 0, 0, 0, 0, -qJD(2) * t62 + t158, -qJD(2) * t64 + t152, t259 * t135 + (t258 * t64 - t23) * t132, -t30 * qJD(2) + t132 * t173 + t135 * t249 + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, -t62 ^ 2 + t64 ^ 2, -t259, -t232, -t179 + (-qJD(5) + t258) * t64, t51, -g(1) * t71 - g(2) * t73 + t115 * t135 - t30 * t64 + t249, g(1) * t72 - g(2) * t74 + t30 * t62 + t261 + (-qJD(5) * t25 - t5) * t135 + (qJD(5) * t24 - t115 - t9) * t132, 0, 0;];
tau_reg = t8;
