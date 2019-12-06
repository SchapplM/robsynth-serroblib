% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:46
% EndTime: 2019-12-05 17:16:56
% DurationCPUTime: 2.93s
% Computational Cost: add. (2361->297), mult. (5460->442), div. (0->0), fcn. (4412->14), ass. (0->178)
t114 = qJD(3) + qJD(4);
t122 = sin(qJ(4));
t126 = cos(qJ(3));
t227 = cos(qJ(4));
t176 = qJD(2) * t227;
t123 = sin(qJ(3));
t199 = qJD(2) * t123;
t238 = -t122 * t199 + t126 * t176;
t239 = t238 * t114;
t71 = qJD(5) - t238;
t237 = t71 - qJD(5);
t231 = pkin(7) + pkin(8);
t117 = qJ(3) + qJ(4);
t111 = sin(t117);
t112 = cos(t117);
t120 = cos(pkin(5));
t119 = sin(pkin(5));
t212 = cos(pkin(10));
t167 = t119 * t212;
t124 = sin(qJ(2));
t207 = t119 * t124;
t118 = sin(pkin(10));
t210 = t118 * t119;
t166 = t212 * t124;
t127 = cos(qJ(2));
t208 = t118 * t127;
t73 = t120 * t166 + t208;
t165 = t212 * t127;
t209 = t118 * t124;
t75 = -t120 * t209 + t165;
t143 = -g(3) * (-t111 * t207 + t120 * t112) - g(2) * (-t73 * t111 - t112 * t167) - g(1) * (-t75 * t111 + t112 * t210);
t113 = qJDD(3) + qJDD(4);
t194 = qJD(1) * qJD(2);
t67 = qJDD(2) * pkin(7) + (qJDD(1) * t124 + t127 * t194) * t119;
t164 = pkin(8) * qJDD(2) + t67;
t201 = qJD(1) * t119;
t181 = t124 * t201;
t171 = t231 * qJD(2) + t181;
t200 = qJD(1) * t120;
t57 = t123 * t200 + t171 * t126;
t192 = t120 * qJDD(1);
t95 = t126 * t192;
t19 = qJDD(3) * pkin(3) - t57 * qJD(3) - t164 * t123 + t95;
t56 = -t171 * t123 + t126 * t200;
t22 = t56 * qJD(3) + t123 * t192 + t164 * t126;
t188 = t227 * t57;
t216 = qJD(3) * pkin(3);
t51 = t56 + t216;
t27 = t122 * t51 + t188;
t233 = t27 * qJD(4) + t122 * t22 - t227 * t19;
t4 = -t113 * pkin(4) + t233;
t140 = t143 - t4;
t72 = -t120 * t165 + t209;
t74 = t120 * t208 + t166;
t159 = g(1) * t74 + g(2) * t72;
t205 = t119 * t127;
t139 = g(3) * t205 - t159;
t137 = t139 * t112;
t168 = qJDD(2) * t227;
t191 = t123 * qJDD(2);
t154 = t122 * t191 - t126 * t168;
t204 = t122 * t126;
t82 = t227 * t123 + t204;
t53 = t114 * t82;
t36 = qJD(2) * t53 + t154;
t33 = qJDD(5) + t36;
t110 = -t126 * pkin(3) - pkin(2);
t144 = -t122 * t123 + t227 * t126;
t44 = -pkin(4) * t144 - t82 * pkin(9) + t110;
t236 = t44 * t33 - t137;
t108 = t122 * pkin(3) + pkin(9);
t80 = -qJD(2) * t204 - t123 * t176;
t43 = -t80 * pkin(4) - pkin(9) * t238;
t235 = (pkin(3) * t199 + qJD(5) * t108 + t43) * t71;
t141 = t123 * t216 - t181;
t158 = g(1) * t75 + g(2) * t73;
t215 = t122 * t57;
t26 = t227 * t51 - t215;
t20 = -t114 * pkin(4) - t26;
t89 = t231 * t123;
t90 = t231 * t126;
t149 = -t122 * t90 - t227 * t89;
t180 = t127 * t201;
t184 = qJD(3) * t231;
t83 = t123 * t184;
t84 = t126 * t184;
t219 = -t149 * qJD(4) + t122 * t84 + t144 * t180 + t227 * t83;
t175 = qJD(4) * t227;
t197 = qJD(4) * t122;
t135 = t122 * t19 + t51 * t175 - t57 * t197 + t227 * t22;
t3 = t113 * pkin(9) + t135;
t70 = t110 * qJD(2) - t180;
t34 = -pkin(4) * t238 + t80 * pkin(9) + t70;
t52 = t114 * t144;
t64 = -t122 * t89 + t227 * t90;
t234 = (qJD(5) * t34 + t3) * t144 + t20 * t52 + t4 * t82 + (-qJD(5) * t44 + t219) * t71 - t64 * t33 - g(3) * t207 - t158;
t121 = sin(qJ(5));
t125 = cos(qJ(5));
t190 = t126 * qJDD(2);
t35 = t122 * t190 + t123 * t168 + t239;
t62 = t121 * t114 - t125 * t80;
t15 = t62 * qJD(5) - t125 * t113 + t121 * t35;
t128 = qJD(3) ^ 2;
t174 = t124 * t194;
t155 = -qJDD(1) * t205 + t119 * t174;
t232 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t128 + t119 * (-g(3) * t127 + t174) - t155 + t159;
t226 = t20 * t238;
t225 = t20 * t82;
t60 = -t125 * t114 - t121 * t80;
t223 = t60 * t71;
t222 = t62 * t71;
t221 = t71 * t80;
t220 = t80 * t238;
t218 = t64 * qJD(4) - t122 * t83 - t82 * t180 + t227 * t84;
t217 = qJD(2) * pkin(2);
t214 = t125 * t33;
t195 = qJD(5) * t125;
t196 = qJD(5) * t121;
t14 = t121 * t113 + t114 * t195 + t125 * t35 + t80 * t196;
t213 = t14 * t121;
t206 = t119 * t126;
t203 = qJDD(1) - g(3);
t115 = t123 ^ 2;
t202 = -t126 ^ 2 + t115;
t198 = qJD(2) * t124;
t193 = qJD(2) * qJD(3);
t187 = t82 * t196;
t186 = t121 * t205;
t185 = t125 * t205;
t183 = t119 * t198;
t182 = qJD(2) * t205;
t21 = t114 * pkin(9) + t27;
t152 = t121 * t21 - t125 * t34;
t177 = -t152 * t80 + t20 * t196;
t173 = t123 * t193;
t172 = t127 * t193;
t162 = t125 * t71;
t28 = t122 * t56 + t188;
t160 = pkin(3) * t197 - t28;
t157 = t53 * pkin(4) - t52 * pkin(9) + t141;
t156 = t33 * t82 + t52 * t71;
t11 = t121 * t34 + t125 * t21;
t129 = qJD(2) ^ 2;
t151 = qJDD(2) * t127 - t124 * t129;
t76 = t120 * t126 - t123 * t207;
t77 = t120 * t123 + t124 * t206;
t150 = -t122 * t77 + t227 * t76;
t40 = t122 * t76 + t227 * t77;
t148 = -g(1) * t118 + t212 * g(2);
t147 = -t121 * t40 - t185;
t146 = -t125 * t40 + t186;
t142 = -t11 * t80 - t140 * t121 + t20 * t195;
t86 = -t180 - t217;
t138 = -t86 * qJD(2) + t158 - t67;
t45 = pkin(3) * t173 + t110 * qJDD(2) + t155;
t134 = -pkin(7) * qJDD(3) + (t180 + t86 - t217) * qJD(3);
t29 = t227 * t56 - t215;
t133 = -t108 * t33 - t226 + (-pkin(3) * t175 + t29) * t71;
t132 = t70 * t80 + t143 - t233;
t47 = -t111 * t167 + t73 * t112;
t49 = t111 * t210 + t75 * t112;
t69 = t120 * t111 + t112 * t207;
t131 = g(1) * t49 + g(2) * t47 + g(3) * t69 - t238 * t70 - t135;
t109 = -t227 * pkin(3) - pkin(4);
t55 = -qJD(3) * t77 - t123 * t182;
t54 = qJD(3) * t76 + t126 * t182;
t37 = -t238 ^ 2 + t80 ^ 2;
t24 = -t154 + (-qJD(2) * t82 - t80) * t114;
t23 = t35 - t239;
t13 = t40 * qJD(4) + t122 * t54 - t227 * t55;
t12 = t150 * qJD(4) + t122 * t55 + t227 * t54;
t9 = t36 * pkin(4) - t35 * pkin(9) + t45;
t8 = t125 * t9;
t7 = t121 * t33 + t71 * t162 + t62 * t80;
t6 = -t71 ^ 2 * t121 - t60 * t80 + t214;
t5 = t62 * t162 + t213;
t1 = (t14 - t223) * t125 + (-t15 - t222) * t121;
t2 = [t203, 0, t151 * t119, (-qJDD(2) * t124 - t127 * t129) * t119, 0, 0, 0, 0, 0, t55 * qJD(3) + t76 * qJDD(3) + (-t123 * t172 + t151 * t126) * t119, -t54 * qJD(3) - t77 * qJDD(3) + (-t151 * t123 - t126 * t172) * t119, 0, 0, 0, 0, 0, t150 * t113 - t13 * t114 + (-t127 * t36 - t198 * t238) * t119, -t40 * t113 - t12 * t114 + (-t127 * t35 - t80 * t198) * t119, 0, 0, 0, 0, 0, (qJD(5) * t146 - t121 * t12 + t125 * t183) * t71 + t147 * t33 + t13 * t60 - t150 * t15, -(qJD(5) * t147 + t125 * t12 + t121 * t183) * t71 + t146 * t33 + t13 * t62 - t150 * t14; 0, qJDD(2), t203 * t205 + t159, -t203 * t207 + t158, t115 * qJDD(2) + 0.2e1 * t126 * t173, 0.2e1 * t123 * t190 - 0.2e1 * t202 * t193, qJDD(3) * t123 + t128 * t126, qJDD(3) * t126 - t128 * t123, 0, t134 * t123 + t232 * t126, -t232 * t123 + t134 * t126, t35 * t82 - t80 * t52, t144 * t35 + t238 * t52 - t82 * t36 + t80 * t53, t82 * t113 + t52 * t114, t113 * t144 - t53 * t114, 0, t110 * t36 + t113 * t149 - t218 * t114 - t141 * t238 - t144 * t45 + t70 * t53 - t137, t110 * t35 + t139 * t111 - t64 * t113 + t219 * t114 - t141 * t80 + t45 * t82 + t70 * t52, -t62 * t187 + (t14 * t82 + t52 * t62) * t125, (-t121 * t62 - t125 * t60) * t52 + (-t213 - t125 * t15 + (t121 * t60 - t125 * t62) * qJD(5)) * t82, t156 * t125 - t14 * t144 - t71 * t187 + t62 * t53, -t82 * t71 * t195 - t156 * t121 + t144 * t15 - t60 * t53, -t144 * t33 + t71 * t53, -t152 * t53 - t149 * t15 - t8 * t144 + t218 * t60 + (t157 * t71 + (t144 * t21 - t64 * t71 + t225) * qJD(5) + t236) * t125 + t234 * t121, -t11 * t53 - t149 * t14 + t218 * t62 + ((-qJD(5) * t21 + t9) * t144 - qJD(5) * t225 + (qJD(5) * t64 - t157) * t71 - t236) * t121 + t234 * t125; 0, 0, 0, 0, -t123 * t129 * t126, t202 * t129, t191, t190, qJDD(3), -g(3) * t76 + t123 * t138 + t148 * t206 + t95, g(3) * t77 + (-t119 * t148 - t192) * t123 + t138 * t126, t220, t37, t23, t24, t113, t28 * t114 + (t227 * t113 - t114 * t197 + t199 * t238) * pkin(3) + t132, t29 * t114 + (-t113 * t122 - t114 * t175 + t80 * t199) * pkin(3) + t131, t5, t1, t7, t6, t221, t109 * t15 + t160 * t60 + t133 * t121 + (t140 - t235) * t125 + t177, t109 * t14 + t121 * t235 + t125 * t133 + t160 * t62 + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, t37, t23, t24, t113, t27 * t114 + t132, t26 * t114 + t131, t5, t1, t7, t6, t221, -pkin(4) * t15 - t27 * t60 + (-pkin(9) * t33 + t26 * t71 - t226) * t121 + ((-pkin(9) * qJD(5) - t43) * t71 + t140) * t125 + t177, -pkin(4) * t14 + (t121 * t43 + t125 * t26) * t71 - t27 * t62 - t125 * t226 + (t71 * t196 - t214) * pkin(9) + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60, -t60 ^ 2 + t62 ^ 2, t14 + t223, -t15 + t222, t33, -t121 * t3 + t8 - t20 * t62 - g(1) * (-t49 * t121 + t74 * t125) - g(2) * (-t47 * t121 + t72 * t125) - g(3) * (-t69 * t121 - t185) + t237 * t11, -t125 * t3 - t121 * t9 + t20 * t60 - g(1) * (-t74 * t121 - t49 * t125) - g(2) * (-t72 * t121 - t47 * t125) - g(3) * (-t69 * t125 + t186) - t237 * t152;];
tau_reg = t2;
