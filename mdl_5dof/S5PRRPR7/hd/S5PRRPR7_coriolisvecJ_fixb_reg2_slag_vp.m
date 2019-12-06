% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:34
% EndTime: 2019-12-05 16:37:47
% DurationCPUTime: 3.26s
% Computational Cost: add. (2883->345), mult. (7619->532), div. (0->0), fcn. (5606->10), ass. (0->177)
t122 = sin(pkin(10));
t128 = sin(qJ(2));
t130 = cos(qJ(3));
t123 = sin(pkin(5));
t184 = qJD(1) * t123;
t124 = cos(pkin(10));
t131 = cos(qJ(2));
t192 = t124 * t131;
t127 = sin(qJ(3));
t142 = pkin(3) * t127 - qJ(4) * t130;
t79 = qJD(3) * t142 - t127 * qJD(4);
t230 = t122 * t79 - (t122 * t128 + t130 * t192) * t184;
t179 = qJD(3) * t127;
t229 = (-pkin(7) * t124 + pkin(8)) * t179 + t230;
t144 = pkin(4) * t122 - pkin(8) * t124;
t137 = pkin(7) + t144;
t162 = t131 * t184;
t150 = t127 * t162;
t177 = qJD(3) * t130;
t228 = -t137 * t177 + t150;
t181 = qJD(2) * t127;
t125 = cos(pkin(5));
t191 = t125 * t130;
t158 = qJD(1) * t191;
t163 = t128 * t184;
t96 = qJD(2) * pkin(7) + t163;
t86 = t127 * t96;
t68 = -t86 + t158;
t92 = t142 * qJD(2);
t34 = t122 * t92 + t124 * t68;
t227 = pkin(8) * t181 - qJD(4) * t124 + t34;
t147 = qJD(2) * t162;
t201 = qJD(3) * t158 + t130 * t147;
t37 = (qJD(4) - t86) * qJD(3) + t201;
t48 = (t79 + t163) * qJD(2);
t14 = t122 * t48 + t124 * t37;
t172 = qJD(2) * qJD(3);
t157 = t127 * t172;
t12 = pkin(8) * t157 + t14;
t126 = sin(qJ(5));
t129 = cos(qJ(5));
t180 = qJD(2) * t130;
t183 = qJD(1) * t127;
t159 = t125 * t183;
t69 = t130 * t96 + t159;
t59 = qJD(3) * qJ(4) + t69;
t101 = -t130 * pkin(3) - t127 * qJ(4) - pkin(2);
t70 = qJD(2) * t101 - t162;
t21 = t122 * t70 + t124 * t59;
t16 = -pkin(8) * t180 + t21;
t154 = -qJD(3) * pkin(3) + qJD(4);
t56 = t154 - t68;
t117 = t124 * qJD(3);
t88 = t122 * t181 - t117;
t161 = t124 * t181;
t173 = t122 * qJD(3);
t90 = t161 + t173;
t18 = t88 * pkin(4) - t90 * pkin(8) + t56;
t140 = t126 * t16 - t129 * t18;
t38 = t159 + (t144 * qJD(2) + t96) * t130;
t19 = qJD(3) * t38 + t127 * t147;
t1 = -qJD(5) * t140 + t129 * t12 + t126 * t19;
t85 = qJD(5) + t88;
t226 = t140 * t85 + t1;
t6 = t126 * t18 + t129 * t16;
t2 = -qJD(5) * t6 - t126 * t12 + t129 * t19;
t225 = t6 * t85 + t2;
t198 = t122 * t130;
t224 = t162 * t198 + (-t163 + t79) * t124;
t156 = t130 * t172;
t146 = t124 * t156;
t62 = t126 * t90 + t129 * t180;
t28 = qJD(5) * t62 - t126 * t157 - t129 * t146;
t223 = -t62 * t85 + t28;
t160 = t126 * t180;
t174 = qJD(5) * t129;
t29 = -qJD(5) * t160 + t126 * t146 - t129 * t157 + t90 * t174;
t64 = t129 * t90 - t160;
t222 = t64 * t85 - t29;
t120 = t127 ^ 2;
t121 = t130 ^ 2;
t221 = qJD(2) * (t120 - 0.2e1 * t121);
t193 = t124 * t130;
t74 = pkin(7) * t193 + t122 * t101;
t61 = -t130 * pkin(8) + t74;
t76 = t137 * t127;
t25 = t126 * t76 + t129 * t61;
t218 = qJD(5) * t25 + t229 * t126 + t228 * t129;
t24 = -t126 * t61 + t129 * t76;
t217 = -qJD(5) * t24 + t228 * t126 - t229 * t129;
t194 = t124 * t129;
t98 = -t124 * pkin(4) - t122 * pkin(8) - pkin(3);
t72 = qJ(4) * t194 + t126 * t98;
t216 = -qJD(5) * t72 + t227 * t126 - t129 * t38;
t182 = qJD(2) * t123;
t165 = t131 * t182;
t40 = t96 * t177 + (qJD(3) * t125 + t165) * t183;
t197 = t123 * t128;
t80 = t127 * t197 - t191;
t215 = t40 * t80;
t213 = t64 * t62;
t71 = -t126 * t124 * qJ(4) + t129 * t98;
t211 = -qJD(5) * t71 + t126 * t38 + t227 * t129;
t168 = pkin(7) * t122 + pkin(4);
t210 = -t168 * t179 - t224;
t171 = pkin(7) * t179;
t152 = t122 * t171;
t209 = t152 + t224;
t151 = t124 * t171;
t208 = -t151 + t230;
t207 = qJD(2) * pkin(2);
t205 = t124 * t88;
t204 = t40 * t122;
t203 = t40 * t124;
t202 = t40 * t127;
t133 = qJD(2) ^ 2;
t199 = t121 * t133;
t196 = t123 * t133;
t195 = t124 * t101;
t190 = t126 * t127;
t189 = t129 * t130;
t132 = qJD(3) ^ 2;
t188 = t132 * t127;
t187 = t132 * t130;
t185 = t120 - t121;
t178 = qJD(3) * t129;
t175 = qJD(5) * t126;
t170 = t128 * t196;
t169 = t124 * t189;
t167 = t122 * t180;
t166 = t128 * t182;
t164 = t130 * t173;
t155 = t85 ^ 2;
t153 = t88 + t117;
t149 = t127 * t165;
t148 = t130 * t165;
t104 = t122 * t156;
t145 = t127 * t156;
t97 = -t162 - t207;
t143 = -t97 - t162;
t13 = -t122 * t37 + t124 * t48;
t20 = -t122 * t59 + t124 * t70;
t33 = -t122 * t68 + t124 * t92;
t81 = t125 * t127 + t130 * t197;
t53 = -t123 * t131 * t122 + t81 * t124;
t23 = t80 * t126 + t53 * t129;
t22 = -t53 * t126 + t80 * t129;
t139 = qJD(2) * (-t90 + t173);
t138 = qJD(2) * t153;
t83 = t124 * t190 + t189;
t136 = qJD(3) * (-t143 - t207);
t135 = -qJ(4) * t179 + (t154 - t56) * t130;
t39 = -t96 * t179 + t201;
t134 = t202 + t39 * t130 + (-t127 * t69 - t130 * t68) * qJD(3);
t119 = t124 ^ 2;
t118 = t122 ^ 2;
t111 = t127 * t133 * t130;
t103 = -0.2e1 * t145;
t93 = t118 * t145;
t84 = -t126 * t130 + t127 * t194;
t78 = (t169 + t190) * qJD(2);
t77 = t124 * t160 - t129 * t181;
t73 = -pkin(7) * t198 + t195;
t60 = t168 * t130 - t195;
t55 = qJD(3) * t81 + t149;
t54 = -qJD(3) * t80 + t148;
t52 = t81 * t122 + t123 * t192;
t42 = -t127 * t178 - t130 * t175 + (t126 * t177 + t127 * t174) * t124;
t41 = -qJD(3) * t169 + qJD(5) * t83 - t126 * t179;
t32 = t122 * t166 + t54 * t124;
t31 = t54 * t122 - t124 * t166;
t26 = -pkin(4) * t181 - t33;
t15 = pkin(4) * t180 - t20;
t11 = -pkin(4) * t157 - t13;
t4 = qJD(5) * t22 + t55 * t126 + t32 * t129;
t3 = -qJD(5) * t23 - t32 * t126 + t55 * t129;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, -t131 * t196, 0, 0, 0, 0, 0, 0, 0, 0, -t130 * t170 + (-t55 - t149) * qJD(3), t127 * t170 + (-t54 - t148) * qJD(3), (t127 * t55 + t130 * t54 + (-t127 * t81 + t130 * t80) * qJD(3)) * qJD(2), t39 * t81 + t215 + t69 * t54 - t68 * t55 + (t97 - t162) * t166, 0, 0, 0, 0, 0, 0, t55 * t88 + (t130 * t31 + (-t127 * t52 + t80 * t198) * qJD(3)) * qJD(2), t55 * t90 + (t130 * t32 + (-t127 * t53 + t80 * t193) * qJD(3)) * qJD(2), t31 * t90 - t32 * t88 + (-t122 * t53 + t124 * t52) * t156, -t13 * t52 + t14 * t53 - t20 * t31 + t21 * t32 + t55 * t56 + t215, 0, 0, 0, 0, 0, 0, t104 * t22 + t52 * t29 + t3 * t85 + t31 * t62, -t104 * t23 - t52 * t28 + t31 * t64 - t4 * t85, t22 * t28 - t23 * t29 - t3 * t64 - t4 * t62, t1 * t23 + t11 * t52 - t140 * t3 + t15 * t31 + t2 * t22 + t4 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t145, -0.2e1 * t185 * t172, t187, t103, -t188, 0, -pkin(7) * t187 + t127 * t136, pkin(7) * t188 + t130 * t136, (-t120 - t121) * t147 + t134, ((t127 * t68 - t130 * t69) * t131 + (-t97 - t207) * t128) * t184 + t134 * pkin(7), (t119 * t181 + t124 * t90) * t177, (-t205 + (-t90 - 0.2e1 * t161) * t122) * t177, (t124 * t221 + t127 * t90) * qJD(3), t88 * t164 + t93, (-t122 * t221 - t127 * t88) * qJD(3), t103, (-t88 * t162 + t204 + (qJD(2) * t73 + t20) * qJD(3)) * t127 + (-t13 + (pkin(7) * t88 + t122 * t56) * qJD(3) + (t152 - t209) * qJD(2)) * t130, (-t90 * t162 + t203 + (-qJD(2) * t74 - t21) * qJD(3)) * t127 + (t14 + (pkin(7) * t90 + t124 * t56) * qJD(3) + (t151 + t208) * qJD(2)) * t130, -t209 * t90 - t208 * t88 + (-t122 * t14 - t124 * t13) * t127 + (-t122 * t21 - t124 * t20 + (-t122 * t74 - t124 * t73) * qJD(2)) * t177, -t56 * t150 + t13 * t73 + t14 * t74 + t208 * t21 + t209 * t20 + (t177 * t56 + t202) * pkin(7), -t28 * t84 - t41 * t64, t28 * t83 - t29 * t84 + t41 * t62 - t42 * t64, -t41 * t85 + (-t127 * t28 + (qJD(2) * t84 + t64) * t177) * t122, t29 * t83 + t42 * t62, -t42 * t85 + (-t127 * t29 + (-qJD(2) * t83 - t62) * t177) * t122, t164 * t85 + t93, t11 * t83 + t15 * t42 + t60 * t29 - t218 * t85 + t210 * t62 + (t127 * t2 + (qJD(2) * t24 - t140) * t177) * t122, t11 * t84 - t15 * t41 - t60 * t28 + t217 * t85 + t210 * t64 + (-t1 * t127 + (-qJD(2) * t25 - t6) * t177) * t122, -t1 * t83 - t140 * t41 - t2 * t84 + t217 * t62 + t218 * t64 + t24 * t28 - t25 * t29 - t42 * t6, t1 * t25 + t11 * t60 + t140 * t218 + t210 * t15 + t2 * t24 - t217 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, t185 * t133, 0, t111, 0, 0, t143 * t181, -t97 * t180 + (t68 + t86) * qJD(3) - t201, 0, 0, t139 * t193, (t122 * t90 + t205 + (-t118 + t119) * qJD(3)) * t180, t124 * t199 + t127 * t139, -t153 * t167, -t122 * t199 + t127 * t138, t111, -t203 - t69 * t88 + (t122 * t135 - t127 * t20 + t130 * t33) * qJD(2), t204 - t69 * t90 + (t124 * t135 + t127 * t21 - t130 * t34) * qJD(2), t33 * t90 + t34 * t88 + (-qJD(4) * t88 + t180 * t20 + t14) * t124 + (qJD(4) * t90 + t180 * t21 - t13) * t122, -t40 * pkin(3) - t20 * t33 - t21 * t34 - t56 * t69 + (-t122 * t20 + t124 * t21) * qJD(4) + (-t13 * t122 + t14 * t124) * qJ(4), -t64 * t78 + (-t129 * t28 - t175 * t64) * t122, t78 * t62 + t64 * t77 + (t126 * t28 - t129 * t29 + (t126 * t62 - t129 * t64) * qJD(5)) * t122, t28 * t124 + (-t122 * t175 - t78) * t85 + (t118 * t178 - t122 * t64) * t180, -t62 * t77 + (t126 * t29 + t174 * t62) * t122, t29 * t124 + (-t122 * t174 + t77) * t85 + (-qJD(3) * t118 * t126 + t122 * t62) * t180, (-t85 - t117) * t167, -t2 * t124 - t15 * t77 - t26 * t62 + t216 * t85 + (t15 * t174 + qJ(4) * t29 + qJD(4) * t62 + t11 * t126 + (qJD(3) * t71 + t140) * t180) * t122, t1 * t124 - t15 * t78 - t26 * t64 + t211 * t85 + (-t15 * t175 - qJ(4) * t28 + qJD(4) * t64 + t11 * t129 + (-qJD(3) * t72 + t6) * t180) * t122, t71 * t28 - t72 * t29 - t140 * t78 + t6 * t77 - t216 * t64 + t211 * t62 + (-t1 * t126 - t2 * t129 + (-t126 * t140 - t129 * t6) * qJD(5)) * t122, t1 * t72 - t15 * t26 + t2 * t71 - t211 * t6 - t216 * t140 + (qJ(4) * t11 + qJD(4) * t15) * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90 * t180 + t104, t130 * t138, -t88 ^ 2 - t90 ^ 2, t20 * t90 + t21 * t88 + t40, 0, 0, 0, 0, 0, 0, t104 * t129 - t126 * t155 - t90 * t62, -t104 * t126 - t129 * t155 - t90 * t64, t222 * t126 + t223 * t129, t226 * t126 + t225 * t129 - t15 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, -t62 ^ 2 + t64 ^ 2, -t223, -t213, t222, t104, -t15 * t64 + t225, t15 * t62 - t226, 0, 0;];
tauc_reg = t5;
