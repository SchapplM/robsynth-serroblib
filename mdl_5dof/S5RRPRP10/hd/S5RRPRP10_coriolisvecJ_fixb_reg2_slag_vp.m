% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:06
% EndTime: 2019-12-31 20:11:13
% DurationCPUTime: 2.41s
% Computational Cost: add. (2246->298), mult. (5205->401), div. (0->0), fcn. (2892->4), ass. (0->173)
t123 = sin(qJ(2));
t172 = qJD(1) * t123;
t108 = qJD(4) + t172;
t124 = cos(qJ(4));
t122 = sin(qJ(4));
t170 = qJD(2) * t122;
t125 = cos(qJ(2));
t171 = qJD(1) * t125;
t77 = t124 * t171 + t170;
t141 = t77 * t108;
t161 = qJD(1) * qJD(2);
t153 = t123 * t161;
t45 = qJD(4) * t77 - t122 * t153;
t209 = t45 - t141;
t208 = t45 + t141;
t198 = pkin(3) + pkin(6);
t165 = qJD(4) * t124;
t166 = qJD(4) * t122;
t107 = pkin(2) * t153;
t144 = pkin(7) * t123 - qJ(3) * t125;
t163 = t123 * qJD(3);
t131 = qJD(2) * t144 - t163;
t42 = qJD(1) * t131 + t107;
t126 = -pkin(2) - pkin(7);
t152 = -qJ(3) * t123 - pkin(1);
t74 = t125 * t126 + t152;
t53 = t74 * qJD(1);
t111 = pkin(6) * t172;
t201 = qJD(3) + t111;
t174 = pkin(3) * t172 + t201;
t56 = qJD(2) * t126 + t174;
t110 = t125 * t161;
t106 = pkin(6) * t110;
t73 = pkin(3) * t110 + t106;
t150 = -t122 * t73 - t124 * t42 - t165 * t56 + t166 * t53;
t26 = -t122 * t53 + t124 * t56;
t207 = -t108 * t26 - t150;
t27 = t122 * t56 + t124 * t53;
t8 = -qJD(4) * t27 - t122 * t42 + t124 * t73;
t130 = t45 * qJ(5) + t8;
t146 = pkin(4) * t110;
t154 = t122 * t171;
t168 = qJD(2) * t124;
t79 = -t154 + t168;
t1 = -t79 * qJD(5) + t130 + t146;
t19 = -qJ(5) * t79 + t26;
t16 = pkin(4) * t108 + t19;
t20 = -qJ(5) * t77 + t27;
t192 = t108 * t20;
t46 = qJD(2) * t165 - qJD(4) * t154 - t124 * t153;
t137 = t46 * qJ(5) + t150;
t3 = -qJD(5) * t77 - t137;
t206 = -(t108 * t16 - t3) * t122 + (t1 + t192) * t124;
t205 = -0.2e1 * t161;
t203 = t27 * t108 + t8;
t190 = t108 * t79;
t30 = -t46 + t190;
t202 = t46 + t190;
t199 = t79 ^ 2;
t197 = t79 * t77;
t196 = t16 - t19;
t162 = t124 * qJD(5);
t175 = qJ(5) - t126;
t181 = t122 * t123;
t115 = pkin(2) * t172;
t62 = qJD(1) * t144 + t115;
t112 = pkin(6) * t171;
t86 = pkin(3) * t171 + t112;
t35 = -t122 * t62 + t124 * t86;
t195 = (pkin(4) * t125 - qJ(5) * t181) * qJD(1) + t35 - t175 * t166 + t162;
t36 = t122 * t86 + t124 * t62;
t90 = t175 * t124;
t194 = qJ(5) * t124 * t172 + qJD(4) * t90 + qJD(5) * t122 + t36;
t97 = t198 * t123;
t40 = t122 * t97 + t124 * t74;
t193 = qJD(2) * pkin(2);
t189 = t124 * t45;
t188 = t125 * t79;
t186 = t46 * t122;
t118 = qJD(2) * qJD(3);
t169 = qJD(2) * t123;
t85 = t198 * t169;
t59 = -qJD(1) * t85 + t118;
t185 = t59 * t122;
t184 = t59 * t124;
t159 = -pkin(4) * t124 - pkin(3);
t183 = pkin(4) * t165 - t159 * t172 + t201;
t92 = -pkin(2) * t125 + t152;
t68 = qJD(1) * t92;
t182 = t108 * t126;
t180 = t123 * t124;
t179 = t124 * t125;
t128 = qJD(1) ^ 2;
t178 = t125 * t128;
t127 = qJD(2) ^ 2;
t177 = t127 * t123;
t176 = t127 * t125;
t98 = t198 * t125;
t120 = t123 ^ 2;
t121 = t125 ^ 2;
t173 = t120 - t121;
t119 = qJD(2) * qJ(3);
t167 = qJD(2) * t125;
t164 = qJD(4) * t125;
t160 = t108 * t180;
t67 = t119 + t86;
t158 = t122 * t164;
t157 = t108 * t165;
t156 = t124 * t164;
t155 = t108 * t171;
t151 = -pkin(4) * t77 - qJD(5);
t149 = qJ(5) * t125 - t74;
t148 = pkin(1) * t205;
t147 = qJD(3) - t193;
t145 = t123 * t110;
t143 = -t122 * t26 + t124 * t27;
t142 = -0.2e1 * qJD(2) * t68;
t140 = -qJD(1) * t121 + t108 * t123;
t139 = t108 * t122;
t133 = -qJ(3) * t167 - t163;
t50 = qJD(1) * t133 + t107;
t114 = pkin(2) * t169;
t64 = t114 + t133;
t138 = pkin(6) * t127 + qJD(1) * t64 + t50;
t31 = t46 * pkin(4) + t59;
t38 = -t151 + t67;
t136 = -t122 * t31 - t165 * t38;
t135 = t31 * t124 - t166 * t38;
t134 = t123 * t67 + t126 * t167;
t48 = t114 + t131;
t87 = t198 * t167;
t12 = t122 * t87 + t124 * t48 + t165 * t97 - t166 * t74;
t88 = pkin(6) * t153 - t118;
t91 = t111 + t147;
t94 = -t112 - t119;
t129 = -t88 * t125 + (t125 * t91 + (t94 + t112) * t123) * qJD(2);
t109 = pkin(4) * t122 + qJ(3);
t105 = t123 * t178;
t100 = t124 * t110;
t96 = -0.2e1 * t145;
t95 = 0.2e1 * t145;
t93 = t173 * t128;
t89 = t175 * t122;
t83 = -qJ(3) * t171 + t115;
t82 = t124 * t97;
t76 = t77 ^ 2;
t72 = t124 * t87;
t66 = t173 * t205;
t65 = pkin(4) * t179 + t98;
t55 = t68 * t172;
t52 = (t108 + t172) * t167;
t41 = -pkin(4) * t158 + (-pkin(6) + t159) * t169;
t39 = -t122 * t74 + t82;
t34 = -qJ(5) * t179 + t40;
t33 = -t76 + t199;
t32 = t123 * pkin(4) + t122 * t149 + t82;
t25 = -t157 - qJD(2) * t79 + (-t122 * t167 - t160) * qJD(1);
t24 = -qJD(2) * t77 - t108 * t139 + t100;
t22 = -t157 + (-t160 + (t77 - t170) * t125) * qJD(1);
t21 = -t108 * t166 + t100 + (-t108 * t181 - t188) * qJD(1);
t18 = t124 * t141 + t186;
t17 = -t139 * t79 - t189;
t15 = -t77 * t158 + (t125 * t46 - t169 * t77) * t124;
t14 = -t79 * t156 + (t125 * t45 + t169 * t79) * t122;
t13 = -qJD(4) * t40 - t122 * t48 + t72;
t11 = -t108 * t156 - t45 * t123 + (t122 * t140 + t188) * qJD(2);
t10 = t108 * t158 - t46 * t123 + (t124 * t140 - t125 * t77) * qJD(2);
t9 = -t125 * t162 + (t123 * t168 + t158) * qJ(5) + t12;
t6 = pkin(4) * t167 + t72 + t149 * t165 + (-qJ(5) * t169 - qJD(4) * t97 + qJD(5) * t125 - t48) * t122;
t5 = t30 * t122 + t124 * t209;
t4 = t122 * t208 - t202 * t124;
t2 = (-t122 * t77 + t124 * t79) * t169 + (t186 + t189 + (t122 * t79 + t124 * t77) * qJD(4)) * t125;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t66, t176, t96, -t177, 0, -pkin(6) * t176 + t123 * t148, pkin(6) * t177 + t125 * t148, 0, 0, 0, -t176, t177, t95, t66, t96, t129, t123 * t142 + t125 * t138, -t123 * t138 + t125 * t142, pkin(6) * t129 + t50 * t92 + t68 * t64, t14, t2, t11, t15, t10, t52, t13 * t108 + t98 * t46 - t85 * t77 + (-t168 * t67 + t8) * t123 + (-t67 * t166 + t184 + (qJD(1) * t39 + t26) * qJD(2)) * t125, -t12 * t108 - t98 * t45 - t85 * t79 + (t170 * t67 + t150) * t123 + (-t67 * t165 - t185 + (-qJD(1) * t40 - t27) * qJD(2)) * t125, -t12 * t77 - t13 * t79 + t39 * t45 - t40 * t46 + t143 * t169 + (t122 * t8 + t124 * t150 + (t122 * t27 + t124 * t26) * qJD(4)) * t125, t12 * t27 + t13 * t26 - t150 * t40 + t39 * t8 + t59 * t98 - t67 * t85, t14, t2, t11, t15, t10, t52, t6 * t108 + t41 * t77 + t65 * t46 + (-t168 * t38 + t1) * t123 + ((qJD(1) * t32 + t16) * qJD(2) + t135) * t125, -t9 * t108 + t41 * t79 - t65 * t45 + (t170 * t38 - t3) * t123 + ((-qJD(1) * t34 - t20) * qJD(2) + t136) * t125, t32 * t45 - t34 * t46 - t6 * t79 - t77 * t9 + (-t122 * t16 + t124 * t20) * t169 + (t1 * t122 - t124 * t3 + (t122 * t20 + t124 * t16) * qJD(4)) * t125, t1 * t32 + t16 * t6 + t20 * t9 + t3 * t34 + t31 * t65 + t38 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t93, 0, t105, 0, 0, t128 * pkin(1) * t123, pkin(1) * t178, 0, 0, 0, 0, 0, -t105, t93, t105, ((-t94 - t119) * t123 + (t147 - t91) * t125) * qJD(1), -t171 * t83 + t55, 0.2e1 * t118 + (t123 * t83 + t125 * t68) * qJD(1), -t88 * qJ(3) - t94 * qJD(3) - t68 * t83 + (-t123 * t94 + (-t91 - t193) * t125) * qJD(1) * pkin(6), t17, t4, t21, t18, t22, -t155, qJ(3) * t46 - t108 * t35 + t185 + t174 * t77 + (-t122 * t182 + t124 * t67) * qJD(4) + (t124 * t134 - t125 * t26) * qJD(1), -qJ(3) * t45 + t108 * t36 + t184 + t174 * t79 + (-t122 * t67 - t124 * t182) * qJD(4) + (-t122 * t134 + t125 * t27) * qJD(1), t35 * t79 + t36 * t77 + (-t27 * t172 + t126 * t45 - t8 + (-t126 * t77 - t27) * qJD(4)) * t124 + (t26 * t172 - t126 * t46 + t150 + (t126 * t79 + t26) * qJD(4)) * t122, t59 * qJ(3) - t26 * t35 - t27 * t36 + t174 * t67 + (qJD(4) * t143 - t122 * t150 + t8 * t124) * t126, t17, t4, t21, t18, t22, -t155, t109 * t46 + t183 * t77 - t195 * t108 + (t38 * t180 + (-qJD(2) * t90 - t16) * t125) * qJD(1) - t136, -t109 * t45 + t183 * t79 + t194 * t108 + (-t38 * t181 + (qJD(2) * t89 + t20) * t125) * qJD(1) + t135, t194 * t77 + t195 * t79 - t45 * t90 + t46 * t89 - t206, -t1 * t90 + t31 * t109 - t16 * t195 + t183 * t38 - t194 * t20 - t3 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t120 * t128 - t127, qJD(2) * t94 + t106 + t55, 0, 0, 0, 0, 0, 0, t24, t25, t5, -t67 * qJD(2) + t122 * t207 + t203 * t124, 0, 0, 0, 0, 0, 0, t24, t25, t5, -t38 * qJD(2) + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t33, -t209, -t197, t30, t110, -t67 * t79 + t203, t67 * t77 - t207, 0, 0, t197, t33, -t209, -t197, t30, t110, 0.2e1 * t146 + t192 + (t151 - t38) * t79 + t130, -t199 * pkin(4) + t19 * t108 + (qJD(5) + t38) * t77 + t137, t45 * pkin(4) - t196 * t77, t196 * t20 + (-t38 * t79 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, -t208, -t76 - t199, t16 * t79 + t20 * t77 + t31;];
tauc_reg = t7;
