% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:51
% EndTime: 2019-12-31 18:19:56
% DurationCPUTime: 1.39s
% Computational Cost: add. (2732->252), mult. (6564->347), div. (0->0), fcn. (4431->8), ass. (0->151)
t147 = cos(pkin(9));
t135 = qJ(4) * qJD(1);
t83 = sin(pkin(8)) * pkin(1) + pkin(6);
t76 = t83 * qJD(1);
t119 = t76 + t135;
t96 = sin(qJ(3));
t139 = t96 * qJD(2);
t98 = cos(qJ(3));
t51 = t119 * t98 + t139;
t43 = t147 * t51;
t89 = t98 * qJD(2);
t50 = -t119 * t96 + t89;
t45 = qJD(3) * pkin(3) + t50;
t92 = sin(pkin(9));
t18 = t92 * t45 + t43;
t16 = qJD(3) * pkin(7) + t18;
t125 = t147 * t96;
t74 = t92 * t98 + t125;
t146 = qJD(1) * t74;
t85 = -cos(pkin(8)) * pkin(1) - pkin(2);
t75 = -pkin(3) * t98 + t85;
t145 = qJD(1) * t75;
t65 = qJD(4) + t145;
t144 = qJD(1) * t96;
t124 = t147 * t98;
t79 = qJD(1) * t124;
t66 = t92 * t144 - t79;
t24 = pkin(4) * t66 - pkin(7) * t146 + t65;
t95 = sin(qJ(5));
t97 = cos(qJ(5));
t114 = t16 * t95 - t24 * t97;
t136 = t98 * qJD(4);
t143 = qJD(3) * t96;
t87 = qJD(3) * t89;
t55 = -t76 * t143 + t87;
t38 = (-qJ(4) * t143 + t136) * qJD(1) + t55;
t126 = t147 * t38;
t138 = t96 * qJD(4);
t181 = -qJD(1) * t138 - t51 * qJD(3);
t12 = t181 * t92 + t126;
t68 = t74 * qJD(3);
t57 = qJD(1) * t68;
t134 = qJD(1) * qJD(3);
t128 = t96 * t134;
t78 = t92 * t128;
t58 = qJD(3) * t79 - t78;
t80 = pkin(3) * t128;
t28 = pkin(4) * t57 - pkin(7) * t58 + t80;
t1 = -t114 * qJD(5) + t97 * t12 + t95 * t28;
t63 = qJD(5) + t66;
t184 = t114 * t63 + t1;
t6 = t16 * t97 + t24 * t95;
t2 = -qJD(5) * t6 - t95 * t12 + t97 * t28;
t183 = t6 * t63 + t2;
t123 = t63 * t95;
t49 = qJD(3) * t95 + t146 * t97;
t182 = t49 * t123;
t142 = qJD(5) * t95;
t158 = t92 * t96;
t107 = t124 - t158;
t71 = t107 * qJD(3);
t162 = t71 * t97;
t108 = t74 * t142 - t162;
t53 = t97 * t57;
t180 = -t108 * t63 + t74 * t53;
t179 = t146 ^ 2;
t178 = pkin(3) * t96;
t11 = -t147 * t181 + t38 * t92;
t148 = qJ(4) + t83;
t72 = t148 * t98;
t35 = t148 * t125 + t72 * t92;
t175 = t11 * t35;
t174 = t11 * t107;
t173 = t11 * t95;
t137 = t97 * qJD(3);
t26 = -qJD(5) * t137 + t142 * t146 - t97 * t58;
t172 = t26 * t95;
t47 = t146 * t95 - t137;
t171 = t47 * t66;
t170 = t47 * t146;
t169 = t47 * t95;
t168 = t49 * t47;
t167 = t49 * t146;
t166 = t49 * t97;
t165 = t57 * t107;
t164 = t146 * t66;
t163 = t71 * t95;
t161 = t74 * t97;
t160 = t76 * t96;
t159 = t92 * t51;
t156 = t95 * t58;
t27 = qJD(5) * t49 + t156;
t21 = t95 * t27;
t157 = t95 * t57;
t99 = qJD(3) ^ 2;
t155 = t99 * t96;
t154 = t99 * t98;
t153 = -t27 * t161 - t47 * t162;
t141 = qJD(5) * t97;
t152 = -t47 * t141 - t21;
t151 = t107 * t26 + t49 * t68;
t150 = -t74 * t57 - t71 * t66;
t149 = t96 ^ 2 - t98 ^ 2;
t77 = qJD(1) * t85;
t140 = t71 * qJD(3);
t133 = t49 * t163;
t131 = pkin(3) * t144;
t130 = pkin(3) * t143;
t100 = qJD(1) ^ 2;
t129 = t96 * t100 * t98;
t122 = t63 * t97;
t121 = 0.2e1 * t146;
t120 = qJD(3) * t148;
t118 = t98 * t128;
t117 = t114 * t97 - t6 * t95;
t116 = -t114 * t95 - t6 * t97;
t113 = t107 * t27 - t47 * t68;
t33 = -pkin(4) * t107 - pkin(7) * t74 + t75;
t36 = t147 * t72 - t148 * t158;
t9 = t33 * t97 - t36 * t95;
t10 = t33 * t95 + t36 * t97;
t112 = t107 * t58 - t146 * t68;
t60 = t76 * t98 + t139;
t111 = 0.2e1 * qJD(3) * t77;
t110 = -t66 * t123 - t63 * t142 + t53;
t109 = t74 * t141 + t163;
t17 = t147 * t45 - t159;
t15 = -qJD(3) * pkin(4) - t17;
t82 = pkin(3) * t92 + pkin(7);
t106 = t63 * t15 - t57 * t82;
t105 = -t98 * t120 - t138;
t103 = -t109 * t63 - t74 * t157;
t102 = t117 * qJD(5) + t1 * t97 - t2 * t95;
t56 = t60 * qJD(3);
t59 = t89 - t160;
t101 = t55 * t98 + t56 * t96 + (-t59 * t98 - t60 * t96) * qJD(3);
t84 = -t147 * pkin(3) - pkin(4);
t64 = t66 ^ 2;
t62 = t68 * qJD(3);
t52 = -t96 * t120 + t136;
t31 = pkin(4) * t68 - pkin(7) * t71 + t130;
t30 = pkin(4) * t146 + pkin(7) * t66 + t131;
t23 = t92 * t105 + t147 * t52;
t22 = -t147 * t105 + t52 * t92;
t20 = t147 * t50 - t159;
t19 = t50 * t92 + t43;
t8 = t20 * t97 + t30 * t95;
t7 = -t20 * t95 + t30 * t97;
t4 = -t10 * qJD(5) - t95 * t23 + t97 * t31;
t3 = t9 * qJD(5) + t97 * t23 + t95 * t31;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t118, -0.2e1 * t149 * t134, t154, -0.2e1 * t118, -t155, 0, t111 * t96 - t83 * t154, t111 * t98 + t83 * t155, t101, t101 * t83, t146 * t71 + t58 * t74, t112 + t150, t140, t66 * t68 - t165, -t62, 0, t75 * t57 + t65 * t68 + (-t22 + (-qJD(1) * t107 + t66) * t178) * qJD(3), t58 * t75 + t65 * t71 + (t121 * t178 - t23) * qJD(3), t107 * t12 + t11 * t74 + t146 * t22 - t17 * t71 - t18 * t68 - t23 * t66 + t35 * t58 - t36 * t57, t175 + t12 * t36 - t17 * t22 + t18 * t23 + (t65 + t145) * t130, -t108 * t49 - t26 * t161, -t133 + (t172 + (-t166 + t169) * qJD(5)) * t74 + t153, t151 + t180, t109 * t47 + t74 * t21, t103 + t113, t63 * t68 - t165, -t107 * t2 + t109 * t15 - t114 * t68 + t74 * t173 + t22 * t47 + t27 * t35 + t4 * t63 + t57 * t9, t1 * t107 - t10 * t57 - t108 * t15 + t11 * t161 + t22 * t49 - t26 * t35 - t3 * t63 - t6 * t68, -t10 * t27 + t26 * t9 - t3 * t47 - t4 * t49 + t117 * t71 + (qJD(5) * t116 - t1 * t95 - t2 * t97) * t74, t1 * t10 - t114 * t4 + t15 * t22 + t2 * t9 + t3 * t6 + t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t154, 0, t55 * t96 - t56 * t98 + (-t59 * t96 + t60 * t98) * qJD(3), 0, 0, 0, 0, 0, 0, -t62, -t140, -t112 + t150, t12 * t74 - t17 * t68 + t18 * t71 - t174, 0, 0, 0, 0, 0, 0, t103 - t113, t151 - t180, t133 + (-t172 + (t166 + t169) * qJD(5)) * t74 + t153, t102 * t74 - t116 * t71 + t15 * t68 - t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t149 * t100, 0, t129, 0, 0, -t77 * t144, -t77 * t98 * qJD(1) - t87 + (t59 + t160) * qJD(3), 0, 0, t164, -t64 + t179, -t78 + (t79 + t66) * qJD(3), -t164, 0, 0, qJD(3) * t19 - t131 * t66 - t146 * t65 - t11, -t126 + t65 * t66 + (-pkin(3) * t146 + qJD(4) * t92) * t144 + (-t92 * (-t135 * t98 - t60) + t20) * qJD(3), (t18 - t19) * t146 + (-t17 + t20) * t66 + (-t147 * t58 - t57 * t92) * pkin(3), t17 * t19 - t18 * t20 + (-t11 * t147 + t12 * t92 - t144 * t65) * pkin(3), t122 * t49 - t172, (-t26 - t171) * t97 - t182 + t152, t122 * t63 + t157 - t167, t123 * t47 - t27 * t97, t110 + t170, -t63 * t146, -t11 * t97 - t19 * t47 + t27 * t84 + t114 * t146 + (-t141 * t82 - t7) * t63 + t106 * t95, t173 - t19 * t49 - t26 * t84 + t6 * t146 + (t142 * t82 + t8) * t63 + t106 * t97, t47 * t8 + t49 * t7 + (-t27 * t82 + t114 * t66 + t1 + (t49 * t82 + t114) * qJD(5)) * t97 + (-t26 * t82 - t6 * t66 - t2 + (t47 * t82 - t6) * qJD(5)) * t95, t102 * t82 + t11 * t84 + t114 * t7 - t15 * t19 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * qJD(3), -t78 + (t79 - t66) * qJD(3), -t64 - t179, t146 * t17 + t18 * t66 + t80, 0, 0, 0, 0, 0, 0, t110 - t170, -t63 ^ 2 * t97 - t157 - t167, (t26 - t171) * t97 + t182 + t152, -t146 * t15 + t183 * t97 + t184 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, -t47 ^ 2 + t49 ^ 2, t47 * t63 - t26, -t168, -t156 + (-qJD(5) + t63) * t49, t57, -t15 * t49 + t183, t15 * t47 - t184, 0, 0;];
tauc_reg = t5;
