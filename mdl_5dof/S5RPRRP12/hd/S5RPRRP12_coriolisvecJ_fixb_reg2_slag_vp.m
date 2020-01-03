% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP12_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:20
% EndTime: 2019-12-31 18:57:26
% DurationCPUTime: 1.84s
% Computational Cost: add. (2015->269), mult. (4348->382), div. (0->0), fcn. (2456->4), ass. (0->163)
t92 = cos(qJ(4));
t140 = t92 * qJD(3);
t93 = cos(qJ(3));
t149 = qJD(1) * t93;
t90 = sin(qJ(4));
t68 = t90 * t149 - t140;
t91 = sin(qJ(3));
t150 = qJD(1) * t91;
t82 = qJD(4) + t150;
t177 = t68 * t82;
t141 = qJD(4) * t93;
t126 = t90 * t141;
t128 = t91 * t140;
t103 = t126 + t128;
t38 = t103 * qJD(1) - qJD(4) * t140;
t192 = -t38 - t177;
t146 = qJD(3) * t93;
t94 = -pkin(1) - pkin(6);
t81 = t94 * qJD(1) + qJD(2);
t130 = t81 * t146;
t142 = qJD(4) * t92;
t143 = qJD(4) * t90;
t112 = pkin(3) * t93 + pkin(7) * t91;
t66 = t112 * qJD(3) + qJD(2);
t50 = t66 * qJD(1);
t75 = pkin(3) * t91 - pkin(7) * t93 + qJ(2);
t56 = t75 * qJD(1);
t74 = t91 * t81;
t59 = qJD(3) * pkin(7) + t74;
t118 = -t92 * t130 - t56 * t142 + t59 * t143 - t90 * t50;
t31 = t92 * t56 - t59 * t90;
t191 = -t31 * t82 - t118;
t133 = 0.2e1 * qJD(1);
t166 = t91 * t94;
t42 = t92 * t166 + t90 * t75;
t32 = t56 * t90 + t59 * t92;
t101 = -qJD(4) * t32 + t92 * t50;
t12 = -t90 * t130 + t101;
t190 = -t32 * t82 - t12;
t131 = t92 * t149;
t148 = qJD(3) * t90;
t70 = t131 + t148;
t175 = t70 * t82;
t147 = qJD(3) * t91;
t129 = t90 * t147;
t144 = qJD(4) * t70;
t39 = -qJD(1) * t129 + t144;
t189 = t39 + t175;
t188 = t70 ^ 2;
t187 = pkin(4) * t68;
t22 = -qJ(5) * t68 + t32;
t186 = t22 * t82;
t29 = pkin(4) * t39 + t81 * t147;
t185 = t29 * t90;
t184 = t29 * t92;
t181 = t38 * t90;
t180 = t38 * t91;
t179 = t39 * t91;
t178 = t39 * t92;
t176 = t70 * t68;
t174 = t70 * t90;
t173 = t81 * t90;
t172 = t81 * t93;
t171 = t82 * t90;
t170 = t82 * t92;
t169 = t90 * t91;
t168 = t90 * t93;
t167 = t91 * t92;
t165 = t92 * t93;
t164 = t93 * t38;
t163 = t93 * t39;
t95 = qJD(3) ^ 2;
t162 = t95 * t91;
t161 = t95 * t93;
t160 = -qJ(5) - pkin(7);
t21 = -qJ(5) * t70 + t31;
t18 = pkin(4) * t82 + t21;
t159 = t18 - t21;
t119 = qJD(4) * t160;
t73 = t112 * qJD(1);
t36 = -t81 * t168 + t92 * t73;
t158 = (pkin(4) * t93 + qJ(5) * t167) * qJD(1) + t36 + t90 * qJD(5) - t92 * t119;
t132 = t90 * t150;
t139 = t92 * qJD(5);
t37 = t81 * t165 + t90 * t73;
t157 = qJ(5) * t132 - t90 * t119 - t139 + t37;
t89 = t93 ^ 2;
t156 = t91 ^ 2 - t89;
t96 = qJD(1) ^ 2;
t155 = -t95 - t96;
t154 = pkin(4) * qJD(1);
t153 = qJ(5) * t93;
t152 = qJD(3) * pkin(3);
t151 = t96 * qJ(2);
t145 = qJD(3) * t94;
t138 = qJ(2) * qJD(3);
t137 = qJD(1) * qJD(3);
t136 = t90 * t166;
t135 = t93 * t96 * t91;
t127 = t93 * t145;
t134 = t92 * t127 + t75 * t142 + t90 * t66;
t125 = t92 * t141;
t124 = t82 * t149;
t123 = qJD(2) * t133;
t122 = -t90 * t94 + pkin(4);
t83 = t93 * t137;
t60 = -t152 - t172;
t121 = -t60 + t172;
t120 = -qJD(5) - t187;
t117 = t82 + t150;
t116 = t68 + t140;
t115 = -t70 + t148;
t114 = qJD(4) * t91 + qJD(1);
t113 = t91 * t83;
t111 = t18 * t92 + t22 * t90;
t110 = t18 * t90 - t22 * t92;
t109 = t31 * t92 + t32 * t90;
t108 = t31 * t90 - t32 * t92;
t107 = t121 * qJD(3);
t106 = qJD(1) * t89 - t82 * t91;
t105 = -pkin(7) * t146 + t60 * t91;
t104 = t39 * qJ(5) + t118;
t102 = t125 - t129;
t100 = -t42 * qJD(4) + t92 * t66;
t98 = t38 * qJ(5) + t101;
t97 = -t109 * qJD(4) - t118 * t92 - t12 * t90;
t86 = qJ(2) * t123;
t85 = -pkin(4) * t92 - pkin(3);
t78 = t160 * t92;
t77 = t160 * t90;
t67 = (pkin(4) * t90 - t94) * t93;
t65 = t68 ^ 2;
t64 = t92 * t75;
t46 = -pkin(4) * t132 + t74;
t45 = t117 * t146;
t41 = t64 - t136;
t40 = t102 * pkin(4) + t91 * t145;
t35 = -t90 * t153 + t42;
t34 = -t120 + t60;
t33 = t122 * t91 - t92 * t153 + t64;
t28 = -t65 + t188;
t26 = t175 - t39;
t25 = -t38 + t177;
t24 = -t90 * t127 + t100;
t23 = -qJD(4) * t136 + t134;
t20 = t82 * t142 + (t115 * t93 + t82 * t167) * qJD(1);
t19 = -t82 * t143 + (t116 * t93 - t82 * t169) * qJD(1);
t17 = t68 * t171 - t178;
t16 = t70 * t170 - t181;
t15 = t102 * t68 + t90 * t163;
t14 = -t103 * t70 - t92 * t164;
t13 = -qJ(5) * t125 + (-qJD(5) * t93 + (qJ(5) * qJD(3) - qJD(4) * t94) * t91) * t90 + t134;
t10 = qJ(5) * t128 + (qJ(5) * t143 + t122 * qJD(3) - t139) * t93 + t100;
t9 = -t82 * t125 - t179 + (-t106 * t90 - t68 * t93) * qJD(3);
t8 = -t82 * t126 - t180 + (t106 * t92 + t70 * t93) * qJD(3);
t7 = -t163 - t114 * t170 + (-t117 * t168 + t68 * t91) * qJD(3);
t6 = t164 + t114 * t171 + (-t82 * t165 + (t70 - t131) * t91) * qJD(3);
t5 = -qJD(5) * t68 - t104;
t4 = -t70 * qJD(5) + (t154 - t173) * t146 + t98;
t3 = -t189 * t90 + t192 * t92;
t2 = (t68 * t92 + t174) * t147 + (t181 - t178 + (t68 * t90 - t70 * t92) * qJD(4)) * t93;
t1 = (t114 * t70 - t68 * t146 - t179) * t92 + (t114 * t68 + t70 * t146 - t180) * t90;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t86, -0.2e1 * t113, 0.2e1 * t156 * t137, -t162, 0.2e1 * t113, -t161, 0, -t94 * t162 + (qJD(2) * t91 + t93 * t138) * t133, -t94 * t161 + (qJD(2) * t93 - t91 * t138) * t133, 0, t86, t14, t2, t8, t15, t9, t45, t12 * t91 + t24 * t82 + (t60 * t142 - t39 * t94) * t93 + ((qJD(1) * t41 + t31) * t93 + (t121 * t90 + t68 * t94) * t91) * qJD(3), t118 * t91 - t23 * t82 + (-t60 * t143 + t38 * t94) * t93 + ((-qJD(1) * t42 - t32) * t93 + (t121 * t92 + t70 * t94) * t91) * qJD(3), -t23 * t68 - t24 * t70 + t38 * t41 - t39 * t42 + t109 * t147 + (qJD(4) * t108 + t118 * t90 - t12 * t92) * t93, -t107 * t166 - t118 * t42 + t12 * t41 + t32 * t23 + t31 * t24, t14, t2, t8, t15, t9, t45, t10 * t82 + t67 * t39 + t40 * t68 + (-t34 * t148 + t4) * t91 + (t34 * t142 + t185 + (qJD(1) * t33 + t18) * qJD(3)) * t93, -t13 * t82 - t67 * t38 + t40 * t70 + (-t140 * t34 - t5) * t91 + (-t34 * t143 + t184 + (-qJD(1) * t35 - t22) * qJD(3)) * t93, -t10 * t70 - t13 * t68 + t33 * t38 - t35 * t39 + t111 * t147 + (qJD(4) * t110 - t4 * t92 - t5 * t90) * t93, t10 * t18 + t13 * t22 + t29 * t67 + t33 * t4 + t34 * t40 + t35 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t151, 0, 0, 0, 0, 0, 0, t155 * t91, t155 * t93, 0, -t151, 0, 0, 0, 0, 0, 0, t7, t6, t1, -t108 * t146 - t109 * qJD(1) + (-t107 + t97) * t91, 0, 0, 0, 0, 0, 0, t7, t6, t1, -t111 * qJD(1) + (-qJD(3) * t110 - t29) * t93 + (qJD(3) * t34 - qJD(4) * t111 - t4 * t90 + t5 * t92) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, -t156 * t96, 0, -t135, 0, 0, -t93 * t151, t91 * t151, 0, 0, t16, t3, t20, t17, t19, -t124, -pkin(3) * t39 - t36 * t82 - t116 * t74 + (-pkin(7) * t170 + t60 * t90) * qJD(4) + (t105 * t90 - t31 * t93) * qJD(1), pkin(3) * t38 + t37 * t82 + t115 * t74 + (pkin(7) * t171 + t60 * t92) * qJD(4) + (t105 * t92 + t32 * t93) * qJD(1), t36 * t70 + t37 * t68 + ((-t39 + t144) * pkin(7) + t191) * t92 + ((qJD(4) * t68 - t38) * pkin(7) + t190) * t90, -t31 * t36 - t32 * t37 + (-t60 - t152) * t74 + t97 * pkin(7), t16, t3, t20, t17, t19, -t124, -t184 + t39 * t85 - t46 * t68 - t158 * t82 + (t34 + t187) * t143 + (t34 * t169 + (qJD(3) * t77 - t18) * t93) * qJD(1), t185 - t38 * t85 - t46 * t70 + t157 * t82 + (pkin(4) * t174 + t34 * t92) * qJD(4) + (t34 * t167 + (qJD(3) * t78 + t22) * t93) * qJD(1), t38 * t77 + t39 * t78 + t158 * t70 + t157 * t68 + (-t18 * t82 + t5) * t92 + (-t4 - t186) * t90, t29 * t85 + t4 * t77 - t5 * t78 + (pkin(4) * t143 - t46) * t34 - t157 * t22 - t158 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t28, t25, -t176, t26, t83, -t60 * t70 - t190, t60 * t68 - t191, 0, 0, t176, t28, t25, -t176, t26, t83, t186 + (0.2e1 * t154 - t173) * t146 + (t120 - t34) * t70 + t98, -t188 * pkin(4) + t21 * t82 + (qJD(5) + t34) * t68 + t104, t38 * pkin(4) - t159 * t68, t159 * t22 + (-t34 * t70 + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t192, -t65 - t188, t18 * t70 + t22 * t68 + t29;];
tauc_reg = t11;
