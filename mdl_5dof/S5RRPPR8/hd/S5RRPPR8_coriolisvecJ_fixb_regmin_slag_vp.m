% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:15
% EndTime: 2019-12-31 19:39:19
% DurationCPUTime: 1.18s
% Computational Cost: add. (1128->216), mult. (2847->311), div. (0->0), fcn. (1870->6), ass. (0->133)
t123 = sin(qJ(5));
t125 = cos(qJ(5));
t153 = qJD(5) * t125;
t154 = qJD(5) * t123;
t120 = sin(pkin(8));
t126 = cos(qJ(2));
t151 = qJD(1) * qJD(2);
t146 = t126 * t151;
t121 = cos(pkin(8));
t124 = sin(qJ(2));
t147 = t124 * t151;
t97 = t121 * t147;
t62 = t120 * t146 - t97;
t80 = t124 * t120 + t126 * t121;
t73 = t80 * qJD(2);
t63 = qJD(1) * t73;
t68 = t80 * qJD(1);
t157 = qJD(1) * t126;
t148 = t120 * t157;
t158 = qJD(1) * t124;
t70 = t121 * t158 - t148;
t1 = -t123 * t62 + t125 * t63 - t68 * t153 - t70 * t154;
t115 = qJD(2) - qJD(5);
t174 = t123 * t70 + t125 * t68;
t166 = t174 * t115;
t181 = t1 - t166;
t28 = -t123 * t68 + t125 * t70;
t165 = t28 * t115;
t2 = t123 * t63 + t125 * t62 + t70 * t153 - t68 * t154;
t180 = -t2 - t165;
t179 = t28 * t174;
t178 = t174 ^ 2 - t28 ^ 2;
t78 = -qJD(1) * pkin(1) - pkin(2) * t157 - qJ(3) * t158;
t45 = pkin(3) * t157 + qJD(4) - t78;
t23 = t68 * pkin(4) + t45;
t116 = qJD(2) * qJD(3);
t156 = qJD(2) * t124;
t170 = pkin(6) - qJ(4);
t65 = -t126 * qJD(4) - t170 * t156;
t41 = t65 * qJD(1) + t116;
t101 = pkin(6) * t146;
t152 = t124 * qJD(4);
t155 = qJD(2) * t126;
t50 = t101 + (-qJ(4) * t155 - t152) * qJD(1);
t145 = t120 * t41 - t121 * t50;
t5 = -t63 * pkin(7) - t145;
t12 = t120 * t50 + t121 * t41;
t6 = -t62 * pkin(7) + t12;
t177 = t123 * t6 - t125 * t5 + t23 * t28;
t176 = t123 * t5 + t125 * t6 - t23 * t174;
t175 = -0.2e1 * t151;
t173 = qJD(5) + t115;
t127 = -pkin(2) - pkin(3);
t172 = t68 * pkin(7);
t171 = t70 * pkin(7);
t96 = t170 * t126;
t66 = qJD(2) * t96 - t152;
t19 = t120 * t66 + t121 * t65;
t149 = t127 * qJD(2);
t107 = pkin(6) * t158;
t86 = qJ(4) * t158 - t107;
t58 = qJD(3) + t149 - t86;
t117 = qJD(2) * qJ(3);
t108 = pkin(6) * t157;
t88 = -qJ(4) * t157 + t108;
t77 = t117 + t88;
t17 = t120 * t58 + t121 * t77;
t33 = t120 * t88 + t121 * t86;
t95 = t170 * t124;
t38 = t120 * t95 + t121 * t96;
t169 = qJD(2) * pkin(2);
t129 = qJD(1) ^ 2;
t164 = t126 * t129;
t128 = qJD(2) ^ 2;
t163 = t128 * t124;
t162 = t128 * t126;
t111 = t124 * qJD(3);
t161 = qJ(3) * t146 + qJD(1) * t111;
t160 = qJ(3) * t155 + t111;
t118 = t124 ^ 2;
t159 = -t126 ^ 2 + t118;
t93 = -t126 * pkin(2) - t124 * qJ(3) - pkin(1);
t150 = t124 * t164;
t16 = -t120 * t77 + t121 * t58;
t18 = -t120 * t65 + t121 * t66;
t32 = -t120 * t86 + t121 * t88;
t37 = -t120 * t96 + t121 * t95;
t144 = qJD(1) * t93 + t78;
t143 = pkin(1) * t175;
t142 = qJD(3) - t169;
t89 = -t120 * qJ(3) + t121 * t127;
t141 = qJD(3) * t120 + t32;
t140 = qJD(3) * t121 - t33;
t79 = t126 * pkin(3) - t93;
t139 = t124 * t149;
t7 = -qJD(2) * pkin(4) + t16 - t171;
t8 = t17 - t172;
t136 = t123 * t8 - t125 * t7;
t135 = -t123 * t7 - t125 * t8;
t134 = t120 * t16 - t121 * t17;
t81 = -t126 * t120 + t124 * t121;
t30 = t123 * t81 + t125 * t80;
t31 = -t123 * t80 + t125 * t81;
t105 = qJ(3) * t157;
t64 = t127 * t158 + t105;
t51 = pkin(2) * t147 - t161;
t67 = pkin(2) * t156 - t160;
t133 = -pkin(6) * t128 - qJD(1) * t67 - t51;
t132 = (-t120 * t125 - t121 * t123) * t115;
t131 = (-t120 * t123 + t121 * t125) * t115;
t44 = t139 + t160;
t36 = qJD(1) * t139 + t161;
t91 = -pkin(6) * t147 + t116;
t92 = t107 + t142;
t94 = t108 + t117;
t130 = t91 * t126 + (t126 * t92 + (-t94 + t108) * t124) * qJD(2);
t90 = t121 * qJ(3) + t120 * t127;
t87 = pkin(2) * t158 - t105;
t85 = -pkin(4) + t89;
t72 = t120 * t155 - t121 * t156;
t34 = t80 * pkin(4) + t79;
t29 = -t70 * pkin(4) + t64;
t22 = t72 * pkin(4) + t44;
t21 = -t80 * pkin(7) + t38;
t20 = -t81 * pkin(7) + t37;
t15 = t62 * pkin(4) + t36;
t14 = t33 + t171;
t13 = t32 - t172;
t10 = -t72 * pkin(7) + t19;
t9 = -t73 * pkin(7) + t18;
t4 = t31 * qJD(5) + t123 * t73 + t125 * t72;
t3 = -t30 * qJD(5) - t123 * t72 + t125 * t73;
t11 = [0, 0, 0, 0.2e1 * t124 * t146, t159 * t175, t162, -t163, 0, -pkin(6) * t162 + t124 * t143, pkin(6) * t163 + t126 * t143, t133 * t126 + t144 * t156, t130, t133 * t124 - t144 * t155, t130 * pkin(6) + t51 * t93 + t78 * t67, -t18 * qJD(2) + t36 * t80 + t44 * t68 + t45 * t72 + t79 * t62, t19 * qJD(2) + t36 * t81 + t44 * t70 + t45 * t73 + t79 * t63, -t12 * t80 + t145 * t81 - t16 * t73 - t17 * t72 - t18 * t70 - t19 * t68 - t37 * t63 - t38 * t62, t12 * t38 - t145 * t37 + t16 * t18 + t17 * t19 + t36 * t79 + t45 * t44, t1 * t31 + t28 * t3, -t1 * t30 - t174 * t3 - t31 * t2 - t28 * t4, -t3 * t115, t4 * t115, 0, t22 * t174 + t34 * t2 + t15 * t30 + t23 * t4 - (-t123 * t10 + t125 * t9 + (-t123 * t20 - t125 * t21) * qJD(5)) * t115, t22 * t28 + t34 * t1 + t15 * t31 + t23 * t3 + (t125 * t10 + t123 * t9 + (-t123 * t21 + t125 * t20) * qJD(5)) * t115; 0, 0, 0, -t150, t159 * t129, 0, 0, 0, t129 * pkin(1) * t124, pkin(1) * t164, (-t124 * t78 + t126 * t87) * qJD(1), ((t94 - t117) * t124 + (t142 - t92) * t126) * qJD(1), 0.2e1 * t116 + (t124 * t87 + t126 * t78) * qJD(1), t91 * qJ(3) + t94 * qJD(3) - t78 * t87 + (t124 * t94 + (-t92 - t169) * t126) * qJD(1) * pkin(6), t141 * qJD(2) + t45 * t70 - t64 * t68 + t145, t140 * qJD(2) - t45 * t68 - t64 * t70 + t12, -t90 * t62 - t89 * t63 + (t141 - t17) * t70 + (-t140 + t16) * t68, -qJD(3) * t134 + t12 * t90 - t145 * t89 - t16 * t32 - t17 * t33 - t45 * t64, -t179, t178, -t181, -t180, 0, -t29 * t174 + (-t123 * t14 + t125 * t13) * t115 - qJD(3) * t132 + (-(-t123 * t85 - t125 * t90) * t115 - t135) * qJD(5) + t177, -t29 * t28 - (t123 * t13 + t125 * t14) * t115 + qJD(3) * t131 + ((-t123 * t90 + t125 * t85) * t115 - t136) * qJD(5) + t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, 0, -t118 * t129 - t128, -t94 * qJD(2) + t78 * t158 + t101, -t120 * t128 - t68 * t158, -t121 * t128 - t70 * t158, -t120 * t62 - t121 * t63 + (-t120 * t70 + t121 * t68) * qJD(2), qJD(2) * t134 + t12 * t120 - t121 * t145 - t158 * t45, 0, 0, 0, 0, 0, t115 * t132 - t158 * t174, -t115 * t131 - t158 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97 + (-t70 + t148) * qJD(2), 0.2e1 * t68 * qJD(2), -t68 ^ 2 - t70 ^ 2, t16 * t70 + t17 * t68 + t36, 0, 0, 0, 0, 0, t2 - t165, t1 + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t178, t181, t180, 0, t173 * t135 - t177, t173 * t136 - t176;];
tauc_reg = t11;
