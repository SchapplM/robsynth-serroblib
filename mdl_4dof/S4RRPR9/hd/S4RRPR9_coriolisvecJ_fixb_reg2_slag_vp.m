% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:04
% EndTime: 2019-12-31 17:10:09
% DurationCPUTime: 1.47s
% Computational Cost: add. (1782->231), mult. (4711->367), div. (0->0), fcn. (3148->6), ass. (0->138)
t108 = sin(qJ(4));
t110 = cos(qJ(4));
t106 = sin(pkin(7));
t109 = sin(qJ(2));
t152 = qJD(1) * t109;
t138 = t106 * t152;
t107 = cos(pkin(7));
t147 = t107 * qJD(2);
t71 = -t138 + t147;
t137 = t107 * t152;
t148 = t106 * qJD(2);
t72 = t137 + t148;
t28 = t108 * t72 - t110 * t71;
t177 = t28 ^ 2;
t31 = t108 * t71 + t110 * t72;
t176 = t31 ^ 2;
t111 = cos(qJ(2));
t146 = t111 * qJD(1);
t97 = -qJD(4) + t146;
t175 = t28 * t97;
t145 = qJD(1) * qJD(2);
t174 = -0.2e1 * t145;
t149 = qJD(4) * t110;
t158 = t108 * t106;
t173 = -qJD(4) * t158 + t107 * t149;
t104 = t109 ^ 2;
t105 = t111 ^ 2;
t172 = qJD(1) * (t104 - 0.2e1 * t105);
t136 = t111 * t145;
t129 = t106 * t136;
t130 = t110 * t111 * t147;
t12 = t108 * (qJD(4) * t72 + t129) - qJD(1) * t130 - t71 * t149;
t82 = -t111 * pkin(2) - t109 * qJ(3) - pkin(1);
t65 = t82 * qJD(1);
t101 = pkin(5) * t146;
t88 = qJD(2) * qJ(3) + t101;
t34 = -t106 * t88 + t107 * t65;
t15 = -pkin(3) * t146 - t72 * pkin(6) + t34;
t35 = t106 * t65 + t107 * t88;
t17 = t71 * pkin(6) + t35;
t124 = t108 * t17 - t110 * t15;
t159 = t107 * t111;
t121 = pkin(3) * t109 - pkin(6) * t159;
t117 = t121 * qJD(2);
t126 = pkin(2) * t109 - qJ(3) * t111;
t59 = qJD(2) * t126 - t109 * qJD(3);
t50 = t59 * qJD(1);
t100 = pkin(5) * t152;
t80 = (qJD(3) - t100) * qJD(2);
t23 = -t106 * t80 + t107 * t50;
t14 = qJD(1) * t117 + t23;
t24 = t106 * t50 + t107 * t80;
t16 = -pkin(6) * t129 + t24;
t1 = -qJD(4) * t124 + t108 * t14 + t110 * t16;
t77 = t126 * qJD(1);
t42 = pkin(5) * t138 + t107 * t77;
t25 = qJD(1) * t121 + t42;
t160 = t107 * t109;
t161 = t106 * t111;
t119 = -pkin(5) * t160 - pkin(6) * t161;
t63 = t106 * t77;
t32 = qJD(1) * t119 + t63;
t168 = pkin(6) + qJ(3);
t85 = t168 * t106;
t86 = t168 * t107;
t39 = -t108 * t86 - t110 * t85;
t74 = -t110 * t107 + t158;
t171 = -qJD(3) * t74 + qJD(4) * t39 - t108 * t25 - t110 * t32;
t40 = -t108 * t85 + t110 * t86;
t75 = t110 * t106 + t108 * t107;
t170 = -qJD(3) * t75 - qJD(4) * t40 + t108 * t32 - t110 * t25;
t169 = t31 * t28;
t118 = t75 * t111;
t62 = t75 * qJD(4);
t167 = qJD(1) * t118 - t62;
t166 = -t74 * t146 - t173;
t151 = qJD(2) * t109;
t143 = pkin(5) * t151;
t36 = t106 * t143 + t107 * t59;
t95 = pkin(5) * t159;
t47 = t106 * t82 + t95;
t96 = pkin(5) * t136;
t58 = pkin(3) * t129 + t96;
t164 = qJD(2) * pkin(2);
t163 = t107 * t71;
t113 = qJD(1) ^ 2;
t162 = t105 * t113;
t157 = t111 * t113;
t112 = qJD(2) ^ 2;
t156 = t112 * t109;
t155 = t112 * t111;
t153 = t104 - t105;
t150 = qJD(2) * t111;
t144 = pkin(5) * t161;
t142 = pkin(3) * t106 + pkin(5);
t141 = t106 * t146;
t135 = t109 * t145;
t134 = pkin(1) * t174;
t133 = qJD(3) - t164;
t132 = -t72 + t148;
t131 = -t71 + t147;
t128 = t111 * t135;
t81 = t100 + t133;
t127 = t133 - t81;
t6 = t108 * t15 + t110 * t17;
t70 = t107 * t82;
t33 = -pkin(6) * t160 + t70 + (-pkin(5) * t106 - pkin(3)) * t111;
t38 = -t106 * t109 * pkin(6) + t47;
t9 = -t108 * t38 + t110 * t33;
t10 = t108 * t33 + t110 * t38;
t123 = qJD(1) * t132;
t122 = qJD(1) * t131;
t120 = t111 * t123;
t115 = qJD(2) * t118;
t2 = -qJD(4) * t6 - t108 * t16 + t110 * t14;
t13 = qJD(1) * t115 + qJD(4) * t31;
t103 = t107 ^ 2;
t102 = t106 ^ 2;
t99 = -t107 * pkin(3) - pkin(2);
t94 = t109 * t157;
t89 = -0.2e1 * t128;
t78 = t142 * t109;
t67 = t142 * t150;
t66 = pkin(3) * t141 + t101;
t56 = t74 * t109;
t55 = t75 * t109;
t51 = t106 * t59;
t46 = t70 - t144;
t43 = -pkin(5) * t137 + t63;
t41 = -t71 * pkin(3) + t81;
t37 = -t107 * t143 + t51;
t26 = qJD(2) * t119 + t51;
t22 = t173 * t109 + t115;
t21 = t108 * t111 * t148 + t109 * t62 - t130;
t20 = t117 + t36;
t4 = -qJD(4) * t10 - t108 * t26 + t110 * t20;
t3 = qJD(4) * t9 + t108 * t20 + t110 * t26;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t128, t153 * t174, t155, t89, -t156, 0, -pkin(5) * t155 + t109 * t134, pkin(5) * t156 + t111 * t134, 0, 0, (t103 * t152 + t107 * t72) * t150, (t163 + (-t72 - 0.2e1 * t137) * t106) * t150, (t107 * t172 + t109 * t72) * qJD(2), (t102 * t152 - t106 * t71) * t150, (-t106 * t172 + t109 * t71) * qJD(2), t89, (-qJD(1) * t36 - t23) * t111 + ((-pkin(5) * t71 + t106 * t81) * t111 + (t34 + (t46 + 0.2e1 * t144) * qJD(1)) * t109) * qJD(2), (qJD(1) * t37 + t24) * t111 + ((pkin(5) * t72 + t107 * t81) * t111 + (-t35 + (-t47 + 0.2e1 * t95) * qJD(1)) * t109) * qJD(2), -t36 * t72 + t37 * t71 + (-t106 * t24 - t107 * t23) * t109 + (-t106 * t35 - t107 * t34 + (-t106 * t47 - t107 * t46) * qJD(1)) * t150, t23 * t46 + t24 * t47 + t34 * t36 + t35 * t37 + (t81 + t100) * pkin(5) * t150, t12 * t56 - t31 * t21, t12 * t55 + t56 * t13 + t21 * t28 - t31 * t22, t12 * t111 + t21 * t97 + (-qJD(1) * t56 + t31) * t151, t13 * t55 + t28 * t22, t13 * t111 + t22 * t97 + (-qJD(1) * t55 - t28) * t151, (-t97 - t146) * t151, -t2 * t111 + t78 * t13 + t41 * t22 + t67 * t28 - t4 * t97 + t58 * t55 + (qJD(1) * t9 - t124) * t151, t1 * t111 - t78 * t12 - t41 * t21 + t3 * t97 + t67 * t31 - t58 * t56 + (-qJD(1) * t10 - t6) * t151, -t1 * t55 - t10 * t13 + t9 * t12 - t124 * t21 + t2 * t56 - t6 * t22 - t3 * t28 - t4 * t31, t1 * t10 - t124 * t4 + t2 * t9 + t6 * t3 + t41 * t67 + t58 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, t153 * t113, 0, t94, 0, 0, t113 * pkin(1) * t109, pkin(1) * t157, 0, 0, t107 * t120, (t106 * t72 - t163 + (-t102 + t103) * qJD(2)) * t146, t107 * t162 + t109 * t123, -t131 * t141, -t106 * t162 + t109 * t122, t94, ((-qJ(3) * t148 - t34) * t109 + (-pkin(5) * t131 + t106 * t127 + t42) * t111) * qJD(1), ((-qJ(3) * t147 + t35) * t109 + (pkin(5) * t132 + t107 * t127 - t43) * t111) * qJD(1), t42 * t72 - t43 * t71 + (qJD(3) * t71 + t146 * t34 + t24) * t107 + (qJD(3) * t72 + t146 * t35 - t23) * t106, -t34 * t42 - t35 * t43 + (-t106 * t34 + t107 * t35) * qJD(3) + (-t23 * t106 + t24 * t107) * qJ(3) + (-t81 - t164) * t101, -t12 * t75 - t166 * t31, t12 * t74 - t75 * t13 + t166 * t28 + t167 * t31, t166 * t97 + (qJD(2) * t75 - t31) * t152, t13 * t74 - t167 * t28, -t167 * t97 + (-qJD(2) * t74 + t28) * t152, t97 * t152, t99 * t13 - t66 * t28 + t58 * t74 - t170 * t97 - t167 * t41 + (qJD(2) * t39 + t124) * t152, -t99 * t12 - t66 * t31 + t58 * t75 + t171 * t97 - t166 * t41 + (-qJD(2) * t40 + t6) * t152, -t1 * t74 + t39 * t12 - t124 * t166 - t40 * t13 + t167 * t6 - t170 * t31 - t171 * t28 - t2 * t75, t1 * t40 - t124 * t170 + t171 * t6 + t2 * t39 - t41 * t66 + t58 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t111 * t122, -t71 ^ 2 - t72 ^ 2, t34 * t72 - t35 * t71 + t96, 0, 0, 0, 0, 0, 0, -t31 * t97 + t13, -t12 + t175, -t176 - t177, -t124 * t31 + t28 * t6 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t176 - t177, -t12 - t175, -t169, -t136 * t75 + (-qJD(4) - t97) * t31, t135, -t41 * t31 - t6 * t97 + t2, t124 * t97 + t41 * t28 - t1, 0, 0;];
tauc_reg = t5;
