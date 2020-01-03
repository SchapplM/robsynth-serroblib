% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR14_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:43
% EndTime: 2019-12-31 19:19:51
% DurationCPUTime: 2.21s
% Computational Cost: add. (3153->223), mult. (9650->450), div. (0->0), fcn. (10248->12), ass. (0->132)
t84 = cos(pkin(6));
t91 = cos(qJ(3));
t157 = t84 * t91;
t81 = sin(pkin(6));
t160 = t81 * t91;
t80 = sin(pkin(11));
t82 = sin(pkin(5));
t83 = cos(pkin(11));
t85 = cos(pkin(5));
t88 = sin(qJ(3));
t171 = t85 * t160 + (t83 * t157 - t80 * t88) * t82;
t87 = sin(qJ(4));
t170 = -0.4e1 * t87;
t166 = pkin(1) * t85;
t73 = t83 * t166;
t48 = t85 * pkin(2) + t73 + (-pkin(8) * t84 - qJ(2)) * t82 * t80;
t55 = (-pkin(8) * t80 * t81 - pkin(2) * t83 - pkin(1)) * t82;
t33 = -t81 * t48 + t84 * t55;
t158 = t83 * t88;
t161 = t81 * t88;
t47 = t85 * t161 + (t84 * t158 + t80 * t91) * t82;
t21 = -pkin(3) * t171 - t47 * pkin(9) + t33;
t159 = t82 * t83;
t59 = t81 * t159 - t85 * t84;
t163 = t48 * t84;
t108 = t55 * t81 + t163;
t150 = qJ(2) * t82;
t153 = t83 * t150 + t80 * t166;
t44 = (t84 * t159 + t81 * t85) * pkin(8) + t153;
t93 = t108 * t88 + t91 * t44;
t24 = -t59 * pkin(9) + t93;
t90 = cos(qJ(4));
t169 = t87 * t21 + t90 * t24;
t89 = cos(qJ(5));
t78 = t89 ^ 2;
t86 = sin(qJ(5));
t152 = t86 ^ 2 - t78;
t122 = t152 * qJD(5);
t149 = qJD(2) * t82;
t131 = t80 * t149;
t118 = t84 * t131;
t147 = qJD(3) * t91;
t126 = t81 * t147;
t130 = t83 * t149;
t17 = -t147 * t163 - t55 * t126 - t91 * t130 + (qJD(3) * t44 + t118) * t88;
t117 = t81 * t131;
t41 = t171 * qJD(3);
t42 = t47 * qJD(3);
t30 = t42 * pkin(3) - t41 * pkin(9) + t117;
t8 = -qJD(4) * t169 + t87 * t17 + t90 * t30;
t167 = 0.2e1 * t42;
t165 = t86 * pkin(9);
t34 = t47 * t87 + t59 * t90;
t26 = -qJD(4) * t34 + t41 * t90;
t35 = t47 * t90 - t59 * t87;
t27 = t171 * t89 + t35 * t86;
t15 = -t27 * qJD(5) + t26 * t89 + t42 * t86;
t164 = t15 * t86;
t156 = t87 * t89;
t155 = t89 * t90;
t77 = t87 ^ 2;
t151 = -t90 ^ 2 + t77;
t148 = qJD(3) * t88;
t146 = qJD(4) * t86;
t145 = qJD(4) * t87;
t144 = qJD(4) * t89;
t143 = qJD(4) * t90;
t142 = qJD(4) * t91;
t141 = qJD(5) * t86;
t140 = qJD(5) * t89;
t139 = qJD(5) * t90;
t138 = t90 * t165;
t137 = pkin(9) * t155;
t136 = -0.2e1 * pkin(3) * qJD(4);
t135 = -0.2e1 * pkin(4) * qJD(5);
t132 = t87 * t161;
t129 = t86 * t139;
t128 = t89 * t139;
t127 = t81 * t148;
t125 = t86 * t140;
t124 = t87 * t143;
t123 = t89 * t143;
t121 = t151 * qJD(4);
t120 = -0.2e1 * t85 * t149;
t119 = 0.2e1 * t124;
t116 = t86 * t123;
t25 = t35 * qJD(4) + t41 * t87;
t115 = t25 * pkin(4) - t26 * pkin(10);
t114 = -t90 * pkin(4) - t87 * pkin(10);
t113 = pkin(4) * t87 - pkin(10) * t90;
t11 = -pkin(10) * t171 + t169;
t23 = t59 * pkin(3) - t108 * t91 + t88 * t44;
t13 = t34 * pkin(4) - t35 * pkin(10) + t23;
t4 = t89 * t11 + t86 * t13;
t111 = t90 * t21 - t87 * t24;
t28 = -t171 * t86 + t35 * t89;
t109 = -t27 * t89 - t28 * t86;
t68 = -pkin(3) + t114;
t57 = t86 * t68 + t137;
t61 = t90 * t161 + t87 * t84;
t107 = t89 * t160 + t86 * t61;
t106 = t86 * t160 - t89 * t61;
t10 = pkin(4) * t171 - t111;
t6 = -t42 * pkin(4) - t8;
t105 = t10 * t140 + t6 * t86;
t104 = t10 * t141 - t6 * t89;
t103 = t113 * t86;
t102 = -t143 * t171 + t87 * t42;
t101 = -t145 * t171 - t90 * t42;
t100 = t34 * t140 + t86 * t25;
t99 = t34 * t141 - t89 * t25;
t51 = t61 * qJD(4) + t87 * t126;
t60 = -t90 * t84 + t132;
t98 = t60 * t140 + t51 * t86;
t97 = t60 * t141 - t51 * t89;
t7 = -t21 * t143 + t24 * t145 + t90 * t17 - t87 * t30;
t96 = -t87 * t141 + t123;
t95 = t87 * t144 + t129;
t94 = -t42 * pkin(10) + t7;
t92 = t93 * qJD(3);
t18 = (t80 * t157 + t158) * t149 + t92;
t56 = t89 * t68 - t138;
t50 = qJD(4) * t132 - t90 * t126 - t84 * t143;
t37 = -t57 * qJD(5) + (t89 * t113 + t87 * t165) * qJD(4);
t36 = t95 * pkin(9) - qJD(4) * t103 - t68 * t140;
t32 = t106 * qJD(5) + t89 * t127 + t86 * t50;
t31 = t107 * qJD(5) - t86 * t127 + t89 * t50;
t14 = t28 * qJD(5) + t26 * t86 - t42 * t89;
t3 = -t86 * t11 + t89 * t13;
t2 = -t4 * qJD(5) + t86 * t94 + (t91 * t118 + t88 * t130 + t115 + t92) * t89;
t1 = t11 * t141 - t13 * t140 + t89 * t94 - t86 * (t115 + t18);
t5 = [0, 0, 0, t80 * t120, t83 * t120, 0.2e1 * (t80 ^ 2 + t83 ^ 2) * t82 ^ 2 * qJD(2), 0.2e1 * (t153 * t83 + (t80 * t150 - t73) * t80) * t149, 0.2e1 * t47 * t41, 0.2e1 * t171 * t41 - 0.2e1 * t47 * t42, -0.2e1 * t41 * t59, t59 * t167, 0, -0.2e1 * t117 * t171 + 0.2e1 * t18 * t59 + 0.2e1 * t33 * t42, 0.2e1 * t47 * t117 - 0.2e1 * t17 * t59 + 0.2e1 * t33 * t41, 0.2e1 * t35 * t26, -0.2e1 * t35 * t25 - 0.2e1 * t26 * t34, -0.2e1 * t171 * t26 + 0.2e1 * t35 * t42, 0.2e1 * t171 * t25 - 0.2e1 * t34 * t42, -t171 * t167, 0.2e1 * t111 * t42 - 0.2e1 * t171 * t8 + 0.2e1 * t18 * t34 + 0.2e1 * t23 * t25, -0.2e1 * t169 * t42 - 0.2e1 * t171 * t7 + 0.2e1 * t18 * t35 + 0.2e1 * t23 * t26, 0.2e1 * t28 * t15, -0.2e1 * t28 * t14 - 0.2e1 * t15 * t27, 0.2e1 * t15 * t34 + 0.2e1 * t28 * t25, -0.2e1 * t14 * t34 - 0.2e1 * t27 * t25, 0.2e1 * t34 * t25, 0.2e1 * t10 * t14 + 0.2e1 * t2 * t34 + 0.2e1 * t3 * t25 + 0.2e1 * t6 * t27, 0.2e1 * t1 * t34 + 0.2e1 * t10 * t15 - 0.2e1 * t4 * t25 + 0.2e1 * t6 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t127 + t84 * t42, t59 * t126 + t84 * t41, 0, 0, 0, 0, 0, -t60 * t42 + t51 * t171 + (t34 * t148 - t25 * t91) * t81, -t61 * t42 - t50 * t171 + (t35 * t148 - t26 * t91) * t81, 0, 0, 0, 0, 0, -t107 * t25 + t60 * t14 + t51 * t27 + t32 * t34, t106 * t25 + t60 * t15 + t51 * t28 + t31 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t42, 0, -t18, t17, t35 * t143 + t26 * t87, -t87 * t25 + t26 * t90 + (-t34 * t90 - t35 * t87) * qJD(4), t102, -t101, 0, -pkin(3) * t25 - t102 * pkin(9) + t23 * t145 - t18 * t90, -pkin(3) * t26 + t101 * pkin(9) + t23 * t143 + t18 * t87, t15 * t156 + t96 * t28, t109 * t143 + (-t14 * t89 - t164 + (t27 * t86 - t28 * t89) * qJD(5)) * t87, (t144 * t34 - t15) * t90 + (qJD(4) * t28 - t99) * t87, (-t146 * t34 + t14) * t90 + (-qJD(4) * t27 - t100) * t87, t145 * t34 - t25 * t90, t56 * t25 + t37 * t34 + (-t2 + (pkin(9) * t27 + t10 * t86) * qJD(4)) * t90 + (pkin(9) * t14 + qJD(4) * t3 + t105) * t87, -t57 * t25 + t36 * t34 + (-t1 + (pkin(9) * t28 + t10 * t89) * qJD(4)) * t90 + (pkin(9) * t15 - qJD(4) * t4 - t104) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t126, 0, 0, 0, 0, 0, (-t87 * t142 - t90 * t148) * t81, (-t90 * t142 + t87 * t148) * t81, 0, 0, 0, 0, 0, (t146 * t60 - t32) * t90 + (-qJD(4) * t107 + t98) * t87, (t144 * t60 - t31) * t90 + (qJD(4) * t106 - t97) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -0.2e1 * t121, 0, 0, 0, t87 * t136, t90 * t136, 0.2e1 * t78 * t124 - 0.2e1 * t77 * t125, t116 * t170 + 0.2e1 * t77 * t122, 0.2e1 * t129 * t87 + 0.2e1 * t144 * t151, -0.2e1 * t121 * t86 + 0.2e1 * t128 * t87, -0.2e1 * t124, 0.2e1 * t56 * t145 - 0.2e1 * t37 * t90 + 0.2e1 * (t119 * t86 + t140 * t77) * pkin(9), -0.2e1 * t57 * t145 - 0.2e1 * t36 * t90 + 0.2e1 * (t119 * t89 - t141 * t77) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t42, t8, t7, t28 * t140 + t164, t109 * qJD(5) - t86 * t14 + t15 * t89, t100, -t99, 0, -pkin(4) * t14 - pkin(10) * t100 + t104, -pkin(4) * t15 + pkin(10) * t99 + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50, 0, 0, 0, 0, 0, t97, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, -t145, 0, -pkin(9) * t143, pkin(9) * t145, -t87 * t122 + t116, t125 * t170 - t143 * t152, t145 * t86 - t128, t95, 0, (pkin(10) * t155 + (-t89 * pkin(4) + t165) * t87) * qJD(5) + (t114 * t86 - t137) * qJD(4), (pkin(9) * t156 + t103) * qJD(5) + (t114 * t89 + t138) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t125, -0.2e1 * t122, 0, 0, 0, t86 * t135, t89 * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, t25, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t140 * t87 - t143 * t86, t145, t37, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t141, 0, -pkin(10) * t140, pkin(10) * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
