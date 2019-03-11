% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:43
% EndTime: 2019-03-08 20:16:49
% DurationCPUTime: 2.00s
% Computational Cost: add. (1272->223), mult. (3282->377), div. (0->0), fcn. (2822->8), ass. (0->145)
t82 = sin(qJ(5));
t76 = t82 ^ 2;
t85 = cos(qJ(5));
t78 = t85 ^ 2;
t144 = t76 - t78;
t163 = 2 * qJD(3);
t162 = 2 * qJD(5);
t161 = pkin(5) * t85;
t86 = cos(qJ(4));
t160 = pkin(9) * t86;
t83 = sin(qJ(4));
t159 = t83 * pkin(4);
t80 = sin(pkin(6));
t87 = cos(qJ(2));
t150 = t80 * t87;
t128 = t83 * t150;
t84 = sin(qJ(2));
t151 = t80 * t84;
t67 = qJD(2) * t151;
t73 = t86 * qJD(4);
t81 = cos(pkin(6));
t21 = -qJD(4) * t128 - t86 * t67 + t81 * t73;
t38 = t86 * t150 + t81 * t83;
t158 = t21 * t38;
t157 = t21 * t86;
t131 = t83 * qJD(4);
t88 = -pkin(2) - pkin(8);
t114 = t88 * t131;
t116 = t82 * t131;
t134 = qJD(5) * t86;
t119 = t85 * t134;
t46 = -t116 + t119;
t28 = t46 * pkin(5) + t114;
t156 = t28 * t82;
t155 = t28 * t85;
t146 = -qJ(6) - pkin(9);
t56 = t146 * t82;
t154 = t56 * t86;
t57 = t146 * t85;
t153 = t57 * t86;
t71 = -pkin(4) - t161;
t152 = t71 * t85;
t149 = t82 * t88;
t148 = t83 * t88;
t147 = t86 * t88;
t139 = qJD(2) * t87;
t117 = t80 * t139;
t145 = qJ(3) * t117 + qJD(3) * t151;
t62 = t85 * t148;
t103 = t159 - t160;
t95 = qJ(3) + t103;
t32 = t82 * t95 + t62;
t143 = t76 + t78;
t77 = t83 ^ 2;
t79 = t86 ^ 2;
t142 = t77 - t79;
t141 = t77 + t79;
t140 = qJ(6) * t86;
t138 = qJD(4) * t38;
t137 = qJD(4) * t82;
t136 = qJD(4) * t85;
t135 = qJD(5) * t82;
t74 = qJD(5) * t85;
t133 = qJD(5) * t88;
t132 = t82 * qJD(6);
t130 = t85 * qJD(6);
t129 = qJ(3) * qJD(4);
t127 = t82 * t148;
t126 = -0.2e1 * t135;
t112 = t88 * t73;
t49 = t85 * t95;
t104 = pkin(4) * t86 + pkin(9) * t83;
t93 = t104 * qJD(4) + qJD(3);
t125 = -qJD(5) * t49 - t85 * t112 - t82 * t93;
t124 = pkin(5) * t135;
t123 = t85 * t140;
t122 = t38 * t135;
t121 = t82 * t134;
t120 = t82 * t133;
t50 = (pkin(5) * t82 - t88) * t86;
t118 = t50 * t135;
t115 = t82 * t74;
t113 = t83 * t73;
t111 = qJ(6) * t131;
t110 = -t71 + t161;
t109 = pkin(5) - t149;
t108 = t143 * t86;
t107 = t142 * qJD(4);
t65 = 0.2e1 * t113;
t106 = t85 * t116;
t105 = t79 * t115;
t102 = pkin(5) * t76 + t152;
t16 = t109 * t83 - t123 + t49;
t19 = -t82 * t140 + t32;
t101 = t16 * t85 + t19 * t82;
t100 = t16 * t82 - t19 * t85;
t39 = t81 * t86 - t128;
t22 = t85 * t151 - t39 * t82;
t23 = t82 * t151 + t39 * t85;
t99 = t22 * t85 + t23 * t82;
t98 = t22 * t82 - t23 * t85;
t31 = t49 - t127;
t97 = t31 * t85 + t32 * t82;
t96 = t31 * t82 - t32 * t85;
t12 = t21 * t82 + t38 * t74;
t13 = -t21 * t85 + t122;
t94 = t85 * t131 + t121;
t45 = t82 * t73 + t83 * t74;
t20 = -t83 * t67 + t138;
t10 = -t23 * qJD(5) + t85 * t117 + t20 * t82;
t11 = -t20 * t85 - t39 * t135 + (t82 * t139 + t84 * t74) * t80;
t4 = -t99 * qJD(5) - t10 * t82 + t11 * t85;
t14 = t83 * t120 + t125;
t90 = -t32 * qJD(5) + t85 * t93;
t15 = -t82 * t112 + t90;
t92 = -t97 * qJD(5) - t14 * t85 - t15 * t82;
t7 = -t20 * t83 - t157 + (t38 * t83 + t39 * t86) * qJD(4);
t34 = -t146 * t135 - t130;
t35 = t146 * t74 - t132;
t91 = -t34 * t85 - t35 * t82 + (-t56 * t85 + t57 * t82) * qJD(5);
t89 = qJ(6) * t121 + t85 * t111 - t86 * t130 + t90;
t75 = qJ(3) * t163;
t64 = -0.2e1 * t115;
t63 = 0.2e1 * t115;
t51 = -0.2e1 * t144 * qJD(5);
t44 = t141 * t74;
t42 = t83 * t135 - t85 * t73;
t41 = t141 * t135;
t40 = qJD(4) * t108;
t30 = -0.2e1 * t78 * t113 - 0.2e1 * t105;
t29 = -0.2e1 * t76 * t113 + 0.2e1 * t105;
t27 = t134 * t144 + t106;
t26 = -0.4e1 * t86 * t115 + t144 * t131;
t25 = 0.2e1 * t82 * t107 - 0.2e1 * t83 * t119;
t24 = -0.2e1 * t83 * t121 - 0.2e1 * t142 * t136;
t18 = t144 * t79 * t162 + 0.4e1 * t86 * t106;
t17 = (-0.1e1 + t143) * t65;
t9 = t86 * t132 - t82 * t111 + (t123 + t127) * qJD(5) + t125;
t8 = t109 * t73 + t89;
t6 = (-t38 * t136 - t11) * t83 + (-qJD(4) * t23 - t13) * t86;
t5 = (-t38 * t137 + t10) * t83 + (qJD(4) * t22 + t12) * t86;
t3 = 0.2e1 * t10 * t22 + 0.2e1 * t11 * t23 + 0.2e1 * t158;
t2 = t99 * t131 + (t98 * qJD(5) - t10 * t85 - t11 * t82) * t86;
t1 = (-t98 * qJD(4) - t21) * t86 + (t4 + t138) * t83;
t33 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t80 ^ 2 * t84 * t139 - 0.2e1 * t20 * t39 + 0.2e1 * t158, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t117, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t117, -pkin(2) * t67 + t145, 0, 0, 0, 0, 0, 0 (t83 * t139 + t84 * t73) * t80 (-t84 * t131 + t86 * t139) * t80, -t7, t7 * t88 + t145, 0, 0, 0, 0, 0, 0, t5, t6, t2, t10 * t31 + t11 * t32 - t23 * t14 + t22 * t15 + (t38 * t131 - t157) * t88, 0, 0, 0, 0, 0, 0, t5, t6, t2, t10 * t16 + t11 * t19 + t21 * t50 + t22 * t8 - t23 * t9 + t28 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, t75, -0.2e1 * t113, 0.2e1 * t107, 0, t65, 0, 0, 0.2e1 * qJD(3) * t83 + 0.2e1 * t86 * t129, 0.2e1 * qJD(3) * t86 - 0.2e1 * t83 * t129, 0, t75, t30, t18, t24, t29, t25, t65, -0.2e1 * t79 * t85 * t133 + 0.2e1 * t15 * t83 + 0.2e1 * (t31 + 0.2e1 * t127) * t73, 0.2e1 * t79 * t120 + 0.2e1 * t14 * t83 + 0.2e1 * (-t32 + 0.2e1 * t62) * t73, 0.2e1 * t97 * t131 + 0.2e1 * (t96 * qJD(5) + t14 * t82 - t15 * t85) * t86, -0.2e1 * t88 ^ 2 * t113 - 0.2e1 * t32 * t14 + 0.2e1 * t31 * t15, t30, t18, t24, t29, t25, t65, 0.2e1 * (-t50 * t137 + t8) * t83 + 0.2e1 * (qJD(4) * t16 + t50 * t74 + t156) * t86, 0.2e1 * (-t50 * t136 + t9) * t83 + 0.2e1 * (-qJD(4) * t19 - t118 + t155) * t86, 0.2e1 * t101 * t131 + 0.2e1 * (t100 * qJD(5) - t8 * t85 + t82 * t9) * t86, 0.2e1 * t16 * t8 - 0.2e1 * t19 * t9 + 0.2e1 * t28 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t41, 0, -t96 * t73 + (t92 - 0.2e1 * t112) * t83, 0, 0, 0, 0, 0, 0, -t44, t41, 0 (-t100 * qJD(4) - t28) * t86 + (qJD(4) * t50 - qJD(5) * t101 - t8 * t82 - t85 * t9) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t20, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12, t4, -t21 * pkin(4) + t4 * pkin(9), 0, 0, 0, 0, 0, 0, t13, t12, t4, pkin(5) * t122 + t10 * t56 - t11 * t57 + t21 * t71 + t22 * t35 - t23 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, 0, -t73, 0, -t114, -t112, 0, 0, -t27, t26, t45, t27, -t42, 0 (-t104 * t85 - t82 * t147) * qJD(5) + (t103 * t82 - t62) * qJD(4) (t104 * t82 - t85 * t147) * qJD(5) + (-t85 * t160 + (pkin(4) * t85 + t149) * t83) * qJD(4), t92, -pkin(4) * t114 + t92 * pkin(9), -t27, t26, t45, t27, -t42, 0, -t155 + t35 * t83 + (-t71 * t82 * t83 + t154) * qJD(4) + (t102 * t86 + t50 * t82) * qJD(5), t156 + t34 * t83 + (-t83 * t152 + t153) * qJD(4) + (t110 * t86 * t82 + t50 * t85) * qJD(5) (t56 * t131 - t35 * t86 - t9 + (-t16 + t153) * qJD(5)) * t85 + (-t57 * t131 + t34 * t86 - t8 + (-t19 + t154) * qJD(5)) * t82, pkin(5) * t118 + t16 * t35 - t19 * t34 + t28 * t71 + t56 * t8 + t57 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t73, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t46, t40 (pkin(9) * t108 - t159) * qJD(4), 0, 0, 0, 0, 0, 0, -t94, -t46, t40 (-t124 + (-t56 * t82 - t57 * t85) * qJD(4)) * t86 + (qJD(4) * t71 + t91) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t51, 0, t64, 0, 0, pkin(4) * t126, -0.2e1 * pkin(4) * t74, 0, 0, t63, t51, 0, t64, 0, 0, t110 * t126, t102 * t162, 0.2e1 * t91, 0.2e1 * t124 * t71 + 0.2e1 * t34 * t57 + 0.2e1 * t35 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, 0, t10 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, 0, -t46, t73, t15, t14, 0, 0, 0, 0, -t94, 0, -t46, t73 (0.2e1 * pkin(5) - t149) * t73 + t89, t9, t94 * pkin(5), t8 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t42, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t42, 0, -t45 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, -t135, 0, -pkin(9) * t74, pkin(9) * t135, 0, 0, 0, 0, t74, 0, -t135, 0, t35, t34, -pkin(5) * t74, t35 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t94, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t74, 0, t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t33;
