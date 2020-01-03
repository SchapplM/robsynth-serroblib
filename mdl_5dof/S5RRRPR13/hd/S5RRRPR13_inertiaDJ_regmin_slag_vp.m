% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR13_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:33
% EndTime: 2019-12-31 21:46:39
% DurationCPUTime: 1.86s
% Computational Cost: add. (1475->229), mult. (4247->447), div. (0->0), fcn. (3805->8), ass. (0->134)
t81 = sin(qJ(3));
t136 = qJD(3) * t81;
t78 = sin(pkin(5));
t85 = cos(qJ(2));
t148 = t78 * t85;
t126 = pkin(7) * t148;
t82 = sin(qJ(2));
t157 = pkin(1) * t82;
t79 = cos(pkin(5));
t45 = t126 + (pkin(8) + t157) * t79;
t46 = (-pkin(2) * t85 - pkin(8) * t82 - pkin(1)) * t78;
t47 = (pkin(2) * t82 - pkin(8) * t85) * t78 * qJD(2);
t138 = qJD(2) * t85;
t139 = qJD(2) * t82;
t66 = t78 * t139;
t48 = -t79 * pkin(1) * t138 + pkin(7) * t66;
t84 = cos(qJ(3));
t72 = qJD(3) * t84;
t11 = t45 * t136 - t46 * t72 - t81 * t47 + t84 * t48;
t95 = (qJ(4) * t139 - qJD(4) * t85) * t78;
t8 = t11 - t95;
t159 = pkin(3) + pkin(9);
t83 = cos(qJ(5));
t145 = t83 * t159;
t140 = t81 * qJ(4);
t163 = t159 * t84 + t140;
t77 = t84 ^ 2;
t113 = qJD(3) * (t81 ^ 2 - t77);
t80 = sin(qJ(5));
t74 = t80 ^ 2;
t143 = -t83 ^ 2 + t74;
t112 = t143 * qJD(5);
t133 = qJD(4) * t84;
t158 = pkin(4) + pkin(8);
t63 = t158 * t84;
t162 = qJD(3) * t163 - qJD(5) * t63 - t133;
t161 = 0.2e1 * t78;
t160 = 0.2e1 * qJD(4);
t156 = pkin(8) * t78;
t118 = t78 * t138;
t149 = t78 * t82;
t52 = t84 * t149 + t79 * t81;
t30 = qJD(3) * t52 + t81 * t118;
t155 = t30 * pkin(3);
t6 = -t30 * pkin(4) - t8;
t154 = t6 * t80;
t153 = t6 * t83;
t147 = t80 * t85;
t124 = t81 * t149;
t51 = -t79 * t84 + t124;
t32 = t78 * t147 + t51 * t83;
t18 = t32 * qJD(5) + t30 * t80 + t83 * t66;
t152 = t18 * t80;
t151 = t18 * t83;
t31 = -qJD(3) * t124 + (qJD(3) * t79 + t118) * t84;
t150 = t31 * t81;
t146 = t80 * t159;
t144 = t84 * t45 + t81 * t46;
t141 = qJ(4) * t84;
t137 = qJD(3) * t80;
t135 = qJD(3) * t83;
t134 = qJD(3) * t85;
t132 = qJD(5) * t80;
t131 = qJD(5) * t83;
t130 = qJD(5) * t84;
t129 = qJD(5) * t159;
t128 = qJ(4) * qJD(5);
t127 = t79 * t157;
t125 = -0.2e1 * pkin(2) * qJD(3);
t123 = t83 * t148;
t122 = pkin(8) * t136;
t73 = t78 ^ 2;
t121 = t73 * t138;
t120 = t80 * t130;
t119 = t83 * t130;
t117 = t81 * t72;
t116 = t81 * t135;
t115 = t80 * t131;
t114 = -t81 * t45 + t84 * t46;
t111 = pkin(3) * t136 - t81 * qJD(4);
t110 = pkin(3) * t66;
t109 = t82 * t121;
t108 = t80 * t116;
t23 = pkin(3) * t148 - t114;
t71 = pkin(8) * t72;
t107 = pkin(4) * t72 + t71;
t106 = -t84 * pkin(3) - t140;
t14 = t52 * pkin(4) + pkin(9) * t148 + t23;
t44 = pkin(7) * t149 + (-pkin(1) * t85 - pkin(2)) * t79;
t93 = -t52 * qJ(4) + t44;
t19 = t159 * t51 + t93;
t4 = t80 * t14 + t83 * t19;
t33 = t51 * t80 - t123;
t105 = t32 * t80 + t33 * t83;
t56 = -pkin(2) - t163;
t62 = t158 * t81;
t29 = t83 * t56 + t80 * t62;
t12 = -t46 * t136 - t45 * t72 + t84 * t47 + t81 * t48;
t103 = -t31 * qJ(4) - t52 * qJD(4);
t22 = qJ(4) * t148 - t144;
t101 = -t52 * t131 - t31 * t80;
t25 = -t52 * t132 + t31 * t83;
t99 = t81 * t134 + t84 * t139;
t98 = -t84 * t134 + t81 * t139;
t97 = t80 * t136 - t119;
t96 = t31 * pkin(4) - t12;
t49 = (t126 + t127) * qJD(2);
t58 = t158 * t136;
t94 = -t58 + (t159 * t81 - t141) * qJD(5);
t92 = t99 * t78;
t91 = t98 * t78;
t90 = t106 * qJD(3) + t133;
t89 = -(pkin(9) * t81 - t141) * qJD(3) - t111;
t9 = -t12 - t110;
t88 = -t8 * t84 + t9 * t81 + (t22 * t81 + t23 * t84) * qJD(3);
t87 = t103 + t49;
t65 = 0.2e1 * t117;
t59 = -pkin(2) + t106;
t54 = -t81 * t131 - t80 * t72;
t53 = -t81 * t132 + t83 * t72;
t50 = -qJ(4) * t72 + t111;
t28 = -t80 * t56 + t83 * t62;
t26 = 0.2e1 * t52 * t31;
t24 = t52 * t72 + t150;
t21 = t51 * pkin(3) + t93;
t20 = -t51 * pkin(4) - t22;
t17 = -t30 * t83 - qJD(5) * t123 + (qJD(5) * t51 + t66) * t80;
t16 = -t29 * qJD(5) + t83 * t107 + t80 * t89;
t15 = -t80 * t107 - t62 * t131 + t56 * t132 + t83 * t89;
t10 = t87 + t155;
t3 = t83 * t14 - t80 * t19;
t2 = -t80 * (t30 * pkin(9) + t103 + t155) + t83 * t96 - t4 * qJD(5) + (-t80 * t127 + (-pkin(7) * t147 - t145 * t82) * t78) * qJD(2);
t1 = t19 * t132 - t14 * t131 - t80 * (-t159 * t66 + t96) - t83 * (t159 * t30 + t87);
t5 = [0, 0, 0, 0.2e1 * t109, 0.2e1 * (-t82 ^ 2 + t85 ^ 2) * t73 * qJD(2), 0.2e1 * t79 * t118, -0.2e1 * t79 * t66, 0, -0.2e1 * t73 * pkin(1) * t139 - 0.2e1 * t49 * t79, -0.2e1 * pkin(1) * t121 + 0.2e1 * t48 * t79, t26, -0.2e1 * t52 * t30 - 0.2e1 * t31 * t51, (t52 * t139 - t31 * t85) * t161, (-t51 * t139 + t30 * t85) * t161, -0.2e1 * t109, 0.2e1 * t44 * t30 + 0.2e1 * t49 * t51 + 0.2e1 * (t114 * t139 - t12 * t85) * t78, 0.2e1 * t44 * t31 + 0.2e1 * t49 * t52 + 0.2e1 * (-t11 * t85 - t144 * t139) * t78, 0.2e1 * t22 * t30 + 0.2e1 * t23 * t31 + 0.2e1 * t8 * t51 + 0.2e1 * t9 * t52, -0.2e1 * t10 * t51 - 0.2e1 * t21 * t30 + 0.2e1 * (t23 * t139 - t85 * t9) * t78, -0.2e1 * t10 * t52 - 0.2e1 * t21 * t31 + 0.2e1 * (-t22 * t139 + t8 * t85) * t78, 0.2e1 * t21 * t10 + 0.2e1 * t22 * t8 + 0.2e1 * t23 * t9, 0.2e1 * t33 * t18, -0.2e1 * t33 * t17 + 0.2e1 * t18 * t32, 0.2e1 * t18 * t52 + 0.2e1 * t33 * t31, -0.2e1 * t17 * t52 + 0.2e1 * t32 * t31, t26, 0.2e1 * t20 * t17 + 0.2e1 * t2 * t52 + 0.2e1 * t3 * t31 - 0.2e1 * t6 * t32, 0.2e1 * t1 * t52 + 0.2e1 * t20 * t18 - 0.2e1 * t4 * t31 + 0.2e1 * t6 * t33; 0, 0, 0, 0, 0, t118, -t66, 0, -t49, t48, t24, -t81 * t30 + t31 * t84 + (-t51 * t84 - t52 * t81) * qJD(3), t91, t92, 0, -pkin(2) * t30 + t44 * t136 - t98 * t156 - t49 * t84, -pkin(2) * t31 - t99 * t156 + t44 * t72 + t49 * t81, (-t30 * t84 + t150 + (t51 * t81 + t52 * t84) * qJD(3)) * pkin(8) + t88, pkin(8) * t91 + t10 * t84 - t21 * t136 - t59 * t30 - t50 * t51, pkin(8) * t92 - t10 * t81 - t21 * t72 - t59 * t31 - t50 * t52, pkin(8) * t88 + t10 * t59 + t21 * t50, -t84 * t152 + t33 * t97, t105 * t136 + (t17 * t80 - t151 + (-t32 * t83 + t33 * t80) * qJD(5)) * t84, (t52 * t137 + t18) * t81 + (qJD(3) * t33 + t101) * t84, (t52 * t135 - t17) * t81 + (qJD(3) * t32 - t25) * t84, t24, t16 * t52 + t63 * t17 + t28 * t31 + t58 * t32 + (-t20 * t135 + t2) * t81 + (qJD(3) * t3 - t20 * t132 + t153) * t84, t15 * t52 + t63 * t18 - t29 * t31 - t58 * t33 + (t20 * t137 + t1) * t81 + (-qJD(3) * t4 - t20 * t131 - t154) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -0.2e1 * t113, 0, 0, 0, t81 * t125, t84 * t125, 0, -0.2e1 * t59 * t136 + 0.2e1 * t50 * t84, -0.2e1 * t50 * t81 - 0.2e1 * t59 * t72, 0.2e1 * t59 * t50, 0.2e1 * t115 * t77 - 0.2e1 * t117 * t74, -0.4e1 * t84 * t108 - 0.2e1 * t77 * t112, 0.2e1 * t113 * t80 - 0.2e1 * t119 * t81, 0.2e1 * t113 * t83 + 0.2e1 * t120 * t81, t65, 0.2e1 * (-t63 * t135 + t16) * t81 + 0.2e1 * (qJD(3) * t28 - t63 * t132 - t58 * t83) * t84, 0.2e1 * (t63 * t137 + t15) * t81 + 0.2e1 * (-qJD(3) * t29 - t63 * t131 + t58 * t80) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, t66, t12, t11, -t31 * pkin(3) - qJ(4) * t30 - qJD(4) * t51, -t12 - 0.2e1 * t110, -t11 + 0.2e1 * t95, -t9 * pkin(3) - t8 * qJ(4) - t22 * qJD(4), -t33 * t132 + t151, -t105 * qJD(5) - t83 * t17 - t152, t25, t101, 0, -t31 * t145 + qJ(4) * t17 - qJD(4) * t32 + t154 + (t52 * t146 + t20 * t83) * qJD(5), t31 * t146 + qJ(4) * t18 + qJD(4) * t33 + t153 + (t52 * t145 - t20 * t80) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t136, 0, -t71, t122, t90, t71, -t122, t90 * pkin(8), t112 * t84 + t108, 0.4e1 * t84 * t115 - t143 * t136, t53, t54, 0, -t162 * t83 + t94 * t80, t162 * t80 + t94 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, qJ(4) * t160, -0.2e1 * t115, 0.2e1 * t112, 0, 0, 0, 0.2e1 * qJD(4) * t80 + 0.2e1 * t128 * t83, 0.2e1 * qJD(4) * t83 - 0.2e1 * t128 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t66, 0, t9, 0, 0, 0, 0, 0, t25, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, t71, 0, 0, 0, 0, 0, t53, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, t31, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t116 + t120, t72, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t131, 0, t80 * t129, t83 * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
