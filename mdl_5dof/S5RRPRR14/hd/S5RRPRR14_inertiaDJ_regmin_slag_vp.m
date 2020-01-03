% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR14_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:38
% EndTime: 2019-12-31 20:38:43
% DurationCPUTime: 1.54s
% Computational Cost: add. (2269->216), mult. (6327->436), div. (0->0), fcn. (6335->10), ass. (0->116)
t90 = sin(pkin(5));
t98 = cos(qJ(2));
t140 = t90 * t98;
t127 = pkin(7) * t140;
t95 = sin(qJ(2));
t146 = pkin(1) * t95;
t92 = cos(pkin(5));
t59 = t127 + (qJ(3) + t146) * t92;
t112 = -pkin(2) * t98 - qJ(3) * t95;
t60 = (-pkin(1) + t112) * t90;
t89 = sin(pkin(10));
t91 = cos(pkin(10));
t37 = -t89 * t59 + t91 * t60;
t141 = t90 * t95;
t66 = t91 * t141 + t89 * t92;
t23 = -pkin(3) * t140 - t66 * pkin(8) + t37;
t38 = t91 * t59 + t89 * t60;
t65 = t89 * t141 - t92 * t91;
t31 = -pkin(8) * t65 + t38;
t94 = sin(qJ(4));
t97 = cos(qJ(4));
t149 = t94 * t23 + t97 * t31;
t96 = cos(qJ(5));
t88 = t96 ^ 2;
t93 = sin(qJ(5));
t133 = t93 ^ 2 - t88;
t119 = t133 * qJD(5);
t139 = t91 * t98;
t54 = (-qJD(3) * t95 + (pkin(2) * t95 - qJ(3) * t98) * qJD(2)) * t90;
t132 = qJD(2) * t95;
t123 = t90 * t132;
t131 = qJD(2) * t98;
t63 = -t92 * pkin(1) * t131 + pkin(7) * t123;
t58 = qJD(3) * t92 - t63;
t32 = t91 * t54 - t89 * t58;
t21 = (pkin(3) * t95 - pkin(8) * t139) * t90 * qJD(2) + t32;
t122 = t90 * t131;
t117 = t89 * t122;
t33 = t89 * t54 + t91 * t58;
t30 = -pkin(8) * t117 + t33;
t6 = -qJD(4) * t149 + t97 * t21 - t94 * t30;
t83 = -pkin(3) * t91 - pkin(2);
t148 = 0.2e1 * t83;
t147 = 0.2e1 * t90;
t41 = t97 * t65 + t66 * t94;
t136 = t97 * t91;
t71 = t89 * t94 - t136;
t27 = -qJD(4) * t41 - t71 * t122;
t42 = -t65 * t94 + t66 * t97;
t34 = t96 * t140 + t93 * t42;
t14 = -qJD(5) * t34 + t93 * t123 + t96 * t27;
t145 = t14 * t93;
t67 = t71 * qJD(4);
t72 = t89 * t97 + t91 * t94;
t144 = t72 * t67;
t143 = t72 * t93;
t142 = t72 * t96;
t138 = t93 * t96;
t137 = t94 * t31;
t135 = pkin(8) + qJ(3);
t130 = qJD(4) * t97;
t129 = qJD(5) * t93;
t128 = qJD(5) * t96;
t126 = -0.2e1 * pkin(4) * qJD(5);
t125 = t93 * t140;
t85 = t90 ^ 2;
t124 = t85 * t131;
t121 = t93 * t128;
t120 = -0.4e1 * t72 * t138;
t118 = t95 * t124;
t116 = 0.2e1 * (t89 ^ 2 + t91 ^ 2) * qJD(3);
t68 = t72 * qJD(4);
t115 = pkin(4) * t68 + pkin(9) * t67;
t114 = pkin(4) * t67 - pkin(9) * t68;
t113 = pkin(4) * t72 + pkin(9) * t71;
t13 = -pkin(9) * t140 + t149;
t62 = pkin(7) * t141 + (-pkin(1) * t98 - pkin(2)) * t92;
t46 = t65 * pkin(3) + t62;
t17 = t41 * pkin(4) - t42 * pkin(9) + t46;
t8 = t13 * t96 + t17 * t93;
t110 = t97 * t23 - t137;
t108 = -t32 * t89 + t33 * t91;
t35 = t96 * t42 - t125;
t107 = -t34 * t96 - t35 * t93;
t47 = pkin(4) * t71 - pkin(9) * t72 + t83;
t76 = t135 * t89;
t77 = t135 * t91;
t52 = -t76 * t94 + t77 * t97;
t26 = t47 * t93 + t52 * t96;
t28 = qJD(4) * t42 + t72 * t122;
t106 = t41 * t128 + t28 * t93;
t105 = t41 * t129 - t28 * t96;
t104 = t72 * t128 - t67 * t93;
t103 = -t72 * t129 - t67 * t96;
t102 = t71 * t128 + t68 * t93;
t5 = qJD(4) * t137 - t23 * t130 - t94 * t21 - t97 * t30;
t64 = (t92 * t146 + t127) * qJD(2);
t50 = pkin(3) * t117 + t64;
t101 = pkin(9) * t123 - t5;
t100 = (t112 * qJD(2) + qJD(3) * t98) * t90;
t99 = t28 * pkin(4) - t27 * pkin(9) + t50;
t70 = t72 ^ 2;
t51 = t76 * t97 + t77 * t94;
t45 = -t71 * t129 + t68 * t96;
t40 = qJD(3) * t72 + qJD(4) * t52;
t39 = t76 * t130 - qJD(3) * t136 + (qJD(3) * t89 + qJD(4) * t77) * t94;
t25 = t47 * t96 - t52 * t93;
t15 = -qJD(5) * t125 - t96 * t123 + t42 * t128 + t27 * t93;
t12 = pkin(4) * t140 - t110;
t11 = -t26 * qJD(5) + t96 * t115 + t93 * t39;
t10 = -t93 * t115 - t47 * t128 + t52 * t129 + t96 * t39;
t7 = -t13 * t93 + t17 * t96;
t4 = -pkin(4) * t123 - t6;
t2 = -t8 * qJD(5) - t93 * t101 + t96 * t99;
t1 = -t101 * t96 - t17 * t128 + t13 * t129 - t93 * t99;
t3 = [0, 0, 0, 0.2e1 * t118, 0.2e1 * (-t95 ^ 2 + t98 ^ 2) * t85 * qJD(2), 0.2e1 * t92 * t122, -0.2e1 * t92 * t123, 0, -0.2e1 * pkin(1) * t85 * t132 - 0.2e1 * t64 * t92, -0.2e1 * pkin(1) * t124 + 0.2e1 * t63 * t92, 0.2e1 * t64 * t65 + 0.2e1 * (-t32 * t98 + (t62 * t89 * t98 + t37 * t95) * qJD(2)) * t90, 0.2e1 * t64 * t66 + 0.2e1 * (t33 * t98 + (t62 * t139 - t38 * t95) * qJD(2)) * t90, -0.2e1 * t32 * t66 - 0.2e1 * t33 * t65 + 0.2e1 * (-t37 * t91 - t38 * t89) * t122, 0.2e1 * t32 * t37 + 0.2e1 * t33 * t38 + 0.2e1 * t62 * t64, 0.2e1 * t42 * t27, -0.2e1 * t27 * t41 - 0.2e1 * t28 * t42, (t42 * t132 - t27 * t98) * t147, (-t41 * t132 + t28 * t98) * t147, -0.2e1 * t118, 0.2e1 * t46 * t28 + 0.2e1 * t50 * t41 + 0.2e1 * (t110 * t132 - t6 * t98) * t90, 0.2e1 * t46 * t27 + 0.2e1 * t50 * t42 + 0.2e1 * (-t132 * t149 - t5 * t98) * t90, 0.2e1 * t35 * t14, -0.2e1 * t14 * t34 - 0.2e1 * t15 * t35, 0.2e1 * t14 * t41 + 0.2e1 * t28 * t35, -0.2e1 * t15 * t41 - 0.2e1 * t28 * t34, 0.2e1 * t41 * t28, 0.2e1 * t12 * t15 + 0.2e1 * t2 * t41 + 0.2e1 * t28 * t7 + 0.2e1 * t34 * t4, 0.2e1 * t1 * t41 + 0.2e1 * t12 * t14 - 0.2e1 * t28 * t8 + 0.2e1 * t35 * t4; 0, 0, 0, 0, 0, t122, -t123, 0, -t64, t63, t100 * t89 - t64 * t91, t100 * t91 + t64 * t89, (-t65 * t91 + t66 * t89) * qJD(3) + t108, -pkin(2) * t64 + (-t37 * t89 + t38 * t91) * qJD(3) + t108 * qJ(3), t27 * t72 - t42 * t67, -t27 * t71 - t28 * t72 + t41 * t67 - t42 * t68, (t72 * t132 + t67 * t98) * t90, (-t71 * t132 + t68 * t98) * t90, 0, t83 * t28 + t46 * t68 + t50 * t71 + (-t51 * t132 + t40 * t98) * t90, t83 * t27 - t46 * t67 + t50 * t72 + (-t52 * t132 - t39 * t98) * t90, t103 * t35 + t14 * t142, -t107 * t67 + (-t145 - t15 * t96 + (t34 * t93 - t35 * t96) * qJD(5)) * t72, t103 * t41 + t14 * t71 + t28 * t142 + t35 * t68, -t104 * t41 - t28 * t143 - t15 * t71 - t34 * t68, t28 * t71 + t41 * t68, t104 * t12 + t11 * t41 + t4 * t143 + t15 * t51 + t2 * t71 + t25 * t28 + t34 * t40 + t68 * t7, t1 * t71 + t10 * t41 + t103 * t12 + t14 * t51 + t4 * t142 - t26 * t28 + t35 * t40 - t68 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, qJ(3) * t116, -0.2e1 * t144, 0.2e1 * t67 * t71 - 0.2e1 * t68 * t72, 0, 0, 0, t68 * t148, -t67 * t148, -0.2e1 * t70 * t121 - 0.2e1 * t88 * t144, 0.2e1 * t70 * t119 - t67 * t120, 0.2e1 * t103 * t71 + 0.2e1 * t68 * t142, -0.2e1 * t104 * t71 - 0.2e1 * t68 * t143, 0.2e1 * t71 * t68, 0.2e1 * t104 * t51 + 0.2e1 * t11 * t71 + 0.2e1 * t40 * t143 + 0.2e1 * t25 * t68, 0.2e1 * t10 * t71 + 0.2e1 * t103 * t51 + 0.2e1 * t40 * t142 - 0.2e1 * t26 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, t91 * t122, 0, t64, 0, 0, 0, 0, 0, t28, t27, 0, 0, 0, 0, 0, -t105, -t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t67, 0, 0, 0, 0, 0, t45, -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t28, t123, t6, t5, t35 * t128 + t145, qJD(5) * t107 + t14 * t96 - t93 * t15, t106, -t105, 0, -pkin(4) * t15 - pkin(9) * t106 + t12 * t129 - t4 * t96, -pkin(4) * t14 + pkin(9) * t105 + t12 * t128 + t4 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t68, 0, -t40, t39, -t72 * t119 - t67 * t138, qJD(5) * t120 + t133 * t67, t102, t45, 0, -t40 * t96 + t114 * t93 + (-t113 * t96 + t51 * t93) * qJD(5), t40 * t93 + t114 * t96 + (t113 * t93 + t51 * t96) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t121, -0.2e1 * t119, 0, 0, 0, t93 * t126, t96 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, t28, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, -t104, t68, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, -t129, 0, -pkin(9) * t128, pkin(9) * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
