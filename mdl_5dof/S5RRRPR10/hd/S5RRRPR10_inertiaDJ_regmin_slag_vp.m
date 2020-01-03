% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:51
% EndTime: 2019-12-31 21:29:56
% DurationCPUTime: 1.50s
% Computational Cost: add. (2512->222), mult. (6932->451), div. (0->0), fcn. (6643->10), ass. (0->124)
t100 = cos(qJ(2));
t93 = sin(pkin(5));
t137 = t100 * t93;
t127 = pkin(7) * t137;
t97 = sin(qJ(2));
t150 = pkin(1) * t97;
t94 = cos(pkin(5));
t62 = t127 + (pkin(8) + t150) * t94;
t63 = (-pkin(2) * t100 - pkin(8) * t97 - pkin(1)) * t93;
t96 = sin(qJ(3));
t99 = cos(qJ(3));
t139 = t99 * t62 + t96 * t63;
t152 = 0.2e1 * t93;
t151 = 2 * qJD(5);
t149 = pkin(8) * t93;
t135 = qJD(2) * t97;
t124 = t93 * t135;
t136 = cos(pkin(10));
t130 = qJD(2) * t100;
t121 = t93 * t130;
t142 = t93 * t97;
t72 = t99 * t142 + t94 * t96;
t49 = t72 * qJD(3) + t96 * t121;
t71 = t96 * t142 - t94 * t99;
t50 = -t71 * qJD(3) + t99 * t121;
t92 = sin(pkin(10));
t33 = t136 * t50 - t92 * t49;
t44 = t136 * t72 - t92 * t71;
t95 = sin(qJ(5));
t98 = cos(qJ(5));
t35 = t98 * t137 + t95 * t44;
t20 = -t35 * qJD(5) + t95 * t124 + t98 * t33;
t148 = t20 * t95;
t118 = t136 * t96;
t75 = t92 * t99 + t118;
t147 = t75 * t95;
t146 = t75 * t98;
t85 = pkin(3) * t92 + pkin(9);
t145 = t85 * t95;
t144 = t85 * t98;
t143 = t92 * t96;
t141 = t95 * t98;
t140 = -qJ(4) - pkin(8);
t65 = (pkin(2) * t97 - pkin(8) * t100) * t93 * qJD(2);
t66 = -t94 * pkin(1) * t130 + pkin(7) * t124;
t26 = -t139 * qJD(3) + t99 * t65 + t66 * t96;
t15 = pkin(3) * t124 - qJ(4) * t50 - qJD(4) * t72 + t26;
t133 = qJD(3) * t99;
t134 = qJD(3) * t96;
t25 = -t63 * t133 + t62 * t134 - t96 * t65 + t99 * t66;
t19 = -qJ(4) * t49 - qJD(4) * t71 - t25;
t6 = t136 * t19 + t92 * t15;
t119 = -t96 * t62 + t99 * t63;
t28 = -pkin(3) * t137 - t72 * qJ(4) + t119;
t34 = -qJ(4) * t71 + t139;
t18 = t136 * t34 + t92 * t28;
t91 = t98 ^ 2;
t138 = t95 ^ 2 - t91;
t132 = qJD(5) * t95;
t131 = qJD(5) * t98;
t129 = qJD(3) * t100;
t128 = -0.2e1 * pkin(2) * qJD(3);
t86 = -t136 * pkin(3) - pkin(4);
t126 = t86 * t151;
t88 = pkin(3) * t134;
t125 = t95 * t137;
t123 = t95 * t131;
t87 = -pkin(3) * t99 - pkin(2);
t89 = t93 ^ 2;
t122 = t89 * t130;
t120 = -0.4e1 * t75 * t141;
t117 = t136 * t99;
t116 = qJD(3) * t140;
t115 = t138 * qJD(5);
t114 = t97 * t122;
t14 = -pkin(9) * t137 + t18;
t43 = t136 * t71 + t72 * t92;
t61 = pkin(7) * t142 + (-pkin(1) * t100 - pkin(2)) * t94;
t46 = t71 * pkin(3) + t61;
t23 = t43 * pkin(4) - t44 * pkin(9) + t46;
t8 = t14 * t98 + t23 * t95;
t36 = t98 * t44 - t125;
t113 = -t35 * t98 - t36 * t95;
t74 = -t117 + t143;
t47 = pkin(4) * t74 - pkin(9) * t75 + t87;
t79 = t140 * t99;
t52 = -t136 * t79 + t140 * t143;
t30 = t95 * t47 + t98 * t52;
t69 = t75 * qJD(3);
t70 = qJD(3) * t117 - t92 * t134;
t112 = -t69 * t85 + t70 * t86;
t111 = t74 * t85 - t75 * t86;
t110 = pkin(4) * t69 - pkin(9) * t70 + t88;
t32 = t136 * t49 + t50 * t92;
t109 = t43 * t131 + t32 * t95;
t108 = t74 * t131 + t69 * t95;
t107 = t75 * t131 + t70 * t95;
t106 = -t75 * t132 + t70 * t98;
t105 = pkin(9) * t124 + t6;
t5 = t136 * t15 - t92 * t19;
t17 = t136 * t28 - t92 * t34;
t104 = t96 * t129 + t99 * t135;
t103 = -t99 * t129 + t96 * t135;
t102 = -qJD(4) * t96 + t99 * t116;
t67 = (t94 * t150 + t127) * qJD(2);
t37 = t49 * pkin(3) + t67;
t101 = t32 * pkin(4) - t33 * pkin(9) + t37;
t73 = t75 ^ 2;
t68 = qJD(4) * t99 + t96 * t116;
t51 = -t140 * t118 - t79 * t92;
t45 = -t74 * t132 + t69 * t98;
t41 = t92 * t102 + t136 * t68;
t40 = -t136 * t102 + t68 * t92;
t29 = t47 * t98 - t52 * t95;
t24 = -t43 * t132 + t32 * t98;
t21 = -qJD(5) * t125 - t98 * t124 + t44 * t131 + t33 * t95;
t13 = pkin(4) * t137 - t17;
t11 = -t30 * qJD(5) + t98 * t110 - t95 * t41;
t10 = -t95 * t110 - t47 * t131 + t52 * t132 - t98 * t41;
t7 = -t14 * t95 + t23 * t98;
t4 = -pkin(4) * t124 - t5;
t2 = -t8 * qJD(5) + t98 * t101 - t95 * t105;
t1 = -t95 * t101 - t98 * t105 - t23 * t131 + t14 * t132;
t3 = [0, 0, 0, 0.2e1 * t114, 0.2e1 * (t100 ^ 2 - t97 ^ 2) * t89 * qJD(2), 0.2e1 * t94 * t121, -0.2e1 * t94 * t124, 0, -0.2e1 * pkin(1) * t89 * t135 - 0.2e1 * t67 * t94, -0.2e1 * pkin(1) * t122 + 0.2e1 * t66 * t94, 0.2e1 * t72 * t50, -0.2e1 * t49 * t72 - 0.2e1 * t50 * t71, (-t100 * t50 + t72 * t135) * t152, (t100 * t49 - t71 * t135) * t152, -0.2e1 * t114, 0.2e1 * t61 * t49 + 0.2e1 * t67 * t71 + 0.2e1 * (-t26 * t100 + t119 * t135) * t93, 0.2e1 * t61 * t50 + 0.2e1 * t67 * t72 + 0.2e1 * (-t25 * t100 - t139 * t135) * t93, -0.2e1 * t17 * t33 - 0.2e1 * t18 * t32 - 0.2e1 * t43 * t6 - 0.2e1 * t44 * t5, 0.2e1 * t17 * t5 + 0.2e1 * t18 * t6 + 0.2e1 * t37 * t46, 0.2e1 * t36 * t20, -0.2e1 * t20 * t35 - 0.2e1 * t21 * t36, 0.2e1 * t20 * t43 + 0.2e1 * t32 * t36, -0.2e1 * t21 * t43 - 0.2e1 * t32 * t35, 0.2e1 * t43 * t32, 0.2e1 * t13 * t21 + 0.2e1 * t2 * t43 + 0.2e1 * t32 * t7 + 0.2e1 * t35 * t4, 0.2e1 * t1 * t43 + 0.2e1 * t13 * t20 - 0.2e1 * t32 * t8 + 0.2e1 * t36 * t4; 0, 0, 0, 0, 0, t121, -t124, 0, -t67, t66, t72 * t133 + t50 * t96, -t49 * t96 + t50 * t99 + (-t71 * t99 - t72 * t96) * qJD(3), t103 * t93, t104 * t93, 0, -pkin(2) * t49 - t103 * t149 + t61 * t134 - t67 * t99, -pkin(2) * t50 - t104 * t149 + t61 * t133 + t67 * t96, -t17 * t70 - t18 * t69 - t32 * t52 + t33 * t51 + t40 * t44 - t41 * t43 - t5 * t75 - t6 * t74, -t17 * t40 + t18 * t41 + t37 * t87 + t46 * t88 - t5 * t51 + t52 * t6, t106 * t36 + t20 * t146, t113 * t70 + (-t148 - t21 * t98 + (t35 * t95 - t36 * t98) * qJD(5)) * t75, t106 * t43 + t32 * t146 + t20 * t74 + t36 * t69, -t107 * t43 - t32 * t147 - t21 * t74 - t35 * t69, t32 * t74 + t43 * t69, t107 * t13 + t11 * t43 + t4 * t147 + t2 * t74 + t21 * t51 + t29 * t32 + t35 * t40 + t69 * t7, t1 * t74 + t10 * t43 + t106 * t13 + t4 * t146 + t20 * t51 - t30 * t32 + t36 * t40 - t69 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t96 * t133, 0.2e1 * (-t96 ^ 2 + t99 ^ 2) * qJD(3), 0, 0, 0, t96 * t128, t99 * t128, 0.2e1 * t40 * t75 - 0.2e1 * t41 * t74 + 0.2e1 * t51 * t70 - 0.2e1 * t52 * t69, 0.2e1 * t40 * t51 + 0.2e1 * t41 * t52 + 0.2e1 * t87 * t88, 0.2e1 * t70 * t75 * t91 - 0.2e1 * t73 * t123, t138 * t73 * t151 + t70 * t120, 0.2e1 * t106 * t74 + 0.2e1 * t69 * t146, -0.2e1 * t107 * t74 - 0.2e1 * t69 * t147, 0.2e1 * t74 * t69, 0.2e1 * t107 * t51 + 0.2e1 * t11 * t74 + 0.2e1 * t40 * t147 + 0.2e1 * t29 * t69, 0.2e1 * t10 * t74 + 0.2e1 * t106 * t51 + 0.2e1 * t40 * t146 - 0.2e1 * t30 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t49, t124, t26, t25, (-t136 * t33 - t32 * t92) * pkin(3), (t136 * t5 + t6 * t92) * pkin(3), t36 * t131 + t148, t113 * qJD(5) + t20 * t98 - t21 * t95, t109, t24, 0, -t32 * t145 + t21 * t86 - t4 * t98 + (t13 * t95 - t43 * t144) * qJD(5), -t32 * t144 + t20 * t86 + t4 * t95 + (t13 * t98 + t43 * t145) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, -t134, 0, -pkin(8) * t133, pkin(8) * t134, (-t136 * t70 - t69 * t92) * pkin(3), (-t136 * t40 + t41 * t92) * pkin(3), -t75 * t115 + t70 * t141, qJD(5) * t120 - t138 * t70, t108, t45, 0, -t40 * t98 + t112 * t95 + (-t111 * t98 + t51 * t95) * qJD(5), t40 * t95 + t112 * t98 + (t111 * t95 + t51 * t98) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t123, -0.2e1 * t115, 0, 0, 0, t95 * t126, t98 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, t24, -t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, t45, -t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t21, t32, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, -t107, t69, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, -t132, 0, -t85 * t131, t85 * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
