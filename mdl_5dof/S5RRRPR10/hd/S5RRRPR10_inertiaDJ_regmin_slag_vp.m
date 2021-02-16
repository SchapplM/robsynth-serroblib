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
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 23:41:17
% EndTime: 2021-01-15 23:41:25
% DurationCPUTime: 1.72s
% Computational Cost: add. (2798->244), mult. (7748->489), div. (0->0), fcn. (7373->10), ass. (0->128)
t100 = cos(qJ(3));
t101 = cos(qJ(2));
t94 = sin(pkin(5));
t139 = t101 * t94;
t129 = pkin(7) * t139;
t98 = sin(qJ(2));
t155 = pkin(1) * t98;
t95 = cos(pkin(5));
t63 = t129 + (pkin(8) + t155) * t95;
t64 = (-pkin(2) * t101 - pkin(8) * t98 - pkin(1)) * t94;
t97 = sin(qJ(3));
t141 = t100 * t63 + t97 * t64;
t157 = 0.2e1 * t94;
t156 = 2 * qJD(5);
t154 = pkin(8) * t94;
t137 = qJD(2) * t98;
t126 = t94 * t137;
t138 = cos(pkin(10));
t133 = qJD(2) * t101;
t123 = t94 * t133;
t148 = t94 * t98;
t73 = t100 * t148 + t95 * t97;
t50 = t73 * qJD(3) + t97 * t123;
t72 = -t95 * t100 + t97 * t148;
t51 = -t72 * qJD(3) + t100 * t123;
t93 = sin(pkin(10));
t34 = t138 * t51 - t93 * t50;
t45 = t138 * t73 - t93 * t72;
t96 = sin(qJ(5));
t99 = cos(qJ(5));
t36 = t99 * t139 + t96 * t45;
t21 = -t36 * qJD(5) + t96 * t126 + t99 * t34;
t153 = t21 * t96;
t44 = t138 * t72 + t93 * t73;
t86 = t93 * pkin(3) + pkin(9);
t152 = t44 * t86;
t120 = t138 * t97;
t76 = t93 * t100 + t120;
t151 = t76 * t96;
t150 = t76 * t99;
t149 = t93 * t97;
t33 = t138 * t50 + t93 * t51;
t147 = t96 * t33;
t70 = t76 * qJD(3);
t146 = t96 * t70;
t145 = t99 * t33;
t144 = t99 * t70;
t117 = t138 * t100;
t136 = qJD(3) * t97;
t71 = qJD(3) * t117 - t93 * t136;
t143 = t99 * t71;
t142 = qJ(4) + pkin(8);
t116 = pkin(3) * t126;
t66 = (pkin(2) * t98 - pkin(8) * t101) * t94 * qJD(2);
t67 = -t95 * pkin(1) * t133 + pkin(7) * t126;
t27 = -t141 * qJD(3) + t100 * t66 + t97 * t67;
t16 = -t51 * qJ(4) - t73 * qJD(4) + t116 + t27;
t132 = qJD(3) * t100;
t26 = t100 * t67 - t64 * t132 + t63 * t136 - t97 * t66;
t20 = -t50 * qJ(4) - t72 * qJD(4) - t26;
t6 = t138 * t20 + t93 * t16;
t121 = t100 * t64 - t97 * t63;
t29 = -pkin(3) * t139 - t73 * qJ(4) + t121;
t35 = -t72 * qJ(4) + t141;
t19 = t138 * t35 + t93 * t29;
t92 = t99 ^ 2;
t140 = t96 ^ 2 - t92;
t135 = qJD(5) * t96;
t134 = qJD(5) * t99;
t131 = qJD(3) * t101;
t130 = -0.2e1 * pkin(2) * qJD(3);
t87 = -t138 * pkin(3) - pkin(4);
t128 = t87 * t156;
t89 = pkin(3) * t136;
t127 = t96 * t139;
t125 = t96 * t134;
t90 = t94 ^ 2;
t124 = t90 * t133;
t88 = -t100 * pkin(3) - pkin(2);
t122 = -0.4e1 * t96 * t150;
t5 = t138 * t16 - t93 * t20;
t119 = qJD(3) * t142;
t118 = t140 * qJD(5);
t115 = t98 * t124;
t15 = -pkin(9) * t139 + t19;
t62 = pkin(7) * t148 + (-pkin(1) * t101 - pkin(2)) * t95;
t47 = t72 * pkin(3) + t62;
t24 = t44 * pkin(4) - t45 * pkin(9) + t47;
t8 = t99 * t15 + t96 * t24;
t37 = t99 * t45 - t127;
t114 = -t36 * t99 - t37 * t96;
t75 = -t117 + t149;
t48 = t75 * pkin(4) - t76 * pkin(9) + t88;
t80 = t142 * t100;
t53 = t138 * t80 - t142 * t149;
t31 = t96 * t48 + t99 * t53;
t113 = -t70 * t86 + t71 * t87;
t112 = t75 * t86 - t76 * t87;
t111 = t70 * pkin(4) - t71 * pkin(9) + t89;
t110 = t44 * t134 + t147;
t109 = t75 * t134 + t146;
t108 = t76 * t134 + t96 * t71;
t107 = -t76 * t135 + t143;
t106 = pkin(9) * t126 + t6;
t18 = t138 * t29 - t93 * t35;
t105 = -t100 * t131 + t97 * t137;
t104 = t100 * t137 + t97 * t131;
t103 = -t97 * qJD(4) - t100 * t119;
t68 = (t95 * t155 + t129) * qJD(2);
t38 = t50 * pkin(3) + t68;
t102 = t33 * pkin(4) - t34 * pkin(9) + t38;
t74 = t76 ^ 2;
t69 = t100 * qJD(4) - t97 * t119;
t52 = t142 * t120 + t93 * t80;
t46 = -t75 * t135 + t144;
t42 = t93 * t103 + t138 * t69;
t41 = -t138 * t103 + t93 * t69;
t30 = t99 * t48 - t96 * t53;
t25 = -t44 * t135 + t145;
t22 = -qJD(5) * t127 - t99 * t126 + t45 * t134 + t96 * t34;
t14 = pkin(4) * t139 - t18;
t11 = -t31 * qJD(5) + t99 * t111 - t96 * t42;
t10 = -t96 * t111 - t48 * t134 + t53 * t135 - t99 * t42;
t7 = -t96 * t15 + t99 * t24;
t4 = -pkin(4) * t126 - t5;
t2 = -t8 * qJD(5) + t99 * t102 - t96 * t106;
t1 = -t96 * t102 - t99 * t106 - t24 * t134 + t15 * t135;
t3 = [0, 0, 0, 0.2e1 * t115, 0.2e1 * (t101 ^ 2 - t98 ^ 2) * t90 * qJD(2), 0.2e1 * t95 * t123, -0.2e1 * t95 * t126, 0, -0.2e1 * t90 * pkin(1) * t137 - 0.2e1 * t68 * t95, -0.2e1 * pkin(1) * t124 + 0.2e1 * t67 * t95, 0.2e1 * t73 * t51, -0.2e1 * t73 * t50 - 0.2e1 * t51 * t72, (-t101 * t51 + t73 * t137) * t157, (t101 * t50 - t72 * t137) * t157, -0.2e1 * t115, 0.2e1 * t62 * t50 + 0.2e1 * t68 * t72 + 0.2e1 * (-t27 * t101 + t121 * t137) * t94, 0.2e1 * t62 * t51 + 0.2e1 * t68 * t73 + 0.2e1 * (-t26 * t101 - t141 * t137) * t94, 0.2e1 * t47 * t33 + 0.2e1 * t38 * t44 + 0.2e1 * (-t101 * t5 + t18 * t137) * t94, 0.2e1 * t47 * t34 + 0.2e1 * t38 * t45 + 0.2e1 * (t101 * t6 - t19 * t137) * t94, -0.2e1 * t18 * t34 - 0.2e1 * t19 * t33 - 0.2e1 * t6 * t44 - 0.2e1 * t5 * t45, 0.2e1 * t18 * t5 + 0.2e1 * t19 * t6 + 0.2e1 * t47 * t38, 0.2e1 * t37 * t21, -0.2e1 * t21 * t36 - 0.2e1 * t37 * t22, 0.2e1 * t21 * t44 + 0.2e1 * t37 * t33, -0.2e1 * t22 * t44 - 0.2e1 * t36 * t33, 0.2e1 * t44 * t33, 0.2e1 * t14 * t22 + 0.2e1 * t2 * t44 + 0.2e1 * t7 * t33 + 0.2e1 * t4 * t36, 0.2e1 * t1 * t44 + 0.2e1 * t14 * t21 - 0.2e1 * t8 * t33 + 0.2e1 * t4 * t37; 0, 0, 0, 0, 0, t123, -t126, 0, -t68, t67, t73 * t132 + t51 * t97, t51 * t100 - t97 * t50 + (-t100 * t72 - t73 * t97) * qJD(3), t105 * t94, t104 * t94, 0, -pkin(2) * t50 - t68 * t100 - t105 * t154 + t62 * t136, -pkin(2) * t51 - t104 * t154 + t62 * t132 + t68 * t97, t44 * t89 + t88 * t33 + t38 * t75 + t47 * t70 + (t101 * t41 - t52 * t137) * t94, t45 * t89 + t88 * t34 + t38 * t76 + t47 * t71 + (t101 * t42 - t53 * t137) * t94, -t18 * t71 - t19 * t70 - t53 * t33 + t52 * t34 + t41 * t45 - t42 * t44 - t5 * t76 - t6 * t75, -t18 * t41 + t19 * t42 + t38 * t88 + t47 * t89 - t5 * t52 + t6 * t53, t107 * t37 + t21 * t150, t114 * t71 + (-t153 - t22 * t99 + (t36 * t96 - t37 * t99) * qJD(5)) * t76, t107 * t44 + t76 * t145 + t21 * t75 + t37 * t70, -t108 * t44 - t76 * t147 - t22 * t75 - t36 * t70, t33 * t75 + t44 * t70, t108 * t14 + t11 * t44 + t4 * t151 + t2 * t75 + t52 * t22 + t30 * t33 + t41 * t36 + t7 * t70, t1 * t75 + t10 * t44 + t107 * t14 + t4 * t150 + t52 * t21 - t31 * t33 + t41 * t37 - t8 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t97 * t132, 0.2e1 * (t100 ^ 2 - t97 ^ 2) * qJD(3), 0, 0, 0, t97 * t130, t100 * t130, 0.2e1 * t88 * t70 + 0.2e1 * t75 * t89, 0.2e1 * t88 * t71 + 0.2e1 * t76 * t89, 0.2e1 * t41 * t76 - 0.2e1 * t42 * t75 + 0.2e1 * t52 * t71 - 0.2e1 * t53 * t70, 0.2e1 * t52 * t41 + 0.2e1 * t53 * t42 + 0.2e1 * t88 * t89, 0.2e1 * t92 * t76 * t71 - 0.2e1 * t74 * t125, t140 * t74 * t156 + t71 * t122, 0.2e1 * t107 * t75 + 0.2e1 * t76 * t144, -0.2e1 * t108 * t75 - 0.2e1 * t76 * t146, 0.2e1 * t75 * t70, 0.2e1 * t108 * t52 + 0.2e1 * t11 * t75 + 0.2e1 * t41 * t151 + 0.2e1 * t30 * t70, 0.2e1 * t10 * t75 + 0.2e1 * t107 * t52 + 0.2e1 * t41 * t150 - 0.2e1 * t31 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t50, t126, t27, t26, t138 * t116 + t5, -t93 * t116 - t6, (-t138 * t34 - t33 * t93) * pkin(3), (t138 * t5 + t6 * t93) * pkin(3), t37 * t134 + t153, qJD(5) * t114 + t21 * t99 - t96 * t22, t110, t25, 0, -t86 * t147 + t87 * t22 - t4 * t99 + (t14 * t96 - t99 * t152) * qJD(5), -t86 * t145 + t87 * t21 + t4 * t96 + (t14 * t99 + t96 * t152) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, -t136, 0, -pkin(8) * t132, pkin(8) * t136, -t41, -t42, (-t138 * t71 - t70 * t93) * pkin(3), (-t138 * t41 + t42 * t93) * pkin(3), -t76 * t118 + t96 * t143, qJD(5) * t122 - t140 * t71, t109, t46, 0, -t41 * t99 + t113 * t96 + (-t112 * t99 + t52 * t96) * qJD(5), t41 * t96 + t113 * t99 + (t112 * t96 + t52 * t99) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t125, -0.2e1 * t118, 0, 0, 0, t96 * t128, t99 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t34, 0, t38, 0, 0, 0, 0, 0, t25, -t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t71, 0, t89, 0, 0, 0, 0, 0, t46, -t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t22, t33, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, -t108, t70, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t135, 0, -t86 * t134, t86 * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
