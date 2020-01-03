% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR15_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:09
% EndTime: 2019-12-31 20:43:15
% DurationCPUTime: 1.83s
% Computational Cost: add. (2026->196), mult. (4308->375), div. (0->0), fcn. (3481->6), ass. (0->125)
t127 = qJD(4) + qJD(5);
t149 = pkin(3) + pkin(6);
t66 = sin(qJ(5));
t67 = sin(qJ(4));
t69 = cos(qJ(5));
t70 = cos(qJ(4));
t40 = t66 * t70 + t67 * t69;
t20 = t127 * t40;
t135 = qJD(5) * t66;
t138 = qJD(4) * t67;
t145 = t69 * t70;
t21 = t127 * t145 - t135 * t67 - t138 * t66;
t41 = -t66 * t67 + t145;
t159 = (-t20 * t69 + t21 * t66 + (t40 * t69 - t41 * t66) * qJD(5)) * pkin(4);
t68 = sin(qJ(2));
t86 = t67 * qJ(3) + t149 * t70;
t83 = pkin(4) + t86;
t150 = pkin(2) + pkin(7);
t128 = pkin(8) + t150;
t71 = cos(qJ(2));
t89 = t128 * t71 + pkin(1);
t156 = t89 * t67 + t83 * t68;
t132 = t68 * qJD(2);
t108 = t70 * t132;
t136 = qJD(4) * t71;
t114 = t67 * t136;
t81 = t108 + t114;
t131 = t68 * qJD(3);
t116 = t150 * t68;
t141 = qJ(3) * t71;
t157 = t116 - t141;
t137 = qJD(4) * t70;
t58 = t71 * qJD(2);
t35 = -t68 * t137 - t67 * t58;
t117 = t150 * t71;
t140 = t68 * qJ(3);
t85 = -t117 - t140;
t82 = -pkin(1) + t85;
t9 = t82 * t138 - t70 * (qJD(2) * t157 - t131) + t35 * t149;
t158 = -pkin(8) * t81 - qJD(5) * t156 + t9;
t103 = t128 * t70;
t155 = t127 * t103;
t111 = t67 * t132;
t112 = t70 * t136;
t154 = t111 - t112;
t65 = t71 ^ 2;
t105 = qJD(2) * (t68 ^ 2 - t65);
t62 = t67 ^ 2;
t64 = t70 ^ 2;
t143 = t62 - t64;
t104 = t143 * qJD(4);
t130 = t71 * qJD(3);
t49 = t149 * t71;
t133 = t49 * qJD(4);
t153 = qJD(2) * t85 + t130 + t133;
t151 = 0.2e1 * qJD(3);
t148 = t40 * t21;
t147 = t41 * t20;
t146 = t67 * t71;
t144 = t70 * t71;
t119 = t67 * t149;
t19 = t119 * t68 + t70 * t82;
t139 = qJD(2) * t49;
t134 = qJD(5) * t69;
t129 = qJ(3) * qJD(4);
t126 = -0.2e1 * pkin(1) * qJD(2);
t125 = t69 * t144;
t124 = pkin(4) * t135;
t123 = pkin(4) * t134;
t122 = pkin(6) * t132;
t121 = pkin(6) * t58;
t120 = t67 * t150;
t118 = t70 * t150;
t110 = t67 * t131;
t109 = t68 * t58;
t107 = t67 * t137;
t106 = qJD(4) * t150;
t44 = t128 * t67;
t102 = t67 * t108;
t101 = t65 * t107;
t100 = pkin(1) + t117;
t99 = -pkin(2) * t71 - t140;
t13 = -qJD(5) * t125 + (t127 * t146 + t108) * t66 + t154 * t69;
t28 = t40 * t71;
t98 = t13 * t41 + t20 * t28;
t14 = t108 * t69 - t111 * t66 + t20 * t71;
t27 = t146 * t66 - t125;
t97 = -t14 * t40 - t21 * t27;
t18 = t100 * t67 + t68 * t86;
t96 = -t18 * t67 + t19 * t70;
t95 = t147 - t148;
t92 = t128 * t138;
t87 = -t21 * t68 - t40 * t58;
t17 = -pkin(8) * t144 + t19;
t80 = (t70 * qJ(3) - t119) * t68;
t72 = t110 + (t70 * t89 + t80) * qJD(4) + (-t44 * t68 + t71 * t83) * qJD(2);
t1 = t17 * t135 + t158 * t69 - t66 * t72;
t2 = -t17 * t134 + t158 * t66 + t69 * t72;
t7 = t156 * t69 - t66 * t17;
t8 = t156 * t66 + t69 * t17;
t78 = -t1 * t40 + t2 * t41 - t20 * t7 + t21 * t8;
t43 = t149 * t132;
t77 = qJD(4) * t157 - t43;
t11 = -t44 * t135 + t155 * t69 - t66 * t92;
t12 = t44 * t134 + t155 * t66 + t69 * t92;
t22 = -t103 * t69 + t66 * t44;
t23 = -t103 * t66 - t69 * t44;
t76 = t11 * t40 - t12 * t41 + t20 * t22 - t21 * t23;
t75 = qJD(2) * t99 + t130;
t61 = qJ(3) * t151;
t56 = pkin(4) * t67 + qJ(3);
t53 = pkin(4) * t137 + qJD(3);
t52 = -0.2e1 * t109;
t51 = 0.2e1 * t109;
t45 = -pkin(1) + t99;
t39 = -0.2e1 * t105;
t34 = -t138 * t68 + t58 * t70;
t33 = pkin(4) * t144 + t49;
t30 = -t131 + (pkin(2) * t68 - t141) * qJD(2);
t26 = -t104 * t71 - t102;
t24 = -pkin(4) * t114 + (-pkin(4) * t70 - t149) * t132;
t15 = -t20 * t68 + t41 * t58;
t10 = t110 + (t100 * t70 + t80) * qJD(4) + (-t116 * t67 + t71 * t86) * qJD(2);
t3 = qJD(4) * t96 + t10 * t70 - t67 * t9;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t39, 0, t52, 0, 0, t68 * t126, t71 * t126, 0, 0, 0, 0, 0, t51, t39, t52, 0, -0.2e1 * t132 * t45 + 0.2e1 * t30 * t71, -0.2e1 * t30 * t68 - 0.2e1 * t45 * t58, 0.2e1 * t45 * t30, -0.2e1 * t109 * t62 + 0.2e1 * t101, -0.4e1 * t102 * t71 - 0.2e1 * t104 * t65, 0.2e1 * t105 * t67 - 0.2e1 * t112 * t68, -0.2e1 * t109 * t64 - 0.2e1 * t101, 0.2e1 * t105 * t70 + 0.2e1 * t114 * t68, t51, 0.2e1 * (-t139 * t70 + t10) * t68 + 0.2e1 * (qJD(2) * t18 - t133 * t67 - t43 * t70) * t71, 0.2e1 * (t139 * t67 + t9) * t68 + 0.2e1 * (-qJD(2) * t19 - t133 * t70 + t43 * t67) * t71, 0.2e1 * t96 * t132 + 0.2e1 * (t10 * t67 + t70 * t9 + (t18 * t70 + t19 * t67) * qJD(4)) * t71, 0.2e1 * t10 * t18 - 0.2e1 * t19 * t9 - 0.2e1 * t43 * t49, -0.2e1 * t28 * t13, 0.2e1 * t13 * t27 - 0.2e1 * t14 * t28, 0.2e1 * t13 * t68 - 0.2e1 * t28 * t58, 0.2e1 * t27 * t14, 0.2e1 * t14 * t68 + 0.2e1 * t27 * t58, t51, -0.2e1 * t14 * t33 + 0.2e1 * t2 * t68 - 0.2e1 * t24 * t27 + 0.2e1 * t58 * t7, 0.2e1 * t1 * t68 + 0.2e1 * t13 * t33 - 0.2e1 * t24 * t28 - 0.2e1 * t58 * t8, -0.2e1 * t1 * t27 - 0.2e1 * t13 * t7 + 0.2e1 * t14 * t8 + 0.2e1 * t2 * t28, -0.2e1 * t1 * t8 + 0.2e1 * t2 * t7 + 0.2e1 * t24 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, -t132, 0, -t121, t122, 0, 0, 0, -t58, t132, 0, 0, 0, t75, t121, -t122, t75 * pkin(6), -t26, 0.4e1 * t107 * t71 - t132 * t143, t34, t26, t35, 0, t153 * t70 + t77 * t67, -t153 * t67 + t77 * t70, -t3, t9 * t120 - t10 * t118 - t43 * qJ(3) + t49 * qJD(3) + (-t118 * t19 + t120 * t18) * qJD(4), t98, -t13 * t40 + t14 * t41 - t20 * t27 + t21 * t28, t15, t97, t87, 0, t12 * t68 - t14 * t56 + t21 * t33 + t22 * t58 + t24 * t40 - t27 * t53, t11 * t68 + t13 * t56 - t20 * t33 - t23 * t58 + t24 * t41 - t28 * t53, -t11 * t27 + t12 * t28 - t13 * t22 + t14 * t23 - t78, -t1 * t23 - t11 * t8 + t12 * t7 + t2 * t22 + t24 * t56 + t33 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t61, -0.2e1 * t107, 0.2e1 * t104, 0, 0.2e1 * t107, 0, 0, 0.2e1 * qJD(3) * t67 + 0.2e1 * t129 * t70, 0.2e1 * qJD(3) * t70 - 0.2e1 * t129 * t67, 0, t61, -0.2e1 * t147, 0.2e1 * t20 * t40 - 0.2e1 * t21 * t41, 0, 0.2e1 * t148, 0, 0, 0.2e1 * t21 * t56 + 0.2e1 * t40 * t53, -0.2e1 * t20 * t56 + 0.2e1 * t41 * t53, 0.2e1 * t76, -0.2e1 * t11 * t23 + 0.2e1 * t12 * t22 + 0.2e1 * t53 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, t121, 0, 0, 0, 0, 0, 0, t34, t35, 0, t3, 0, 0, 0, 0, 0, 0, t15, t87, -t97 - t98, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t95, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, 0, t81, t58, t10, t9, 0, 0, 0, 0, t13, 0, t14, t58, pkin(4) * t58 * t69 - t124 * t68 + t2, (-t134 * t68 - t58 * t66) * pkin(4) + t1, (-t13 * t69 + t14 * t66 + (t27 * t69 - t28 * t66) * qJD(5)) * pkin(4), (-t1 * t66 + t2 * t69 + (-t66 * t7 + t69 * t8) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, 0, -t137, 0, t67 * t106, t70 * t106, 0, 0, 0, 0, -t20, 0, -t21, 0, t12, t11, -t159, (-t11 * t66 + t12 * t69 + (-t22 * t66 + t23 * t69) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, -t137, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t124, -0.2e1 * t123, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t14, t58, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, -t21, 0, t12, t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t123, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
