% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:24
% EndTime: 2019-12-31 21:09:30
% DurationCPUTime: 1.82s
% Computational Cost: add. (959->218), mult. (2525->365), div. (0->0), fcn. (1738->4), ass. (0->129)
t84 = cos(qJ(2));
t132 = t84 * qJD(2);
t83 = cos(qJ(3));
t73 = qJD(3) * t83;
t81 = sin(qJ(3));
t82 = sin(qJ(2));
t31 = t81 * t132 + t82 * t73;
t135 = qJD(3) * t84;
t118 = t83 * t135;
t70 = t82 * qJD(2);
t32 = -t81 * t70 + t118;
t80 = -pkin(3) - qJ(5);
t157 = t80 * t82;
t71 = qJD(4) * t83;
t139 = t81 * qJ(4);
t99 = -t83 * pkin(3) - t139;
t156 = t99 * qJD(3) + t71;
t76 = t81 ^ 2;
t78 = t83 ^ 2;
t142 = t76 - t78;
t96 = t80 * t83 - t139;
t45 = t142 * qJD(3);
t153 = pkin(7) * t84;
t101 = pkin(2) * t82 - t153;
t42 = t101 * qJD(2);
t150 = t82 * pkin(7);
t102 = -t84 * pkin(2) - t150;
t46 = -pkin(1) + t102;
t144 = t81 * t42 + t46 * t73;
t120 = t81 * t135;
t30 = t83 * t70 + t120;
t7 = t30 * pkin(6) - t144;
t44 = -pkin(2) + t99;
t128 = qJ(4) * qJD(3);
t110 = t82 * t128;
t66 = pkin(6) * t132;
t106 = t31 * pkin(3) + t81 * t110 + t66;
t129 = qJ(4) * qJD(2);
t111 = t84 * t129;
t6 = (-qJD(4) * t82 - t111) * t83 + t106;
t155 = (t44 * t82 + t153) * qJD(3) - t6;
t77 = t82 ^ 2;
t141 = -t84 ^ 2 + t77;
t108 = t141 * qJD(2);
t14 = -0.2e1 * t81 * t108 + 0.2e1 * t82 * t118;
t154 = pkin(4) + pkin(7);
t134 = qJD(5) * t81;
t140 = qJ(4) * t83;
t95 = qJ(5) * t81 - t140;
t3 = t95 * t132 + (t134 + (qJ(5) * qJD(3) - qJD(4)) * t83) * t82 + t106;
t152 = t3 * t81;
t151 = t3 * t83;
t148 = t81 * t82;
t147 = t81 * t84;
t146 = t82 * t83;
t145 = t83 * t84;
t143 = pkin(3) * t148 + t82 * pkin(6);
t64 = pkin(6) * t145;
t26 = t81 * t46 + t64;
t138 = qJD(2) * t81;
t137 = qJD(2) * t83;
t72 = qJD(3) * t81;
t136 = qJD(3) * t82;
t133 = t81 * qJD(4);
t131 = t84 * qJD(4);
t130 = t84 * qJD(5);
t127 = qJ(4) * qJD(4);
t126 = pkin(4) * t145;
t63 = pkin(6) * t147;
t125 = -0.2e1 * pkin(1) * qJD(2);
t124 = -0.2e1 * pkin(2) * qJD(3);
t123 = pkin(3) * t70;
t122 = pkin(7) * t72;
t121 = t81 * t136;
t35 = -pkin(2) + t96;
t117 = t35 * t73;
t114 = t81 * t73;
t113 = t82 * t132;
t112 = t83 * t132;
t109 = t82 * t129;
t25 = t83 * t46 - t63;
t107 = 0.2e1 * t113;
t105 = t32 * pkin(6) - t83 * t42 + t46 * t72;
t104 = t81 * t112;
t103 = t77 * t114;
t20 = t84 * qJ(4) - t26;
t75 = t84 * pkin(3);
t22 = -t25 + t75;
t98 = t20 * t81 + t22 * t83;
t97 = -t25 * t83 - t26 * t81;
t94 = -t83 * qJD(5) - t133;
t67 = pkin(3) * t72;
t15 = t95 * qJD(3) + t67 + t94;
t93 = -t15 * t83 + t35 * t72;
t16 = t95 * t82 + t143;
t92 = -qJD(3) * t16 - t35 * t132;
t90 = -pkin(4) * t121 + t105;
t4 = -t109 + t7 + t131;
t5 = t105 - t123;
t89 = t98 * qJD(3) - t4 * t83 + t5 * t81;
t88 = t97 * qJD(3) + t105 * t81 - t7 * t83;
t27 = -t82 * t140 + t143;
t28 = -t83 * t128 - t133 + t67;
t87 = -qJD(3) * t27 - t28 * t82 + (-t44 * t84 + t150) * qJD(2);
t86 = 0.2e1 * t109 - 0.2e1 * t131 - t7;
t85 = 0.2e1 * qJD(4);
t68 = pkin(7) * t73;
t57 = -0.2e1 * t113;
t56 = -0.2e1 * t114;
t55 = 0.2e1 * t114;
t49 = t154 * t83;
t48 = t154 * t81;
t41 = pkin(4) * t73 + t68;
t40 = t154 * t72;
t39 = -0.2e1 * t45;
t29 = -t112 + t121;
t24 = 0.2e1 * t78 * t113 - 0.2e1 * t103;
t23 = 0.2e1 * t76 * t113 + 0.2e1 * t103;
t21 = t142 * t136 - t104;
t18 = t82 * t120 + t141 * t137;
t17 = 0.4e1 * t82 * t114 + t142 * t132;
t13 = 0.2e1 * t18;
t12 = -0.2e1 * t82 * t104 + t77 * t45;
t11 = 0.2e1 * t12;
t10 = -pkin(4) * t148 - t20;
t9 = t84 * qJ(5) + t63 + t75 + (pkin(4) * t82 - t46) * t83;
t2 = -t131 + (-pkin(4) * t146 - t63) * qJD(3) + (-pkin(4) * t147 + (-pkin(6) * t83 + qJ(4)) * t82) * qJD(2) + t144;
t1 = t130 + (t126 + t157) * qJD(2) + t90;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, -0.2e1 * t108, 0, t57, 0, 0, t82 * t125, t84 * t125, 0, 0, t24, t11, t13, t23, t14, t57, 0.2e1 * t25 * t70 + 0.2e1 * t105 * t84 + 0.2e1 * (t81 * t107 + t77 * t73) * pkin(6), -0.2e1 * t26 * t70 - 0.2e1 * t7 * t84 + 0.2e1 * (t83 * t107 - t77 * t72) * pkin(6), 0.2e1 * t97 * t132 + 0.2e1 * (t7 * t81 + t105 * t83 + (t25 * t81 - t26 * t83) * qJD(3)) * t82, 0.2e1 * pkin(6) ^ 2 * t113 - 0.2e1 * t105 * t25 - 0.2e1 * t26 * t7, t57, -0.2e1 * t18, -t14, t24, t11, t23, 0.2e1 * t98 * t132 + 0.2e1 * (t4 * t81 + t5 * t83 + (t20 * t83 - t22 * t81) * qJD(3)) * t82, 0.2e1 * (-t27 * t138 - t5) * t84 + 0.2e1 * (qJD(2) * t22 - t27 * t73 - t6 * t81) * t82, 0.2e1 * (-t27 * t137 + t4) * t84 + 0.2e1 * (-qJD(2) * t20 + t27 * t72 - t6 * t83) * t82, 0.2e1 * t20 * t4 + 0.2e1 * t22 * t5 + 0.2e1 * t27 * t6, t57, -t14, t13, t23, -0.2e1 * t12, t24, 0.2e1 * (-t10 * t81 + t83 * t9) * t132 + 0.2e1 * (t1 * t83 - t2 * t81 + (-t10 * t83 - t81 * t9) * qJD(3)) * t82, 0.2e1 * (-t16 * t137 - t2) * t84 + 0.2e1 * (qJD(2) * t10 + t16 * t72 - t151) * t82, 0.2e1 * (t16 * t138 + t1) * t84 + 0.2e1 * (-qJD(2) * t9 + t16 * t73 + t152) * t82, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t16 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, 0, -t70, 0, -t66, pkin(6) * t70, 0, 0, -t21, -t17, -t32, t21, t30, 0, (pkin(7) * t145 + (-pkin(2) * t83 + pkin(6) * t81) * t82) * qJD(3) + (t102 * t81 - t64) * qJD(2), (pkin(6) * t146 + t101 * t81) * qJD(3) + (t102 * t83 + t63) * qJD(2), t88, -pkin(2) * t66 + t88 * pkin(7), 0, t32, -t30, -t21, -t17, t21, t89, -t155 * t83 + t87 * t81, t155 * t81 + t87 * t83, t89 * pkin(7) + t27 * t28 + t6 * t44, 0, -t30, -t32, t21, t17, -t21, (t48 * t132 + t41 * t82 + t2 + (-t49 * t82 + t9) * qJD(3)) * t83 + (-t49 * t132 + t40 * t82 + t1 + (-t48 * t82 - t10) * qJD(3)) * t81, -t152 + t40 * t84 + t92 * t83 + (qJD(2) * t49 + t93) * t82, -t151 + t41 * t84 + (-qJD(2) * t48 + t117) * t82 + (t15 * t82 - t92) * t81, t1 * t48 - t10 * t40 + t16 * t15 + t2 * t49 + t3 * t35 + t9 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t39, 0, t56, 0, 0, t81 * t124, t83 * t124, 0, 0, 0, 0, 0, t55, t39, t56, 0, 0.2e1 * t28 * t83 - 0.2e1 * t44 * t72, -0.2e1 * t28 * t81 - 0.2e1 * t44 * t73, 0.2e1 * t44 * t28, 0, 0, 0, t56, 0.2e1 * t45, t55, -0.2e1 * t40 * t83 + 0.2e1 * t41 * t81 + 0.2e1 * (t48 * t83 - t49 * t81) * qJD(3), -0.2e1 * t15 * t81 - 0.2e1 * t117, 0.2e1 * t93, 0.2e1 * t35 * t15 - 0.2e1 * t49 * t40 + 0.2e1 * t48 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, -t31, t70, -t105, t7, 0, 0, t70, t29, t31, 0, 0, 0, (-pkin(3) * t132 - t110) * t83 + (-t111 + (pkin(3) * qJD(3) - qJD(4)) * t82) * t81, t105 - 0.2e1 * t123, t86, -t5 * pkin(3) - t4 * qJ(4) - t20 * qJD(4), t70, t31, -t29, 0, 0, 0, t96 * t132 + ((-t80 * t81 - t140) * qJD(3) + t94) * t82, -pkin(4) * t31 + t86, -0.2e1 * t130 + (-t126 - 0.2e1 * t157) * qJD(2) - t90, t2 * qJ(4) + t10 * qJD(4) - t9 * qJD(5) + t1 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, -t72, 0, -t68, t122, 0, 0, 0, -t73, t72, 0, 0, 0, t156, t68, -t122, t156 * pkin(7), 0, t72, t73, 0, 0, 0, qJD(3) * t96 - t134 + t71, -t40, -t41, -t40 * qJ(4) + t49 * qJD(4) - t48 * qJD(5) + t41 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0.2e1 * t127, 0, 0, 0, 0, 0, 0, 0, t85, 0.2e1 * qJD(5), -0.2e1 * t80 * qJD(5) + 0.2e1 * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t70, 0, t5, 0, 0, 0, 0, 0, 0, -t29, 0, -t70, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, t68, 0, 0, 0, 0, 0, 0, t73, 0, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t70, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
