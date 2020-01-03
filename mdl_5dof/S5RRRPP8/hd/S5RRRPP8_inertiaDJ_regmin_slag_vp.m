% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:24
% EndTime: 2019-12-31 21:09:29
% DurationCPUTime: 1.41s
% Computational Cost: add. (825->191), mult. (2038->334), div. (0->0), fcn. (1404->4), ass. (0->112)
t58 = cos(qJ(2));
t106 = t58 * qJD(2);
t57 = cos(qJ(3));
t47 = qJD(3) * t57;
t55 = sin(qJ(3));
t56 = sin(qJ(2));
t135 = t55 * t106 + t56 * t47;
t45 = t56 * qJD(2);
t109 = qJD(3) * t58;
t92 = t57 * t109;
t134 = t55 * t45 - t92;
t133 = -0.4e1 * t56;
t54 = -pkin(3) - qJ(5);
t132 = t54 * t56;
t46 = qJD(4) * t57;
t113 = t55 * qJ(4);
t72 = -t57 * pkin(3) - t113;
t131 = t72 * qJD(3) + t46;
t70 = t54 * t57 - t113;
t52 = t57 ^ 2;
t116 = t55 ^ 2 - t52;
t81 = t116 * qJD(3);
t128 = pkin(7) * t58;
t74 = pkin(2) * t56 - t128;
t27 = t74 * qJD(2);
t125 = t56 * pkin(7);
t75 = -t58 * pkin(2) - t125;
t29 = -pkin(1) + t75;
t119 = t55 * t27 + t29 * t47;
t94 = t55 * t109;
t64 = t57 * t45 + t94;
t7 = t64 * pkin(6) - t119;
t28 = -pkin(2) + t72;
t41 = pkin(6) * t106;
t102 = qJ(4) * qJD(3);
t84 = t56 * t102;
t78 = t135 * pkin(3) + t55 * t84 + t41;
t103 = qJ(4) * qJD(2);
t85 = t58 * t103;
t6 = (-qJD(4) * t56 - t85) * t57 + t78;
t130 = qJD(3) * (t28 * t56 + t128) - t6;
t129 = pkin(4) + pkin(7);
t108 = qJD(5) * t55;
t114 = qJ(4) * t57;
t69 = qJ(5) * t55 - t114;
t3 = t69 * t106 + (t108 + (qJ(5) * qJD(3) - qJD(4)) * t57) * t56 + t78;
t127 = t3 * t55;
t126 = t3 * t57;
t123 = t55 * t56;
t122 = t55 * t58;
t121 = t56 * t57;
t120 = t57 * t58;
t118 = pkin(3) * t123 + t56 * pkin(6);
t39 = pkin(6) * t120;
t117 = t55 * t29 + t39;
t51 = t56 ^ 2;
t115 = -t58 ^ 2 + t51;
t112 = qJD(2) * t55;
t111 = qJD(2) * t57;
t110 = qJD(3) * t55;
t107 = t55 * qJD(4);
t105 = t58 * qJD(4);
t104 = t58 * qJD(5);
t101 = qJ(4) * qJD(4);
t100 = pkin(4) * t120;
t38 = pkin(6) * t122;
t99 = -0.2e1 * pkin(1) * qJD(2);
t98 = -0.2e1 * pkin(2) * qJD(3);
t97 = pkin(3) * t45;
t96 = pkin(7) * t110;
t95 = t56 * t110;
t21 = -pkin(2) + t70;
t91 = t21 * t47;
t88 = t55 * t47;
t87 = t56 * t106;
t86 = t57 * t106;
t83 = t56 * t103;
t82 = t57 * t29 - t38;
t80 = t115 * qJD(2);
t79 = 0.2e1 * t87;
t77 = -t134 * pkin(6) + t29 * t110 - t57 * t27;
t76 = t55 * t86;
t13 = t58 * qJ(4) - t117;
t49 = t58 * pkin(3);
t14 = t49 - t82;
t71 = t13 * t55 + t14 * t57;
t68 = -t57 * qJD(5) - t107;
t42 = pkin(3) * t110;
t11 = t69 * qJD(3) + t42 + t68;
t67 = -t11 * t57 + t21 * t110;
t12 = t69 * t56 + t118;
t66 = -qJD(3) * t12 - t21 * t106;
t63 = -pkin(4) * t95 + t77;
t4 = t7 - t83 + t105;
t5 = t77 - t97;
t62 = t71 * qJD(3) - t4 * t57 + t5 * t55;
t15 = -t56 * t114 + t118;
t16 = -t57 * t102 - t107 + t42;
t61 = -qJD(3) * t15 - t16 * t56 + (-t28 * t58 + t125) * qJD(2);
t60 = -0.2e1 * t105 - t7 + 0.2e1 * t83;
t59 = 0.2e1 * qJD(4);
t43 = pkin(7) * t47;
t31 = t129 * t57;
t30 = t129 * t55;
t26 = pkin(4) * t47 + t43;
t25 = t129 * t110;
t17 = t86 - t95;
t10 = -pkin(4) * t123 - t13;
t9 = t58 * qJ(5) + t38 + t49 + (pkin(4) * t56 - t29) * t57;
t2 = -t105 + (-pkin(4) * t121 - t38) * qJD(3) + (-pkin(4) * t122 + (-pkin(6) * t57 + qJ(4)) * t56) * qJD(2) + t119;
t1 = t104 + (t100 + t132) * qJD(2) + t63;
t8 = [0, 0, 0, t79, -0.2e1 * t80, 0, 0, 0, t56 * t99, t58 * t99, -0.2e1 * t51 * t88 + 0.2e1 * t52 * t87, t76 * t133 + 0.2e1 * t51 * t81, 0.2e1 * t115 * t111 + 0.2e1 * t56 * t94, -0.2e1 * t55 * t80 + 0.2e1 * t56 * t92, -0.2e1 * t87, 0.2e1 * t77 * t58 + 0.2e1 * t82 * t45 + 0.2e1 * (t51 * t47 + t55 * t79) * pkin(6), -0.2e1 * t7 * t58 - 0.2e1 * t117 * t45 + 0.2e1 * (-t51 * t110 + t57 * t79) * pkin(6), 0.2e1 * t71 * t106 + 0.2e1 * (t4 * t55 + t5 * t57 + (t13 * t57 - t14 * t55) * qJD(3)) * t56, 0.2e1 * (-t15 * t112 - t5) * t58 + 0.2e1 * (qJD(2) * t14 - t15 * t47 - t6 * t55) * t56, 0.2e1 * (-t15 * t111 + t4) * t58 + 0.2e1 * (-qJD(2) * t13 + t15 * t110 - t6 * t57) * t56, 0.2e1 * t13 * t4 + 0.2e1 * t14 * t5 + 0.2e1 * t15 * t6, 0.2e1 * (-t10 * t55 + t57 * t9) * t106 + 0.2e1 * (t1 * t57 - t2 * t55 + (-t10 * t57 - t55 * t9) * qJD(3)) * t56, 0.2e1 * (-t12 * t111 - t2) * t58 + 0.2e1 * (qJD(2) * t10 + t12 * t110 - t126) * t56, 0.2e1 * (t12 * t112 + t1) * t58 + 0.2e1 * (-qJD(2) * t9 + t12 * t47 + t127) * t56, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t12 * t3; 0, 0, 0, 0, 0, t106, -t45, 0, -t41, pkin(6) * t45, -t56 * t81 + t76, -t116 * t106 + t88 * t133, t134, t64, 0, (pkin(7) * t120 + (-pkin(2) * t57 + pkin(6) * t55) * t56) * qJD(3) + (t75 * t55 - t39) * qJD(2), (pkin(6) * t121 + t55 * t74) * qJD(3) + (t57 * t75 + t38) * qJD(2), t62, -t130 * t57 + t61 * t55, t130 * t55 + t61 * t57, pkin(7) * t62 + t15 * t16 + t6 * t28, (t30 * t106 + t26 * t56 + t2 + (-t31 * t56 + t9) * qJD(3)) * t57 + (-t31 * t106 + t25 * t56 + t1 + (-t30 * t56 - t10) * qJD(3)) * t55, t25 * t58 - t127 + t66 * t57 + (qJD(2) * t31 + t67) * t56, t26 * t58 - t126 + (-qJD(2) * t30 + t91) * t56 + (t11 * t56 - t66) * t55, t1 * t30 - t10 * t25 + t12 * t11 + t2 * t31 + t3 * t21 + t9 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t88, -0.2e1 * t81, 0, 0, 0, t55 * t98, t57 * t98, 0, -0.2e1 * t28 * t110 + 0.2e1 * t16 * t57, -0.2e1 * t16 * t55 - 0.2e1 * t28 * t47, 0.2e1 * t28 * t16, -0.2e1 * t25 * t57 + 0.2e1 * t26 * t55 + 0.2e1 * (t30 * t57 - t31 * t55) * qJD(3), -0.2e1 * t11 * t55 - 0.2e1 * t91, 0.2e1 * t67, 0.2e1 * t21 * t11 - 0.2e1 * t31 * t25 + 0.2e1 * t30 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t135, t45, -t77, t7, (-pkin(3) * t106 - t84) * t57 + (-t85 + (pkin(3) * qJD(3) - qJD(4)) * t56) * t55, t77 - 0.2e1 * t97, t60, -t5 * pkin(3) - t4 * qJ(4) - t13 * qJD(4), t70 * t106 + ((-t54 * t55 - t114) * qJD(3) + t68) * t56, -pkin(4) * t135 + t60, -0.2e1 * t104 + (-t100 - 0.2e1 * t132) * qJD(2) - t63, t2 * qJ(4) + t10 * qJD(4) - t9 * qJD(5) + t1 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t110, 0, -t43, t96, t131, t43, -t96, t131 * pkin(7), qJD(3) * t70 - t108 + t46, -t25, -t26, -t25 * qJ(4) + t31 * qJD(4) - t30 * qJD(5) + t26 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0.2e1 * t101, 0, t59, 0.2e1 * qJD(5), -0.2e1 * t54 * qJD(5) + 0.2e1 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t45, 0, t5, t17, 0, -t45, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, t43, t47, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, t45, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, 0, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
