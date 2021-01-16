% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:51:50
% EndTime: 2021-01-16 01:51:57
% DurationCPUTime: 1.37s
% Computational Cost: add. (865->187), mult. (2207->322), div. (0->0), fcn. (1867->8), ass. (0->117)
t55 = sin(qJ(4));
t101 = t55 * qJD(4);
t58 = cos(qJ(4));
t130 = qJ(6) * t101 - t58 * qJD(6);
t129 = 2 * qJD(3);
t128 = 2 * qJD(5);
t57 = cos(qJ(5));
t127 = t57 * pkin(5);
t126 = t58 * pkin(9);
t54 = sin(qJ(5));
t103 = qJD(5) * t58;
t91 = t57 * t103;
t30 = -t54 * t101 + t91;
t60 = -pkin(2) - pkin(8);
t86 = t60 * t101;
t17 = t30 * pkin(5) + t86;
t125 = t17 * t54;
t124 = t17 * t57;
t115 = -qJ(6) - pkin(9);
t37 = t115 * t54;
t123 = t37 * t58;
t38 = t115 * t57;
t122 = t38 * t58;
t45 = -pkin(4) - t127;
t121 = t45 * t57;
t52 = sin(pkin(6));
t56 = sin(qJ(2));
t120 = t52 * t56;
t59 = cos(qJ(2));
t119 = t52 * t59;
t118 = t54 * t60;
t117 = t55 * t60;
t116 = t58 * t60;
t73 = t55 * pkin(4) - t126;
t67 = qJ(3) + t73;
t32 = t54 * t67;
t42 = t57 * t117;
t114 = t42 + t32;
t48 = t54 ^ 2;
t50 = t57 ^ 2;
t113 = t48 - t50;
t112 = t48 + t50;
t49 = t55 ^ 2;
t51 = t58 ^ 2;
t111 = t49 - t51;
t110 = t49 + t51;
t109 = qJ(6) * t58;
t108 = qJD(2) * t59;
t53 = cos(pkin(6));
t23 = t58 * t119 + t53 * t55;
t107 = qJD(4) * t23;
t106 = qJD(4) * t54;
t105 = qJD(4) * t57;
t104 = qJD(5) * t54;
t47 = qJD(5) * t57;
t102 = qJD(5) * t60;
t100 = t58 * qJD(4);
t98 = qJ(3) * qJD(4);
t97 = -0.2e1 * pkin(4) * qJD(5);
t33 = t57 * t67;
t74 = pkin(4) * t58 + pkin(9) * t55;
t65 = t74 * qJD(4) + qJD(3);
t83 = t60 * t100;
t96 = -qJD(5) * t33 - t54 * t65 - t57 * t83;
t95 = pkin(5) * t104;
t94 = t57 * t109;
t93 = t54 * t103;
t92 = t54 * t102;
t90 = t23 * t104;
t34 = (pkin(5) * t54 - t60) * t58;
t89 = t34 * t104;
t43 = qJD(2) * t120;
t88 = t52 * t108;
t87 = t54 * t47;
t85 = t57 * t101;
t84 = t55 * t100;
t81 = -t45 + t127;
t80 = pkin(5) - t118;
t79 = t113 * qJD(5);
t78 = t111 * qJD(4);
t77 = 0.2e1 * t84;
t76 = t54 * t83;
t75 = t54 * t85;
t72 = pkin(5) * t48 + t121;
t11 = t80 * t55 + t33 - t94;
t12 = -t54 * t109 + t114;
t71 = t11 * t57 + t12 * t54;
t70 = t11 * t54 - t12 * t57;
t24 = -t55 * t119 + t53 * t58;
t15 = t57 * t120 - t24 * t54;
t16 = t54 * t120 + t24 * t57;
t69 = t15 * t57 + t16 * t54;
t68 = t15 * t54 - t16 * t57;
t14 = t24 * qJD(4) - t58 * t43;
t7 = t14 * t54 + t23 * t47;
t8 = -t14 * t57 + t90;
t66 = t85 + t93;
t29 = t54 * t100 + t55 * t47;
t13 = -t55 * t43 + t107;
t5 = -t16 * qJD(5) + t13 * t54 + t57 * t88;
t6 = -t13 * t57 - t24 * t104 + (t54 * t108 + t56 * t47) * t52;
t64 = -t69 * qJD(5) - t5 * t54 + t6 * t57;
t19 = -t57 * qJD(6) - t115 * t104;
t20 = -t54 * qJD(6) + t115 * t47;
t63 = -t19 * t57 - t20 * t54 + (-t37 * t57 + t38 * t54) * qJD(5);
t62 = -t114 * qJD(5) + t57 * t65;
t61 = qJ(6) * t93 + t130 * t57 + t62;
t28 = t110 * t47;
t26 = -t57 * t100 + t55 * t104;
t25 = t110 * t104;
t10 = t62 - t76;
t9 = t55 * t92 + t96;
t4 = t94 * qJD(5) + t96 + (t117 * qJD(5) - t130) * t54;
t3 = t80 * t100 + t61;
t2 = (-t23 * t105 - t6) * t55 + (-qJD(4) * t16 - t8) * t58;
t1 = (-t23 * t106 + t5) * t55 + (qJD(4) * t15 + t7) * t58;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t23 * t14 + 0.2e1 * t15 * t5 + 0.2e1 * t16 * t6; 0, 0, -t43, -t88, t43, t88, (qJD(3) * t56 + (-pkin(2) * t56 + qJ(3) * t59) * qJD(2)) * t52, 0, 0, 0, 0, 0, (t56 * t100 + t55 * t108) * t52, (-t56 * t101 + t58 * t108) * t52, 0, 0, 0, 0, 0, t1, t2, t1, t2, t69 * t101 + (t68 * qJD(5) - t5 * t57 - t54 * t6) * t58, t5 * t11 + t6 * t12 + t14 * t34 + t15 * t3 - t16 * t4 + t23 * t17; 0, 0, 0, 0, 0, t129, qJ(3) * t129, -0.2e1 * t84, 0.2e1 * t78, 0, 0, 0, 0.2e1 * qJD(3) * t55 + 0.2e1 * t58 * t98, 0.2e1 * qJD(3) * t58 - 0.2e1 * t55 * t98, -0.2e1 * t50 * t84 - 0.2e1 * t51 * t87, t113 * t51 * t128 + 0.4e1 * t58 * t75, -0.2e1 * t111 * t105 - 0.2e1 * t55 * t93, 0.2e1 * t54 * t78 - 0.2e1 * t55 * t91, t77, -0.2e1 * t51 * t57 * t102 + 0.2e1 * t33 * t100 + 0.2e1 * (t10 + t76) * t55, 0.2e1 * t51 * t92 + 0.2e1 * t9 * t55 + 0.2e1 * (-t32 + t42) * t100, 0.2e1 * (-t34 * t106 + t3) * t55 + 0.2e1 * (qJD(4) * t11 + t34 * t47 + t125) * t58, 0.2e1 * (-t34 * t105 + t4) * t55 + 0.2e1 * (-qJD(4) * t12 + t124 - t89) * t58, 0.2e1 * t71 * t101 + 0.2e1 * (t70 * qJD(5) - t3 * t57 + t4 * t54) * t58, 0.2e1 * t11 * t3 - 0.2e1 * t12 * t4 + 0.2e1 * t34 * t17; 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t68 * qJD(4) - t14) * t58 + (t64 + t107) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t25, -t28, t25, 0, (-t70 * qJD(4) - t17) * t58 + (qJD(4) * t34 - t71 * qJD(5) - t3 * t54 - t4 * t57) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t112) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t13, 0, 0, 0, 0, 0, t8, t7, t8, t7, t64, pkin(5) * t90 + t14 * t45 + t15 * t20 - t16 * t19 + t5 * t37 - t6 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t100, 0, -t86, -t83, -t58 * t79 - t75, t113 * t101 - 0.4e1 * t58 * t87, t29, -t26, 0, (-t54 * t116 - t74 * t57) * qJD(5) + (t73 * t54 - t42) * qJD(4), (-t57 * t116 + t74 * t54) * qJD(5) + (-t57 * t126 + (pkin(4) * t57 + t118) * t55) * qJD(4), -t124 + t20 * t55 + (-t45 * t54 * t55 + t123) * qJD(4) + (t34 * t54 + t72 * t58) * qJD(5), t125 + t19 * t55 + (-t55 * t121 + t122) * qJD(4) + (t81 * t58 * t54 + t34 * t57) * qJD(5), (t37 * t101 - t20 * t58 - t4 + (-t11 + t122) * qJD(5)) * t57 + (-t38 * t101 + t19 * t58 - t3 + (-t12 + t123) * qJD(5)) * t54, pkin(5) * t89 + t11 * t20 - t12 * t19 + t17 * t45 + t3 * t37 + t4 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t100, 0, 0, 0, 0, 0, -t66, -t30, -t66, -t30, t112 * t100, (-t95 + (-t37 * t54 - t38 * t57) * qJD(4)) * t58 + (qJD(4) * t45 + t63) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t87, -0.2e1 * t79, 0, 0, 0, t54 * t97, t57 * t97, -0.2e1 * t81 * t104, t72 * t128, 0.2e1 * t63, 0.2e1 * t38 * t19 + 0.2e1 * t37 * t20 + 0.2e1 * t45 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, t5, -t6, 0, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t30, t100, t10, t9, (0.2e1 * pkin(5) - t118) * t100 + t61, t4, t66 * pkin(5), t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t26, -t29, t26, 0, -t29 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t104, 0, -pkin(9) * t47, pkin(9) * t104, t20, t19, -pkin(5) * t47, t20 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t66, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t47, 0, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t18;
