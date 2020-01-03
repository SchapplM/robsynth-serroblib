% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR13_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:52
% EndTime: 2019-12-31 20:33:56
% DurationCPUTime: 1.08s
% Computational Cost: add. (1517->167), mult. (3916->321), div. (0->0), fcn. (3669->8), ass. (0->105)
t85 = cos(pkin(9));
t88 = sin(qJ(2));
t127 = t85 * t88;
t91 = cos(qJ(2));
t98 = -t91 * pkin(2) - t88 * qJ(3);
t67 = -pkin(1) + t98;
t60 = t85 * t67;
t84 = sin(pkin(9));
t38 = -pkin(7) * t127 + t60 + (-pkin(6) * t84 - pkin(3)) * t91;
t129 = t84 * t88;
t126 = t85 * t91;
t76 = pkin(6) * t126;
t45 = t84 * t67 + t76;
t43 = -pkin(7) * t129 + t45;
t87 = sin(qJ(4));
t90 = cos(qJ(4));
t119 = t87 * t38 + t90 * t43;
t120 = pkin(7) + qJ(3);
t69 = t120 * t84;
t70 = t120 * t85;
t117 = -t87 * t69 + t90 * t70;
t113 = qJD(4) * t90;
t121 = t90 * t85;
t21 = (qJD(3) * t84 + qJD(4) * t70) * t87 - qJD(3) * t121 + t69 * t113;
t78 = -t85 * pkin(3) - pkin(2);
t133 = 0.2e1 * t78;
t65 = t90 * t84 + t87 * t85;
t56 = t65 * qJD(4);
t132 = pkin(4) * t56;
t112 = t91 * qJD(2);
t114 = qJD(4) * t88;
t123 = t87 * t84;
t30 = t65 * t112 + t113 * t127 - t114 * t123;
t124 = t87 * t43;
t80 = t88 * qJD(2);
t106 = pkin(6) * t80;
t53 = -t88 * qJD(3) + (pkin(2) * t88 - qJ(3) * t91) * qJD(2);
t41 = t84 * t106 + t85 * t53;
t26 = (pkin(3) * t88 - pkin(7) * t126) * qJD(2) + t41;
t128 = t84 * t91;
t48 = t84 * t53;
t32 = t48 + (-pkin(6) * t127 - pkin(7) * t128) * qJD(2);
t8 = qJD(4) * t124 - t38 * t113 - t87 * t26 - t90 * t32;
t7 = -t30 * pkin(8) - t8;
t89 = cos(qJ(5));
t131 = t89 * t7;
t130 = t91 * pkin(4);
t51 = t65 * t88;
t13 = -t51 * pkin(8) + t119;
t86 = sin(qJ(5));
t125 = t86 * t13;
t122 = t89 * t13;
t105 = t84 * t112;
t79 = pkin(6) * t112;
t57 = pkin(3) * t105 + t79;
t66 = pkin(3) * t129 + t88 * pkin(6);
t116 = pkin(4) * qJD(5);
t115 = qJD(3) * t91;
t111 = pkin(6) * t128;
t110 = -0.2e1 * pkin(1) * qJD(2);
t109 = pkin(4) * t80;
t108 = t86 * t116;
t107 = t89 * t116;
t104 = t88 * t112;
t64 = -t121 + t123;
t29 = -t64 * t112 - t65 * t114;
t9 = -t119 * qJD(4) + t90 * t26 - t87 * t32;
t6 = -t29 * pkin(8) + t109 + t9;
t103 = t89 * t6 - t86 * t7;
t101 = t90 * t38 - t124;
t52 = t64 * t88;
t12 = t52 * pkin(8) + t101 - t130;
t102 = -t12 + t130;
t100 = -t90 * t69 - t87 * t70;
t99 = 0.2e1 * (t84 ^ 2 + t85 ^ 2) * qJD(3);
t97 = t89 * t12 - t125;
t96 = t86 * t12 + t122;
t27 = -t65 * pkin(8) + t100;
t28 = -t64 * pkin(8) + t117;
t95 = t89 * t27 - t86 * t28;
t94 = t86 * t27 + t89 * t28;
t42 = -t85 * t106 + t48;
t93 = -t41 * t84 + t42 * t85;
t24 = t89 * t51 - t86 * t52;
t25 = -t86 * t51 - t89 * t52;
t36 = t89 * t64 + t86 * t65;
t37 = -t86 * t64 + t89 * t65;
t22 = -t65 * qJD(3) - t117 * qJD(4);
t73 = -0.2e1 * t104;
t55 = t64 * qJD(4);
t47 = t64 * pkin(4) + t78;
t44 = t60 - t111;
t40 = t51 * pkin(4) + t66;
t18 = t30 * pkin(4) + t57;
t17 = t55 * pkin(8) + t22;
t16 = -t56 * pkin(8) - t21;
t15 = t37 * qJD(5) - t86 * t55 + t89 * t56;
t14 = -t36 * qJD(5) - t89 * t55 - t86 * t56;
t11 = t25 * qJD(5) + t86 * t29 + t89 * t30;
t10 = -t24 * qJD(5) + t89 * t29 - t86 * t30;
t4 = -t94 * qJD(5) - t86 * t16 + t89 * t17;
t3 = -t95 * qJD(5) - t89 * t16 - t86 * t17;
t2 = -t96 * qJD(5) + t103;
t1 = -t97 * qJD(5) - t86 * t6 - t131;
t5 = [0, 0, 0, 0.2e1 * t104, 0.2e1 * (-t88 ^ 2 + t91 ^ 2) * qJD(2), 0, 0, 0, t88 * t110, t91 * t110, -0.2e1 * t41 * t91 + 0.2e1 * (t44 + 0.2e1 * t111) * t80, 0.2e1 * t42 * t91 + 0.2e1 * (-t45 + 0.2e1 * t76) * t80, 0.2e1 * (-t41 * t85 - t42 * t84) * t88 + 0.2e1 * (-t44 * t85 - t45 * t84) * t112, 0.2e1 * pkin(6) ^ 2 * t104 + 0.2e1 * t44 * t41 + 0.2e1 * t45 * t42, -0.2e1 * t52 * t29, -0.2e1 * t29 * t51 + 0.2e1 * t52 * t30, -0.2e1 * t29 * t91 - 0.2e1 * t52 * t80, 0.2e1 * t30 * t91 - 0.2e1 * t51 * t80, t73, 0.2e1 * t101 * t80 + 0.2e1 * t66 * t30 + 0.2e1 * t57 * t51 - 0.2e1 * t9 * t91, -0.2e1 * t119 * t80 + 0.2e1 * t66 * t29 - 0.2e1 * t57 * t52 - 0.2e1 * t8 * t91, 0.2e1 * t25 * t10, -0.2e1 * t10 * t24 - 0.2e1 * t25 * t11, -0.2e1 * t10 * t91 + 0.2e1 * t25 * t80, 0.2e1 * t11 * t91 - 0.2e1 * t24 * t80, t73, 0.2e1 * t40 * t11 + 0.2e1 * t18 * t24 - 0.2e1 * t2 * t91 + 0.2e1 * t97 * t80, -0.2e1 * t1 * t91 + 0.2e1 * t40 * t10 + 0.2e1 * t18 * t25 - 0.2e1 * t96 * t80; 0, 0, 0, 0, 0, t112, -t80, 0, -t79, t106, t84 * t115 + (t98 * t84 - t76) * qJD(2), t85 * t115 + (t98 * t85 + t111) * qJD(2), t93, -pkin(2) * t79 + (-t44 * t84 + t45 * t85) * qJD(3) + t93 * qJ(3), t29 * t65 + t52 * t55, -t29 * t64 - t65 * t30 + t55 * t51 + t52 * t56, t55 * t91 + t65 * t80, t56 * t91 - t64 * t80, 0, t100 * t80 - t22 * t91 + t78 * t30 + t66 * t56 + t57 * t64, -t117 * t80 - t21 * t91 + t78 * t29 - t66 * t55 + t57 * t65, t10 * t37 + t25 * t14, -t10 * t36 - t37 * t11 - t14 * t24 - t25 * t15, -t14 * t91 + t37 * t80, t15 * t91 - t36 * t80, 0, t47 * t11 + t24 * t132 + t40 * t15 + t18 * t36 - t4 * t91 + t95 * t80, t47 * t10 + t25 * t132 + t40 * t14 + t18 * t37 - t3 * t91 - t94 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, qJ(3) * t99, -0.2e1 * t65 * t55, 0.2e1 * t55 * t64 - 0.2e1 * t65 * t56, 0, 0, 0, t56 * t133, -t55 * t133, 0.2e1 * t37 * t14, -0.2e1 * t14 * t36 - 0.2e1 * t37 * t15, 0, 0, 0, 0.2e1 * t36 * t132 + 0.2e1 * t47 * t15, 0.2e1 * t37 * t132 + 0.2e1 * t47 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t85 * t112, 0, t79, 0, 0, 0, 0, 0, t30, t29, 0, 0, 0, 0, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t55, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t30, t80, t9, t8, 0, 0, t10, -t11, t80, t89 * t109 + (t102 * t86 - t122) * qJD(5) + t103, -t131 + (-t6 - t109) * t86 + (t102 * t89 + t125) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t56, 0, t22, t21, 0, 0, t14, -t15, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t108, -0.2e1 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t80, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
