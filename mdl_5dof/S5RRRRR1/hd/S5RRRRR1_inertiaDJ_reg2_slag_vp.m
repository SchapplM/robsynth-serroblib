% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:51:00
% EndTime: 2019-12-05 18:51:07
% DurationCPUTime: 1.45s
% Computational Cost: add. (1703->136), mult. (4188->270), div. (0->0), fcn. (4131->8), ass. (0->101)
t54 = sin(qJ(5));
t52 = t54 ^ 2;
t57 = cos(qJ(5));
t53 = t57 ^ 2;
t111 = t52 + t53;
t115 = sin(qJ(3));
t116 = cos(qJ(4));
t55 = sin(qJ(4));
t117 = cos(qJ(3));
t98 = t117 * pkin(2);
t84 = t98 + pkin(3);
t73 = t116 * t84;
t87 = qJD(3) * t98;
t90 = t115 * qJD(3);
t16 = -qJD(4) * t73 - t116 * t87 + (t115 * qJD(4) + t90) * t55 * pkin(2);
t128 = t111 * t16;
t83 = t116 * t115;
t127 = (qJD(3) + qJD(4)) * pkin(2) * (t117 * t55 + t83);
t112 = t52 - t53;
t125 = t112 * qJD(5);
t124 = qJD(2) + qJD(3);
t56 = sin(qJ(2));
t58 = cos(qJ(2));
t35 = t115 * t58 + t117 * t56;
t71 = t115 * t56 - t117 * t58;
t67 = t55 * t71;
t23 = -t116 * t35 + t67;
t109 = qJD(4) * t55;
t65 = t71 * qJD(3);
t59 = qJD(2) * t71 + t65;
t66 = t35 * qJD(3);
t60 = qJD(2) * t35 + t66;
t63 = t116 * t71;
t7 = -qJD(4) * t63 - t35 * t109 - t116 * t59 - t55 * t60;
t122 = t23 * t7;
t91 = qJD(4) * t116;
t8 = qJD(4) * t67 - t116 * t60 - t35 * t91 + t55 * t59;
t121 = t54 * t8;
t120 = t56 * pkin(2);
t119 = t57 * t7;
t118 = t57 * t8;
t49 = t58 * pkin(2) + pkin(1);
t99 = pkin(3) * t109;
t17 = t99 + t127;
t97 = t115 * pkin(2);
t31 = -t55 * t97 + t73;
t27 = -pkin(4) - t31;
t50 = qJD(5) * t57;
t114 = t17 * t54 + t27 * t50;
t85 = pkin(3) * t91;
t30 = t111 * t85;
t48 = -t116 * pkin(3) - pkin(4);
t113 = t48 * t50 + t54 * t99;
t32 = pkin(2) * t83 + t55 * t84;
t110 = qJD(4) * pkin(3);
t107 = qJD(5) * t54;
t106 = t56 * qJD(2);
t105 = t58 * qJD(2);
t22 = -t35 * t55 - t63;
t104 = 0.2e1 * t22 * t8;
t103 = t54 * t119;
t102 = -0.2e1 * t106;
t101 = pkin(4) * t107;
t100 = pkin(4) * t50;
t96 = t54 * t50;
t95 = t56 * t105;
t26 = -pkin(3) * t71 + t49;
t9 = t22 * pkin(4) - t23 * pkin(6) + t26;
t94 = t111 * t9;
t24 = t27 * t107;
t92 = -t17 * t57 + t24;
t89 = -0.2e1 * t99;
t21 = t23 ^ 2;
t88 = t21 * t96;
t86 = pkin(2) * t90;
t82 = t16 * t22 + t17 * t23;
t28 = pkin(6) + t32;
t81 = t22 * t28 - t23 * t27;
t47 = t55 * pkin(3) + pkin(6);
t80 = t22 * t47 - t23 * t48;
t36 = t48 * t107;
t79 = -t57 * t99 + t36;
t14 = -pkin(3) * t66 + (-pkin(3) * t35 - t120) * qJD(2);
t2 = t8 * pkin(4) + t7 * pkin(6) + t14;
t78 = -t2 * t54 - t50 * t9;
t77 = -t107 * t9 + t2 * t57;
t5 = t22 * t50 + t121;
t76 = t107 * t22 - t118;
t75 = -t23 * t50 + t54 * t7;
t74 = t107 * t23 + t119;
t72 = t111 * t116;
t70 = (-t116 * t22 + t23 * t55) * qJD(4);
t68 = t49 * t35;
t64 = -t27 * t7 - t28 * t8 + t82;
t62 = pkin(3) * t70 - t47 * t8 - t48 * t7;
t45 = -0.2e1 * t96;
t44 = 0.2e1 * t96;
t34 = -0.2e1 * t125;
t3 = t23 * t125 + t103;
t1 = t112 * t7 - 0.4e1 * t23 * t96;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t95, 0.2e1 * (-t56 ^ 2 + t58 ^ 2) * qJD(2), 0, -0.2e1 * t95, 0, 0, pkin(1) * t102, -0.2e1 * pkin(1) * t105, 0, 0, -0.2e1 * t35 * t59, 0.2e1 * t124 * (-t35 ^ 2 + t71 ^ 2), 0, 0.2e1 * t71 * t60, 0, 0, -0.2e1 * qJD(3) * t68 + 0.2e1 * (t71 * t120 - t68) * qJD(2), 0.2e1 * t49 * t65 + 0.2e1 * (t35 * t120 + t49 * t71) * qJD(2), 0, t49 * pkin(2) * t102, -0.2e1 * t122, 0.2e1 * t22 * t7 - 0.2e1 * t23 * t8, 0, t104, 0, 0, 0.2e1 * t14 * t22 + 0.2e1 * t26 * t8, 0.2e1 * t14 * t23 - 0.2e1 * t26 * t7, 0, 0.2e1 * t26 * t14, -0.2e1 * t53 * t122 - 0.2e1 * t88, 0.4e1 * t103 * t23 + 0.2e1 * t21 * t125, 0.2e1 * t23 * t118 - 0.2e1 * t22 * t74, -0.2e1 * t52 * t122 + 0.2e1 * t88, -0.2e1 * t23 * t121 + 0.2e1 * t22 * t75, t104, 0.2e1 * t9 * t118 + 0.2e1 * t22 * t77, -0.2e1 * t9 * t121 + 0.2e1 * t22 * t78, -0.2e1 * t111 * t2 * t23 + 0.2e1 * t7 * t94, 0.2e1 * t2 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, 0, t106, 0, 0, 0, 0, 0, 0, 0, t59, 0, t60, 0, 0, 0, -t35 * t86 + t60 * t97 + (-t124 * t71 + t65) * t98, 0, 0, 0, -t7, 0, -t8, 0, 0, 0, t31 * t7 - t32 * t8 + t82, 0, -t3, t1, t5, t3, -t76, 0, -t50 * t81 + t54 * t64, t107 * t81 + t57 * t64, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t86, -0.2e1 * t87, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17, 0.2e1 * t16, 0, -0.2e1 * t16 * t32 - 0.2e1 * t17 * t31, t44, t34, 0, t45, 0, 0, 0.2e1 * t92, 0.2e1 * t114, -0.2e1 * t128, -0.2e1 * t128 * t28 + 0.2e1 * t27 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t60, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, 0, 0, (t116 * t7 - t55 * t8 + t70) * pkin(3), 0, -t3, t1, t5, t3, -t76, 0, -t50 * t80 + t54 * t62, t107 * t80 + t57 * t62, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t87, 0, 0, 0, 0, 0, 0, 0, 0, t89 - t127, -t85 + t16, 0, (-t116 * t17 - t16 * t55 + (t116 * t32 - t31 * t55) * qJD(4)) * pkin(3), t44, t34, 0, t45, 0, 0, t24 + t36 + (-t17 - t99) * t57, t113 + t114, t30 - t128, t17 * t48 - t47 * t128 + (t27 * t55 + t28 * t72) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -0.2e1 * t85, 0, 0, t44, t34, 0, t45, 0, 0, 0.2e1 * t79, 0.2e1 * t113, 0.2e1 * t30, 0.2e1 * (t47 * t72 + t48 * t55) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, 0, 0, 0, 0, -t3, t1, t5, t3, -t76, 0, pkin(4) * t75 - pkin(6) * t5, pkin(4) * t74 + pkin(6) * t76, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, 0, t44, t34, 0, t45, 0, 0, t92 - t101, -t100 + t114, -t128, -t17 * pkin(4) - pkin(6) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t85, 0, 0, t44, t34, 0, t45, 0, 0, t79 - t101, -t100 + t113, t30, (-pkin(4) * t55 + pkin(6) * t72) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t34, 0, t45, 0, 0, -0.2e1 * t101, -0.2e1 * t100, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, 0, t75, t8, t77, t78, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t107, 0, t16 * t54 - t28 * t50, t107 * t28 + t16 * t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t107, 0, -t47 * t50 - t54 * t85, t107 * t47 - t57 * t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t107, 0, -pkin(6) * t50, pkin(6) * t107, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
