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
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
% StartTime: 2018-11-16 14:52:10
% EndTime: 2018-11-16 14:52:14
% DurationCPUTime: 1.42s
% Computational Cost: add. (1703->136), mult. (4188->270), div. (0->0), fcn. (4131->8), ass. (0->101)
t54 = sin(qJ(5));
t52 = t54 ^ 2;
t57 = cos(qJ(5));
t53 = t57 ^ 2;
t112 = t52 + t53;
t116 = sin(qJ(3));
t117 = cos(qJ(4));
t55 = sin(qJ(4));
t118 = cos(qJ(3));
t99 = t118 * pkin(2);
t85 = t99 + pkin(3);
t74 = t117 * t85;
t88 = qJD(3) * t99;
t91 = t116 * qJD(3);
t16 = -qJD(4) * t74 - t117 * t88 + (t116 * qJD(4) + t91) * t55 * pkin(2);
t129 = t112 * t16;
t84 = t117 * t116;
t128 = (qJD(3) + qJD(4)) * pkin(2) * (t118 * t55 + t84);
t113 = t52 - t53;
t126 = t113 * qJD(5);
t125 = qJD(2) + qJD(3);
t56 = sin(qJ(2));
t58 = cos(qJ(2));
t72 = t116 * t56 - t118 * t58;
t67 = t55 * t72;
t71 = -t116 * t58 - t118 * t56;
t23 = t117 * t71 + t67;
t110 = qJD(4) * t55;
t65 = t72 * qJD(3);
t59 = qJD(2) * t72 + t65;
t66 = t71 * qJD(3);
t60 = -qJD(2) * t71 - t66;
t63 = t117 * t72;
t7 = -qJD(4) * t63 + t110 * t71 - t117 * t59 - t55 * t60;
t123 = t23 * t7;
t92 = qJD(4) * t117;
t8 = qJD(4) * t67 - t117 * t60 + t55 * t59 + t71 * t92;
t122 = t54 * t8;
t121 = t56 * pkin(2);
t120 = t57 * t7;
t119 = t57 * t8;
t49 = t58 * pkin(2) + pkin(1);
t100 = pkin(3) * t110;
t17 = t100 + t128;
t98 = t116 * pkin(2);
t31 = -t55 * t98 + t74;
t27 = -pkin(4) - t31;
t50 = qJD(5) * t57;
t115 = t17 * t54 + t27 * t50;
t86 = pkin(3) * t92;
t30 = t112 * t86;
t48 = -t117 * pkin(3) - pkin(4);
t114 = t54 * t100 + t48 * t50;
t32 = pkin(2) * t84 + t55 * t85;
t111 = qJD(4) * pkin(3);
t108 = qJD(5) * t54;
t107 = t56 * qJD(2);
t106 = t58 * qJD(2);
t22 = t55 * t71 - t63;
t105 = 0.2e1 * t22 * t8;
t104 = t54 * t120;
t103 = -0.2e1 * t107;
t102 = pkin(4) * t108;
t101 = pkin(4) * t50;
t97 = t54 * t50;
t96 = t56 * t106;
t26 = -pkin(3) * t72 + t49;
t9 = t22 * pkin(4) - t23 * pkin(6) + t26;
t95 = t112 * t9;
t24 = t27 * t108;
t93 = -t17 * t57 + t24;
t90 = -0.2e1 * t100;
t21 = t23 ^ 2;
t89 = t21 * t97;
t87 = pkin(2) * t91;
t83 = t16 * t22 + t17 * t23;
t28 = pkin(6) + t32;
t82 = t22 * t28 - t23 * t27;
t47 = t55 * pkin(3) + pkin(6);
t81 = t22 * t47 - t23 * t48;
t36 = t48 * t108;
t80 = -t100 * t57 + t36;
t14 = pkin(3) * t66 + (pkin(3) * t71 - t121) * qJD(2);
t2 = t8 * pkin(4) + t7 * pkin(6) + t14;
t79 = -t2 * t54 - t50 * t9;
t78 = -t108 * t9 + t2 * t57;
t5 = t22 * t50 + t122;
t77 = t108 * t22 - t119;
t76 = -t23 * t50 + t54 * t7;
t75 = t108 * t23 + t120;
t73 = t112 * t117;
t70 = (-t117 * t22 + t23 * t55) * qJD(4);
t68 = t49 * t71;
t64 = -t27 * t7 - t28 * t8 + t83;
t62 = pkin(3) * t70 - t47 * t8 - t48 * t7;
t45 = -0.2e1 * t97;
t44 = 0.2e1 * t97;
t34 = -0.2e1 * t126;
t3 = t126 * t23 + t104;
t1 = t113 * t7 - 0.4e1 * t23 * t97;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t96, 0.2e1 * (-t56 ^ 2 + t58 ^ 2) * qJD(2), 0, -0.2e1 * t96, 0, 0, pkin(1) * t103, -0.2e1 * pkin(1) * t106, 0, 0, 0.2e1 * t71 * t59, 0.2e1 * t125 * (-t71 ^ 2 + t72 ^ 2) 0, 0.2e1 * t72 * t60, 0, 0, 0.2e1 * qJD(3) * t68 + 0.2e1 * (t72 * t121 + t68) * qJD(2), 0.2e1 * t49 * t65 + 0.2e1 * (-t121 * t71 + t49 * t72) * qJD(2), 0, t49 * pkin(2) * t103, -0.2e1 * t123, 0.2e1 * t22 * t7 - 0.2e1 * t23 * t8, 0, t105, 0, 0, 0.2e1 * t14 * t22 + 0.2e1 * t26 * t8, 0.2e1 * t14 * t23 - 0.2e1 * t26 * t7, 0, 0.2e1 * t26 * t14, -0.2e1 * t53 * t123 - 0.2e1 * t89, 0.4e1 * t104 * t23 + 0.2e1 * t21 * t126, 0.2e1 * t23 * t119 - 0.2e1 * t22 * t75, -0.2e1 * t52 * t123 + 0.2e1 * t89, -0.2e1 * t23 * t122 + 0.2e1 * t22 * t76, t105, 0.2e1 * t9 * t119 + 0.2e1 * t22 * t78, -0.2e1 * t9 * t122 + 0.2e1 * t22 * t79, -0.2e1 * t112 * t2 * t23 + 0.2e1 * t7 * t95, 0.2e1 * t2 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, 0, t107, 0, 0, 0, 0, 0, 0, 0, t59, 0, t60, 0, 0, 0, t71 * t87 + t60 * t98 + (-t125 * t72 + t65) * t99, 0, 0, 0, -t7, 0, -t8, 0, 0, 0, t31 * t7 - t32 * t8 + t83, 0, -t3, t1, t5, t3, -t77, 0, -t50 * t82 + t54 * t64, t108 * t82 + t57 * t64, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t87, -0.2e1 * t88, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17, 0.2e1 * t16, 0, -0.2e1 * t16 * t32 - 0.2e1 * t17 * t31, t44, t34, 0, t45, 0, 0, 0.2e1 * t93, 0.2e1 * t115, -0.2e1 * t129, -0.2e1 * t129 * t28 + 0.2e1 * t27 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t60, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, 0, 0 (t117 * t7 - t55 * t8 + t70) * pkin(3), 0, -t3, t1, t5, t3, -t77, 0, -t50 * t81 + t54 * t62, t108 * t81 + t57 * t62, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t88, 0, 0, 0, 0, 0, 0, 0, 0, t90 - t128, -t86 + t16, 0 (-t117 * t17 - t16 * t55 + (t117 * t32 - t31 * t55) * qJD(4)) * pkin(3), t44, t34, 0, t45, 0, 0, t24 + t36 + (-t17 - t100) * t57, t114 + t115, t30 - t129, t17 * t48 - t47 * t129 + (t27 * t55 + t28 * t73) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -0.2e1 * t86, 0, 0, t44, t34, 0, t45, 0, 0, 0.2e1 * t80, 0.2e1 * t114, 0.2e1 * t30, 0.2e1 * (t47 * t73 + t48 * t55) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, 0, 0, 0, 0, -t3, t1, t5, t3, -t77, 0, pkin(4) * t76 - pkin(6) * t5, pkin(4) * t75 + pkin(6) * t77, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, 0, t44, t34, 0, t45, 0, 0, t93 - t102, -t101 + t115, -t129, -t17 * pkin(4) - pkin(6) * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t86, 0, 0, t44, t34, 0, t45, 0, 0, t80 - t102, -t101 + t114, t30 (-pkin(4) * t55 + pkin(6) * t73) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t34, 0, t45, 0, 0, -0.2e1 * t102, -0.2e1 * t101, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, 0, t76, t8, t78, t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t108, 0, t16 * t54 - t28 * t50, t108 * t28 + t16 * t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t108, 0, -t47 * t50 - t54 * t86, t108 * t47 - t57 * t86, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t108, 0, -pkin(6) * t50, pkin(6) * t108, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t4;
