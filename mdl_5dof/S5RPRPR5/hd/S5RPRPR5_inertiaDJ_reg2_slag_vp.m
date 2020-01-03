% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:02
% EndTime: 2020-01-03 11:43:12
% DurationCPUTime: 1.55s
% Computational Cost: add. (1780->150), mult. (4305->292), div. (0->0), fcn. (4009->8), ass. (0->96)
t70 = sin(pkin(9));
t71 = sin(pkin(8));
t113 = qJ(4) * t71;
t73 = cos(pkin(8));
t58 = -t73 * pkin(2) - t71 * pkin(6) - pkin(1);
t76 = cos(qJ(3));
t63 = t76 * t73 * qJ(2);
t75 = sin(qJ(3));
t38 = t75 * t58 + t63;
t80 = t75 * t113 - t38;
t126 = t70 * t80;
t72 = cos(pkin(9));
t125 = t72 * t80;
t108 = qJD(3) * t76;
t99 = t71 * t108;
t59 = t72 * t99;
t109 = qJD(3) * t75;
t98 = t70 * t109;
t89 = t71 * t98 - t59;
t124 = -0.2e1 * t89;
t119 = cos(qJ(5));
t121 = pkin(3) * t70;
t74 = sin(qJ(5));
t96 = t72 * pkin(3) + pkin(4);
t46 = t119 * t96 - t74 * t121;
t111 = qJD(2) * t76;
t115 = t58 * t108 + t73 * t111;
t104 = qJ(2) * qJD(3);
t95 = t75 * t104;
t29 = t73 * t95 - t115;
t112 = qJD(2) * t75;
t97 = t73 * t112;
t30 = -t38 * qJD(3) - t97;
t114 = qJ(2) * t75;
t101 = t73 * t114;
t53 = t76 * t58;
t37 = t53 - t101;
t123 = t29 * t75 - t30 * t76 + (t37 * t75 - t38 * t76) * qJD(3);
t122 = 0.2e1 * t73;
t55 = t70 * t76 + t72 * t75;
t107 = qJD(4) * t71;
t83 = -t76 * t107 - t97;
t88 = -t75 * t107 + t115;
t79 = t70 * t83 + t72 * t88;
t100 = t76 * t113;
t84 = -t100 - t101;
t82 = t72 * t84;
t9 = (t82 + t126) * qJD(3) + t79;
t120 = t9 * t55;
t116 = t72 * t76;
t54 = -t70 * t75 + t116;
t49 = t54 * qJD(3);
t118 = t55 * t49;
t117 = t71 * t75;
t24 = -t100 + t53 + (-pkin(3) - t114) * t73;
t18 = t70 * t24 - t125;
t51 = pkin(3) * t99 + t71 * qJD(2);
t56 = pkin(3) * t117 + t71 * qJ(2);
t110 = qJD(3) * t71;
t106 = qJD(5) * t74;
t105 = qJ(2) * qJD(2);
t17 = t72 * t24 + t126;
t94 = qJD(5) * t119;
t93 = t110 * t122;
t68 = t71 ^ 2;
t92 = t68 * t75 * t108;
t44 = t71 * t116 - t70 * t117;
t13 = -t73 * pkin(4) - t44 * pkin(7) + t17;
t85 = t55 * t71;
t14 = -pkin(7) * t85 + t18;
t77 = -t59 * pkin(7) + (t82 + (-t63 + (-t58 + (pkin(7) + qJ(4)) * t71) * t75) * t70) * qJD(3) + t79;
t39 = t55 * t110;
t8 = -t70 * t88 + t72 * t83 + (-t70 * t84 + t125) * qJD(3);
t78 = t39 * pkin(7) + t8;
t1 = t14 * t106 - t119 * t77 - t13 * t94 - t74 * t78;
t26 = t119 * t55 + t74 * t54;
t81 = t119 * t85;
t47 = t119 * t121 + t74 * t96;
t20 = t119 * t44 - t74 * t85;
t2 = -t13 * t106 + t119 * t78 - t14 * t94 - t74 * t77;
t69 = t73 ^ 2;
t65 = t68 * t105;
t48 = t55 * qJD(3);
t41 = t47 * qJD(5);
t40 = t46 * qJD(5);
t27 = pkin(4) * t85 + t56;
t25 = t119 * t54 - t74 * t55;
t21 = -t89 * pkin(4) + t51;
t19 = t74 * t44 + t81;
t16 = -t26 * qJD(5) - t119 * t48 - t74 * t49;
t15 = t55 * t106 - t119 * t49 + t74 * t48 - t54 * t94;
t12 = t20 * qJD(5) - t119 * t89 - t74 * t39;
t11 = qJD(5) * t81 + t44 * t106 + t119 * t39 - t74 * t89;
t4 = t119 * t14 + t74 * t13;
t3 = t119 * t13 - t74 * t14;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t68 + t69) * qJD(2), 0.2e1 * t69 * t105 + 0.2e1 * t65, -0.2e1 * t92, 0.2e1 * (t75 ^ 2 - t76 ^ 2) * t68 * qJD(3), t75 * t93, 0.2e1 * t92, t76 * t93, 0, -0.2e1 * t30 * t73 + 0.2e1 * (t76 * t104 + t112) * t68, -0.2e1 * t29 * t73 + 0.2e1 * (-t95 + t111) * t68, 0.2e1 * t123 * t71, -0.2e1 * t38 * t29 + 0.2e1 * t37 * t30 + 0.2e1 * t65, -0.2e1 * t44 * t39, -0.2e1 * t44 * t59 + 0.2e1 * (t39 * t55 + t44 * t98) * t71, t39 * t122, t85 * t124, t73 * t124, 0, 0.2e1 * t56 * t59 - 0.2e1 * t8 * t73 + 0.2e1 * (t51 * t55 - t56 * t98) * t71, -0.2e1 * t56 * t39 + 0.2e1 * t51 * t44 + 0.2e1 * t9 * t73, 0.2e1 * t17 * t39 - 0.2e1 * t18 * t59 - 0.2e1 * t8 * t44 + 0.2e1 * (t18 * t98 - t120) * t71, 0.2e1 * t17 * t8 + 0.2e1 * t18 * t9 + 0.2e1 * t56 * t51, -0.2e1 * t20 * t11, 0.2e1 * t11 * t19 - 0.2e1 * t20 * t12, t11 * t122, 0.2e1 * t19 * t12, t12 * t122, 0, 0.2e1 * t27 * t12 + 0.2e1 * t21 * t19 - 0.2e1 * t2 * t73, -0.2e1 * t1 * t73 - 0.2e1 * t27 * t11 + 0.2e1 * t21 * t20, 0.2e1 * t1 * t19 + 0.2e1 * t3 * t11 - 0.2e1 * t4 * t12 - 0.2e1 * t2 * t20, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t27 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t109, t73 * t108, 0, -t123, 0, 0, 0, 0, 0, 0, t48 * t73, t49 * t73, t54 * t39 + t48 * t44 - t55 * t59 + (t55 * t98 - t118) * t71, -t17 * t48 + t18 * t49 + t8 * t54 + t120, 0, 0, 0, 0, 0, 0, -t16 * t73, -t15 * t73, t25 * t11 - t26 * t12 + t15 * t19 - t16 * t20, -t1 * t26 - t4 * t15 + t3 * t16 + t2 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t54 * t48 + 0.2e1 * t118, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t26 * t15 + 0.2e1 * t25 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71 * t109, 0, -t99, 0, t30, t29, 0, 0, 0, 0, -t39, 0, t89, 0, t8, -t9, (t39 * t72 + t70 * t89) * pkin(3), (t70 * t9 + t72 * t8) * pkin(3), 0, 0, -t11, 0, -t12, 0, t41 * t73 + t2, t40 * t73 + t1, t46 * t11 - t47 * t12 - t40 * t19 + t41 * t20, -t1 * t47 + t2 * t46 - t3 * t41 + t4 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t108, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t49, 0, (-t48 * t72 + t49 * t70) * pkin(3), 0, 0, 0, 0, 0, 0, t16, t15, 0, -t15 * t47 + t16 * t46 - t25 * t41 + t26 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t41, -0.2e1 * t40, 0, 0.2e1 * t47 * t40 - 0.2e1 * t46 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t39, 0, t51, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t12, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
