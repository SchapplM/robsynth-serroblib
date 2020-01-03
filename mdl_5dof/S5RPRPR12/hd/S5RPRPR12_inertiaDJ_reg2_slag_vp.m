% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR12
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR12_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:15
% EndTime: 2019-12-31 18:30:19
% DurationCPUTime: 1.17s
% Computational Cost: add. (2021->144), mult. (4691->281), div. (0->0), fcn. (4637->8), ass. (0->95)
t64 = cos(pkin(9));
t56 = -pkin(4) * t64 - pkin(3);
t126 = 0.2e1 * t56;
t121 = cos(qJ(3));
t63 = sin(pkin(8));
t65 = cos(pkin(8));
t67 = sin(qJ(3));
t80 = t121 * t65 - t63 * t67;
t40 = t80 * qJD(3);
t47 = t121 * t63 + t65 * t67;
t41 = t47 * qJD(3);
t85 = t40 * t80 - t41 * t47;
t125 = 0.2e1 * t85;
t66 = sin(qJ(5));
t102 = qJD(5) * t66;
t62 = sin(pkin(9));
t120 = cos(qJ(5));
t95 = qJD(5) * t120;
t124 = -t102 * t62 + t64 * t95;
t107 = pkin(6) + qJ(2);
t48 = t107 * t63;
t50 = t107 * t65;
t81 = -t121 * t48 - t50 * t67;
t44 = -t120 * t64 + t62 * t66;
t74 = qJD(2) * t80 + qJD(3) * t81;
t123 = 0.2e1 * t40;
t122 = 0.2e1 * t41;
t33 = t121 * t50 - t48 * t67;
t21 = qJD(2) * t47 + qJD(3) * t33;
t119 = t81 * t21;
t110 = t66 * t64;
t46 = t120 * t62 + t110;
t39 = t46 * qJD(5);
t118 = t44 * t39;
t117 = t46 * t124;
t116 = t47 * t62;
t115 = t62 * t33;
t114 = t62 * t40;
t113 = t62 * t41;
t111 = t64 * t40;
t106 = pkin(7) + qJ(4);
t13 = t124 * t47 + t46 * t40;
t23 = t46 * t47;
t105 = -t124 * t23 - t13 * t46;
t57 = -t65 * pkin(2) - pkin(1);
t83 = -pkin(3) * t80 + t57;
t77 = -qJ(4) * t47 + t83;
t15 = t33 * t64 + t62 * t77;
t58 = t62 ^ 2;
t60 = t64 ^ 2;
t103 = t58 + t60;
t27 = t80 * t122;
t101 = t47 * t123;
t100 = t62 * t111;
t96 = t106 * t62;
t94 = t120 * qJD(4);
t93 = 0.2e1 * t103 * qJD(4);
t92 = 0.2e1 * (t63 ^ 2 + t65 ^ 2) * qJD(2);
t72 = t62 * t74;
t84 = t41 * pkin(3) - t47 * qJD(4);
t78 = -qJ(4) * t40 + t84;
t8 = t64 * t78 - t72;
t71 = t64 * t74;
t9 = t62 * t78 + t71;
t90 = t62 * t9 + t64 * t8;
t89 = -t62 * t8 + t64 * t9;
t12 = t47 * t39 + t40 * t44;
t24 = t44 * t47;
t88 = -t12 * t44 - t24 * t39;
t87 = t21 * t47 - t40 * t81;
t86 = t124 * t80 - t41 * t46;
t82 = t120 * t96;
t79 = -pkin(3) * t40 - qJ(4) * t41 + qJD(4) * t80;
t76 = -t106 * t40 + t84;
t75 = -t80 * pkin(4) - t115 + (-t106 * t47 + t83) * t64;
t73 = t66 * t75;
t70 = t120 * t75;
t69 = t62 * t76 + t71;
t68 = t41 * pkin(4) + t64 * t76 - t72;
t49 = t106 * t64;
t36 = t64 * t41;
t32 = t120 * t49 - t66 * t96;
t30 = -t66 * t49 - t82;
t22 = pkin(4) * t116 - t81;
t20 = -t49 * t95 - qJD(4) * t110 + (t102 * t106 - t94) * t62;
t19 = qJD(5) * t82 - t64 * t94 + (qJD(4) * t62 + qJD(5) * t49) * t66;
t17 = t39 * t80 - t41 * t44;
t16 = pkin(4) * t114 + t21;
t14 = t64 * t77 - t115;
t11 = -pkin(7) * t116 + t15;
t4 = t11 * t120 + t73;
t3 = -t66 * t11 + t70;
t2 = -qJD(5) * t73 - t11 * t95 + t120 * t68 - t66 * t69;
t1 = -qJD(5) * t70 + t102 * t11 - t120 * t69 - t66 * t68;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, qJ(2) * t92, t101, t125, 0, -t27, 0, 0, t57 * t122, t57 * t123, -0.2e1 * t33 * t41 + 0.2e1 * t74 * t80 + 0.2e1 * t87, 0.2e1 * t33 * t74 - 0.2e1 * t119, t60 * t101, -0.4e1 * t47 * t100, -0.2e1 * t85 * t64, t58 * t101, t62 * t125, -t27, 0.2e1 * t14 * t41 + 0.2e1 * t62 * t87 - 0.2e1 * t8 * t80, -0.2e1 * t15 * t41 + 0.2e1 * t64 * t87 + 0.2e1 * t80 * t9, -0.2e1 * t90 * t47 - 0.2e1 * (t14 * t64 + t15 * t62) * t40, 0.2e1 * t14 * t8 + 0.2e1 * t15 * t9 - 0.2e1 * t119, 0.2e1 * t24 * t12, 0.2e1 * t12 * t23 + 0.2e1 * t13 * t24, 0.2e1 * t12 * t80 - 0.2e1 * t24 * t41, 0.2e1 * t23 * t13, 0.2e1 * t13 * t80 - 0.2e1 * t23 * t41, -t27, 0.2e1 * t13 * t22 + 0.2e1 * t16 * t23 - 0.2e1 * t2 * t80 + 0.2e1 * t3 * t41, -0.2e1 * t1 * t80 - 0.2e1 * t12 * t22 - 0.2e1 * t16 * t24 - 0.2e1 * t4 * t41, 0.2e1 * t1 * t23 + 0.2e1 * t12 * t3 - 0.2e1 * t13 * t4 + 0.2e1 * t2 * t24, -0.2e1 * t1 * t4 + 0.2e1 * t16 * t22 + 0.2e1 * t2 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t113, -t103 * t40, t90, 0, 0, 0, 0, 0, 0, t17, t86, t88 + t105, -t1 * t46 + t124 * t4 - t2 * t44 - t3 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t117 + 0.2e1 * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t41, 0, -t21, -t74, 0, 0, t100, -(t58 - t60) * t40, t113, -t100, t36, 0, -t21 * t64 + t62 * t79, t21 * t62 + t64 * t79, t89, -pkin(3) * t21 + (-t14 * t62 + t15 * t64) * qJD(4) + t89 * qJ(4), -t12 * t46 - t124 * t24, -t88 + t105, -t86, t13 * t44 + t23 * t39, t17, 0, t13 * t56 + t16 * t44 - t20 * t80 + t22 * t39 + t30 * t41, -t12 * t56 + t124 * t22 + t16 * t46 - t19 * t80 - t32 * t41, t1 * t44 + t12 * t30 - t124 * t3 - t13 * t32 + t19 * t23 - t2 * t46 + t20 * t24 - t39 * t4, -t1 * t32 + t16 * t56 - t19 * t4 + t2 * t30 + t20 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124 * t32 - t19 * t46 - t20 * t44 - t30 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, qJ(4) * t93, 0.2e1 * t117, -0.2e1 * t124 * t44 - 0.2e1 * t39 * t46, 0, 0.2e1 * t118, 0, 0, t39 * t126, t124 * t126, -0.2e1 * t124 * t30 + 0.2e1 * t19 * t44 - 0.2e1 * t20 * t46 - 0.2e1 * t32 * t39, -0.2e1 * t19 * t32 + 0.2e1 * t20 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t111, 0, t21, 0, 0, 0, 0, 0, 0, t13, -t12, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t124, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, -t13, t41, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t124, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, 0, -t39, 0, t20, t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
