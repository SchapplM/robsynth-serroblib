% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRR12
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR12_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:37
% EndTime: 2019-12-31 20:30:42
% DurationCPUTime: 1.40s
% Computational Cost: add. (1665->165), mult. (3593->309), div. (0->0), fcn. (3019->6), ass. (0->111)
t62 = sin(qJ(4));
t106 = qJD(4) * t62;
t63 = sin(qJ(2));
t65 = cos(qJ(2));
t122 = cos(qJ(4));
t88 = qJD(4) * t122;
t89 = t63 * t122;
t57 = t65 * qJD(2);
t94 = t62 * t57;
t17 = -qJD(2) * t89 - t65 * t106 + t63 * t88 + t94;
t113 = t63 * t62;
t33 = t65 * t122 + t113;
t18 = t33 * qJD(2) - t63 * t106 - t65 * t88;
t101 = t63 * qJD(3);
t109 = qJ(3) * t65;
t128 = pkin(2) + pkin(3);
t21 = (-t128 * t63 + t109) * qJD(2) + t101;
t135 = t17 * pkin(4) - t18 * pkin(8) + t21;
t61 = sin(qJ(5));
t59 = t61 ^ 2;
t64 = cos(qJ(5));
t60 = t64 ^ 2;
t110 = t59 - t60;
t42 = t110 * qJD(5);
t134 = 0.2e1 * t42;
t133 = t59 + t60;
t127 = pkin(6) - pkin(7);
t44 = t127 * t65;
t20 = t127 * t113 + t122 * t44;
t108 = t63 * qJ(3);
t30 = t128 * t65 + pkin(1) + t108;
t34 = -t65 * t62 + t89;
t69 = t33 * pkin(4) - t34 * pkin(8) + t30;
t67 = t64 * t69;
t6 = -t61 * t20 + t67;
t7 = t64 * t20 + t61 * t69;
t132 = -t6 * t61 + t64 * t7;
t102 = t63 * qJD(2);
t38 = t127 * t102;
t84 = t127 * t122;
t11 = t20 * qJD(4) - t62 * t38 - t84 * t57;
t83 = t122 * t128;
t39 = -t62 * qJ(3) - t83;
t36 = pkin(4) - t39;
t40 = t122 * qJ(3) - t62 * t128;
t37 = -pkin(8) + t40;
t130 = -t11 + (t33 * t37 - t34 * t36) * qJD(5);
t129 = 0.2e1 * qJD(3);
t121 = t11 * t61;
t120 = t11 * t64;
t19 = t62 * t44 - t63 * t84;
t119 = t19 * t11;
t23 = t62 * qJD(3) + t40 * qJD(4);
t118 = t19 * t23;
t117 = t23 * t61;
t116 = t23 * t64;
t115 = t34 * t18;
t114 = t61 * t17;
t112 = t64 * t17;
t111 = t64 * t18;
t107 = qJD(4) * t19;
t104 = qJD(5) * t61;
t103 = qJD(5) * t64;
t100 = 0.2e1 * t33 * t17;
t99 = -0.2e1 * pkin(1) * qJD(2);
t98 = -0.2e1 * pkin(4) * qJD(5);
t97 = t61 * t111;
t96 = pkin(6) * t102;
t95 = pkin(6) * t57;
t93 = t61 * t103;
t92 = t63 * t57;
t91 = t11 * t122;
t90 = t23 * t122;
t22 = qJ(3) * t106 - t122 * qJD(3) + qJD(4) * t83;
t15 = t133 * t22;
t87 = qJD(5) * (pkin(4) + t36);
t86 = qJD(5) * t122;
t31 = t34 ^ 2;
t85 = t31 * t93;
t82 = -pkin(4) * t18 - pkin(8) * t17;
t81 = pkin(4) * t34 + pkin(8) * t33;
t80 = -t6 * t64 - t61 * t7;
t78 = -t65 * pkin(2) - t108;
t77 = t22 * t33 + t23 * t34;
t75 = t34 * t103 + t61 * t18;
t74 = -t34 * t104 + t111;
t73 = t122 * t34 + t33 * t62;
t72 = t133 * t122;
t70 = t78 * qJD(2) + t65 * qJD(3);
t10 = t122 * t38 - t127 * t94 + t107;
t2 = -qJD(5) * t67 + t64 * t10 + t20 * t104 - t61 * t135;
t3 = -t7 * qJD(5) + t61 * t10 + t135 * t64;
t1 = t80 * qJD(5) - t2 * t64 - t3 * t61;
t68 = -qJD(5) * t19 - t17 * t37 + t18 * t36 + t77;
t66 = -t122 * t18 - t62 * t17 + (-t122 * t33 + t34 * t62) * qJD(4);
t48 = -0.2e1 * t92;
t47 = 0.2e1 * t92;
t46 = -0.2e1 * t93;
t45 = 0.2e1 * t93;
t43 = (-t63 ^ 2 + t65 ^ 2) * qJD(2);
t41 = -pkin(1) + t78;
t32 = -0.2e1 * t42;
t28 = t64 * t106 + t61 * t86;
t27 = t61 * t106 - t64 * t86;
t26 = t72 * qJD(4);
t25 = -t101 + (pkin(2) * t63 - t109) * qJD(2);
t13 = t33 * t103 + t114;
t12 = t33 * t104 - t112;
t8 = t34 * t42 - t97;
t5 = t110 * t18 + 0.4e1 * t34 * t93;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0.2e1 * t43, 0, t48, 0, 0, t63 * t99, t65 * t99, 0, 0, t47, 0, -0.2e1 * t43, 0, 0, t48, 0.2e1 * t41 * t102 - 0.2e1 * t25 * t65, 0, -0.2e1 * t25 * t63 - 0.2e1 * t41 * t57, 0.2e1 * t41 * t25, 0.2e1 * t115, -0.2e1 * t34 * t17 - 0.2e1 * t18 * t33, 0, t100, 0, 0, 0.2e1 * t30 * t17 + 0.2e1 * t21 * t33, 0.2e1 * t30 * t18 + 0.2e1 * t21 * t34, 0.2e1 * t10 * t33 + 0.2e1 * t11 * t34 - 0.2e1 * t20 * t17 + 0.2e1 * t19 * t18, -0.2e1 * t20 * t10 + 0.2e1 * t30 * t21 + 0.2e1 * t119, 0.2e1 * t60 * t115 - 0.2e1 * t85, t31 * t134 - 0.4e1 * t34 * t97, 0.2e1 * t34 * t112 + 0.2e1 * t33 * t74, 0.2e1 * t59 * t115 + 0.2e1 * t85, -0.2e1 * t34 * t114 - 0.2e1 * t33 * t75, t100, 0.2e1 * t34 * t121 + 0.2e1 * t6 * t17 + 0.2e1 * t19 * t75 + 0.2e1 * t3 * t33, 0.2e1 * t34 * t120 - 0.2e1 * t7 * t17 + 0.2e1 * t19 * t74 + 0.2e1 * t2 * t33, 0.2e1 * t80 * t18 + 0.2e1 * (-qJD(5) * t132 + t2 * t61 - t3 * t64) * t34, -0.2e1 * t7 * t2 + 0.2e1 * t6 * t3 + 0.2e1 * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, -t102, 0, -t95, t96, 0, 0, 0, t57, 0, 0, t102, 0, -t95, t70, -t96, t70 * pkin(6), 0, 0, -t18, 0, t17, 0, t11, -t10, -t40 * t17 - t39 * t18 + t77, -t10 * t40 - t11 * t39 - t20 * t22 + t118, t8, t5, -t13, -t8, t12, 0, -t130 * t64 + t68 * t61, t130 * t61 + t68 * t64, -t1, t1 * t37 + t11 * t36 - t132 * t22 + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, qJ(3) * t129, 0, 0, 0, 0, 0, 0, 0.2e1 * t23, -0.2e1 * t22, 0, -0.2e1 * t40 * t22 - 0.2e1 * t39 * t23, t45, t32, 0, t46, 0, 0, -0.2e1 * t36 * t104 + 0.2e1 * t116, -0.2e1 * t36 * t103 - 0.2e1 * t117, 0.2e1 * t15, -0.2e1 * t15 * t37 + 0.2e1 * t36 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, t95, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t91 - t10 * t62 + (t122 * t20 + t19 * t62) * qJD(4), 0, 0, 0, 0, 0, 0, -t73 * t103 + t61 * t66, t73 * t104 + t64 * t66, 0, -t91 + (t1 + t107) * t62 + t132 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t88, 0, -t90 - t22 * t62 + (t122 * t40 - t39 * t62) * qJD(4), 0, 0, 0, 0, 0, 0, t28, -t27, -t26, -t90 - t62 * t15 + (t36 * t62 + t37 * t72) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-t122 + t72) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t17, 0, -t11, t10, 0, 0, -t8, -t5, t13, t8, -t12, 0, -t120 + t82 * t61 + (t19 * t61 - t64 * t81) * qJD(5), t121 + t82 * t64 + (t19 * t64 + t61 * t81) * qJD(5), t1, -t11 * pkin(4) + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t22, 0, 0, t46, t134, 0, t45, 0, 0, t61 * t87 - t116, t64 * t87 + t117, -t15, -t23 * pkin(4) - pkin(8) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t88, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t27, t26, (-pkin(4) * t62 + pkin(8) * t72) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t32, 0, t46, 0, 0, t61 * t98, t64 * t98, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, -t75, t17, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, 0, t104, 0, -t37 * t103 + t61 * t22, t37 * t104 + t64 * t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t103 - t61 * t88, t62 * t104 - t64 * t88, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, -t104, 0, -pkin(8) * t103, pkin(8) * t104, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
