% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:13
% EndTime: 2019-12-31 17:28:17
% DurationCPUTime: 1.19s
% Computational Cost: add. (1219->176), mult. (3132->348), div. (0->0), fcn. (2507->6), ass. (0->103)
t119 = -pkin(7) - pkin(6);
t55 = cos(qJ(3));
t121 = t119 * t55;
t52 = sin(qJ(4));
t53 = sin(qJ(3));
t115 = cos(qJ(4));
t71 = t119 * t115;
t66 = t53 * t71;
t18 = t121 * t52 + t66;
t54 = sin(qJ(2));
t123 = -0.4e1 * t54;
t80 = t115 * qJD(4);
t122 = t115 * qJD(3) + t80;
t48 = t53 ^ 2;
t50 = t55 ^ 2;
t108 = t48 - t50;
t79 = qJD(3) * t108;
t49 = t54 ^ 2;
t56 = cos(qJ(2));
t78 = (-t56 ^ 2 + t49) * qJD(2);
t120 = qJD(3) + qJD(4);
t116 = t56 * pkin(2);
t69 = -t54 * pkin(6) - t116;
t65 = -pkin(1) + t69;
t32 = t55 * t65;
t68 = pkin(2) * t54 - pkin(6) * t56;
t64 = t68 * t53;
t109 = qJD(2) * t64 + qJD(3) * t32;
t47 = t54 * qJD(2);
t104 = qJD(3) * t56;
t91 = t53 * t104;
t63 = t55 * t47 + t91;
t12 = t63 * pkin(5) - t109;
t118 = pkin(5) * t53;
t117 = t55 * pkin(2);
t114 = t52 * t53;
t113 = t53 * t54;
t112 = t53 * t56;
t111 = t54 * t55;
t110 = t55 * t56;
t44 = pkin(5) * t110;
t23 = t53 * t65 + t44;
t106 = qJD(3) * t53;
t105 = qJD(3) * t55;
t103 = qJD(4) * t52;
t102 = t56 * qJD(2);
t101 = pkin(5) * t111;
t100 = pkin(5) * t112;
t99 = -0.2e1 * pkin(1) * qJD(2);
t98 = -0.2e1 * pkin(2) * qJD(3);
t97 = t52 * t113;
t96 = pkin(3) * t106;
t95 = pkin(3) * t103;
t94 = pkin(5) * t102;
t92 = t54 * t106;
t90 = t54 * t105;
t89 = t55 * t104;
t88 = t53 * t102;
t87 = t53 * t105;
t86 = t54 * t102;
t85 = t55 * t102;
t84 = -pkin(3) - t118;
t17 = -pkin(7) * t113 + t23;
t83 = t115 * t17;
t82 = t115 * t55;
t77 = 0.2e1 * t86;
t76 = t119 * t114;
t75 = t53 * t85;
t74 = t49 * t87;
t73 = -pkin(7) * t111 + t32;
t72 = pkin(3) * t80;
t70 = t115 * t102;
t22 = t32 - t100;
t67 = -t22 * t55 - t23 * t53;
t57 = (-t44 + (-t119 * t54 + pkin(1) + t116) * t53) * qJD(3) + (t56 * t121 + (-t84 + t117) * t54) * qJD(2);
t62 = t88 + t90;
t58 = -t62 * pkin(7) - t12;
t61 = t84 * t56 + t73;
t60 = t115 * t61;
t1 = -qJD(4) * t60 + t17 * t103 - t115 * t58 - t52 * t57;
t34 = t115 * t53 + t52 * t55;
t19 = -t115 * t121 + t76;
t13 = -t23 * qJD(3) + (pkin(5) * t113 + t55 * t68) * qJD(2);
t59 = t67 * qJD(3) - t12 * t55 - t13 * t53;
t6 = t52 * t61 + t83;
t16 = t120 * t34;
t46 = -t55 * pkin(3) - pkin(2);
t41 = -0.2e1 * t86;
t35 = (pkin(3) * t53 + pkin(5)) * t54;
t33 = -t82 + t114;
t25 = t54 * t82 - t97;
t24 = t34 * t54;
t21 = t62 * pkin(3) + t94;
t20 = t54 * t79 - t75;
t15 = t120 * t114 - t122 * t55;
t11 = -t19 * qJD(4) + (t55 * t71 - t76) * qJD(3);
t10 = -t18 * qJD(3) - qJD(4) * t66 - t121 * t103;
t9 = t53 * t70 - t52 * t92 - qJD(4) * t97 + (t52 * t102 + t122 * t54) * t55;
t8 = t16 * t54 + t52 * t88 - t55 * t70;
t5 = -t52 * t17 + t60;
t4 = t115 * t57;
t2 = -t6 * qJD(4) - t52 * t58 + t4;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -0.2e1 * t78, 0, t41, 0, 0, t54 * t99, t56 * t99, 0, 0, 0.2e1 * t50 * t86 - 0.2e1 * t74, t75 * t123 + 0.2e1 * t49 * t79, 0.2e1 * t54 * t91 + 0.2e1 * t55 * t78, 0.2e1 * t48 * t86 + 0.2e1 * t74, -0.2e1 * t53 * t78 + 0.2e1 * t54 * t89, t41, 0.2e1 * t22 * t47 - 0.2e1 * t13 * t56 + 0.2e1 * (t49 * t105 + t53 * t77) * pkin(5), -0.2e1 * t23 * t47 - 0.2e1 * t12 * t56 + 0.2e1 * (-t49 * t106 + t55 * t77) * pkin(5), 0.2e1 * t67 * t102 + 0.2e1 * (t12 * t53 - t13 * t55 + (t22 * t53 - t23 * t55) * qJD(3)) * t54, 0.2e1 * pkin(5) ^ 2 * t86 - 0.2e1 * t23 * t12 + 0.2e1 * t22 * t13, -0.2e1 * t25 * t8, 0.2e1 * t8 * t24 - 0.2e1 * t25 * t9, 0.2e1 * t25 * t47 + 0.2e1 * t8 * t56, 0.2e1 * t24 * t9, -0.2e1 * t24 * t47 + 0.2e1 * t9 * t56, t41, -0.2e1 * t2 * t56 + 0.2e1 * t21 * t24 + 0.2e1 * t35 * t9 + 0.2e1 * t5 * t47, -0.2e1 * t1 * t56 + 0.2e1 * t21 * t25 - 0.2e1 * t35 * t8 - 0.2e1 * t6 * t47, 0.2e1 * t1 * t24 - 0.2e1 * t2 * t25 + 0.2e1 * t5 * t8 - 0.2e1 * t6 * t9, -0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * t35 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, -t47, 0, -t94, pkin(5) * t47, 0, 0, -t20, -t108 * t102 + t87 * t123, t53 * t47 - t89, t20, t63, 0, (pkin(6) * t110 + (-t117 + t118) * t54) * qJD(3) + (t69 * t53 - t44) * qJD(2), (t64 + t101) * qJD(3) + (t69 * t55 + t100) * qJD(2), t59, -pkin(2) * t94 + t59 * pkin(6), -t25 * t15 - t8 * t34, t15 * t24 - t25 * t16 + t8 * t33 - t34 * t9, t15 * t56 + t34 * t47, t24 * t16 + t9 * t33, t16 * t56 - t33 * t47, 0, -t11 * t56 + t35 * t16 + t18 * t47 + t21 * t33 + t24 * t96 + t46 * t9, -t10 * t56 - t35 * t15 - t19 * t47 + t21 * t34 + t25 * t96 - t46 * t8, t1 * t33 + t10 * t24 - t11 * t25 + t5 * t15 - t6 * t16 + t18 * t8 - t19 * t9 - t2 * t34, -t1 * t19 - t6 * t10 + t5 * t11 + t2 * t18 + t21 * t46 + t35 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t87, -0.2e1 * t79, 0, -0.2e1 * t87, 0, 0, t53 * t98, t55 * t98, 0, 0, -0.2e1 * t34 * t15, 0.2e1 * t15 * t33 - 0.2e1 * t34 * t16, 0, 0.2e1 * t33 * t16, 0, 0, 0.2e1 * t46 * t16 + 0.2e1 * t33 * t96, -0.2e1 * t46 * t15 + 0.2e1 * t34 * t96, 0.2e1 * t10 * t33 - 0.2e1 * t11 * t34 + 0.2e1 * t18 * t15 - 0.2e1 * t19 * t16, -0.2e1 * t19 * t10 + 0.2e1 * t18 * t11 + 0.2e1 * t46 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 - t92, 0, -t62, t47, t13, t12, 0, 0, 0, 0, -t8, 0, -t9, t47, -t52 * (-pkin(5) * t91 - pkin(7) * t90 + t109) + t4 + (t115 * pkin(3) * t54 - t52 * (-pkin(7) * t112 - t101)) * qJD(2) + (-t83 + ((0.2e1 * pkin(3) + t118) * t56 - t73) * t52) * qJD(4), (-t52 * t47 + t56 * t80) * pkin(3) + t1, (t115 * t8 - t52 * t9 + (-t115 * t24 + t25 * t52) * qJD(4)) * pkin(3), (t115 * t2 - t1 * t52 + (t115 * t6 - t5 * t52) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, 0, -t106, 0, -pkin(6) * t105, pkin(6) * t106, 0, 0, 0, 0, -t15, 0, -t16, 0, t11, t10, (t115 * t15 - t16 * t52 + (-t115 * t33 + t34 * t52) * qJD(4)) * pkin(3), (t115 * t11 - t10 * t52 + (t115 * t19 - t18 * t52) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t95, -0.2e1 * t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, -t9, t47, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, -t16, 0, t11, t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
