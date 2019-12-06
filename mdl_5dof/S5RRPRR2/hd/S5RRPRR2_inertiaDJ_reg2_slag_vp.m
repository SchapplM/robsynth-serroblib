% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRR2
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
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:27
% EndTime: 2019-12-05 18:28:31
% DurationCPUTime: 1.24s
% Computational Cost: add. (3171->137), mult. (6850->264), div. (0->0), fcn. (6781->8), ass. (0->80)
t103 = -qJ(3) - pkin(6);
t75 = cos(qJ(2));
t108 = t103 * t75;
t74 = sin(qJ(2));
t63 = t103 * t74;
t70 = sin(pkin(9));
t71 = cos(pkin(9));
t109 = t108 * t71 - t70 * t63;
t39 = t108 * t70 + t71 * t63;
t106 = cos(qJ(4));
t107 = pkin(2) * t70;
t73 = sin(qJ(4));
t93 = pkin(2) * t71 + pkin(3);
t53 = t106 * t93 - t73 * t107;
t105 = cos(qJ(5));
t54 = t106 * t107 + t73 * t93;
t72 = sin(qJ(5));
t104 = t72 * t54;
t58 = t70 * t75 + t71 * t74;
t86 = -pkin(7) * t58 + t39;
t28 = t73 * t86;
t57 = t70 * t74 - t71 * t75;
t35 = -pkin(7) * t57 - t109;
t17 = t106 * t35 + t28;
t102 = qJD(4) * t73;
t101 = qJD(5) * t72;
t100 = t74 * qJD(2);
t99 = t75 * qJD(2);
t97 = -0.2e1 * pkin(1) * qJD(2);
t47 = t53 * qJD(4);
t48 = t54 * qJD(4);
t50 = pkin(4) + t53;
t88 = qJD(5) * t105;
t96 = -t105 * t47 + t48 * t72 - t50 * t88;
t69 = pkin(2) * t100;
t95 = pkin(4) * t101;
t94 = t74 * t99;
t68 = -t75 * pkin(2) - pkin(1);
t92 = t105 * t54;
t55 = t58 * qJD(2);
t41 = t55 * pkin(3) + t69;
t29 = t106 * t86;
t16 = -t35 * t73 + t29;
t90 = -t105 * t48 - t72 * t47;
t89 = qJD(4) * t106;
t87 = pkin(4) * t88;
t46 = t57 * pkin(3) + t68;
t38 = t106 * t58 - t57 * t73;
t85 = -pkin(8) * t38 + t16;
t84 = t106 * t57 + t58 * t73;
t32 = t50 * t72 + t92;
t33 = qJD(2) * t109 - t58 * qJD(3);
t56 = -t100 * t70 + t71 * t99;
t77 = -t56 * pkin(7) + t33;
t34 = qJD(2) * t39 - t57 * qJD(3);
t78 = -t55 * pkin(7) + t34;
t9 = -qJD(4) * t29 + t102 * t35 - t106 * t78 - t73 * t77;
t83 = -t102 * t57 + t106 * t55 + t56 * t73 + t58 * t89;
t82 = t72 * t85;
t81 = t105 * t85;
t80 = t105 * t84;
t23 = t105 * t38 - t72 * t84;
t79 = pkin(8) * t83 + t9;
t10 = -qJD(4) * t28 + t106 * t77 - t35 * t89 - t73 * t78;
t21 = qJD(4) * t84 - t106 * t56 + t73 * t55;
t76 = t21 * pkin(8) + t10;
t31 = t105 * t50 - t104;
t26 = pkin(4) * t84 + t46;
t22 = t38 * t72 + t80;
t15 = -qJD(5) * t32 + t90;
t14 = t101 * t54 + t96;
t13 = pkin(4) * t83 + t41;
t12 = -pkin(8) * t84 + t17;
t8 = qJD(5) * t23 + t105 * t83 - t72 * t21;
t7 = qJD(5) * t80 + t101 * t38 + t105 * t21 + t72 * t83;
t6 = t105 * t12 + t82;
t5 = -t12 * t72 + t81;
t2 = -qJD(5) * t82 + t105 * t76 - t12 * t88 + t72 * t79;
t1 = -qJD(5) * t81 + t101 * t12 + t105 * t79 - t72 * t76;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t94, 0.2e1 * (-t74 ^ 2 + t75 ^ 2) * qJD(2), 0, -0.2e1 * t94, 0, 0, t74 * t97, t75 * t97, 0, 0, 0.2e1 * t58 * t56, -0.2e1 * t55 * t58 - 0.2e1 * t56 * t57, 0, 0.2e1 * t57 * t55, 0, 0, 0.2e1 * t55 * t68 + 0.2e1 * t57 * t69, 0.2e1 * t56 * t68 + 0.2e1 * t58 * t69, 0.2e1 * t109 * t55 - 0.2e1 * t33 * t58 - 0.2e1 * t34 * t57 - 0.2e1 * t39 * t56, -0.2e1 * t109 * t34 + 0.2e1 * t33 * t39 + 0.2e1 * t68 * t69, -0.2e1 * t38 * t21, 0.2e1 * t21 * t84 - 0.2e1 * t38 * t83, 0, 0.2e1 * t84 * t83, 0, 0, 0.2e1 * t41 * t84 + 0.2e1 * t46 * t83, -0.2e1 * t21 * t46 + 0.2e1 * t38 * t41, -0.2e1 * t10 * t38 + 0.2e1 * t16 * t21 - 0.2e1 * t17 * t83 + 0.2e1 * t84 * t9, 0.2e1 * t10 * t16 - 0.2e1 * t17 * t9 + 0.2e1 * t41 * t46, -0.2e1 * t23 * t7, 0.2e1 * t22 * t7 - 0.2e1 * t23 * t8, 0, 0.2e1 * t22 * t8, 0, 0, 0.2e1 * t13 * t22 + 0.2e1 * t26 * t8, 0.2e1 * t13 * t23 - 0.2e1 * t26 * t7, 0.2e1 * t1 * t22 - 0.2e1 * t2 * t23 + 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t1 * t6 + 0.2e1 * t13 * t26 + 0.2e1 * t2 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, -t100, 0, -pkin(6) * t99, pkin(6) * t100, 0, 0, 0, 0, t56, 0, -t55, 0, t33, -t34, (-t55 * t70 - t56 * t71) * pkin(2), (t33 * t71 + t34 * t70) * pkin(2), 0, 0, -t21, 0, -t83, 0, t10, t9, t21 * t53 + t38 * t48 - t47 * t84 - t54 * t83, t10 * t53 - t16 * t48 + t17 * t47 - t54 * t9, 0, 0, -t7, 0, -t8, 0, t2, t1, t14 * t22 - t15 * t23 + t31 * t7 - t32 * t8, -t1 * t32 - t14 * t6 + t15 * t5 + t2 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t48, -0.2e1 * t47, 0, 0.2e1 * t47 * t54 - 0.2e1 * t48 * t53, 0, 0, 0, 0, 0, 0, 0.2e1 * t15, 0.2e1 * t14, 0, -0.2e1 * t14 * t32 + 0.2e1 * t15 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t56, 0, t69, 0, 0, 0, 0, 0, 0, t83, -t21, 0, t41, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, -t83, 0, t10, t9, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, (t105 * t7 - t72 * t8 + (-t105 * t22 + t23 * t72) * qJD(5)) * pkin(4), (t105 * t2 - t1 * t72 + (t105 * t6 - t5 * t72) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, 0, 0, 0, 0, 0, 0, 0, 0, (-t92 + (-pkin(4) - t50) * t72) * qJD(5) + t90, (-pkin(4) * t105 + t104) * qJD(5) + t96, 0, (t105 * t15 - t14 * t72 + (t105 * t32 - t31 * t72) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t95, -0.2e1 * t87, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t87, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
