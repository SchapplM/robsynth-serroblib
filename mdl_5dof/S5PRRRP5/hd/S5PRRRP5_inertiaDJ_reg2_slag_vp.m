% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:18
% EndTime: 2019-12-05 16:49:22
% DurationCPUTime: 0.78s
% Computational Cost: add. (721->101), mult. (1967->187), div. (0->0), fcn. (1728->6), ass. (0->70)
t50 = sin(qJ(4));
t51 = sin(qJ(3));
t80 = t50 * t51;
t85 = qJD(3) + qJD(4);
t86 = t85 * t80;
t52 = sin(qJ(2));
t53 = cos(qJ(3));
t81 = cos(qJ(4));
t61 = t81 * t53;
t55 = t61 - t80;
t28 = t55 * t52;
t84 = -pkin(7) - pkin(6);
t83 = pkin(3) * t50;
t59 = t81 * qJD(4);
t58 = pkin(3) * t59;
t35 = t50 * t53 + t81 * t51;
t21 = t85 * t35;
t54 = cos(qJ(2));
t73 = t54 * qJD(2);
t62 = t51 * t73;
t7 = t21 * t52 + t50 * t62 - t61 * t73;
t82 = t28 * t58 - t7 * t83;
t79 = -t21 * t83 + t55 * t58;
t40 = t84 * t53;
t69 = t50 * t84;
t23 = -t81 * t40 + t51 * t69;
t48 = t51 ^ 2;
t49 = t53 ^ 2;
t78 = t48 + t49;
t77 = qJD(4) * t50;
t76 = t51 * qJD(3);
t75 = t52 * qJD(2);
t74 = t53 * qJD(3);
t72 = -0.2e1 * pkin(2) * qJD(3);
t71 = pkin(3) * t76;
t70 = pkin(3) * t77;
t68 = t81 * pkin(3);
t27 = t35 * t52;
t67 = t27 * t77;
t66 = t35 * t77;
t65 = t51 * t74;
t64 = t52 * t73;
t63 = t54 * t76;
t47 = -t53 * pkin(3) - pkin(2);
t60 = t78 * t54;
t57 = t84 * t81;
t37 = t51 * t57;
t22 = t50 * t40 + t37;
t56 = qJD(3) * t57;
t9 = -qJD(4) * t37 - t40 * t77 - t51 * t56 - t69 * t74;
t10 = t40 * t59 + t53 * t56 - t84 * t86;
t46 = t68 + pkin(4);
t44 = -0.2e1 * t58;
t43 = -0.2e1 * t70;
t25 = -pkin(4) * t55 + t47;
t20 = -qJD(3) * t61 - t53 * t59 + t86;
t17 = t21 * pkin(4) + t71;
t16 = qJ(5) * t55 + t23;
t15 = -t35 * qJ(5) + t22;
t14 = -0.2e1 * t35 * t20;
t13 = -0.2e1 * t55 * t21;
t12 = t54 * t20 + t35 * t75;
t11 = -t54 * t21 - t55 * t75;
t8 = -t85 * t28 - t35 * t73;
t5 = -0.2e1 * t20 * t55 - 0.2e1 * t35 * t21;
t4 = t20 * qJ(5) - t35 * qJD(5) + t10;
t3 = t21 * qJ(5) - qJD(5) * t55 + t9;
t2 = -0.2e1 * t27 * t8 - 0.2e1 * t28 * t7 - 0.2e1 * t64;
t1 = -t27 * t20 - t28 * t21 - t8 * t35 - t55 * t7;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t78) * t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t73, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t75 - t63, t51 * t75 - t54 * t74, qJD(2) * t60, (-pkin(2) * t52 + pkin(6) * t60) * qJD(2), 0, 0, 0, 0, 0, 0, t11, t12, t1, -pkin(3) * t63 - t27 * t10 + t8 * t22 - t7 * t23 - t28 * t9 + t47 * t75, 0, 0, 0, 0, 0, 0, t11, t12, t1, t8 * t15 - t7 * t16 - t54 * t17 + t25 * t75 - t27 * t4 - t28 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t65, 0.2e1 * (-t48 + t49) * qJD(3), 0, -0.2e1 * t65, 0, 0, t51 * t72, t53 * t72, 0, 0, t14, t5, 0, t13, 0, 0, 0.2e1 * t47 * t21 - 0.2e1 * t55 * t71, -0.2e1 * t47 * t20 + 0.2e1 * t35 * t71, -0.2e1 * t10 * t35 + 0.2e1 * t22 * t20 - 0.2e1 * t23 * t21 - 0.2e1 * t55 * t9, 0.2e1 * t22 * t10 - 0.2e1 * t23 * t9 + 0.2e1 * t47 * t71, t14, t5, 0, t13, 0, 0, -0.2e1 * t17 * t55 + 0.2e1 * t25 * t21, 0.2e1 * t17 * t35 - 0.2e1 * t25 * t20, 0.2e1 * t15 * t20 - 0.2e1 * t16 * t21 - 0.2e1 * t3 * t55 - 0.2e1 * t4 * t35, 0.2e1 * t15 * t4 - 0.2e1 * t16 * t3 + 0.2e1 * t25 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * t74 - t62, t52 * t76 - t53 * t73, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, (t81 * t8 + t67) * pkin(3) + t82, 0, 0, 0, 0, 0, 0, t8, t7, 0, pkin(3) * t67 + t8 * t46 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, -t76, 0, -pkin(6) * t74, pkin(6) * t76, 0, 0, 0, 0, -t20, 0, -t21, 0, t10, t9, (t81 * t20 + t66) * pkin(3) + t79, (t81 * t10 - t50 * t9 + (-t22 * t50 + t81 * t23) * qJD(4)) * pkin(3), 0, 0, -t20, 0, -t21, 0, t4, t3, pkin(3) * t66 + t46 * t20 + t79, t4 * t46 + (-t3 * t50 + (-t15 * t50 + t81 * t16) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t44, 0, 0, 0, 0, 0, 0, 0, 0, t43, t44, 0, 0.2e1 * (t68 - t46) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, t8 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, -t21, 0, t10, t9, 0, 0, 0, 0, -t20, 0, -t21, 0, t4, t3, t20 * pkin(4), t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t58, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t58, 0, -pkin(4) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
