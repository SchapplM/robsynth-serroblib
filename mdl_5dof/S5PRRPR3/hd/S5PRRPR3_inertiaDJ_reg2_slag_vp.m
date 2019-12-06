% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:41
% EndTime: 2019-12-05 16:19:44
% DurationCPUTime: 0.63s
% Computational Cost: add. (900->90), mult. (2105->177), div. (0->0), fcn. (1978->6), ass. (0->56)
t45 = sin(qJ(5));
t70 = cos(pkin(9));
t58 = t70 * pkin(3) + pkin(4);
t69 = sin(pkin(9));
t63 = t69 * pkin(3);
t74 = cos(qJ(5));
t29 = -t45 * t63 + t74 * t58;
t46 = sin(qJ(3));
t47 = cos(qJ(3));
t33 = t70 * t46 + t69 * t47;
t53 = t69 * t46 - t70 * t47;
t51 = t74 * t53;
t15 = t45 * t33 + t51;
t16 = t74 * t33 - t45 * t53;
t60 = qJD(3) * t69;
t61 = qJD(3) * t70;
t32 = -t46 * t60 + t47 * t61;
t71 = t46 * t61 + t47 * t60;
t8 = t16 * qJD(5) + t45 * t32 + t74 * t71;
t76 = t15 * t8;
t68 = qJD(5) * t45;
t7 = qJD(5) * t51 - t74 * t32 + t33 * t68 + t45 * t71;
t75 = t16 * t7;
t73 = t33 * t32;
t72 = -qJ(4) - pkin(6);
t38 = t72 * t46;
t39 = t72 * t47;
t20 = t69 * t38 - t70 * t39;
t67 = t46 * qJD(3);
t66 = t47 * qJD(3);
t65 = -0.2e1 * pkin(2) * qJD(3);
t44 = pkin(3) * t67;
t64 = t46 * t66;
t43 = -t47 * pkin(3) - pkin(2);
t62 = qJD(3) * t72;
t19 = t70 * t38 + t69 * t39;
t57 = t33 * pkin(7) - t19;
t55 = -t46 * qJD(4) + t47 * t62;
t54 = t47 * qJD(4) + t46 * t62;
t52 = t74 * t57;
t50 = t53 * t71;
t30 = t45 * t58 + t74 * t63;
t14 = -t53 * pkin(7) + t20;
t6 = t74 * t14 - t45 * t57;
t13 = t70 * t54 + t69 * t55;
t12 = -t69 * t54 + t70 * t55;
t49 = -t32 * pkin(7) + t12;
t48 = -t71 * pkin(7) + t13;
t24 = t30 * qJD(5);
t23 = t29 * qJD(5);
t22 = t53 * pkin(4) + t43;
t21 = t71 * pkin(4) + t44;
t5 = -t45 * t14 - t52;
t2 = -t6 * qJD(5) - t45 * t48 + t74 * t49;
t1 = qJD(5) * t52 + t14 * t68 - t45 * t49 - t74 * t48;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t50 + 0.2e1 * t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t75 + 0.2e1 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t12 + t33 * t13 - t71 * t19 + t32 * t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t1 - t15 * t2 - t8 * t5 - t7 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t64, 0.2e1 * (-t46 ^ 2 + t47 ^ 2) * qJD(3), 0, -0.2e1 * t64, 0, 0, t46 * t65, t47 * t65, 0, 0, 0.2e1 * t73, -0.2e1 * t32 * t53 - 0.2e1 * t33 * t71, 0, 0.2e1 * t50, 0, 0, 0.2e1 * t43 * t71 + 0.2e1 * t53 * t44, 0.2e1 * t43 * t32 + 0.2e1 * t33 * t44, -0.2e1 * t12 * t33 - 0.2e1 * t13 * t53 - 0.2e1 * t19 * t32 - 0.2e1 * t20 * t71, 0.2e1 * t19 * t12 + 0.2e1 * t20 * t13 + 0.2e1 * t43 * t44, -0.2e1 * t75, 0.2e1 * t15 * t7 - 0.2e1 * t16 * t8, 0, 0.2e1 * t76, 0, 0, 0.2e1 * t21 * t15 + 0.2e1 * t22 * t8, 0.2e1 * t21 * t16 - 0.2e1 * t22 * t7, 0.2e1 * t1 * t15 - 0.2e1 * t2 * t16 + 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * t22 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t66, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t32, 0, (t32 * t69 - t71 * t70) * pkin(3), 0, 0, 0, 0, 0, 0, -t8, t7, 0, t15 * t24 + t16 * t23 - t8 * t29 - t7 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, -t67, 0, -pkin(6) * t66, pkin(6) * t67, 0, 0, 0, 0, t32, 0, -t71, 0, t12, -t13, (-t70 * t32 - t69 * t71) * pkin(3), (t12 * t70 + t13 * t69) * pkin(3), 0, 0, -t7, 0, -t8, 0, t2, t1, -t23 * t15 + t24 * t16 + t29 * t7 - t30 * t8, -t1 * t30 + t2 * t29 + t6 * t23 - t5 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t24, -0.2e1 * t23, 0, 0.2e1 * t30 * t23 - 0.2e1 * t29 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t32, 0, t44, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
