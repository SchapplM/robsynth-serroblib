% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:58
% EndTime: 2019-12-05 17:08:00
% DurationCPUTime: 0.28s
% Computational Cost: add. (239->61), mult. (635->93), div. (0->0), fcn. (508->6), ass. (0->59)
t71 = qJD(4) + qJD(5);
t70 = -pkin(7) - pkin(8);
t46 = cos(qJ(3));
t69 = t46 * pkin(2);
t43 = sin(qJ(3));
t36 = t43 * pkin(2) + pkin(7);
t68 = -pkin(8) - t36;
t41 = sin(qJ(5));
t42 = sin(qJ(4));
t44 = cos(qJ(5));
t45 = cos(qJ(4));
t23 = t41 * t45 + t44 * t42;
t12 = t71 * t23;
t64 = t41 * t42;
t22 = -t44 * t45 + t64;
t59 = t42 * qJD(4);
t54 = pkin(4) * t59;
t61 = pkin(2) * qJD(3);
t56 = t43 * t61;
t26 = t54 + t56;
t38 = -t45 * pkin(4) - pkin(3);
t29 = t38 - t69;
t67 = t29 * t12 + t26 * t22;
t39 = t45 * qJD(4);
t60 = qJD(5) * t44;
t11 = -t44 * t39 - t45 * t60 + t71 * t64;
t66 = -t29 * t11 + t26 * t23;
t65 = -t38 * t11 + t23 * t54;
t63 = t38 * t12 + t22 * t54;
t37 = -pkin(3) - t69;
t62 = t37 * t39 + t42 * t56;
t58 = pkin(3) * t59;
t57 = pkin(3) * t39;
t55 = t46 * t61;
t53 = qJD(5) * t41 * pkin(4);
t52 = pkin(4) * t60;
t51 = qJD(4) * t70;
t50 = qJD(4) * t68;
t49 = t42 * t55;
t48 = t45 * t55;
t47 = t37 * t59 - t45 * t56;
t40 = t45 * pkin(8);
t33 = 0.2e1 * t42 * t39;
t31 = t45 * pkin(7) + t40;
t30 = t70 * t42;
t25 = t45 * t51;
t24 = t42 * t51;
t21 = 0.2e1 * (-t42 ^ 2 + t45 ^ 2) * qJD(4);
t20 = t45 * t36 + t40;
t19 = t68 * t42;
t16 = t45 * t50 - t49;
t15 = t42 * t50 + t48;
t6 = -0.2e1 * t23 * t11;
t5 = -t41 * t24 + t44 * t25 + (-t30 * t41 - t31 * t44) * qJD(5);
t4 = -t44 * t24 - t41 * t25 + (-t30 * t44 + t31 * t41) * qJD(5);
t3 = 0.2e1 * t11 * t22 - 0.2e1 * t23 * t12;
t2 = -t41 * t15 + t44 * t16 + (-t19 * t41 - t20 * t44) * qJD(5);
t1 = -t44 * t15 - t41 * t16 + (-t19 * t44 + t20 * t41) * qJD(5);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -0.2e1 * t56, -0.2e1 * t55, t33, t21, 0, 0, 0, 0.2e1 * t47, 0.2e1 * t62, t6, t3, 0, 0, 0, 0.2e1 * t67, 0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t56, -t55, t33, t21, 0, 0, 0, t47 - t58, -t57 + t62, t6, t3, 0, 0, 0, t63 + t67, t65 + t66; 0, 0, 0, 0, 0, 0, 0, t33, t21, 0, 0, 0, -0.2e1 * t58, -0.2e1 * t57, t6, t3, 0, 0, 0, 0.2e1 * t63, 0.2e1 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t39, 0, 0, 0, 0, 0, -t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t59, 0, -t36 * t39 - t49, t36 * t59 - t48, 0, 0, -t11, -t12, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t59, 0, -pkin(7) * t39, pkin(7) * t59, 0, 0, -t11, -t12, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t53, -0.2e1 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
