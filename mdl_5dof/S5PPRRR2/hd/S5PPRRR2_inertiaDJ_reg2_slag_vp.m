% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:49
% EndTime: 2019-12-05 15:14:52
% DurationCPUTime: 0.60s
% Computational Cost: add. (536->79), mult. (1461->157), div. (0->0), fcn. (1416->8), ass. (0->52)
t32 = sin(pkin(9));
t35 = sin(qJ(3));
t55 = cos(pkin(9));
t61 = cos(qJ(3));
t21 = t61 * t32 + t35 * t55;
t36 = cos(qJ(4));
t60 = cos(qJ(5));
t45 = t60 * t36;
t33 = sin(qJ(5));
t34 = sin(qJ(4));
t58 = t33 * t34;
t37 = t45 - t58;
t8 = t37 * t21;
t20 = t35 * t32 - t61 * t55;
t64 = qJD(4) + qJD(5);
t63 = -pkin(7) - pkin(6);
t23 = t33 * t36 + t60 * t34;
t43 = t60 * qJD(5);
t9 = -qJD(4) * t45 - t36 * t43 + t64 * t58;
t62 = t23 * t9;
t16 = t21 * qJD(3);
t11 = t20 * t16;
t10 = t64 * t23;
t59 = t37 * t10;
t15 = t20 * qJD(3);
t57 = t34 * t15;
t54 = qJD(5) * t33;
t53 = t34 * qJD(4);
t52 = t36 * qJD(4);
t51 = -0.2e1 * pkin(3) * qJD(4);
t50 = pkin(4) * t53;
t49 = pkin(4) * t54;
t48 = t33 * t63;
t47 = t20 * t53;
t46 = t34 * t52;
t30 = t34 ^ 2;
t31 = t36 ^ 2;
t44 = (-t30 - t31) * t15;
t42 = t34 * t48;
t41 = pkin(4) * t43;
t40 = t63 * t60;
t38 = t34 * t40;
t24 = t63 * t36;
t13 = -t60 * t24 + t42;
t29 = -pkin(4) * t36 - pkin(3);
t12 = t33 * t24 + t38;
t7 = t23 * t21;
t4 = -t13 * qJD(5) + (t36 * t40 - t42) * qJD(4);
t3 = -t24 * t54 - t64 * t38 - t48 * t52;
t2 = t23 * t15 - t64 * t8;
t1 = t10 * t21 + t15 * t45 - t33 * t57;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t15 * t21 + 0.2e1 * t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t21 * t44 + 0.2e1 * t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t1 * t8 - 0.2e1 * t2 * t7 + 0.2e1 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t23 + t10 * t7 + t2 * t37 - t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t59 - 0.2e1 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t36 + t47, t16 * t34 + t20 * t52, t44, -pkin(3) * t16 + pkin(6) * t44, 0, 0, 0, 0, 0, 0, t10 * t20 - t16 * t37, t16 * t23 - t20 * t9, -t1 * t37 - t10 * t8 - t2 * t23 - t7 * t9, pkin(4) * t47 - t1 * t13 + t12 * t2 + t16 * t29 - t3 * t8 - t4 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t12 - t13 * t9 - t23 * t3 + t37 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t46, 0.2e1 * (-t30 + t31) * qJD(4), 0, -0.2e1 * t46, 0, 0, t34 * t51, t36 * t51, 0, 0, -0.2e1 * t62, -0.2e1 * t23 * t10 - 0.2e1 * t37 * t9, 0, -0.2e1 * t59, 0, 0, 0.2e1 * t10 * t29 - 0.2e1 * t37 * t50, 0.2e1 * t23 * t50 - 0.2e1 * t29 * t9, -0.2e1 * t10 * t13 + 0.2e1 * t12 * t9 - 0.2e1 * t23 * t4 - 0.2e1 * t3 * t37, 0.2e1 * t12 * t4 - 0.2e1 * t13 * t3 + 0.2e1 * t29 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t52 + t57, t15 * t36 + t21 * t53, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, (t60 * t2 - t1 * t33 + (t33 * t7 + t60 * t8) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, (-t60 * t10 - t33 * t9 + (t60 * t23 - t33 * t37) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, -t53, 0, -pkin(6) * t52, pkin(6) * t53, 0, 0, 0, 0, -t9, 0, -t10, 0, t4, t3, (t60 * t9 - t10 * t33 + (t23 * t33 + t37 * t60) * qJD(5)) * pkin(4), (t60 * t4 - t3 * t33 + (-t12 * t33 + t60 * t13) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t49, -0.2e1 * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, -t10, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
