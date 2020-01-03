% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:13
% EndTime: 2019-12-31 18:09:15
% DurationCPUTime: 0.43s
% Computational Cost: add. (403->57), mult. (925->113), div. (0->0), fcn. (789->6), ass. (0->45)
t33 = sin(pkin(8));
t35 = sin(qJ(3));
t36 = cos(qJ(3));
t54 = cos(pkin(8));
t25 = t33 * t36 + t54 * t35;
t21 = t25 * qJD(3);
t46 = t54 * t36;
t52 = t35 * qJD(3);
t22 = qJD(3) * t46 - t33 * t52;
t24 = t33 * t35 - t46;
t56 = t25 * t21 + t22 * t24;
t61 = 0.2e1 * qJD(3);
t60 = 2 * qJD(5);
t58 = t24 * t21;
t53 = t25 * qJD(5);
t51 = t36 * qJD(3);
t50 = 0.2e1 * t58;
t31 = -cos(pkin(7)) * pkin(1) - pkin(2);
t49 = t31 * t61;
t32 = pkin(3) * t52;
t48 = t35 * t51;
t45 = sin(pkin(7)) * pkin(1) + pkin(6);
t43 = qJ(4) + t45;
t23 = t43 * t36;
t40 = t43 * t35;
t10 = t33 * t23 + t54 * t40;
t11 = t54 * t23 - t33 * t40;
t39 = qJD(3) * t43;
t19 = t36 * qJD(4) - t35 * t39;
t37 = -t35 * qJD(4) - t36 * t39;
t7 = t33 * t19 - t54 * t37;
t8 = t54 * t19 + t33 * t37;
t47 = t10 * t7 + t11 * t8;
t26 = -t36 * pkin(3) + t31;
t17 = t25 * t22;
t44 = 0.2e1 * t17 + 0.2e1 * t58;
t42 = t10 * t21 + t11 * t22 + t7 * t24 + t8 * t25;
t41 = t45 * qJD(3);
t38 = 0.2e1 * t10 * t22 - 0.2e1 * t11 * t21 - 0.2e1 * t8 * t24 + 0.2e1 * t7 * t25;
t30 = -t54 * pkin(3) - pkin(4);
t28 = t33 * pkin(3) + qJ(5);
t12 = 0.2e1 * t17;
t9 = t24 * pkin(4) - t25 * qJ(5) + t26;
t4 = t21 * pkin(4) - t22 * qJ(5) + t32 - t53;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t48, (-t35 ^ 2 + t36 ^ 2) * t61, 0, -0.2e1 * t48, 0, 0, t35 * t49, t36 * t49, 0, 0, t12, -0.2e1 * t56, 0, t50, 0, 0, 0.2e1 * t26 * t21 + 0.2e1 * t24 * t32, 0.2e1 * t26 * t22 + 0.2e1 * t25 * t32, t38, 0.2e1 * t26 * t32 + 0.2e1 * t47, t12, 0, 0.2e1 * t56, 0, 0, t50, 0.2e1 * t9 * t21 + 0.2e1 * t4 * t24, t38, -0.2e1 * t9 * t22 - 0.2e1 * t4 * t25, 0.2e1 * t9 * t4 + 0.2e1 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, -t52, 0, -t36 * t41, t35 * t41, 0, 0, 0, 0, t22, 0, -t21, 0, -t7, -t8, (-t21 * t33 - t54 * t22) * pkin(3), (t33 * t8 - t54 * t7) * pkin(3), 0, t22, 0, 0, t21, 0, -t7, -qJD(5) * t24 - t28 * t21 + t30 * t22, t8, t11 * qJD(5) + t8 * t28 + t7 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t51, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, 0, (-t54 * t21 + t22 * t33) * pkin(3), 0, 0, 0, 0, 0, 0, -t21, 0, t22, t21 * t30 + t22 * t28 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t28 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22, 0, t32, 0, 0, 0, 0, 0, 0, t21, 0, -t22, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
