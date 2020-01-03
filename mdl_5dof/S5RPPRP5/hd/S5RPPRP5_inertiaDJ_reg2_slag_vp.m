% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:43
% EndTime: 2019-12-31 17:53:44
% DurationCPUTime: 0.39s
% Computational Cost: add. (286->60), mult. (712->104), div. (0->0), fcn. (614->4), ass. (0->38)
t36 = sin(pkin(7));
t34 = t36 ^ 2;
t37 = cos(pkin(7));
t58 = (t37 ^ 2 + t34) * qJD(2);
t57 = 2 * qJD(5);
t56 = -pkin(6) + qJ(2);
t55 = qJ(2) * t58;
t54 = t36 * qJD(2);
t53 = t36 * qJD(3);
t52 = t37 * qJD(2);
t38 = sin(qJ(4));
t51 = t38 * qJD(4);
t39 = cos(qJ(4));
t50 = t39 * qJD(4);
t16 = t36 * t50 - t37 * t51;
t21 = t36 * t38 + t37 * t39;
t48 = 0.2e1 * t21 * t16;
t47 = -t37 * pkin(2) - t36 * qJ(3) - pkin(1);
t23 = t56 * t37;
t45 = t56 * t36;
t44 = t39 * t45;
t11 = t38 * t23 - t44;
t12 = t39 * t23 + t38 * t45;
t5 = -qJD(4) * t44 + t23 * t51 - t38 * t54 - t39 * t52;
t6 = t12 * qJD(4) + t38 * t52 - t39 * t54;
t46 = t11 * t6 - t12 * t5;
t17 = t37 * pkin(3) - t47;
t15 = t21 * qJD(4);
t22 = t36 * t39 - t37 * t38;
t43 = t15 * t21 - t22 * t16;
t42 = t11 * t51 + t12 * t50 - t5 * t38 - t6 * t39;
t41 = t39 * t15 - t38 * t16 - t21 * t50 + t22 * t51;
t40 = -0.2e1 * t11 * t15 - 0.2e1 * t12 * t16 + 0.2e1 * t5 * t21 + 0.2e1 * t6 * t22;
t20 = 0.2e1 * t58;
t10 = -0.2e1 * t22 * t15;
t7 = t21 * pkin(4) - t22 * qJ(5) + t17;
t3 = t16 * pkin(4) + t15 * qJ(5) - t22 * qJD(5) + t53;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0.2e1 * t55, 0, 0, 0, 0, 0, 0, 0.2e1 * t37 * t53, t20, 0.2e1 * t34 * qJD(3), -0.2e1 * t47 * t53 + 0.2e1 * t55, t10, 0.2e1 * t43, 0, t48, 0, 0, 0.2e1 * t17 * t16 + 0.2e1 * t21 * t53, -0.2e1 * t17 * t15 + 0.2e1 * t22 * t53, t40, 0.2e1 * t17 * t53 + 0.2e1 * t46, t10, 0, -0.2e1 * t43, 0, 0, t48, 0.2e1 * t7 * t16 + 0.2e1 * t3 * t21, t40, 0.2e1 * t7 * t15 - 0.2e1 * t3 * t22, 0.2e1 * t7 * t3 + 0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, 0, 0, 0, 0, 0, 0, -t16, t15, 0, -t53, 0, 0, 0, 0, 0, 0, -t16, 0, -t15, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, t41, t42, 0, 0, 0, 0, 0, 0, 0, t41, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, -t16, 0, -t6, t5, 0, 0, 0, -t15, 0, 0, t16, 0, -t6, t15 * pkin(4) - t16 * qJ(5) - t21 * qJD(5), -t5, -t6 * pkin(4) - t5 * qJ(5) + t12 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50, 0, 0, 0, 0, 0, 0, 0, 0, -t51, 0, t50, t38 * qJD(5) + (-pkin(4) * t38 + qJ(5) * t39) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, qJ(5) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
