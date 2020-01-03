% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPR8
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:42
% EndTime: 2019-12-31 17:42:43
% DurationCPUTime: 0.46s
% Computational Cost: add. (418->65), mult. (1152->117), div. (0->0), fcn. (1042->8), ass. (0->54)
t63 = qJD(2) + qJD(3);
t44 = cos(qJ(3));
t60 = cos(qJ(2));
t45 = t60 * qJD(2);
t48 = t44 * t60;
t41 = sin(qJ(3));
t42 = sin(qJ(2));
t58 = t41 * t42;
t10 = -qJD(3) * t48 - t44 * t45 + t63 * t58;
t24 = t41 * t60 + t44 * t42;
t11 = t63 * t24;
t38 = sin(pkin(9));
t39 = cos(pkin(9));
t4 = -t38 * t10 + t39 * t11;
t23 = t48 - t58;
t8 = -t39 * t23 + t38 * t24;
t62 = t8 * t4;
t55 = pkin(2) * qJD(3);
t59 = t39 * t41;
t17 = (t38 * t44 + t59) * t55;
t61 = t8 * t17;
t34 = t44 * pkin(2) + pkin(3);
t19 = -t38 * t41 * pkin(2) + t39 * t34;
t15 = -pkin(4) - t19;
t43 = cos(qJ(5));
t35 = t43 * qJD(5);
t40 = sin(qJ(5));
t57 = t15 * t35 + t17 * t40;
t20 = pkin(2) * t59 + t38 * t34;
t36 = t40 ^ 2;
t37 = t43 ^ 2;
t56 = t36 + t37;
t54 = t40 * qJD(5);
t53 = t41 * t55;
t52 = t44 * t55;
t33 = -t39 * pkin(3) - pkin(4);
t51 = t33 * t54;
t50 = t33 * t35;
t49 = t40 * t35;
t5 = -t39 * t10 - t38 * t11;
t1 = t56 * t5;
t18 = -t38 * t53 + t39 * t52;
t6 = t56 * t18;
t32 = t38 * pkin(3) + pkin(7);
t47 = t56 * t32;
t46 = t15 * t54 - t17 * t43;
t28 = -0.2e1 * t49;
t27 = 0.2e1 * t49;
t22 = 0.2e1 * (-t36 + t37) * qJD(5);
t16 = pkin(7) + t20;
t9 = t38 * t23 + t39 * t24;
t3 = -t4 * t43 + t8 * t54;
t2 = t8 * t35 + t4 * t40;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t24 * t10 - 0.2e1 * t23 * t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t9 * t5 + 0.2e1 * t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t9 * t1 + 0.2e1 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42 * qJD(2), -t45, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, (-t10 * t41 - t11 * t44 + (-t23 * t41 + t24 * t44) * qJD(3)) * pkin(2), 0, 0, 0, 0, 0, 0, -t4, -t5, 0, t9 * t18 - t4 * t19 + t5 * t20 + t61, 0, 0, 0, 0, 0, 0, t3, t2, t1, t16 * t1 + t4 * t15 + t9 * t6 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t53, -0.2e1 * t52, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17, -0.2e1 * t18, 0, -0.2e1 * t19 * t17 + 0.2e1 * t20 * t18, t27, t22, 0, t28, 0, 0, 0.2e1 * t46, 0.2e1 * t57, 0.2e1 * t6, 0.2e1 * t15 * t17 + 0.2e1 * t16 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t5, 0, (t38 * t5 - t39 * t4) * pkin(3), 0, 0, 0, 0, 0, 0, t3, t2, t1, t4 * t33 + t5 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, (-t17 * t39 + t18 * t38) * pkin(3), t27, t22, 0, t28, 0, 0, t46 + t51, t50 + t57, t6, t17 * t33 + t18 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t22, 0, t28, 0, 0, 0.2e1 * t51, 0.2e1 * t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t35 - t40 * t5, -t43 * t5 + t9 * t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, -t54, 0, -t16 * t35 - t40 * t18, t16 * t54 - t43 * t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, -t54, 0, -t32 * t35, t32 * t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
