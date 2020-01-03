% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:53
% EndTime: 2019-12-31 17:59:55
% DurationCPUTime: 0.59s
% Computational Cost: add. (445->81), mult. (895->163), div. (0->0), fcn. (644->6), ass. (0->65)
t24 = cos(qJ(5));
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t64 = pkin(4) * t25;
t32 = pkin(7) * t23 + t64;
t67 = t24 * t32;
t22 = sin(qJ(5));
t18 = t22 ^ 2;
t20 = t24 ^ 2;
t59 = t18 + t20;
t19 = t23 ^ 2;
t21 = t25 ^ 2;
t58 = t19 - t21;
t66 = qJD(4) * t58;
t60 = t18 - t20;
t35 = qJD(5) * t60;
t65 = 2 * qJD(3);
t63 = t23 * pkin(4);
t15 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(6);
t62 = t22 * t15;
t61 = t24 * t25;
t57 = t19 + t21;
t56 = qJD(5) * t22;
t55 = qJD(5) * t24;
t54 = qJD(5) * t25;
t53 = t23 * qJD(4);
t52 = t25 * qJD(4);
t51 = -0.2e1 * pkin(4) * qJD(5);
t50 = t23 * t62;
t49 = t25 * t62;
t48 = t24 * t23 * t15;
t16 = sin(pkin(8)) * pkin(1) + qJ(3);
t47 = t16 * t65;
t46 = qJD(5) * t15 * t21;
t45 = t23 * t56;
t44 = t22 * t54;
t43 = t24 * t54;
t42 = t22 * t55;
t41 = t15 * t53;
t40 = t24 * t53;
t39 = t23 * t52;
t38 = t15 * t52;
t37 = t59 * t23;
t36 = t59 * t25;
t34 = t22 * t40;
t33 = t21 * t42;
t31 = -t25 * pkin(7) + t63;
t28 = t16 + t31;
t27 = t24 * t28;
t3 = t27 - t50;
t4 = t22 * t28 + t48;
t30 = t22 * t4 + t24 * t3;
t29 = t22 * t3 - t24 * t4;
t1 = t15 * t45 - t24 * t38 - qJD(5) * t27 - t22 * (t32 * qJD(4) + qJD(3));
t2 = t24 * qJD(3) - qJD(5) * t4 + (-t49 + t67) * qJD(4);
t26 = -t30 * qJD(5) - t1 * t24 - t2 * t22;
t17 = qJD(4) * t19;
t14 = 0.2e1 * t39;
t12 = t22 * t53 - t43;
t11 = t22 * t52 + t23 * t55;
t10 = t40 + t44;
t9 = -t24 * t52 + t45;
t6 = t25 * t35 + t34;
t5 = (-0.1e1 + t59) * t39;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t47, -0.2e1 * t39, -0.2e1 * t21 * qJD(4) + 0.2e1 * t17, 0, t14, 0, 0, 0.2e1 * qJD(3) * t23 + 0.2e1 * t16 * t52, 0.2e1 * qJD(3) * t25 - 0.2e1 * t16 * t53, 0, t47, -0.2e1 * t20 * t39 - 0.2e1 * t33, 0.2e1 * t21 * t35 + 0.4e1 * t25 * t34, -0.2e1 * t23 * t44 - 0.2e1 * t24 * t66, -0.2e1 * t18 * t39 + 0.2e1 * t33, 0.2e1 * t22 * t66 - 0.2e1 * t23 * t43, t14, -0.2e1 * t24 * t46 + 0.2e1 * t2 * t23 + 0.2e1 * (t3 + 0.2e1 * t50) * t52, 0.2e1 * t22 * t46 + 0.2e1 * t1 * t23 + 0.2e1 * (-t4 + 0.2e1 * t48) * t52, 0.2e1 * t30 * t53 + 0.2e1 * (qJD(5) * t29 + t1 * t22 - t2 * t24) * t25, -0.2e1 * t15 ^ 2 * t39 - 0.2e1 * t4 * t1 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t25 + (t15 * t58 + t23 * t29) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 * t55, t57 * t56, 0, -t29 * t52 + (t26 - 0.2e1 * t38) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 + (-t59 * t58 - t21) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, 0, -t52, 0, -t41, -t38, 0, 0, -t6, -0.4e1 * t25 * t42 + t53 * t60, t11, t6, -t9, 0, (-t49 - t67) * qJD(5) + (t22 * t31 - t48) * qJD(4), (-t15 * t61 + t22 * t32) * qJD(5) + (-pkin(7) * t61 + (t24 * pkin(4) + t62) * t23) * qJD(4), t26, -pkin(4) * t41 + pkin(7) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t53, 0, 0, 0, 0, 0, 0, 0, 0, t9, t11, -qJD(4) * t37, (-pkin(7) * t37 - t64) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t12, qJD(4) * t36, (pkin(7) * t36 - t63) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t42, -0.2e1 * t35, 0, -0.2e1 * t42, 0, 0, t22 * t51, t24 * t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, t12, t52, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, -t56, 0, -pkin(7) * t55, pkin(7) * t56, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
