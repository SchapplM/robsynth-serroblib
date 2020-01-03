% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR10_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:55
% EndTime: 2019-12-31 17:11:57
% DurationCPUTime: 0.50s
% Computational Cost: add. (339->72), mult. (805->159), div. (0->0), fcn. (547->4), ass. (0->64)
t27 = sin(qJ(2));
t67 = pkin(2) + pkin(6);
t49 = t67 * t27;
t29 = cos(qJ(2));
t63 = qJ(3) * t29;
t71 = t49 - t63;
t62 = t27 * qJ(3);
t70 = t67 * t29 + t62;
t25 = t29 ^ 2;
t41 = qJD(2) * (t27 ^ 2 - t25);
t26 = sin(qJ(4));
t22 = t26 ^ 2;
t28 = cos(qJ(4));
t24 = t28 ^ 2;
t65 = t22 - t24;
t40 = t65 * qJD(4);
t66 = pkin(3) + pkin(5);
t17 = t66 * t29;
t54 = t29 * qJD(3);
t69 = t70 * qJD(2) - qJD(4) * t17 - t54;
t68 = 0.2e1 * qJD(3);
t61 = qJD(2) * t17;
t60 = qJD(4) * t26;
t59 = qJD(4) * t28;
t58 = qJD(4) * t29;
t57 = qJD(4) * t67;
t56 = t27 * qJD(2);
t55 = t27 * qJD(3);
t20 = t29 * qJD(2);
t53 = qJ(3) * qJD(4);
t52 = -0.2e1 * pkin(1) * qJD(2);
t51 = pkin(5) * t56;
t50 = pkin(5) * t20;
t48 = t66 * t27;
t47 = t26 * t58;
t46 = t28 * t58;
t45 = t26 * t20;
t44 = t27 * t20;
t43 = t28 * t56;
t42 = t26 * t59;
t39 = t28 * t48;
t38 = t26 * t43;
t37 = t25 * t42;
t12 = -pkin(1) - t70;
t4 = -t26 * t12 + t39;
t5 = t28 * t12 + t26 * t48;
t36 = -t26 * t4 + t28 * t5;
t35 = -t29 * pkin(2) - t62;
t15 = t66 * t56;
t32 = t71 * qJD(4) - t15;
t31 = t35 * qJD(2) + t54;
t2 = t12 * t60 - t66 * t45 - qJD(4) * t39 - t28 * (t71 * qJD(2) - t55);
t3 = -t12 * t59 + (-qJD(4) * t66 + qJD(3)) * t26 * t27 + (-t26 * t49 + (t26 * qJ(3) + t28 * t66) * t29) * qJD(2);
t1 = t36 * qJD(4) - t2 * t26 + t3 * t28;
t21 = qJ(3) * t68;
t19 = -0.2e1 * t44;
t18 = 0.2e1 * t44;
t16 = -pkin(1) + t35;
t14 = -0.2e1 * t41;
t10 = -t27 * t59 - t45;
t9 = t28 * t20 - t27 * t60;
t8 = -t55 + (pkin(2) * t27 - t63) * qJD(2);
t7 = -t29 * t40 - t38;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t14, 0, t19, 0, 0, t27 * t52, t29 * t52, 0, 0, 0, 0, 0, t18, t14, t19, 0, -0.2e1 * t16 * t56 + 0.2e1 * t8 * t29, -0.2e1 * t16 * t20 - 0.2e1 * t8 * t27, 0.2e1 * t16 * t8, -0.2e1 * t22 * t44 + 0.2e1 * t37, -0.2e1 * t25 * t40 - 0.4e1 * t29 * t38, 0.2e1 * t26 * t41 - 0.2e1 * t27 * t46, -0.2e1 * t24 * t44 - 0.2e1 * t37, 0.2e1 * t27 * t47 + 0.2e1 * t28 * t41, t18, 0.2e1 * (-t28 * t61 + t3) * t27 + 0.2e1 * (qJD(2) * t4 - t15 * t28 - t17 * t60) * t29, 0.2e1 * (t26 * t61 + t2) * t27 + 0.2e1 * (-qJD(2) * t5 + t15 * t26 - t17 * t59) * t29, 0.2e1 * t36 * t56 + 0.2e1 * (t2 * t28 + t26 * t3 + (t26 * t5 + t28 * t4) * qJD(4)) * t29, -0.2e1 * t17 * t15 - 0.2e1 * t5 * t2 + 0.2e1 * t4 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t56, 0, -t50, t51, 0, 0, 0, -t20, t56, 0, 0, 0, t31, t50, -t51, t31 * pkin(5), -t7, 0.4e1 * t29 * t42 - t65 * t56, t9, t7, t10, 0, t32 * t26 - t28 * t69, t69 * t26 + t32 * t28, -t1, -t15 * qJ(3) + t17 * qJD(3) - t1 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t21, -0.2e1 * t42, 0.2e1 * t40, 0, 0.2e1 * t42, 0, 0, 0.2e1 * qJD(3) * t26 + 0.2e1 * t28 * t53, 0.2e1 * qJD(3) * t28 - 0.2e1 * t26 * t53, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, t50, 0, 0, 0, 0, 0, 0, t9, t10, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t56 - t46, 0, t43 + t47, t20, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, -t59, 0, t26 * t57, t28 * t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
