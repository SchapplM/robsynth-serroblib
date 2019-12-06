% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:35
% EndTime: 2019-12-05 15:51:37
% DurationCPUTime: 0.37s
% Computational Cost: add. (215->73), mult. (708->157), div. (0->0), fcn. (653->10), ass. (0->67)
t32 = cos(qJ(5));
t23 = t32 ^ 2;
t29 = sin(qJ(5));
t66 = t29 ^ 2 - t23;
t44 = t66 * qJD(5);
t25 = sin(pkin(10));
t19 = t25 * pkin(2) + pkin(7);
t68 = t19 * t29;
t33 = cos(qJ(4));
t67 = t32 * t33;
t30 = sin(qJ(4));
t22 = t30 ^ 2;
t65 = -t33 ^ 2 + t22;
t26 = sin(pkin(5));
t27 = cos(pkin(10));
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t12 = (t25 * t34 + t27 * t31) * t26;
t28 = cos(pkin(5));
t7 = t12 * t30 - t28 * t33;
t64 = qJD(4) * t7;
t63 = qJD(2) * t26;
t62 = qJD(4) * t32;
t61 = qJD(5) * t29;
t60 = qJD(5) * t32;
t59 = qJD(5) * t33;
t58 = t30 * qJD(4);
t57 = t33 * qJD(4);
t56 = -0.2e1 * pkin(4) * qJD(5);
t55 = t33 * t68;
t54 = t19 * t67;
t20 = -t27 * pkin(2) - pkin(3);
t53 = 0.2e1 * qJD(4) * t20;
t52 = t29 * t59;
t51 = t32 * t59;
t50 = t29 * t58;
t49 = t29 * t60;
t48 = t30 * t57;
t47 = t19 * t58;
t46 = t32 * t58;
t45 = t32 * t57;
t43 = t65 * qJD(4);
t42 = t30 * t45;
t41 = -t33 * pkin(4) - t30 * pkin(8);
t40 = pkin(4) * t30 - pkin(8) * t33;
t37 = t25 * t31 - t27 * t34;
t11 = t37 * t26;
t8 = t12 * t33 + t28 * t30;
t39 = t11 * t32 - t8 * t29;
t38 = t11 * t29 + t8 * t32;
t10 = t37 * t63;
t5 = t8 * qJD(4) - t10 * t30;
t36 = t5 * t29 + t7 * t60;
t35 = -t5 * t32 + t7 * t61;
t18 = t40 * qJD(4);
t17 = t20 + t41;
t16 = t50 - t51;
t15 = -t29 * t57 - t30 * t60;
t14 = t46 + t52;
t13 = t30 * t61 - t45;
t9 = qJD(2) * t12;
t6 = -t10 * t33 - t64;
t4 = t29 * t47 + t32 * t18 + (-t29 * t17 - t54) * qJD(5);
t3 = t19 * t46 - t29 * t18 + (-t32 * t17 + t55) * qJD(5);
t2 = t39 * qJD(5) + t9 * t29 + t6 * t32;
t1 = -t38 * qJD(5) - t6 * t29 + t9 * t32;
t21 = [0, 0, 0, 0, -0.2e1 * t12 * t10 + 0.2e1 * t11 * t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t31 * t63, -t34 * t63, (-t10 * t25 - t27 * t9) * pkin(2), 0, 0, 0, 0, 0, t11 * t58 - t9 * t33, t11 * t57 + t9 * t30, 0, 0, 0, 0, 0, (t29 * t64 - t1) * t33 + (t39 * qJD(4) + t36) * t30, (t7 * t62 + t2) * t33 + (-t38 * qJD(4) - t35) * t30; 0, 0, 0, 0, 0, 0.2e1 * t48, -0.2e1 * t43, 0, 0, 0, t30 * t53, t33 * t53, -0.2e1 * t22 * t49 + 0.2e1 * t23 * t48, 0.2e1 * t22 * t44 - 0.4e1 * t29 * t42, 0.2e1 * t30 * t52 + 0.2e1 * t65 * t62, -0.2e1 * t29 * t43 + 0.2e1 * t30 * t51, -0.2e1 * t48, 0.2e1 * t17 * t46 - 0.2e1 * t4 * t33 + 0.2e1 * (t22 * t60 + t29 * t48) * t19, -0.2e1 * t17 * t50 - 0.2e1 * t3 * t33 + 0.2e1 * (-t22 * t61 + t42) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, t35, t36; 0, 0, 0, 0, 0, 0, 0, t57, -t58, 0, -t19 * t57, t47, t29 * t45 - t30 * t44, -0.4e1 * t30 * t49 - t66 * t57, t16, t14, 0, (pkin(8) * t67 + (-pkin(4) * t32 + t68) * t30) * qJD(5) + (t41 * t29 - t54) * qJD(4), (t19 * t30 * t32 + t40 * t29) * qJD(5) + (t41 * t32 + t55) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t57, 0, 0, 0, 0, 0, -t14, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t49, -0.2e1 * t44, 0, 0, 0, t29 * t56, t32 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t15, t58, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t61, 0, -pkin(8) * t60, pkin(8) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t21;
