% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:16
% EndTime: 2019-12-05 15:11:18
% DurationCPUTime: 0.52s
% Computational Cost: add. (165->53), mult. (599->99), div. (0->0), fcn. (519->6), ass. (0->52)
t33 = cos(pkin(8));
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t32 = sin(pkin(8));
t37 = cos(qJ(3));
t63 = t32 * t37;
t12 = t33 * t36 + t34 * t63;
t13 = -t33 * t34 + t36 * t63;
t35 = sin(qJ(3));
t58 = t35 * qJD(3);
t48 = t32 * t58;
t7 = t12 * qJD(4) + t36 * t48;
t4 = t7 * t36;
t28 = t36 * qJD(4);
t50 = t37 * t28;
t6 = -t32 * t50 + (qJD(4) * t33 + t48) * t34;
t40 = -t6 * t34 - t4 + (t12 * t36 - t13 * t34) * qJD(4);
t57 = t37 * qJD(3);
t29 = t34 ^ 2;
t30 = t36 ^ 2;
t61 = t29 + t30;
t15 = t61 * t57;
t16 = t35 * t28 + t34 * t57;
t64 = 2 * qJD(5);
t62 = pkin(6) * t15;
t60 = qJD(3) * t32;
t59 = t34 * qJD(4);
t56 = -0.2e1 * pkin(3) * qJD(4);
t54 = pkin(6) * t59;
t53 = pkin(6) * t28;
t52 = t35 * t57;
t51 = t36 * t57;
t49 = t34 * t28;
t46 = t32 * t57;
t45 = -pkin(4) * t36 - qJ(5) * t34;
t44 = pkin(4) * t34 - qJ(5) * t36;
t43 = 0.2e1 * t32 ^ 2 * t52 - 0.2e1 * t12 * t6 - 0.2e1 * t13 * t7;
t14 = t35 * t59 - t51;
t41 = t45 * qJD(4) + t36 * qJD(5);
t39 = t40 * pkin(6);
t38 = -t37 ^ 2 * t60 + t13 * t51 + t16 * t12 + (-t4 + (-qJD(4) * t13 - t6) * t34 + t60 * t35) * t35;
t26 = -0.2e1 * t49;
t25 = 0.2e1 * t49;
t22 = (-t29 + t30) * qJD(4);
t19 = -pkin(3) + t45;
t18 = -t36 * t58 - t37 * t59;
t17 = t34 * t58 - t50;
t11 = t44 * qJD(4) - t34 * qJD(5);
t10 = t14 * t32;
t9 = t16 * t32;
t5 = 0.2e1 * (-0.1e1 + t61) * t52;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t48, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, t40, -pkin(3) * t46 + t39, 0, 0, 0, 0, 0, 0, t10, t40, -t9, (t11 * t35 + t19 * t57) * t32 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t57, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, t15, -pkin(3) * t58 + t62, 0, 0, 0, 0, 0, 0, t18, t15, -t17, -t11 * t37 + t19 * t58 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0.2e1 * t22, 0, t26, 0, 0, t34 * t56, t36 * t56, 0, 0, t25, 0, -0.2e1 * t22, 0, 0, t26, -0.2e1 * t11 * t36 + 0.2e1 * t19 * t59, 0, -0.2e1 * t11 * t34 - 0.2e1 * t19 * t28, 0.2e1 * t19 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, -t7, pkin(4) * t6 - qJ(5) * t7 + qJD(5) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t14, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, -t14, t41 * t35 - t44 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t59, 0, -t53, t54, 0, 0, 0, t28, 0, 0, t59, 0, -t53, t41, -t54, t41 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, qJ(5) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
