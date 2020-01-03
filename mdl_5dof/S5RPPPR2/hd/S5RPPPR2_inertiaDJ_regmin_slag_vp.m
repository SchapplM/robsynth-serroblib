% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:55
% EndTime: 2020-01-03 11:22:57
% DurationCPUTime: 0.31s
% Computational Cost: add. (280->63), mult. (752->150), div. (0->0), fcn. (739->8), ass. (0->56)
t33 = sin(pkin(9));
t36 = cos(pkin(9));
t38 = cos(pkin(7));
t35 = sin(pkin(7));
t37 = cos(pkin(8));
t58 = t35 * t37;
t21 = -t38 * t33 + t36 * t58;
t39 = sin(qJ(5));
t34 = sin(pkin(8));
t40 = cos(qJ(5));
t59 = t34 * t40;
t41 = -t21 * t39 + t35 * t59;
t8 = t41 * qJD(5);
t64 = 0.2e1 * t8;
t53 = qJD(3) * t35;
t54 = qJD(2) * t38;
t22 = t34 * t54 + t37 * t53;
t63 = t22 * t37;
t62 = t33 * t34;
t61 = t34 * t35;
t60 = t34 * t39;
t25 = -t38 * pkin(2) - t35 * qJ(3) - pkin(1);
t55 = qJ(2) * t38;
t56 = t34 * t25 + t37 * t55;
t13 = -t38 * qJ(4) + t56;
t16 = (pkin(3) * t34 - qJ(4) * t37 + qJ(2)) * t35;
t57 = t36 * t13 + t33 * t16;
t52 = qJD(5) * t39;
t51 = qJD(5) * t40;
t31 = t35 ^ 2;
t50 = t31 * qJD(2);
t49 = qJ(2) * qJD(2);
t48 = t33 * t52;
t47 = t33 * t51;
t46 = t37 * t25 - t34 * t55;
t45 = t38 * pkin(3) - t46;
t23 = -t34 * t53 + t37 * t54;
t44 = -t33 * t13 + t36 * t16;
t43 = -t23 * t34 + t63;
t42 = (-qJD(4) * t37 + qJD(2)) * t35;
t12 = t21 * t40 + t35 * t60;
t32 = t38 ^ 2;
t29 = t31 * t49;
t20 = t33 * t58 + t38 * t36;
t18 = (-t36 * t59 + t37 * t39) * qJD(5);
t17 = (t36 * t60 + t37 * t40) * qJD(5);
t15 = -t38 * qJD(4) + t23;
t9 = t12 * qJD(5);
t7 = t36 * t15 + t33 * t42;
t6 = t33 * t15 - t36 * t42;
t5 = pkin(6) * t61 + t57;
t4 = -pkin(4) * t61 - t44;
t3 = t20 * pkin(4) - t21 * pkin(6) + t45;
t2 = t40 * t22 - t39 * t7 + (-t3 * t39 - t40 * t5) * qJD(5);
t1 = -t39 * t22 - t40 * t7 + (-t3 * t40 + t39 * t5) * qJD(5);
t10 = [0, 0, 0, 0, 0, 0.2e1 * (t31 + t32) * qJD(2), 0.2e1 * t32 * t49 + 0.2e1 * t29, 0.2e1 * t22 * t38 + 0.2e1 * t34 * t50, 0.2e1 * t23 * t38 + 0.2e1 * t37 * t50, 0.2e1 * t43 * t35, -0.2e1 * t46 * t22 + 0.2e1 * t56 * t23 + 0.2e1 * t29, 0.2e1 * t22 * t20 - 0.2e1 * t6 * t61, 0.2e1 * t22 * t21 - 0.2e1 * t7 * t61, -0.2e1 * t7 * t20 + 0.2e1 * t6 * t21, 0.2e1 * t45 * t22 - 0.2e1 * t44 * t6 + 0.2e1 * t57 * t7, t12 * t64, -0.2e1 * t12 * t9 + 0.2e1 * t41 * t8, t20 * t64, -0.2e1 * t9 * t20, 0, 0.2e1 * t2 * t20 + 0.2e1 * t4 * t9 - 0.2e1 * t41 * t6, 0.2e1 * t1 * t20 + 0.2e1 * t6 * t12 + 0.2e1 * t4 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, 0, 0, -t63 + (t33 * t6 + t36 * t7) * t34, 0, 0, 0, 0, 0, t18 * t20 + t9 * t62, t17 * t20 + t8 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * qJD(2), 0, 0, 0, t7 * t33 - t6 * t36, 0, 0, 0, 0, 0, -t20 * t47 - t36 * t9, t20 * t48 - t36 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, -t20 * t52, -t20 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
