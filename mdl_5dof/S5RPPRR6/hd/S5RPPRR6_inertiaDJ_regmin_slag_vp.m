% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:07
% EndTime: 2019-12-31 17:58:09
% DurationCPUTime: 0.31s
% Computational Cost: add. (342->66), mult. (836->128), div. (0->0), fcn. (783->8), ass. (0->55)
t38 = cos(qJ(5));
t32 = t38 ^ 2;
t36 = sin(qJ(5));
t55 = t36 ^ 2 - t32;
t46 = t55 * qJD(5);
t33 = sin(pkin(9));
t34 = cos(pkin(9));
t37 = sin(qJ(4));
t62 = cos(qJ(4));
t23 = t62 * t33 + t37 * t34;
t18 = t23 * qJD(4);
t66 = 0.2e1 * t18;
t27 = sin(pkin(8)) * pkin(1) + qJ(3);
t63 = pkin(6) + t27;
t19 = t63 * t33;
t20 = t63 * t34;
t11 = -t37 * t19 + t62 * t20;
t4 = t23 * qJD(3) + t11 * qJD(4);
t65 = t4 * t36;
t64 = t4 * t38;
t47 = qJD(4) * t62;
t59 = t37 * t33;
t17 = qJD(4) * t59 - t34 * t47;
t61 = t23 * t17;
t60 = t36 * t18;
t58 = t38 * t17;
t57 = t38 * t18;
t48 = t62 * t34;
t22 = -t48 + t59;
t56 = -t22 * t58 + t23 * t57;
t54 = qJD(5) * t36;
t53 = qJD(5) * t38;
t52 = -0.2e1 * pkin(4) * qJD(5);
t51 = t36 * t58;
t50 = t23 * t54;
t49 = t36 * t53;
t45 = 0.2e1 * (t33 ^ 2 + t34 ^ 2) * qJD(3);
t44 = pkin(4) * t17 - pkin(7) * t18;
t43 = pkin(4) * t23 + pkin(7) * t22;
t24 = -cos(pkin(8)) * pkin(1) - pkin(2) - t34 * pkin(3);
t5 = t22 * pkin(4) - t23 * pkin(7) + t24;
t42 = t38 * t11 + t36 * t5;
t41 = t36 * t11 - t38 * t5;
t40 = t17 * t22 - t23 * t18;
t39 = -t36 * t17 + t23 * t53;
t7 = t50 + t58;
t8 = t22 * t53 + t60;
t21 = t23 ^ 2;
t12 = t18 * pkin(4) + t17 * pkin(7);
t10 = t62 * t19 + t37 * t20;
t6 = t22 * t54 - t57;
t3 = t19 * t47 - qJD(3) * t48 + (qJD(3) * t33 + qJD(4) * t20) * t37;
t2 = -t42 * qJD(5) + t38 * t12 + t36 * t3;
t1 = t41 * qJD(5) - t36 * t12 + t38 * t3;
t9 = [0, 0, 0, 0, 0, 0, t45, t27 * t45, -0.2e1 * t61, 0.2e1 * t40, 0, 0, 0, t24 * t66, -0.2e1 * t24 * t17, -0.2e1 * t21 * t49 - 0.2e1 * t32 * t61, 0.2e1 * t21 * t46 + 0.4e1 * t23 * t51, -0.2e1 * t22 * t50 + 0.2e1 * t56, -0.2e1 * t39 * t22 - 0.2e1 * t23 * t60, t22 * t66, 0.2e1 * t39 * t10 - 0.2e1 * t41 * t18 + 0.2e1 * t2 * t22 + 0.2e1 * t23 * t65, 0.2e1 * t1 * t22 - 0.2e1 * t7 * t10 - 0.2e1 * t42 * t18 + 0.2e1 * t23 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t38 + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, 0, 0, 0, 0, -t6, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, -t4, t3, -t23 * t46 - t51, t55 * t17 - 0.4e1 * t23 * t49, t8, -t6, 0, -t64 + t44 * t36 + (t10 * t36 - t43 * t38) * qJD(5), t65 + t44 * t38 + (t10 * t38 + t43 * t36) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0, 0, 0, 0, 0, t6, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t49, -0.2e1 * t46, 0, 0, 0, t36 * t52, t38 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t39, t18, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t54, 0, -pkin(7) * t53, pkin(7) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
