% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:23
% DurationCPUTime: 0.43s
% Computational Cost: add. (357->70), mult. (840->146), div. (0->0), fcn. (782->6), ass. (0->55)
t36 = cos(qJ(5));
t31 = t36 ^ 2;
t35 = sin(qJ(5));
t55 = t35 ^ 2 - t31;
t45 = qJD(5) * t55;
t32 = sin(pkin(8));
t33 = cos(pkin(8));
t21 = (t32 ^ 2 + t33 ^ 2) * qJD(3);
t62 = sin(qJ(4));
t63 = cos(qJ(4));
t37 = -t62 * t32 + t63 * t33;
t15 = t37 ^ 2;
t67 = 2 * qJD(2);
t34 = -pkin(1) - qJ(3);
t64 = -pkin(6) + t34;
t19 = t64 * t32;
t20 = t64 * t33;
t10 = t63 * t19 + t62 * t20;
t4 = t37 * qJD(3) + t10 * qJD(4);
t66 = t36 * t4;
t65 = t4 * t35;
t46 = qJD(4) * t62;
t47 = qJD(4) * t63;
t13 = -t32 * t47 - t33 * t46;
t61 = t13 * t36;
t60 = t37 * t13;
t14 = -t32 * t46 + t33 * t47;
t17 = t63 * t32 + t62 * t33;
t59 = t17 * t14;
t58 = t35 * t13;
t57 = t35 * t36;
t25 = t32 * pkin(3) + qJ(2);
t54 = qJD(5) * t35;
t53 = qJD(5) * t36;
t52 = qJ(2) * qJD(2);
t51 = -0.2e1 * pkin(4) * qJD(5);
t50 = t35 * t53;
t49 = t17 ^ 2 + t15;
t48 = -0.4e1 * t37 * t57;
t44 = -pkin(4) * t13 - pkin(7) * t14;
t43 = -pkin(4) * t37 - pkin(7) * t17;
t8 = pkin(4) * t17 - pkin(7) * t37 + t25;
t42 = t10 * t36 + t35 * t8;
t41 = t35 * t10 - t36 * t8;
t40 = -t37 * t53 - t58;
t39 = -t37 * t54 + t61;
t6 = t35 * t14 + t17 * t53;
t5 = -t36 * t14 + t17 * t54;
t38 = -0.2e1 * t59 - 0.2e1 * t60;
t9 = t62 * t19 - t63 * t20;
t7 = pkin(4) * t14 - pkin(7) * t13 + qJD(2);
t3 = t17 * qJD(3) + t19 * t46 - t20 * t47;
t2 = -t42 * qJD(5) + t35 * t3 + t36 * t7;
t1 = t41 * qJD(5) + t3 * t36 - t35 * t7;
t11 = [0, 0, 0, 0, t67, 0.2e1 * t52, t32 * t67, t33 * t67, 0.2e1 * t21, -0.2e1 * t34 * t21 + 0.2e1 * t52, 0.2e1 * t60, -0.2e1 * t13 * t17 - 0.2e1 * t14 * t37, 0, 0, 0, 0.2e1 * qJD(2) * t17 + 0.2e1 * t14 * t25, 0.2e1 * qJD(2) * t37 + 0.2e1 * t13 * t25, -0.2e1 * t15 * t50 + 0.2e1 * t31 * t60, t13 * t48 + 0.2e1 * t15 * t45, 0.2e1 * t17 * t61 - 0.2e1 * t37 * t5, -0.2e1 * t17 * t58 - 0.2e1 * t37 * t6, 0.2e1 * t59, 0.2e1 * t2 * t17 - 0.2e1 * t41 * t14 + 0.2e1 * t9 * t58 - 0.2e1 * (-t9 * t53 - t65) * t37, 0.2e1 * t1 * t17 - 0.2e1 * t42 * t14 + 0.2e1 * t9 * t61 - 0.2e1 * (t9 * t54 - t66) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t38 - t49 * t53, t36 * t38 + t49 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, t14, t13, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, 0, -t4, t3, t13 * t57 - t37 * t45, qJD(5) * t48 - t55 * t13, t6, -t5, 0, -t66 + t44 * t35 + (t35 * t9 + t43 * t36) * qJD(5), t65 + t44 * t36 + (-t43 * t35 + t36 * t9) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, 0, 0, 0, 0, 0, t39, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t50, -0.2e1 * t45, 0, 0, 0, t35 * t51, t36 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t40, t14, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t54, 0, -pkin(7) * t53, pkin(7) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
