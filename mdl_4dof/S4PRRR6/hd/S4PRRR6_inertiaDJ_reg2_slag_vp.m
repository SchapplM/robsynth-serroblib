% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:06
% EndTime: 2019-12-31 16:35:08
% DurationCPUTime: 0.46s
% Computational Cost: add. (293->61), mult. (861->132), div. (0->0), fcn. (735->6), ass. (0->45)
t53 = qJD(3) + qJD(4);
t24 = sin(qJ(4));
t27 = cos(qJ(4));
t25 = sin(qJ(3));
t52 = -pkin(6) - pkin(5);
t39 = t52 * t25;
t30 = t27 * t39;
t28 = cos(qJ(3));
t54 = t52 * t28;
t7 = t24 * t54 + t30;
t55 = -t24 * t39 + t27 * t54;
t51 = t24 * t25;
t14 = -t27 * t28 + t51;
t26 = sin(qJ(2));
t10 = t14 * t26;
t22 = t25 ^ 2;
t23 = t28 ^ 2;
t50 = t22 + t23;
t49 = qJD(4) * t24;
t48 = qJD(4) * t27;
t47 = t25 * qJD(3);
t46 = t26 * qJD(2);
t45 = t28 * qJD(3);
t29 = cos(qJ(2));
t44 = t29 * qJD(2);
t43 = -0.2e1 * pkin(2) * qJD(3);
t42 = pkin(3) * t47;
t41 = pkin(3) * t49;
t40 = pkin(3) * t48;
t37 = t25 * t45;
t36 = t26 * t44;
t35 = t29 * t47;
t34 = t25 * t44;
t33 = t28 * t44;
t32 = t50 * t29;
t15 = t24 * t28 + t27 * t25;
t6 = t53 * t15;
t21 = -t28 * pkin(3) - pkin(2);
t9 = t15 * t26;
t5 = -t27 * t45 - t28 * t48 + t53 * t51;
t4 = t53 * t55;
t3 = -qJD(3) * t7 - qJD(4) * t30 - t54 * t49;
t2 = t53 * t10 - t15 * t44;
t1 = t24 * t34 + t6 * t26 - t27 * t33;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t50) * t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t1 - 0.2e1 * t9 * t2 - 0.2e1 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t44, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t46 - t35, t25 * t46 - t29 * t45, qJD(2) * t32, (-pkin(2) * t26 + pkin(5) * t32) * qJD(2), 0, 0, 0, 0, 0, 0, t14 * t46 - t29 * t6, t15 * t46 + t29 * t5, t1 * t14 + t10 * t6 - t2 * t15 - t9 * t5, -pkin(3) * t35 + t1 * t55 + t10 * t3 + t2 * t7 + t21 * t46 - t9 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t37, 0.2e1 * (-t22 + t23) * qJD(3), 0, -0.2e1 * t37, 0, 0, t25 * t43, t28 * t43, 0, 0, -0.2e1 * t15 * t5, 0.2e1 * t5 * t14 - 0.2e1 * t15 * t6, 0, 0.2e1 * t14 * t6, 0, 0, 0.2e1 * t14 * t42 + 0.2e1 * t21 * t6, 0.2e1 * t15 * t42 - 0.2e1 * t21 * t5, 0.2e1 * t3 * t14 - 0.2e1 * t4 * t15 + 0.2e1 * t7 * t5 + 0.2e1 * t55 * t6, 0.2e1 * t21 * t42 + 0.2e1 * t3 * t55 + 0.2e1 * t7 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26 * t45 - t34, t26 * t47 - t33, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, (-t1 * t24 + t2 * t27 + (-t10 * t27 + t24 * t9) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, -t47, 0, -pkin(5) * t45, pkin(5) * t47, 0, 0, 0, 0, -t5, 0, -t6, 0, t4, t3, (-t24 * t6 + t27 * t5 + (-t14 * t27 + t15 * t24) * qJD(4)) * pkin(3), (-t24 * t3 + t27 * t4 + (-t24 * t7 - t27 * t55) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t41, -0.2e1 * t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
