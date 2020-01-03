% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:54
% EndTime: 2019-12-31 17:04:55
% DurationCPUTime: 0.36s
% Computational Cost: add. (632->72), mult. (1483->158), div. (0->0), fcn. (1330->6), ass. (0->50)
t46 = sin(pkin(7));
t64 = pkin(2) * t46;
t63 = cos(qJ(4));
t62 = -qJ(3) - pkin(5);
t49 = sin(qJ(2));
t38 = t62 * t49;
t50 = cos(qJ(2));
t39 = t62 * t50;
t47 = cos(pkin(7));
t17 = t46 * t38 - t47 * t39;
t58 = t50 * qJD(2);
t59 = t49 * qJD(2);
t61 = t46 * t58 + t47 * t59;
t48 = sin(qJ(4));
t60 = qJD(4) * t48;
t57 = t48 * t64;
t56 = -0.2e1 * pkin(1) * qJD(2);
t45 = pkin(2) * t59;
t55 = t49 * t58;
t44 = -t50 * pkin(2) - pkin(1);
t16 = t47 * t38 + t46 * t39;
t54 = qJD(4) * t63;
t53 = qJD(2) * t62;
t28 = t50 * qJD(3) + t49 * t53;
t29 = -t49 * qJD(3) + t50 * t53;
t11 = t47 * t28 + t46 * t29;
t10 = -t46 * t28 + t47 * t29;
t33 = t46 * t50 + t47 * t49;
t12 = -t33 * pkin(6) + t16;
t32 = t46 * t49 - t47 * t50;
t13 = -t32 * pkin(6) + t17;
t4 = t48 * t12 + t63 * t13;
t15 = -t48 * t32 + t63 * t33;
t43 = t47 * pkin(2) + pkin(3);
t27 = t48 * t43 + t63 * t64;
t31 = -t46 * t59 + t47 * t58;
t52 = -t31 * pkin(6) + t10;
t51 = -t61 * pkin(6) + t11;
t26 = t63 * t43 - t57;
t21 = t27 * qJD(4);
t20 = qJD(4) * t57 - t43 * t54;
t19 = t32 * pkin(3) + t44;
t18 = t61 * pkin(3) + t45;
t14 = t63 * t32 + t48 * t33;
t6 = t15 * qJD(4) + t48 * t31 + t63 * t61;
t5 = -t63 * t31 + t32 * t54 + t33 * t60 + t48 * t61;
t3 = t63 * t12 - t48 * t13;
t2 = -t4 * qJD(4) - t48 * t51 + t63 * t52;
t1 = -t12 * t54 + t13 * t60 - t48 * t52 - t63 * t51;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t55, 0.2e1 * (-t49 ^ 2 + t50 ^ 2) * qJD(2), 0, -0.2e1 * t55, 0, 0, t49 * t56, t50 * t56, 0, 0, 0.2e1 * t33 * t31, -0.2e1 * t31 * t32 - 0.2e1 * t33 * t61, 0, 0.2e1 * t32 * t61, 0, 0, 0.2e1 * t32 * t45 + 0.2e1 * t44 * t61, 0.2e1 * t44 * t31 + 0.2e1 * t33 * t45, -0.2e1 * t10 * t33 - 0.2e1 * t11 * t32 - 0.2e1 * t16 * t31 - 0.2e1 * t17 * t61, 0.2e1 * t16 * t10 + 0.2e1 * t17 * t11 + 0.2e1 * t44 * t45, -0.2e1 * t15 * t5, 0.2e1 * t5 * t14 - 0.2e1 * t15 * t6, 0, 0.2e1 * t14 * t6, 0, 0, 0.2e1 * t18 * t14 + 0.2e1 * t19 * t6, 0.2e1 * t18 * t15 - 0.2e1 * t19 * t5, 0.2e1 * t1 * t14 - 0.2e1 * t2 * t15 + 0.2e1 * t3 * t5 - 0.2e1 * t4 * t6, -0.2e1 * t4 * t1 + 0.2e1 * t19 * t18 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, -t59, 0, -pkin(5) * t58, pkin(5) * t59, 0, 0, 0, 0, t31, 0, -t61, 0, t10, -t11, (-t47 * t31 - t46 * t61) * pkin(2), (t10 * t47 + t11 * t46) * pkin(2), 0, 0, -t5, 0, -t6, 0, t2, t1, t20 * t14 + t21 * t15 + t26 * t5 - t27 * t6, -t1 * t27 + t2 * t26 - t4 * t20 - t3 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t21, 0.2e1 * t20, 0, -0.2e1 * t27 * t20 - 0.2e1 * t26 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t31, 0, t45, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t20, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
