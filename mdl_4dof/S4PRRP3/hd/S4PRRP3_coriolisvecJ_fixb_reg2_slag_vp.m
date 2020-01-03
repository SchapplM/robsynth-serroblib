% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:55
% EndTime: 2019-12-31 16:26:56
% DurationCPUTime: 0.26s
% Computational Cost: add. (188->70), mult. (554->113), div. (0->0), fcn. (260->2), ass. (0->56)
t39 = (qJD(2) * qJD(3));
t58 = -2 * t39;
t27 = sin(qJ(3));
t51 = -qJ(4) - pkin(5);
t15 = t51 * t27;
t28 = cos(qJ(3));
t23 = t28 * qJD(1);
t4 = qJD(2) * t15 + t23;
t48 = qJD(3) * pkin(3);
t3 = t4 + t48;
t57 = t3 - t4;
t25 = t27 ^ 2;
t56 = pkin(3) * t25;
t55 = t28 * pkin(3);
t16 = t51 * t28;
t43 = t27 * qJD(1);
t6 = -qJD(2) * t16 + t43;
t54 = t28 * t6;
t30 = qJD(2) ^ 2;
t53 = t28 * t30;
t29 = qJD(3) ^ 2;
t52 = t29 * t27;
t24 = t29 * t28;
t37 = t27 * t39;
t40 = qJD(1) * qJD(3);
t50 = pkin(5) * t37 - t28 * t40;
t26 = t28 ^ 2;
t49 = t25 - t26;
t22 = -pkin(2) - t55;
t47 = qJD(2) * t22;
t46 = qJD(2) * t27;
t45 = qJD(2) * t28;
t44 = qJD(3) * t27;
t42 = t28 * qJD(4);
t13 = qJD(4) + t47;
t41 = -qJD(4) - t13;
t19 = t27 * t53;
t38 = qJ(4) * t44;
t36 = t28 * t39;
t35 = qJD(3) * t51;
t34 = pkin(2) * t58;
t33 = t27 * t36;
t32 = t28 * t35;
t12 = pkin(5) * t45 + t43;
t7 = -t27 * qJD(4) + t32;
t10 = t12 * qJD(3);
t11 = -pkin(5) * t46 + t23;
t31 = t10 * t27 - t50 * t28 + (-t11 * t28 - t12 * t27) * qJD(3);
t18 = -0.2e1 * t33;
t17 = 0.2e1 * t33;
t14 = t49 * t30;
t8 = t49 * t58;
t5 = t27 * t35 + t42;
t2 = qJD(2) * t7 - t27 * t40;
t1 = (-t38 + t42) * qJD(2) - t50;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t24, 0, -t10 * t28 - t50 * t27 + (-t11 * t27 + t12 * t28) * qJD(3), 0, 0, 0, 0, 0, 0, -t52, -t24, 0, t1 * t27 + t2 * t28 + (-t27 * t3 + t54) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t8, t24, t18, -t52, 0, -pkin(5) * t24 + t27 * t34, pkin(5) * t52 + t28 * t34, t31, t31 * pkin(5), t17, t8, t24, t18, -t52, 0, (t7 + (t13 + (t22 - 0.2e1 * t55) * qJD(2)) * t27) * qJD(3), (t13 * t28 - t5 + (t22 * t28 + 0.2e1 * t56) * qJD(2)) * qJD(3), t1 * t28 - t2 * t27 + (-t27 * t6 - t28 * t3) * qJD(3) + (-t27 * t7 + t28 * t5 + (-t15 * t28 + t16 * t27) * qJD(3)) * qJD(2), -t1 * t16 + t2 * t15 + t3 * t7 + t6 * t5 + (t13 + t47) * pkin(3) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t14, 0, t19, 0, 0, t30 * pkin(2) * t27, pkin(2) * t53 + t11 * qJD(3) + t50, 0, 0, -t19, t14, 0, t19, 0, 0, pkin(3) * t19 + (t6 - t43) * qJD(3) + (t41 * t27 + t32) * qJD(2), -t30 * t56 + t4 * qJD(3) + (t41 * t28 + t38) * qJD(2) + t50, (-t48 + t57) * t45, t57 * t6 + (-t13 * t46 + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t37, 0.2e1 * t36, (-t25 - t26) * t30, (-t54 + (t3 + t48) * t27) * qJD(2);];
tauc_reg = t9;
