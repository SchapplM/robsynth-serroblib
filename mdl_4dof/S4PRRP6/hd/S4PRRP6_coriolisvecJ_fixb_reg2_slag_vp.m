% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRRP6
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:43
% EndTime: 2019-12-31 16:30:44
% DurationCPUTime: 0.27s
% Computational Cost: add. (214->75), mult. (617->109), div. (0->0), fcn. (313->4), ass. (0->60)
t28 = sin(qJ(2));
t31 = qJD(3) ^ 2;
t32 = qJD(2) ^ 2;
t61 = (t31 + t32) * t28;
t47 = t28 * qJD(1);
t52 = qJD(2) * pkin(5);
t21 = t47 + t52;
t27 = sin(qJ(3));
t60 = t27 * t21;
t29 = cos(qJ(3));
t59 = t27 * t29;
t58 = t31 * t27;
t24 = t31 * t29;
t30 = cos(qJ(2));
t57 = t32 * t30;
t46 = t30 * qJD(1);
t43 = qJD(2) * t46;
t19 = t27 * t43;
t48 = qJD(3) * t29;
t5 = t21 * t48 + t19;
t25 = t27 ^ 2;
t26 = t29 ^ 2;
t56 = -t25 + t26;
t55 = t25 + t26;
t53 = qJD(2) * pkin(2);
t14 = -t29 * pkin(3) - t27 * qJ(4) - pkin(2);
t51 = qJD(2) * t14;
t50 = qJD(2) * t27;
t49 = qJD(3) * t27;
t45 = qJD(3) * qJ(4);
t44 = qJD(2) * qJD(3);
t22 = -t46 - t53;
t42 = t22 - t53;
t6 = -t46 + t51;
t41 = t6 + t51;
t40 = -qJD(3) * pkin(3) + qJD(4);
t39 = -0.2e1 * t30 * t44;
t38 = t44 * t59;
t8 = t40 + t60;
t9 = t29 * t21 + t45;
t37 = t27 * t9 - t29 * t8;
t36 = t27 * t8 + t29 * t9;
t35 = pkin(3) * t27 - qJ(4) * t29;
t7 = t35 * qJD(3) - t27 * qJD(4);
t1 = (t7 + t47) * qJD(2);
t34 = -pkin(5) * t31 - t1 + (-t7 + t47) * qJD(2);
t20 = t29 * t43;
t2 = t20 + (qJD(4) - t60) * qJD(3);
t33 = -t37 * qJD(3) + t2 * t29 + t5 * t27;
t23 = t32 * t59;
t18 = t46 * t49;
t17 = -0.2e1 * t38;
t16 = 0.2e1 * t38;
t15 = t56 * t32;
t13 = t35 * qJD(2);
t12 = t56 * t44;
t10 = t55 * t57;
t4 = t27 * t39 - t29 * t61;
t3 = t27 * t61 + t29 * t39;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32 * t28, -t57, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, t10, (t22 * t28 + (-t47 + (t21 + t47) * t55) * t30) * qJD(2), 0, 0, 0, 0, 0, 0, t4, t10, -t3, (t36 * qJD(2) - t1) * t30 + (qJD(2) * t6 + t33) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0.2e1 * t12, t24, t17, -t58, 0, -pkin(5) * t24 + t42 * t49 + t18, pkin(5) * t58 + (t42 + t46) * t48, 0, ((-t22 - t53) * t28 + (-t21 + t52) * t30 * t55) * qJD(1), t16, t24, -0.2e1 * t12, 0, t58, t17, t34 * t29 + t41 * t49 + t18, -t55 * t43 + t33, (-t41 - t46) * t48 + t34 * t27, t1 * t14 + t6 * t7 + (-t28 * t6 - t36 * t30) * qJD(1) + t33 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t15, 0, t23, 0, 0, -t22 * t50 - t19, -t22 * t29 * qJD(2) - t20, 0, 0, -t23, 0, t15, 0, 0, t23, -t19 + (t13 * t29 - t27 * t6) * qJD(2), ((t9 - t45) * t27 + (t40 - t8) * t29) * qJD(2), 0.2e1 * qJD(3) * qJD(4) + t20 + (t13 * t27 + t29 * t6) * qJD(2), -t5 * pkin(3) + t2 * qJ(4) + t9 * qJD(4) - t6 * t13 + t37 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t25 * t32 - t31, -t9 * qJD(3) + t6 * t50 + t5;];
tauc_reg = t11;
