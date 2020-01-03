% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:53
% EndTime: 2019-12-31 16:43:54
% DurationCPUTime: 0.24s
% Computational Cost: add. (256->58), mult. (656->96), div. (0->0), fcn. (316->4), ass. (0->48)
t29 = cos(qJ(3));
t44 = t29 * qJD(2);
t20 = sin(pkin(6)) * pkin(1) + pkin(5);
t15 = t20 * qJD(1);
t28 = sin(qJ(3));
t50 = t28 * t15;
t8 = t44 - t50;
t53 = qJD(4) - t8;
t2 = -qJD(3) * pkin(3) + t53;
t9 = qJD(2) * t28 + t15 * t29;
t3 = qJD(3) * qJ(4) + t9;
t7 = t9 * qJD(3);
t52 = t7 * t28;
t51 = t7 * t29;
t49 = t28 * t29;
t30 = qJD(3) ^ 2;
t48 = t30 * t28;
t23 = t30 * t29;
t24 = t28 ^ 2;
t47 = t29 ^ 2 - t24;
t37 = pkin(3) * t28 - qJ(4) * t29;
t10 = qJD(3) * t37 - t28 * qJD(4);
t5 = qJD(1) * t10;
t21 = -cos(pkin(6)) * pkin(1) - pkin(2);
t11 = -pkin(3) * t29 - qJ(4) * t28 + t21;
t4 = qJD(1) * t11;
t16 = qJD(1) * t21;
t46 = qJD(1) * t28;
t45 = qJD(3) * t28;
t42 = qJD(1) * qJD(3);
t41 = t8 + t50;
t40 = 0.2e1 * t4;
t38 = t42 * t49;
t36 = 0.2e1 * qJD(3) * t16;
t35 = -t20 * t30 - 0.2e1 * t5;
t22 = qJD(3) * t44;
t1 = t22 + (qJD(4) - t50) * qJD(3);
t33 = t1 * t29 + t52 + (t2 * t29 - t28 * t3) * qJD(3);
t6 = -t15 * t45 + t22;
t32 = t52 + t6 * t29 + (-t28 * t9 - t29 * t8) * qJD(3);
t31 = qJD(1) ^ 2;
t19 = t31 * t49;
t18 = -0.2e1 * t38;
t17 = 0.2e1 * t38;
t14 = t47 * t31;
t13 = t37 * qJD(1);
t12 = t47 * t42;
t25 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0.2e1 * t12, t23, t18, -t48, 0, -t20 * t23 + t28 * t36, t20 * t48 + t29 * t36, t32, t32 * t20, t17, t23, -0.2e1 * t12, 0, t48, t18, t29 * t35 + t40 * t45, t33, -qJD(3) * t29 * t40 + t28 * t35, t4 * t10 + t5 * t11 + t20 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t23, 0, t6 * t28 - t51 + (-t28 * t8 + t29 * t9) * qJD(3), 0, 0, 0, 0, 0, 0, -t48, 0, t23, t1 * t28 - t51 + (t2 * t28 + t29 * t3) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t14, 0, t19, 0, 0, -t16 * t46, -t16 * t29 * qJD(1) + qJD(3) * t41 - t22, 0, 0, -t19, 0, t14, 0, 0, t19, (t13 * t29 - t28 * t4) * qJD(1), 0, t22 + (t13 * t28 + t29 * t4) * qJD(1) + (0.2e1 * qJD(4) - t41) * qJD(3), -t7 * pkin(3) + t1 * qJ(4) - t4 * t13 - t2 * t9 + t3 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, -t24 * t31 - t30, t4 * t46 + (-t3 + t9) * qJD(3);];
tauc_reg = t25;
