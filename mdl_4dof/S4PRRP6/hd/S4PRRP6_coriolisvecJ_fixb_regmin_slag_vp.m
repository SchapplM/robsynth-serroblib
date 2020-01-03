% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:42
% EndTime: 2019-12-31 16:30:43
% DurationCPUTime: 0.21s
% Computational Cost: add. (194->68), mult. (512->98), div. (0->0), fcn. (262->4), ass. (0->54)
t21 = sin(qJ(2));
t24 = qJD(3) ^ 2;
t25 = qJD(2) ^ 2;
t55 = (t24 + t25) * t21;
t41 = t21 * qJD(1);
t16 = qJD(2) * pkin(5) + t41;
t20 = sin(qJ(3));
t54 = t20 * t16;
t22 = cos(qJ(3));
t53 = t20 * t22;
t52 = t24 * t20;
t51 = t24 * t22;
t23 = cos(qJ(2));
t50 = t25 * t23;
t40 = t23 * qJD(1);
t35 = qJD(2) * t40;
t14 = t20 * t35;
t42 = qJD(3) * t22;
t5 = t16 * t42 + t14;
t18 = t20 ^ 2;
t19 = t22 ^ 2;
t49 = t18 - t19;
t48 = t18 + t19;
t46 = qJD(2) * pkin(2);
t12 = -t22 * pkin(3) - t20 * qJ(4) - pkin(2);
t45 = qJD(2) * t12;
t44 = qJD(2) * t20;
t43 = qJD(3) * t20;
t39 = qJD(3) * qJ(4);
t38 = qJD(2) * qJD(3);
t37 = t25 * t53;
t36 = 0.2e1 * t38;
t17 = -t40 - t46;
t34 = t17 - t46;
t6 = -t40 + t45;
t33 = t6 + t45;
t32 = -qJD(3) * pkin(3) + qJD(4);
t31 = -0.2e1 * t23 * t38;
t8 = t32 + t54;
t9 = t22 * t16 + t39;
t30 = t20 * t9 - t22 * t8;
t29 = t20 * t8 + t22 * t9;
t28 = pkin(3) * t20 - qJ(4) * t22;
t7 = t28 * qJD(3) - t20 * qJD(4);
t1 = (t7 + t41) * qJD(2);
t27 = -pkin(5) * t24 - t1 + (-t7 + t41) * qJD(2);
t15 = t22 * t35;
t2 = t15 + (qJD(4) - t54) * qJD(3);
t26 = -t30 * qJD(3) + t2 * t22 + t5 * t20;
t13 = t40 * t43;
t11 = t28 * qJD(2);
t4 = t20 * t31 - t22 * t55;
t3 = t20 * t55 + t22 * t31;
t10 = [0, 0, -t25 * t21, -t50, 0, 0, 0, 0, 0, t4, t3, t4, t48 * t50, -t3, (t29 * qJD(2) - t1) * t23 + (qJD(2) * t6 + t26) * t21; 0, 0, 0, 0, t36 * t53, -t49 * t36, t51, -t52, 0, -pkin(5) * t51 + t34 * t43 + t13, pkin(5) * t52 + (t34 + t40) * t42, t27 * t22 + t33 * t43 + t13, -t48 * t35 + t26, (-t33 - t40) * t42 + t27 * t20, t1 * t12 + t6 * t7 + (-t21 * t6 - t29 * t23) * qJD(1) + t26 * pkin(5); 0, 0, 0, 0, -t37, t49 * t25, 0, 0, 0, -t17 * t44 - t14, -t17 * t22 * qJD(2) - t15, -t14 + (t11 * t22 - t20 * t6) * qJD(2), ((t9 - t39) * t20 + (t32 - t8) * t22) * qJD(2), 0.2e1 * qJD(3) * qJD(4) + t15 + (t11 * t20 + t22 * t6) * qJD(2), -t5 * pkin(3) + t2 * qJ(4) + t9 * qJD(4) - t6 * t11 + t30 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, -t18 * t25 - t24, -t9 * qJD(3) + t6 * t44 + t5;];
tauc_reg = t10;
