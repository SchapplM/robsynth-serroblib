% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:14
% EndTime: 2019-12-31 16:46:15
% DurationCPUTime: 0.26s
% Computational Cost: add. (210->69), mult. (474->104), div. (0->0), fcn. (187->2), ass. (0->62)
t32 = -pkin(1) - pkin(5);
t19 = t32 * qJD(1) + qJD(2);
t30 = sin(qJ(3));
t44 = qJD(1) * qJD(3);
t41 = t30 * t44;
t21 = qJ(4) * t41;
t31 = cos(qJ(3));
t48 = t31 * qJD(4);
t52 = qJD(3) * t30;
t1 = -qJD(1) * t48 - t19 * t52 + t21;
t51 = qJD(3) * t31;
t42 = qJ(4) * t51;
t49 = t30 * qJD(4);
t2 = t19 * t51 + (-t42 - t49) * qJD(1);
t56 = qJD(3) * pkin(3);
t45 = qJ(4) * qJD(1);
t9 = t31 * t19;
t7 = -t31 * t45 + t9;
t3 = t7 + t56;
t62 = t30 * t19;
t6 = -t30 * t45 + t62;
t65 = -t1 * t31 - t2 * t30 + (t3 * t30 - t31 * t6) * qJD(3);
t27 = qJD(1) * qJD(2);
t64 = 0.2e1 * t27;
t63 = t3 - t7;
t33 = qJD(3) ^ 2;
t61 = t33 * t30;
t60 = t33 * t31;
t34 = qJD(1) ^ 2;
t59 = t34 * t30;
t40 = t31 * t44;
t10 = pkin(3) * t40 + t27;
t28 = t30 ^ 2;
t29 = t31 ^ 2;
t58 = t28 - t29;
t57 = -t33 - t34;
t55 = t34 * qJ(2);
t54 = qJ(4) - t32;
t24 = t30 * pkin(3) + qJ(2);
t53 = qJD(1) * t24;
t11 = qJD(4) + t53;
t50 = t11 * qJD(1);
t47 = qJD(4) + t11;
t46 = qJ(2) * qJD(3);
t43 = 0.2e1 * qJD(1);
t13 = t54 * t31;
t20 = pkin(3) * t51 + qJD(2);
t39 = qJD(1) * t20 + t10;
t38 = t11 + t53;
t37 = t30 * t40;
t26 = qJ(2) * t64;
t22 = t31 * t59;
t18 = -0.2e1 * t37;
t17 = 0.2e1 * t37;
t16 = t57 * t31;
t15 = t57 * t30;
t14 = t58 * t34;
t12 = t54 * t30;
t8 = 0.2e1 * t58 * t44;
t5 = -qJD(3) * t13 - t49;
t4 = t54 * t52 - t48;
t23 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t26, t18, t8, -t61, t17, -t60, 0, -t32 * t61 + (qJD(2) * t30 + t31 * t46) * t43, -t32 * t60 + (qJD(2) * t31 - t30 * t46) * t43, 0, t26, t18, t8, -t61, t17, -t60, 0, t39 * t30 + (t38 * t31 + t4) * qJD(3), t39 * t31 + (-t38 * t30 - t5) * qJD(3), (-t30 * t5 - t31 * t4 + (t12 * t31 - t13 * t30) * qJD(3)) * qJD(1) + t65, -t1 * t13 + t10 * t24 + t11 * t20 - t2 * t12 + t3 * t4 + t6 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t55, 0, 0, 0, 0, 0, 0, t15, t16, 0, -t55, 0, 0, 0, 0, 0, 0, t15, t16, 0, -t50 - t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t14, 0, -t22, 0, 0, -t31 * t55, t30 * t55, 0, 0, t22, -t14, 0, -t22, 0, 0, t21 + (t6 - t62) * qJD(3) + (-pkin(3) * t59 - t47 * qJD(1)) * t31, -t29 * t34 * pkin(3) + (t7 - t9) * qJD(3) + (t47 * t30 + t42) * qJD(1), (t56 - t63) * t30 * qJD(1), t63 * t6 + (-t31 * t50 + t1) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t40, -0.2e1 * t41, (-t28 - t29) * t34, (t3 * t31 + t30 * t6) * qJD(1) + t10;];
tauc_reg = t23;
