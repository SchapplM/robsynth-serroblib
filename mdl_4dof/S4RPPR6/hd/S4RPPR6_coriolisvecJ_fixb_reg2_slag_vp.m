% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:44
% EndTime: 2019-12-31 16:40:45
% DurationCPUTime: 0.34s
% Computational Cost: add. (339->82), mult. (944->131), div. (0->0), fcn. (623->4), ass. (0->58)
t39 = cos(qJ(4));
t36 = sin(pkin(6));
t54 = qJD(1) * t36;
t49 = t39 * t54;
t37 = cos(pkin(6));
t38 = sin(qJ(4));
t60 = t37 * t38;
t50 = qJD(1) * t60;
t18 = t49 - t50;
t22 = t36 * t39 - t60;
t21 = t36 * t38 + t37 * t39;
t55 = qJD(1) * t21;
t64 = t55 ^ 2;
t63 = t18 ^ 2;
t62 = t18 * t55;
t59 = -pkin(5) + qJ(2);
t34 = t36 ^ 2;
t35 = t37 ^ 2;
t58 = t34 + t35;
t57 = qJ(2) * t35;
t47 = t36 * qJ(3) + pkin(1);
t19 = (pkin(2) + pkin(3)) * t37 + t47;
t56 = qJD(1) * t19;
t53 = qJD(3) * t36;
t29 = qJ(2) * t54 + qJD(3);
t52 = qJD(1) * qJD(2);
t51 = qJD(1) * qJD(3);
t26 = t59 * t37;
t41 = qJD(1) ^ 2;
t24 = t58 * t41;
t48 = t36 * t51;
t46 = qJ(2) * t52;
t45 = 0.2e1 * t55;
t20 = -pkin(5) * t54 + t29;
t23 = qJD(1) * t26;
t5 = t20 * t39 - t23 * t38;
t6 = t20 * t38 + t23 * t39;
t25 = t59 * t36;
t7 = t25 * t39 - t26 * t38;
t8 = t25 * t38 + t26 * t39;
t44 = -pkin(2) * t37 - t47;
t43 = t18 * qJD(2);
t42 = t21 * qJD(2);
t14 = t21 * qJD(4);
t40 = qJD(4) ^ 2;
t30 = t34 * t46;
t28 = qJD(4) * t50;
t15 = t22 * qJD(4);
t13 = 0.2e1 * t58 * t52;
t12 = qJD(1) * t44 + qJD(2);
t11 = qJD(4) * t49 - t28;
t10 = qJD(1) * t14;
t9 = -qJD(2) + t56;
t4 = qJD(2) * t22 - qJD(4) * t8;
t3 = qJD(4) * t7 + t42;
t2 = -t6 * qJD(4) + t43;
t1 = qJD(1) * t42 + t5 * qJD(4);
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0.2e1 * t35 * t46 + 0.2e1 * t30, 0, 0, 0, 0, 0, 0, 0.2e1 * t37 * t48, t13, 0.2e1 * t34 * t51, t30 + (t29 * qJD(2) - t12 * qJD(3)) * t36 + (0.2e1 * qJD(2) * t57 - t44 * t53) * qJD(1), -t10 * t22 - t14 * t18, t10 * t21 - t11 * t22 + t14 * t55 - t15 * t18, -t14 * qJD(4), t11 * t21 + t15 * t55, -t15 * qJD(4), 0, t4 * qJD(4) + t19 * t11 + t9 * t15 + t45 * t53, -t3 * qJD(4) - t19 * t10 - t9 * t14 + (qJD(1) * t22 + t18) * t53, -t1 * t21 + t10 * t7 - t11 * t8 + t14 * t5 - t15 * t6 - t18 * t4 - t2 * t22 - t3 * t55, t1 * t8 + t2 * t7 + t6 * t3 + t5 * t4 + (t9 + t56) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -qJ(2) * t24, 0, 0, 0, 0, 0, 0, 0, -t24, 0, -t41 * t57 + (-qJD(3) - t29) * t54, 0, 0, 0, 0, 0, 0, t28 + (-t18 - t49) * qJD(4), t45 * qJD(4), t63 + t64, -t18 * t5 - t55 * t6 - t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 * t41 * t37, 0, -t34 * t41, (qJD(2) + t12) * t54, 0, 0, 0, 0, 0, 0, -t38 * t40 - t54 * t55, -t18 * t54 - t39 * t40, t39 * t10 - t38 * t11 + (t18 * t38 - t39 * t55) * qJD(4), -t9 * t54 + t1 * t38 + t2 * t39 + (-t38 * t5 + t39 * t6) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t63 - t64, 0, -t62, t28 + (t18 - t49) * qJD(4), 0, -t9 * t18 + t43, -(qJD(2) - t9) * t55, 0, 0;];
tauc_reg = t16;
