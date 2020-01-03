% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:41
% EndTime: 2019-12-31 16:41:42
% DurationCPUTime: 0.28s
% Computational Cost: add. (429->71), mult. (943->106), div. (0->0), fcn. (597->4), ass. (0->55)
t38 = cos(pkin(6));
t41 = cos(qJ(4));
t58 = t41 * t38;
t52 = qJD(1) * t58;
t40 = sin(qJ(4));
t37 = sin(pkin(6));
t55 = qJD(1) * t37;
t53 = t40 * t55;
t17 = t52 - t53;
t65 = -t40 * t37 + t58;
t39 = -pkin(1) - qJ(3);
t64 = qJD(1) * t39;
t57 = t37 ^ 2 + t38 ^ 2;
t63 = t57 * qJD(3);
t62 = t17 ^ 2;
t35 = qJD(1) * qJD(2);
t51 = 0.2e1 * t35;
t61 = -pkin(5) + t39;
t21 = t41 * t37 + t40 * t38;
t56 = qJD(1) * t21;
t60 = t17 * t56;
t20 = t65 * qJD(4);
t54 = t20 * qJD(4);
t36 = qJD(1) * qJ(2);
t32 = qJD(3) + t36;
t28 = qJD(2) + t64;
t50 = -pkin(5) * qJD(1) + t28;
t49 = qJD(1) * t57;
t19 = t21 * qJD(4);
t9 = qJD(1) * t19;
t48 = -t17 * t19 - t65 * t9;
t27 = qJD(4) * t53;
t10 = qJD(4) * t52 - t27;
t47 = t21 * t10 + t20 * t56;
t11 = t50 * t37;
t12 = t50 * t38;
t6 = t41 * t11 + t40 * t12;
t5 = -t40 * t11 + t41 * t12;
t23 = t61 * t37;
t24 = t61 * t38;
t8 = t41 * t23 + t40 * t24;
t7 = -t40 * t23 + t41 * t24;
t45 = t17 * qJD(3);
t44 = t21 * qJD(3);
t1 = -qJD(1) * t44 + t5 * qJD(4);
t2 = -t6 * qJD(4) - t45;
t43 = t1 * t21 - t5 * t19 + t2 * t65 + t6 * t20;
t42 = qJD(1) ^ 2;
t31 = t37 * pkin(3) + qJ(2);
t25 = pkin(3) * t55 + t32;
t14 = t56 ^ 2;
t13 = t19 * qJD(4);
t4 = -qJD(3) * t65 - t8 * qJD(4);
t3 = t7 * qJD(4) - t44;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, qJ(2) * t51, 0, 0, 0, 0, 0, 0, t37 * t51, t38 * t51, 0.2e1 * qJD(3) * t49, (t32 + t36) * qJD(2) + (-t28 - t64) * t63, t48, -t10 * t65 - t17 * t20 + t19 * t56 + t9 * t21, -t13, t47, -t54, 0, 0.2e1 * t56 * qJD(2) + t4 * qJD(4) + t31 * t10 + t25 * t20, -t3 * qJD(4) - t25 * t19 - t31 * t9 + (qJD(1) * t65 + t17) * qJD(2), -t8 * t10 - t4 * t17 - t3 * t56 + t7 * t9 - t43, t1 * t8 + t2 * t7 + t6 * t3 + t5 * t4 + (qJD(1) * t31 + t25) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t42 * qJ(2), 0, 0, 0, 0, 0, 0, -t42 * t37, -t42 * t38, 0, (-t32 - t63) * qJD(1), 0, 0, 0, 0, 0, 0, -qJD(1) * t56 - t13, -qJD(1) * t17 - t54, -t47 - t48, -t25 * qJD(1) + t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 * t42, t28 * t49 + t35, 0, 0, 0, 0, 0, 0, -t27 + (t17 + t52) * qJD(4), -0.2e1 * t56 * qJD(4), -t14 - t62, t5 * t17 + t56 * t6 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t14 + t62, 0, -t60, t27 + (t17 - t52) * qJD(4), 0, -t25 * t17 - t45, (qJD(3) + t25) * t56, 0, 0;];
tauc_reg = t15;
