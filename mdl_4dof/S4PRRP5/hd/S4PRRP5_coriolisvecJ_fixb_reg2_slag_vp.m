% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRRP5
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:22
% EndTime: 2019-12-31 16:29:23
% DurationCPUTime: 0.35s
% Computational Cost: add. (236->83), mult. (686->129), div. (0->0), fcn. (352->4), ass. (0->72)
t52 = (qJD(2) * qJD(3));
t73 = -2 * t52;
t31 = sin(qJ(2));
t34 = qJD(3) ^ 2;
t35 = qJD(2) ^ 2;
t72 = (t34 + t35) * t31;
t55 = t31 * qJD(1);
t25 = qJD(2) * t55;
t30 = sin(qJ(3));
t48 = t30 * t52;
t13 = pkin(3) * t48 + t25;
t32 = cos(qJ(3));
t58 = qJD(3) * pkin(3);
t59 = qJD(2) * pkin(5);
t21 = t55 + t59;
t45 = qJ(4) * qJD(2) + t21;
t6 = t45 * t30;
t5 = -t6 + t58;
t7 = t45 * t32;
t41 = t30 * t5 - t32 * t7;
t71 = t41 * qJD(2) + t13;
t70 = t5 + t6;
t28 = t30 ^ 2;
t69 = pkin(3) * t28;
t68 = t32 * pkin(3);
t67 = t32 * t21;
t66 = t34 * t30;
t27 = t34 * t32;
t33 = cos(qJ(2));
t65 = t35 * t33;
t64 = -qJ(4) - pkin(5);
t29 = t32 ^ 2;
t63 = t28 - t29;
t62 = t28 + t29;
t60 = qJD(2) * pkin(2);
t26 = -pkin(2) - t68;
t54 = t33 * qJD(1);
t11 = t26 * qJD(2) + qJD(4) - t54;
t57 = qJD(2) * t11;
t56 = qJD(3) * t30;
t53 = qJ(4) * qJD(3);
t23 = t30 * t35 * t32;
t51 = t30 * t53;
t50 = t32 * t53;
t49 = qJD(3) * t54;
t47 = t32 * t52;
t46 = qJD(3) * t64;
t44 = t33 * t73;
t43 = qJD(4) + t54;
t42 = t30 * t47;
t40 = -t13 + t25;
t39 = -t11 - t43;
t22 = -t54 - t60;
t38 = qJD(3) * (t22 - t60);
t37 = qJD(2) * (-t22 - t54);
t1 = -t21 * t56 + (t43 * t32 - t51) * qJD(2);
t2 = -qJD(3) * t67 + (-t43 * t30 - t50) * qJD(2);
t36 = t1 * t32 - t2 * t30 + (-t30 * t7 - t32 * t5) * qJD(3);
t20 = t32 * t49;
t19 = t30 * t49;
t18 = -0.2e1 * t42;
t17 = 0.2e1 * t42;
t16 = t64 * t32;
t15 = t64 * t30;
t14 = t63 * t35;
t12 = t62 * t65;
t10 = t63 * t73;
t9 = -t30 * qJD(4) + t32 * t46;
t8 = t32 * qJD(4) + t30 * t46;
t4 = t30 * t44 - t32 * t72;
t3 = t30 * t72 + t32 * t44;
t24 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35 * t31, -t65, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, t12, (t22 * t31 + (-t55 + (t21 + t55) * t62) * t33) * qJD(2), 0, 0, 0, 0, 0, 0, t4, t3, t12, -t71 * t33 + (t36 + t57) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t10, t27, t18, -t66, 0, -pkin(5) * t27 + t30 * t38 + t19, pkin(5) * t66 + t32 * t38 + t20, 0, ((-t22 - t60) * t31 + (-t21 + t59) * t33 * t62) * qJD(1), t17, t10, t27, t18, -t66, 0, t19 + t40 * t32 + (t9 + (t11 + (t26 - t68) * qJD(2)) * t30) * qJD(3), t20 - t40 * t30 + (t11 * t32 - t8 + (t26 * t32 + t69) * qJD(2)) * qJD(3), (-t30 * t9 + t32 * t8 + (-t15 * t32 + t16 * t30) * qJD(3) - t62 * t54) * qJD(2) + t36, t11 * pkin(3) * t56 - t1 * t16 + t13 * t26 + t2 * t15 + t5 * t9 + t7 * t8 + (-t11 * t31 + t41 * t33) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t14, 0, t23, 0, 0, t30 * t37, t32 * t37, 0, 0, -t23, t14, 0, t23, 0, 0, pkin(3) * t23 + (t7 - t67) * qJD(3) + (t39 * t30 - t50) * qJD(2), -t35 * t69 + (t30 * t21 - t6) * qJD(3) + (t39 * t32 + t51) * qJD(2), (-t58 + t70) * t32 * qJD(2), t70 * t7 + (-t30 * t57 + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t48, 0.2e1 * t47, -t62 * t35, t71;];
tauc_reg = t24;
