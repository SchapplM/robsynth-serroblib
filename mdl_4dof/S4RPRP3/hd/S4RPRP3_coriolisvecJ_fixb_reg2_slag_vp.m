% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRP3
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:48
% EndTime: 2019-12-31 16:42:49
% DurationCPUTime: 0.30s
% Computational Cost: add. (270->73), mult. (697->114), div. (0->0), fcn. (342->4), ass. (0->59)
t34 = cos(qJ(3));
t24 = sin(pkin(6)) * pkin(1) + pkin(5);
t19 = t24 * qJD(1);
t41 = qJ(4) * qJD(1) + t19;
t33 = sin(qJ(3));
t50 = t33 * qJD(2);
t5 = t41 * t34 + t50;
t27 = t34 * qJD(2);
t4 = -t41 * t33 + t27;
t56 = qJD(3) * pkin(3);
t3 = t4 + t56;
t64 = t3 - t4;
t29 = t33 ^ 2;
t63 = pkin(3) * t29;
t62 = pkin(3) * t34;
t61 = t34 * t5;
t36 = qJD(1) ^ 2;
t60 = t34 * t36;
t35 = qJD(3) ^ 2;
t59 = t35 * t33;
t28 = t35 * t34;
t51 = qJD(3) * t33;
t58 = -qJD(3) * t27 + t19 * t51;
t30 = t34 ^ 2;
t57 = t29 - t30;
t55 = qJ(4) + t24;
t25 = -cos(pkin(6)) * pkin(1) - pkin(2);
t17 = t25 - t62;
t54 = qJD(1) * t17;
t20 = qJD(1) * t25;
t53 = qJD(1) * t33;
t52 = qJD(1) * t34;
t49 = t33 * qJD(4);
t48 = t34 * qJD(4);
t12 = qJD(4) + t54;
t47 = -qJD(4) - t12;
t46 = qJD(1) * qJD(3);
t45 = 0.2e1 * t46;
t44 = qJ(4) * t51;
t43 = t34 * t46;
t42 = qJD(3) * t55;
t40 = t33 * t43;
t11 = t19 * t34 + t50;
t39 = 0.2e1 * qJD(3) * t20;
t10 = -t19 * t33 + t27;
t9 = qJD(3) * t11;
t37 = t33 * t9 - t34 * t58 + (-t10 * t34 - t11 * t33) * qJD(3);
t23 = t33 * t60;
t22 = -0.2e1 * t40;
t21 = 0.2e1 * t40;
t18 = t57 * t36;
t15 = t55 * t34;
t14 = t55 * t33;
t13 = t57 * t45;
t7 = -t34 * t42 - t49;
t6 = -t33 * t42 + t48;
t2 = -qJD(1) * t49 - qJD(3) * t5;
t1 = (-t44 + t48) * qJD(1) - t58;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t13, t28, t22, -t59, 0, -t24 * t28 + t33 * t39, t24 * t59 + t34 * t39, t37, t37 * t24, t21, -t13, t28, t22, -t59, 0, (t7 + (t12 + (t17 - 0.2e1 * t62) * qJD(1)) * t33) * qJD(3), (t12 * t34 - t6 + (t17 * t34 + 0.2e1 * t63) * qJD(1)) * qJD(3), t1 * t34 - t2 * t33 + (-t3 * t34 - t33 * t5) * qJD(3) + (-t33 * t7 + t34 * t6 + (t14 * t34 - t15 * t33) * qJD(3)) * qJD(1), t1 * t15 - t14 * t2 + t3 * t7 + t5 * t6 + (t12 + t54) * pkin(3) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t28, 0, -t33 * t58 - t34 * t9 + (-t10 * t33 + t11 * t34) * qJD(3), 0, 0, 0, 0, 0, 0, -t59, -t28, 0, t1 * t33 + t2 * t34 + (-t3 * t33 + t61) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t18, 0, t23, 0, 0, -t20 * t53, qJD(3) * t10 - t20 * t52 + t58, 0, 0, -t23, t18, 0, t23, 0, 0, (pkin(3) * t60 + t47 * qJD(1)) * t33, -t36 * t63 + qJD(3) * t4 + (t47 * t34 + t44) * qJD(1) + t58, (-t56 + t64) * t52, t64 * t5 + (-t12 * t53 + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 * t45, 0.2e1 * t43, (-t29 - t30) * t36, (-t61 + (t3 + t56) * t33) * qJD(1);];
tauc_reg = t8;
