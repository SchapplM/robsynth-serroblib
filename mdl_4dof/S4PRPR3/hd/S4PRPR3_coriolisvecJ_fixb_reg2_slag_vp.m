% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:55
% EndTime: 2019-12-31 16:20:56
% DurationCPUTime: 0.24s
% Computational Cost: add. (320->59), mult. (872->95), div. (0->0), fcn. (624->4), ass. (0->50)
t40 = cos(pkin(7));
t42 = cos(qJ(4));
t39 = sin(pkin(7));
t41 = sin(qJ(4));
t59 = t41 * t39;
t47 = -t42 * t40 + t59;
t45 = t47 * qJD(3);
t63 = qJD(2) * t45;
t24 = t42 * t39 + t41 * t40;
t19 = t24 * qJD(2);
t62 = t19 ^ 2;
t22 = t24 * qJD(4);
t12 = qJD(2) * t22;
t55 = qJD(2) * t40;
t51 = t42 * t55;
t52 = qJD(2) * t59;
t17 = -t51 + t52;
t21 = t47 * qJD(4);
t61 = -t24 * t12 + t21 * t17;
t60 = t19 * t17;
t57 = pkin(5) + qJ(3);
t53 = qJ(3) * qJD(2);
t26 = t39 * qJD(1) + t40 * t53;
t56 = t39 ^ 2 + t40 ^ 2;
t54 = t21 * qJD(4);
t34 = -t40 * pkin(3) - pkin(2);
t28 = t57 * t39;
t50 = t56 * qJD(2);
t30 = qJD(4) * t51;
t11 = qJD(4) * t52 - t30;
t49 = t11 * t47 - t19 * t22;
t36 = t40 * qJD(1);
t14 = -qJD(2) * t28 + t36;
t15 = pkin(5) * t55 + t26;
t5 = t42 * t14 - t41 * t15;
t6 = t41 * t14 + t42 * t15;
t48 = (-t39 * t53 + t36) * t39 - t26 * t40;
t29 = t57 * t40;
t9 = -t42 * t28 - t41 * t29;
t10 = -t41 * t28 + t42 * t29;
t46 = t24 * qJD(3);
t44 = qJD(2) * t46;
t27 = t34 * qJD(2) + qJD(3);
t16 = t17 ^ 2;
t13 = t22 * qJD(4);
t4 = -t10 * qJD(4) - t46;
t3 = t9 * qJD(4) - t45;
t2 = -t6 * qJD(4) - t44;
t1 = t5 * qJD(4) - t63;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t54, -t49 + t61, t1 * t24 - t2 * t47 - t6 * t21 - t5 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t50, (qJ(3) * t50 - t48) * qJD(3), -t11 * t24 - t19 * t21, t49 + t61, -t54, t12 * t47 + t17 * t22, -t13, 0, t4 * qJD(4) + t34 * t12 + t27 * t22, -t3 * qJD(4) - t34 * t11 - t27 * t21, -t1 * t47 - t10 * t12 + t9 * t11 - t3 * t17 - t4 * t19 - t2 * t24 + t5 * t21 - t6 * t22, t1 * t10 + t2 * t9 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56 * qJD(2) ^ 2, t48 * qJD(2), 0, 0, 0, 0, 0, 0, 0.2e1 * t19 * qJD(4), t30 + (-t17 - t52) * qJD(4), -t16 - t62, t6 * t17 + t5 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t16 + t62, t30 + (t17 - t52) * qJD(4), -t60, 0, 0, -t27 * t19 - t44, t27 * t17 + t63, 0, 0;];
tauc_reg = t7;
