% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:57
% DurationCPUTime: 0.24s
% Computational Cost: add. (404->62), mult. (1016->98), div. (0->0), fcn. (708->6), ass. (0->51)
t43 = cos(pkin(7));
t46 = cos(qJ(4));
t41 = sin(pkin(7));
t45 = sin(qJ(4));
t61 = t45 * t41;
t51 = -t46 * t43 + t61;
t49 = t51 * qJD(3);
t66 = qJD(1) * t49;
t30 = t46 * t41 + t45 * t43;
t22 = t30 * qJD(1);
t65 = t22 ^ 2;
t36 = sin(pkin(6)) * pkin(1) + qJ(3);
t64 = pkin(5) + t36;
t25 = t30 * qJD(4);
t16 = qJD(1) * t25;
t58 = qJD(1) * t43;
t55 = t46 * t58;
t56 = qJD(1) * t61;
t20 = -t55 + t56;
t24 = t51 * qJD(4);
t63 = -t30 * t16 + t24 * t20;
t62 = t22 * t20;
t32 = t36 * qJD(1);
t14 = t41 * qJD(2) + t43 * t32;
t59 = t41 ^ 2 + t43 ^ 2;
t57 = t24 * qJD(4);
t54 = qJD(1) * t59;
t38 = t43 * qJD(2);
t11 = t38 + (-pkin(5) * qJD(1) - t32) * t41;
t12 = pkin(5) * t58 + t14;
t3 = t46 * t11 - t45 * t12;
t4 = t45 * t11 + t46 * t12;
t53 = (-t41 * t32 + t38) * t41 - t14 * t43;
t33 = qJD(4) * t55;
t15 = qJD(4) * t56 - t33;
t52 = t15 * t51 - t22 * t25;
t26 = t64 * t41;
t27 = t64 * t43;
t7 = -t46 * t26 - t45 * t27;
t8 = -t45 * t26 + t46 * t27;
t31 = -cos(pkin(6)) * pkin(1) - pkin(2) - t43 * pkin(3);
t50 = t30 * qJD(3);
t48 = qJD(1) * t50;
t19 = t20 ^ 2;
t18 = t31 * qJD(1) + qJD(3);
t17 = t25 * qJD(4);
t6 = -t8 * qJD(4) - t50;
t5 = t7 * qJD(4) - t49;
t2 = -t4 * qJD(4) - t48;
t1 = t3 * qJD(4) - t66;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t54, (t36 * t54 - t53) * qJD(3), -t15 * t30 - t22 * t24, t52 + t63, -t57, t16 * t51 + t20 * t25, -t17, 0, t6 * qJD(4) + t31 * t16 + t18 * t25, -t5 * qJD(4) - t31 * t15 - t18 * t24, -t1 * t51 + t7 * t15 - t8 * t16 - t2 * t30 - t5 * t20 - t6 * t22 + t3 * t24 - t4 * t25, t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t57, -t52 + t63, t1 * t30 - t2 * t51 - t4 * t24 - t3 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59 * qJD(1) ^ 2, t53 * qJD(1), 0, 0, 0, 0, 0, 0, 0.2e1 * t22 * qJD(4), t33 + (-t20 - t56) * qJD(4), -t19 - t65, t4 * t20 + t3 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t19 + t65, t33 + (t20 - t56) * qJD(4), -t62, 0, 0, -t18 * t22 - t48, t18 * t20 + t66, 0, 0;];
tauc_reg = t9;
