% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:39
% EndTime: 2019-12-31 16:31:40
% DurationCPUTime: 0.23s
% Computational Cost: add. (236->49), mult. (535->89), div. (0->0), fcn. (255->4), ass. (0->51)
t23 = cos(qJ(4));
t24 = cos(qJ(3));
t47 = pkin(2) * qJD(3);
t36 = qJD(2) * t47;
t34 = t24 * t36;
t18 = qJD(2) + qJD(3);
t22 = sin(qJ(3));
t48 = pkin(2) * qJD(2);
t11 = t18 * pkin(6) + t22 * t48;
t21 = sin(qJ(4));
t4 = t23 * qJD(1) - t21 * t11;
t2 = t4 * qJD(4) + t23 * t34;
t5 = t21 * qJD(1) + t23 * t11;
t3 = -t5 * qJD(4) - t21 * t34;
t57 = t2 * t23 - t3 * t21 + (-t21 * t5 - t23 * t4) * qJD(4);
t19 = t21 ^ 2;
t20 = t23 ^ 2;
t56 = t18 * (t19 + t20);
t33 = t21 * t4 - t23 * t5;
t55 = t33 * t24;
t40 = t24 * t48;
t12 = -t18 * pkin(3) - t40;
t45 = qJD(4) * t23;
t53 = t21 * t22;
t54 = t12 * t45 + t36 * t53;
t52 = t22 * t23;
t25 = qJD(4) ^ 2;
t51 = t25 * t21;
t16 = t25 * t23;
t50 = t19 - t20;
t46 = qJD(4) * t21;
t44 = qJD(4) * t24;
t43 = -qJD(2) - t18;
t42 = -qJD(3) + t18;
t17 = t18 ^ 2;
t41 = t21 * t17 * t23;
t39 = t18 * t46;
t38 = t18 * t45;
t37 = t21 * t44;
t35 = t21 * t38;
t32 = t42 * t48;
t31 = t43 * t47;
t30 = -t18 * t53 + t23 * t44;
t29 = -t12 * t18 - t34;
t15 = -t24 * pkin(2) - pkin(3);
t14 = t22 * pkin(2) + pkin(6);
t10 = -0.2e1 * t35;
t9 = 0.2e1 * t35;
t7 = t12 * t46;
t6 = -0.2e1 * t50 * t18 * qJD(4);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t16, 0, -t33 * qJD(4) + t2 * t21 + t3 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t31, t24 * t31, 0, 0, t9, t6, t16, t10, -t51, 0, t15 * t39 - t14 * t16 + t7 + (t43 * t52 - t37) * t47, t14 * t51 + t15 * t38 - t30 * t47 + t54, t24 * t47 * t56 + t57, t57 * t14 + (-t55 + (qJD(2) * t15 + t12) * t22) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t32, t24 * t32, 0, 0, t9, t6, t16, t10, -t51, 0, -pkin(3) * t39 - pkin(6) * t16 + t7 + (t42 * t52 + t37) * t48, -pkin(3) * t38 + pkin(6) * t51 + t30 * t48 + t54, -t40 * t56 + t57, t57 * pkin(6) + (t55 + (-pkin(3) * qJD(3) - t12) * t22) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t50 * t17, 0, t41, 0, 0, t29 * t21, t29 * t23, 0, 0;];
tauc_reg = t1;
