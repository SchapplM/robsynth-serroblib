% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:05
% EndTime: 2019-03-09 00:11:06
% DurationCPUTime: 0.15s
% Computational Cost: add. (514->57), mult. (1170->132), div. (0->0), fcn. (836->10), ass. (0->50)
t55 = qJD(2) ^ 2;
t65 = t55 / 0.2e1;
t51 = sin(qJ(2));
t46 = sin(pkin(6));
t63 = qJD(1) * t46;
t36 = qJD(2) * pkin(8) + t51 * t63;
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t47 = cos(pkin(6));
t62 = qJD(1) * t47;
t28 = t53 * t36 + t50 * t62;
t24 = qJD(3) * pkin(9) + t28;
t54 = cos(qJ(2));
t58 = t54 * t63;
t29 = -t58 + (-pkin(3) * t53 - pkin(9) * t50 - pkin(2)) * qJD(2);
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t13 = t52 * t24 + t49 * t29;
t61 = qJD(2) * t50;
t33 = -t52 * qJD(3) + t49 * t61;
t10 = -t33 * pkin(10) + t13;
t48 = sin(qJ(5));
t64 = cos(qJ(5));
t12 = -t49 * t24 + t52 * t29;
t35 = t49 * qJD(3) + t52 * t61;
t60 = t53 * qJD(2);
t42 = -qJD(4) + t60;
t7 = -t42 * pkin(4) - t35 * pkin(10) + t12;
t4 = t64 * t10 + t48 * t7;
t59 = qJD(2) * qJD(3);
t3 = -t48 * t10 + t64 * t7;
t57 = qJD(2) * t63;
t27 = -t50 * t36 + t53 * t62;
t23 = -qJD(3) * pkin(3) - t27;
t16 = t33 * pkin(4) + t23;
t56 = qJD(1) ^ 2;
t39 = -qJD(5) + t42;
t38 = t39 ^ 2 / 0.2e1;
t37 = -qJD(2) * pkin(2) - t58;
t21 = -t48 * t33 + t64 * t35;
t19 = t64 * t33 + t48 * t35;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t15 = t21 * t39;
t14 = t19 * t39;
t11 = t21 * t19;
t8 = t19 * pkin(5) + qJD(6) + t16;
t2 = -t19 * qJ(6) + t4;
t1 = -t39 * pkin(5) - t21 * qJ(6) + t3;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t56 / 0.2e1, 0, 0, 0, 0, 0, t65, t54 * t57, -t51 * t57, 0 (t47 ^ 2 / 0.2e1 + (t51 ^ 2 / 0.2e1 + t54 ^ 2 / 0.2e1) * t46 ^ 2) * t56, t50 ^ 2 * t65, t50 * t55 * t53, t50 * t59, t53 ^ 2 * t65, t53 * t59, qJD(3) ^ 2 / 0.2e1, t27 * qJD(3) - t37 * t60, -t28 * qJD(3) + t37 * t61 (-t27 * t50 + t28 * t53) * qJD(2), t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, -t35 * t42, t33 ^ 2 / 0.2e1, t33 * t42, t42 ^ 2 / 0.2e1, -t12 * t42 + t23 * t33, t13 * t42 + t23 * t35, -t12 * t35 - t13 * t33, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t18, -t11, -t15, t17, t14, t38, t16 * t19 - t3 * t39, t16 * t21 + t4 * t39, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t18, -t11, -t15, t17, t14, t38, -t1 * t39 + t8 * t19, t2 * t39 + t8 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg  = t5;
