% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:53
% EndTime: 2019-12-05 17:44:53
% DurationCPUTime: 0.19s
% Computational Cost: add. (382->51), mult. (1068->128), div. (0->0), fcn. (709->8), ass. (0->45)
t45 = qJD(1) ^ 2;
t57 = t45 / 0.2e1;
t56 = cos(qJ(5));
t39 = sin(pkin(8));
t36 = t39 ^ 2;
t55 = t36 * t45;
t40 = cos(pkin(9));
t54 = t39 * t40;
t41 = cos(pkin(8));
t24 = qJD(2) + (-pkin(2) * t41 - qJ(3) * t39 - pkin(1)) * qJD(1);
t23 = t40 * t24;
t38 = sin(pkin(9));
t12 = t23 + (-pkin(6) * t54 + (-qJ(2) * t38 - pkin(3)) * t41) * qJD(1);
t51 = qJ(2) * qJD(1);
t47 = t41 * t51;
t17 = t38 * t24 + t40 * t47;
t53 = qJD(1) * t39;
t49 = t38 * t53;
t15 = -pkin(6) * t49 + t17;
t43 = sin(qJ(4));
t44 = cos(qJ(4));
t6 = t43 * t12 + t44 * t15;
t52 = t41 * qJD(1);
t30 = t39 * t51 + qJD(3);
t50 = t39 * t45 * t41;
t48 = t55 / 0.2e1;
t25 = pkin(3) * t49 + t30;
t5 = t44 * t12 - t43 * t15;
t31 = -qJD(4) + t52;
t42 = sin(qJ(5));
t37 = t41 ^ 2;
t35 = -qJD(1) * pkin(1) + qJD(2);
t33 = t37 * t57;
t28 = -qJD(5) + t31;
t21 = (-t38 * t43 + t40 * t44) * t53;
t19 = (-t38 * t44 - t40 * t43) * t53;
t16 = -t38 * t47 + t23;
t14 = -t19 * pkin(4) + t25;
t9 = t42 * t19 + t56 * t21;
t7 = -t56 * t19 + t42 * t21;
t4 = t19 * pkin(7) + t6;
t3 = -t31 * pkin(4) - t21 * pkin(7) + t5;
t2 = t42 * t3 + t56 * t4;
t1 = t56 * t3 - t42 * t4;
t8 = [0, 0, 0, 0, 0, t57, 0, 0, 0, 0, t48, t50, 0, t33, 0, 0, -t35 * t52, t35 * t53, (t36 + t37) * t45 * qJ(2), t35 ^ 2 / 0.2e1 + (t37 / 0.2e1 + t36 / 0.2e1) * qJ(2) ^ 2 * t45, t40 ^ 2 * t48, -t40 * t38 * t55, -t40 * t50, t38 ^ 2 * t48, t38 * t50, t33, (t30 * t38 * t39 - t16 * t41) * qJD(1), (t17 * t41 + t30 * t54) * qJD(1), (-t16 * t40 - t17 * t38) * t53, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, t21 * t19, -t21 * t31, t19 ^ 2 / 0.2e1, -t19 * t31, t31 ^ 2 / 0.2e1, -t25 * t19 - t5 * t31, t25 * t21 + t6 * t31, t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, -t9 * t28, t7 ^ 2 / 0.2e1, t7 * t28, t28 ^ 2 / 0.2e1, -t1 * t28 + t14 * t7, t14 * t9 + t2 * t28, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1;];
T_reg = t8;
