% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:50
% EndTime: 2019-12-05 18:38:50
% DurationCPUTime: 0.18s
% Computational Cost: add. (549->52), mult. (1414->127), div. (0->0), fcn. (1001->8), ass. (0->44)
t44 = qJD(1) ^ 2;
t56 = t44 / 0.2e1;
t55 = -pkin(7) - pkin(6);
t54 = cos(qJ(3));
t53 = cos(qJ(5));
t43 = cos(qJ(2));
t52 = t43 * t44;
t42 = sin(qJ(2));
t50 = qJD(1) * t42;
t30 = qJD(2) * pkin(2) + t55 * t50;
t49 = qJD(1) * t43;
t31 = t55 * t49;
t41 = sin(qJ(3));
t20 = t54 * t30 + t41 * t31;
t28 = (t41 * t43 + t54 * t42) * qJD(1);
t36 = qJD(2) + qJD(3);
t13 = t36 * pkin(3) - t28 * qJ(4) + t20;
t21 = t41 * t30 - t54 * t31;
t26 = t41 * t50 - t54 * t49;
t15 = -t26 * qJ(4) + t21;
t39 = sin(pkin(9));
t51 = cos(pkin(9));
t6 = t39 * t13 + t51 * t15;
t48 = qJD(1) * qJD(2);
t47 = t42 * t48;
t46 = t43 * t48;
t5 = t51 * t13 - t39 * t15;
t32 = (-pkin(2) * t43 - pkin(1)) * qJD(1);
t22 = t26 * pkin(3) + qJD(4) + t32;
t40 = sin(qJ(5));
t38 = t43 ^ 2;
t37 = t42 ^ 2;
t35 = qJD(5) + t36;
t34 = t36 ^ 2 / 0.2e1;
t19 = -t39 * t26 + t51 * t28;
t17 = t51 * t26 + t39 * t28;
t10 = t17 * pkin(4) + t22;
t9 = -t40 * t17 + t53 * t19;
t7 = t53 * t17 + t40 * t19;
t4 = -t17 * pkin(8) + t6;
t3 = t36 * pkin(4) - t19 * pkin(8) + t5;
t2 = t40 * t3 + t53 * t4;
t1 = t53 * t3 - t40 * t4;
t8 = [0, 0, 0, 0, 0, t56, 0, 0, 0, 0, t37 * t56, t42 * t52, t47, t38 * t56, t46, qJD(2) ^ 2 / 0.2e1, pkin(1) * t52 - pkin(6) * t47, -t44 * pkin(1) * t42 - pkin(6) * t46, (t37 + t38) * t44 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t38 / 0.2e1 + t37 / 0.2e1) * pkin(6) ^ 2) * t44, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * t36, t26 ^ 2 / 0.2e1, -t26 * t36, t34, t20 * t36 + t32 * t26, -t21 * t36 + t32 * t28, -t20 * t28 - t21 * t26, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * t36, t17 ^ 2 / 0.2e1, -t17 * t36, t34, t22 * t17 + t5 * t36, t22 * t19 - t6 * t36, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t35, t7 ^ 2 / 0.2e1, -t7 * t35, t35 ^ 2 / 0.2e1, t1 * t35 + t10 * t7, t10 * t9 - t2 * t35, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
