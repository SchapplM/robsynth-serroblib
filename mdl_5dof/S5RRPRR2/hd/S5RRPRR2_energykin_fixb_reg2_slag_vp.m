% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:26
% EndTime: 2019-12-05 18:28:26
% DurationCPUTime: 0.18s
% Computational Cost: add. (529->52), mult. (1414->127), div. (0->0), fcn. (1001->8), ass. (0->44)
t44 = qJD(1) ^ 2;
t56 = t44 / 0.2e1;
t55 = cos(qJ(4));
t54 = cos(qJ(5));
t43 = cos(qJ(2));
t53 = t43 * t44;
t52 = pkin(6) + qJ(3);
t42 = sin(qJ(2));
t50 = qJD(1) * t42;
t30 = qJD(2) * pkin(2) - t52 * t50;
t49 = qJD(1) * t43;
t31 = t52 * t49;
t39 = sin(pkin(9));
t51 = cos(pkin(9));
t20 = t51 * t30 - t39 * t31;
t28 = (t39 * t43 + t51 * t42) * qJD(1);
t14 = qJD(2) * pkin(3) - t28 * pkin(7) + t20;
t21 = t39 * t30 + t51 * t31;
t26 = t39 * t50 - t51 * t49;
t15 = -t26 * pkin(7) + t21;
t41 = sin(qJ(4));
t6 = t41 * t14 + t55 * t15;
t48 = qJD(1) * qJD(2);
t35 = qJD(2) + qJD(4);
t47 = t42 * t48;
t46 = t43 * t48;
t5 = t55 * t14 - t41 * t15;
t32 = qJD(3) + (-pkin(2) * t43 - pkin(1)) * qJD(1);
t22 = t26 * pkin(3) + t32;
t40 = sin(qJ(5));
t38 = t43 ^ 2;
t37 = t42 ^ 2;
t36 = qJD(2) ^ 2 / 0.2e1;
t34 = qJD(5) + t35;
t19 = -t41 * t26 + t55 * t28;
t17 = t55 * t26 + t41 * t28;
t10 = t17 * pkin(4) + t22;
t9 = -t40 * t17 + t54 * t19;
t7 = t54 * t17 + t40 * t19;
t4 = -t17 * pkin(8) + t6;
t3 = t35 * pkin(4) - t19 * pkin(8) + t5;
t2 = t40 * t3 + t54 * t4;
t1 = t54 * t3 - t40 * t4;
t8 = [0, 0, 0, 0, 0, t56, 0, 0, 0, 0, t37 * t56, t42 * t53, t47, t38 * t56, t46, t36, pkin(1) * t53 - pkin(6) * t47, -t44 * pkin(1) * t42 - pkin(6) * t46, (t37 + t38) * t44 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t38 / 0.2e1 + t37 / 0.2e1) * pkin(6) ^ 2) * t44, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * qJD(2), t26 ^ 2 / 0.2e1, -t26 * qJD(2), t36, t20 * qJD(2) + t32 * t26, -t21 * qJD(2) + t32 * t28, -t20 * t28 - t21 * t26, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * t35, t17 ^ 2 / 0.2e1, -t17 * t35, t35 ^ 2 / 0.2e1, t22 * t17 + t5 * t35, t22 * t19 - t6 * t35, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t34, t7 ^ 2 / 0.2e1, -t7 * t34, t34 ^ 2 / 0.2e1, t1 * t34 + t10 * t7, t10 * t9 - t2 * t34, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
