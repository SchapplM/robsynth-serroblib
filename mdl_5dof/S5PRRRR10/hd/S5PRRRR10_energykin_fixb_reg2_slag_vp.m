% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:04
% EndTime: 2019-12-05 17:26:05
% DurationCPUTime: 0.19s
% Computational Cost: add. (412->48), mult. (1070->127), div. (0->0), fcn. (832->12), ass. (0->47)
t46 = cos(qJ(2));
t38 = sin(pkin(5));
t57 = qJD(1) * t38;
t28 = qJD(2) * pkin(2) + t46 * t57;
t37 = sin(pkin(6));
t39 = cos(pkin(6));
t40 = cos(pkin(5));
t56 = qJD(1) * t40;
t62 = t28 * t39 + t37 * t56;
t44 = sin(qJ(2));
t55 = qJD(2) * t37;
t26 = pkin(8) * t55 + t44 * t57;
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t11 = -t43 * t26 + t62 * t45;
t12 = t45 * t26 + t62 * t43;
t34 = t39 * qJD(2) + qJD(3);
t10 = t34 * pkin(9) + t12;
t33 = t39 * t56;
t14 = t33 + (-t28 + (-pkin(3) * t45 - pkin(9) * t43) * qJD(2)) * t37;
t42 = sin(qJ(4));
t61 = cos(qJ(4));
t6 = t61 * t10 + t42 * t14;
t60 = cos(qJ(5));
t47 = qJD(2) ^ 2;
t58 = t37 ^ 2 * t47;
t53 = t43 * t55;
t52 = t45 * t55;
t51 = t58 / 0.2e1;
t50 = qJD(2) * t57;
t20 = -t61 * t34 + t42 * t53;
t5 = -t42 * t10 + t61 * t14;
t9 = -t34 * pkin(3) - t11;
t48 = qJD(1) ^ 2;
t41 = sin(qJ(5));
t31 = -qJD(4) + t52;
t22 = t42 * t34 + t61 * t53;
t19 = qJD(5) + t20;
t18 = -t37 * t28 + t33;
t17 = t60 * t22 - t41 * t31;
t15 = t41 * t22 + t60 * t31;
t7 = t20 * pkin(4) - t22 * pkin(10) + t9;
t4 = -t31 * pkin(10) + t6;
t3 = t31 * pkin(4) - t5;
t2 = t60 * t4 + t41 * t7;
t1 = -t41 * t4 + t60 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t48 / 0.2e1, 0, 0, 0, 0, 0, t47 / 0.2e1, t46 * t50, -t44 * t50, 0, (t40 ^ 2 / 0.2e1 + (t44 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1) * t38 ^ 2) * t48, t43 ^ 2 * t51, t43 * t45 * t58, t34 * t53, t45 ^ 2 * t51, t34 * t52, t34 ^ 2 / 0.2e1, t11 * t34 - t18 * t52, -t12 * t34 + t18 * t53, (-t11 * t43 + t12 * t45) * t55, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, -t22 * t31, t20 ^ 2 / 0.2e1, t20 * t31, t31 ^ 2 / 0.2e1, t9 * t20 - t5 * t31, t9 * t22 + t6 * t31, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t19, t15 ^ 2 / 0.2e1, -t15 * t19, t19 ^ 2 / 0.2e1, t1 * t19 + t3 * t15, t3 * t17 - t2 * t19, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
