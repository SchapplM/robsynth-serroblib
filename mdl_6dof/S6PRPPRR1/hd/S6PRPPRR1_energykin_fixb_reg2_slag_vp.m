% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:13
% EndTime: 2019-03-08 19:16:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (380->51), mult. (911->126), div. (0->0), fcn. (669->12), ass. (0->46)
t48 = qJD(2) ^ 2;
t37 = t48 / 0.2e1;
t47 = cos(qJ(2));
t40 = sin(pkin(6));
t55 = qJD(1) * t40;
t28 = qJD(2) * pkin(2) + t47 * t55;
t39 = sin(pkin(11));
t42 = cos(pkin(11));
t46 = sin(qJ(2));
t52 = t46 * t55;
t19 = t39 * t28 + t42 * t52;
t17 = qJD(2) * qJ(4) + t19;
t43 = cos(pkin(6));
t34 = t43 * qJD(1) + qJD(3);
t41 = cos(pkin(12));
t30 = t41 * t34;
t38 = sin(pkin(12));
t10 = t30 + (-pkin(8) * qJD(2) - t17) * t38;
t13 = t41 * t17 + t38 * t34;
t53 = qJD(2) * t41;
t11 = pkin(8) * t53 + t13;
t45 = sin(qJ(5));
t57 = cos(qJ(5));
t6 = t45 * t10 + t57 * t11;
t56 = cos(qJ(6));
t54 = qJD(2) * t38;
t51 = qJD(2) * t55;
t18 = t42 * t28 - t39 * t52;
t24 = t45 * t54 - t57 * t53;
t50 = qJD(4) - t18;
t5 = t57 * t10 - t45 * t11;
t14 = (-pkin(4) * t41 - pkin(3)) * qJD(2) + t50;
t49 = qJD(1) ^ 2;
t44 = sin(qJ(6));
t26 = (t57 * t38 + t41 * t45) * qJD(2);
t23 = qJD(6) + t24;
t22 = t44 * qJD(5) + t56 * t26;
t20 = -t56 * qJD(5) + t44 * t26;
t16 = -qJD(2) * pkin(3) + t50;
t12 = -t38 * t17 + t30;
t7 = t24 * pkin(5) - t26 * pkin(9) + t14;
t4 = qJD(5) * pkin(9) + t6;
t3 = -qJD(5) * pkin(5) - t5;
t2 = t56 * t4 + t44 * t7;
t1 = -t44 * t4 + t56 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t49 / 0.2e1, 0, 0, 0, 0, 0, t37, t47 * t51, -t46 * t51, 0 (t43 ^ 2 / 0.2e1 + (t46 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1) * t40 ^ 2) * t49, 0, 0, 0, 0, 0, t37, t18 * qJD(2), -t19 * qJD(2), 0, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t38 ^ 2 * t37, t38 * t48 * t41, 0, t41 ^ 2 * t37, 0, 0, -t16 * t53, t16 * t54 (-t12 * t38 + t13 * t41) * qJD(2), t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, t26 * qJD(5), t24 ^ 2 / 0.2e1, -t24 * qJD(5), qJD(5) ^ 2 / 0.2e1, t5 * qJD(5) + t14 * t24, -t6 * qJD(5) + t14 * t26, -t6 * t24 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t23, t20 ^ 2 / 0.2e1, -t20 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t3 * t20, -t2 * t23 + t3 * t22, -t1 * t22 - t2 * t20, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
