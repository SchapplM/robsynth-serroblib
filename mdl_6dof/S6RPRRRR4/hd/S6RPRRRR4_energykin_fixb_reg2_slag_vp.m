% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:05:52
% EndTime: 2019-03-09 07:05:53
% DurationCPUTime: 0.22s
% Computational Cost: add. (1126->67), mult. (2960->155), div. (0->0), fcn. (2307->10), ass. (0->50)
t57 = qJD(1) ^ 2;
t66 = t57 / 0.2e1;
t51 = sin(pkin(11));
t60 = qJD(1) * t51;
t61 = pkin(7) + qJ(2);
t41 = t61 * t60;
t52 = cos(pkin(11));
t59 = qJD(1) * t52;
t42 = t61 * t59;
t56 = sin(qJ(3));
t65 = cos(qJ(3));
t31 = -t65 * t41 - t56 * t42;
t40 = (t65 * t51 + t52 * t56) * qJD(1);
t25 = qJD(3) * pkin(3) - t40 * pkin(8) + t31;
t32 = -t56 * t41 + t65 * t42;
t38 = t56 * t60 - t65 * t59;
t26 = -t38 * pkin(8) + t32;
t55 = sin(qJ(4));
t64 = cos(qJ(4));
t13 = t55 * t25 + t64 * t26;
t28 = t64 * t38 + t55 * t40;
t11 = -t28 * pkin(9) + t13;
t54 = sin(qJ(5));
t63 = cos(qJ(5));
t12 = t64 * t25 - t55 * t26;
t30 = -t55 * t38 + t64 * t40;
t50 = qJD(3) + qJD(4);
t9 = t50 * pkin(4) - t30 * pkin(9) + t12;
t6 = t63 * t11 + t54 * t9;
t62 = cos(qJ(6));
t18 = t63 * t28 + t54 * t30;
t5 = -t54 * t11 + t63 * t9;
t43 = qJD(2) + (-pkin(2) * t52 - pkin(1)) * qJD(1);
t33 = t38 * pkin(3) + t43;
t21 = t28 * pkin(4) + t33;
t53 = sin(qJ(6));
t49 = t52 ^ 2;
t48 = t51 ^ 2;
t47 = qJD(5) + t50;
t46 = -qJD(1) * pkin(1) + qJD(2);
t20 = -t54 * t28 + t63 * t30;
t17 = qJD(6) + t18;
t16 = t62 * t20 + t53 * t47;
t14 = t53 * t20 - t62 * t47;
t7 = t18 * pkin(5) - t20 * pkin(10) + t21;
t4 = t47 * pkin(10) + t6;
t3 = -t47 * pkin(5) - t5;
t2 = t62 * t4 + t53 * t7;
t1 = -t53 * t4 + t62 * t7;
t8 = [0, 0, 0, 0, 0, t66, 0, 0, 0, 0, t48 * t66, t51 * t57 * t52, 0, t49 * t66, 0, 0, -t46 * t59, t46 * t60 (t48 + t49) * t57 * qJ(2), t46 ^ 2 / 0.2e1 + (t49 / 0.2e1 + t48 / 0.2e1) * qJ(2) ^ 2 * t57, t40 ^ 2 / 0.2e1, -t40 * t38, t40 * qJD(3), t38 ^ 2 / 0.2e1, -t38 * qJD(3), qJD(3) ^ 2 / 0.2e1, t31 * qJD(3) + t43 * t38, -t32 * qJD(3) + t43 * t40, -t31 * t40 - t32 * t38, t32 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t50 * t30, t28 ^ 2 / 0.2e1, -t50 * t28, t50 ^ 2 / 0.2e1, t12 * t50 + t33 * t28, -t13 * t50 + t33 * t30, -t12 * t30 - t13 * t28, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t47, t18 ^ 2 / 0.2e1, -t18 * t47, t47 ^ 2 / 0.2e1, t21 * t18 + t5 * t47, t21 * t20 - t6 * t47, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t17, t14 ^ 2 / 0.2e1, -t14 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t14, t3 * t16 - t2 * t17, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
