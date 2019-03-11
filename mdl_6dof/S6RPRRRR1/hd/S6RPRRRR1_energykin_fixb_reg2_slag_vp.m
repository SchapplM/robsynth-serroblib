% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRR1
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:15
% EndTime: 2019-03-09 06:55:15
% DurationCPUTime: 0.16s
% Computational Cost: add. (689->61), mult. (1534->150), div. (0->0), fcn. (1047->10), ass. (0->47)
t51 = qJD(1) ^ 2;
t43 = t51 / 0.2e1;
t60 = pkin(1) * t51;
t44 = sin(pkin(11));
t35 = (pkin(1) * t44 + pkin(7)) * qJD(1);
t50 = cos(qJ(3));
t41 = t50 * qJD(2);
t49 = sin(qJ(3));
t25 = qJD(3) * pkin(3) + t41 + (-pkin(8) * qJD(1) - t35) * t49;
t29 = t49 * qJD(2) + t50 * t35;
t55 = qJD(1) * t50;
t26 = pkin(8) * t55 + t29;
t48 = sin(qJ(4));
t59 = cos(qJ(4));
t13 = t48 * t25 + t59 * t26;
t56 = qJD(1) * t49;
t30 = t48 * t56 - t59 * t55;
t11 = -t30 * pkin(9) + t13;
t47 = sin(qJ(5));
t58 = cos(qJ(5));
t12 = t59 * t25 - t48 * t26;
t32 = (t48 * t50 + t59 * t49) * qJD(1);
t42 = qJD(3) + qJD(4);
t9 = t42 * pkin(4) - t32 * pkin(9) + t12;
t6 = t58 * t11 + t47 * t9;
t57 = cos(qJ(6));
t54 = qJD(1) * qJD(3);
t45 = cos(pkin(11));
t53 = -pkin(1) * t45 - pkin(2);
t18 = t58 * t30 + t47 * t32;
t5 = -t47 * t11 + t58 * t9;
t33 = (-pkin(3) * t50 + t53) * qJD(1);
t21 = t30 * pkin(4) + t33;
t46 = sin(qJ(6));
t39 = qJD(5) + t42;
t36 = t53 * qJD(1);
t28 = -t49 * t35 + t41;
t20 = -t47 * t30 + t58 * t32;
t17 = qJD(6) + t18;
t16 = t57 * t20 + t46 * t39;
t14 = t46 * t20 - t57 * t39;
t7 = t18 * pkin(5) - t20 * pkin(10) + t21;
t4 = t39 * pkin(10) + t6;
t3 = -t39 * pkin(5) - t5;
t2 = t57 * t4 + t46 * t7;
t1 = -t46 * t4 + t57 * t7;
t8 = [0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t45 * t60, -t44 * t60, 0, qJD(2) ^ 2 / 0.2e1 + (t44 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t51, t49 ^ 2 * t43, t49 * t51 * t50, t49 * t54, t50 ^ 2 * t43, t50 * t54, qJD(3) ^ 2 / 0.2e1, t28 * qJD(3) - t36 * t55, -t29 * qJD(3) + t36 * t56 (-t28 * t49 + t29 * t50) * qJD(1), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t30, t32 * t42, t30 ^ 2 / 0.2e1, -t30 * t42, t42 ^ 2 / 0.2e1, t12 * t42 + t33 * t30, -t13 * t42 + t33 * t32, -t12 * t32 - t13 * t30, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t39, t18 ^ 2 / 0.2e1, -t18 * t39, t39 ^ 2 / 0.2e1, t21 * t18 + t5 * t39, t21 * t20 - t6 * t39, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t17, t14 ^ 2 / 0.2e1, -t14 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t14, t3 * t16 - t2 * t17, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
