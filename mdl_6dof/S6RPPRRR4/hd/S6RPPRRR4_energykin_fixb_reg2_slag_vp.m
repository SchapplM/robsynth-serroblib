% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:46
% EndTime: 2019-03-09 02:26:46
% DurationCPUTime: 0.14s
% Computational Cost: add. (451->53), mult. (822->122), div. (0->0), fcn. (435->8), ass. (0->39)
t46 = qJD(1) ^ 2;
t38 = t46 / 0.2e1;
t31 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t39 = sin(pkin(10));
t40 = cos(pkin(10));
t48 = qJ(2) * qJD(1);
t24 = t39 * t31 + t40 * t48;
t22 = -qJD(1) * pkin(7) + t24;
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t18 = t43 * qJD(3) + t45 * t22;
t12 = qJD(4) * pkin(8) + t18;
t23 = t40 * t31 - t39 * t48;
t21 = qJD(1) * pkin(3) - t23;
t13 = (pkin(4) * t45 + pkin(8) * t43) * qJD(1) + t21;
t44 = cos(qJ(5));
t51 = sin(qJ(5));
t6 = t44 * t12 + t51 * t13;
t52 = cos(qJ(6));
t50 = qJD(1) * t43;
t49 = t45 * qJD(1);
t47 = qJD(1) * qJD(4);
t17 = t45 * qJD(3) - t43 * t22;
t32 = qJD(5) + t49;
t5 = -t51 * t12 + t44 * t13;
t11 = -qJD(4) * pkin(4) - t17;
t42 = sin(qJ(6));
t35 = -pkin(1) * qJD(1) + qJD(2);
t30 = qJD(6) + t32;
t27 = -t51 * qJD(4) + t44 * t50;
t26 = t44 * qJD(4) + t51 * t50;
t16 = t42 * t26 - t52 * t27;
t14 = -t52 * t26 - t42 * t27;
t7 = -t26 * pkin(5) + t11;
t4 = t26 * pkin(9) + t6;
t3 = t32 * pkin(5) + t27 * pkin(9) + t5;
t2 = t42 * t3 + t52 * t4;
t1 = t52 * t3 - t42 * t4;
t8 = [0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, -t35 * qJD(1), 0, t46 * qJ(2), qJ(2) ^ 2 * t38 + t35 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t38, -t23 * qJD(1), t24 * qJD(1), 0, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t43 ^ 2 * t38, t43 * t46 * t45, -t43 * t47, t45 ^ 2 * t38, -t45 * t47, qJD(4) ^ 2 / 0.2e1, t17 * qJD(4) + t21 * t49, -t18 * qJD(4) - t21 * t50 (t17 * t43 - t18 * t45) * qJD(1), t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t26, -t27 * t32, t26 ^ 2 / 0.2e1, t26 * t32, t32 ^ 2 / 0.2e1, -t11 * t26 + t5 * t32, -t11 * t27 - t6 * t32, t6 * t26 + t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t30, t14 ^ 2 / 0.2e1, -t14 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t14, t7 * t16 - t2 * t30, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
