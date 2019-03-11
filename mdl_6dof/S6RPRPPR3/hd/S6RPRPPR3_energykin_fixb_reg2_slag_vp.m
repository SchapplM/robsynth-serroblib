% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:23
% EndTime: 2019-03-09 02:45:23
% DurationCPUTime: 0.15s
% Computational Cost: add. (235->52), mult. (530->113), div. (0->0), fcn. (237->6), ass. (0->41)
t42 = qJD(1) ^ 2;
t33 = t42 / 0.2e1;
t51 = -pkin(3) - pkin(4);
t50 = pkin(1) * t42;
t35 = sin(pkin(9));
t19 = (pkin(1) * t35 + pkin(7)) * qJD(1);
t39 = sin(qJ(3));
t41 = cos(qJ(3));
t13 = t39 * qJD(2) + t41 * t19;
t36 = cos(pkin(9));
t20 = (-pkin(1) * t36 - pkin(2)) * qJD(1);
t49 = qJD(1) * t41;
t48 = t39 * qJD(1);
t47 = qJ(5) * qJD(1);
t46 = qJD(1) * qJD(3);
t10 = qJD(3) * qJ(4) + t13;
t11 = -pkin(3) * t49 - qJ(4) * t48 + t20;
t12 = t41 * qJD(2) - t39 * t19;
t45 = qJD(4) - t12;
t7 = pkin(4) * t49 + qJD(5) - t11;
t8 = t41 * t47 - t10;
t44 = -t39 * t47 + t45;
t40 = cos(qJ(6));
t38 = sin(qJ(6));
t32 = qJD(3) ^ 2 / 0.2e1;
t26 = t41 * t46;
t25 = t39 * t46;
t24 = t41 ^ 2 * t33;
t23 = t39 ^ 2 * t33;
t22 = qJD(6) + t48;
t21 = t39 * t42 * t41;
t17 = t38 * qJD(3) + t40 * t49;
t16 = -t40 * qJD(3) + t38 * t49;
t9 = -qJD(3) * pkin(3) + t45;
t6 = qJD(3) * pkin(5) - t8;
t5 = t51 * qJD(3) + t44;
t4 = (-pkin(8) + t51) * qJD(3) + t44;
t3 = (pkin(5) * t39 + pkin(8) * t41) * qJD(1) + t7;
t2 = t38 * t3 + t40 * t4;
t1 = t40 * t3 - t38 * t4;
t14 = [0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t36 * t50, -t35 * t50, 0, qJD(2) ^ 2 / 0.2e1 + (t35 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t42, t23, t21, t25, t24, t26, t32, t12 * qJD(3) - t20 * t49, -t13 * qJD(3) + t20 * t48 (-t12 * t39 + t13 * t41) * qJD(1), t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t23, t25, -t21, t32, -t26, t24, -t9 * qJD(3) - t11 * t49 (t10 * t41 + t39 * t9) * qJD(1), t10 * qJD(3) - t11 * t48, t10 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t24, t21, t26, t23, t25, t32, -t8 * qJD(3) + t7 * t48, t5 * qJD(3) - t7 * t49 (-t39 * t5 + t41 * t8) * qJD(1), t5 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t16, -t17 * t22, t16 ^ 2 / 0.2e1, t16 * t22, t22 ^ 2 / 0.2e1, t1 * t22 - t6 * t16, -t6 * t17 - t2 * t22, t1 * t17 + t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg  = t14;
