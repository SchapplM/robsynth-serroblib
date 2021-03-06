% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPPRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:58
% EndTime: 2019-03-09 01:33:58
% DurationCPUTime: 0.15s
% Computational Cost: add. (396->48), mult. (765->115), div. (0->0), fcn. (425->8), ass. (0->38)
t47 = qJD(1) ^ 2;
t38 = t47 / 0.2e1;
t29 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t40 = sin(pkin(9));
t42 = cos(pkin(9));
t48 = qJ(2) * qJD(1);
t22 = t40 * t29 + t42 * t48;
t20 = -qJD(1) * qJ(4) + t22;
t41 = cos(pkin(10));
t36 = t41 * qJD(3);
t39 = sin(pkin(10));
t10 = t36 + (pkin(7) * qJD(1) - t20) * t39;
t13 = t39 * qJD(3) + t41 * t20;
t49 = qJD(1) * t41;
t11 = -pkin(7) * t49 + t13;
t45 = sin(qJ(5));
t46 = cos(qJ(5));
t6 = t45 * t10 + t46 * t11;
t51 = cos(qJ(6));
t50 = qJD(1) * t39;
t21 = t42 * t29 - t40 * t48;
t25 = t45 * t50 - t46 * t49;
t5 = t46 * t10 - t45 * t11;
t19 = qJD(1) * pkin(3) + qJD(4) - t21;
t14 = pkin(4) * t49 + t19;
t44 = sin(qJ(6));
t34 = -pkin(1) * qJD(1) + qJD(2);
t26 = (-t39 * t46 - t41 * t45) * qJD(1);
t23 = -qJD(6) + t25;
t18 = t44 * qJD(5) + t51 * t26;
t16 = -t51 * qJD(5) + t44 * t26;
t12 = -t39 * t20 + t36;
t7 = -t25 * pkin(5) - t26 * pkin(8) + t14;
t4 = qJD(5) * pkin(8) + t6;
t3 = -qJD(5) * pkin(5) - t5;
t2 = t51 * t4 + t44 * t7;
t1 = -t44 * t4 + t51 * t7;
t8 = [0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, -t34 * qJD(1), 0, t47 * qJ(2), qJ(2) ^ 2 * t38 + t34 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t38, -t21 * qJD(1), t22 * qJD(1), 0, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t39 ^ 2 * t38, t39 * t47 * t41, 0, t41 ^ 2 * t38, 0, 0, t19 * t49, -t19 * t50 (t12 * t39 - t13 * t41) * qJD(1), t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, t26 * t25, t26 * qJD(5), t25 ^ 2 / 0.2e1, t25 * qJD(5), qJD(5) ^ 2 / 0.2e1, t5 * qJD(5) - t14 * t25, -t6 * qJD(5) + t14 * t26, t6 * t25 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, -t18 * t23, t16 ^ 2 / 0.2e1, t16 * t23, t23 ^ 2 / 0.2e1, -t1 * t23 + t3 * t16, t3 * t18 + t2 * t23, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
