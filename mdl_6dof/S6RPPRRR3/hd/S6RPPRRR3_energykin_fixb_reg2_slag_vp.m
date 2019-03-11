% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:47
% EndTime: 2019-03-09 02:23:47
% DurationCPUTime: 0.16s
% Computational Cost: add. (351->52), mult. (700->118), div. (0->0), fcn. (384->8), ass. (0->41)
t41 = qJD(1) ^ 2;
t34 = t41 / 0.2e1;
t52 = pkin(1) * t41;
t51 = cos(qJ(5));
t50 = cos(qJ(6));
t36 = cos(pkin(10));
t44 = -pkin(1) * t36 - pkin(2);
t21 = qJD(3) + (-pkin(7) + t44) * qJD(1);
t39 = sin(qJ(4));
t40 = cos(qJ(4));
t17 = t40 * qJD(2) + t39 * t21;
t13 = qJD(4) * pkin(8) + t17;
t35 = sin(pkin(10));
t43 = -pkin(1) * t35 - qJ(3);
t18 = (pkin(4) * t39 - pkin(8) * t40 - t43) * qJD(1);
t38 = sin(qJ(5));
t6 = t51 * t13 + t38 * t18;
t49 = qJD(1) * t40;
t26 = t43 * qJD(1);
t48 = t26 * qJD(1);
t47 = t39 * qJD(1);
t46 = t26 ^ 2 / 0.2e1;
t45 = qJD(1) * qJD(4);
t5 = -t38 * t13 + t51 * t18;
t16 = -t39 * qJD(2) + t40 * t21;
t29 = qJD(5) + t47;
t12 = -qJD(4) * pkin(4) - t16;
t37 = sin(qJ(6));
t33 = qJD(2) ^ 2 / 0.2e1;
t28 = qJD(6) + t29;
t25 = t44 * qJD(1) + qJD(3);
t24 = t38 * qJD(4) + t51 * t49;
t22 = -t51 * qJD(4) + t38 * t49;
t10 = -t37 * t22 + t50 * t24;
t8 = t50 * t22 + t37 * t24;
t7 = t22 * pkin(5) + t12;
t4 = -t22 * pkin(9) + t6;
t3 = t29 * pkin(5) - t24 * pkin(9) + t5;
t2 = t37 * t3 + t50 * t4;
t1 = t50 * t3 - t37 * t4;
t9 = [0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t36 * t52, -t35 * t52, 0, t33 + (t35 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t41, t34, 0, 0, 0, 0, 0, 0, t25 * qJD(1), -t48, t33 + t46 + t25 ^ 2 / 0.2e1, t40 ^ 2 * t34, -t40 * t41 * t39, t40 * t45, t39 ^ 2 * t34, -t39 * t45, qJD(4) ^ 2 / 0.2e1, t16 * qJD(4) - t26 * t47, -t17 * qJD(4) - t40 * t48 (-t16 * t40 - t17 * t39) * qJD(1), t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t46, t24 ^ 2 / 0.2e1, -t24 * t22, t24 * t29, t22 ^ 2 / 0.2e1, -t22 * t29, t29 ^ 2 / 0.2e1, t12 * t22 + t5 * t29, t12 * t24 - t6 * t29, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t28, t8 ^ 2 / 0.2e1, -t8 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t8, t7 * t10 - t2 * t28, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t9;
