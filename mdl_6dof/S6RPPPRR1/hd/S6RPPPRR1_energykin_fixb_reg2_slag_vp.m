% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:13
% EndTime: 2019-03-09 01:30:13
% DurationCPUTime: 0.12s
% Computational Cost: add. (166->41), mult. (330->92), div. (0->0), fcn. (135->6), ass. (0->33)
t30 = qJD(1) ^ 2;
t23 = t30 / 0.2e1;
t40 = pkin(1) * t30;
t28 = sin(qJ(5));
t29 = cos(qJ(5));
t25 = sin(pkin(9));
t17 = (-pkin(1) * t25 - qJ(3)) * qJD(1);
t15 = qJD(4) - t17;
t9 = -qJD(1) * pkin(7) + t15;
t7 = t29 * qJD(2) + t28 * t9;
t39 = cos(qJ(6));
t38 = qJD(1) * t29;
t26 = cos(pkin(9));
t33 = -pkin(1) * t26 - pkin(2);
t32 = -qJ(4) + t33;
t10 = -t32 * qJD(1) - qJD(3);
t37 = t10 * qJD(1);
t36 = t28 * qJD(1);
t35 = t10 ^ 2 / 0.2e1;
t34 = qJD(1) * qJD(5);
t6 = -t28 * qJD(2) + t29 * t9;
t27 = sin(qJ(6));
t22 = qJD(2) ^ 2 / 0.2e1;
t18 = qJD(6) + t36;
t16 = t33 * qJD(1) + qJD(3);
t14 = t27 * qJD(5) + t39 * t38;
t12 = -t39 * qJD(5) + t27 * t38;
t5 = -qJD(3) + (pkin(5) * t28 - pkin(8) * t29 - t32) * qJD(1);
t4 = qJD(5) * pkin(8) + t7;
t3 = -qJD(5) * pkin(5) - t6;
t2 = t27 * t5 + t39 * t4;
t1 = -t27 * t4 + t39 * t5;
t8 = [0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t26 * t40, -t25 * t40, 0, t22 + (t25 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t30, t23, 0, 0, 0, 0, 0, 0, t16 * qJD(1), -t17 * qJD(1), t22 + t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t23, 0, 0, 0, 0, 0, 0, t15 * qJD(1), t37, t22 + t35 + t15 ^ 2 / 0.2e1, t29 ^ 2 * t23, -t29 * t30 * t28, t29 * t34, t28 ^ 2 * t23, -t28 * t34, qJD(5) ^ 2 / 0.2e1, t6 * qJD(5) + t10 * t36, -t7 * qJD(5) + t29 * t37 (-t28 * t7 - t29 * t6) * qJD(1), t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t35, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t18, t12 ^ 2 / 0.2e1, -t12 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t3 * t12, t3 * t14 - t2 * t18, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
