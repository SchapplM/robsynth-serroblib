% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPPRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:41
% EndTime: 2019-03-09 01:37:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (245->39), mult. (395->96), div. (0->0), fcn. (155->6), ass. (0->30)
t34 = qJD(1) ^ 2;
t27 = t34 / 0.2e1;
t23 = qJ(2) * qJD(1) + qJD(3);
t19 = pkin(3) * qJD(1) + t23;
t20 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t29 = sin(pkin(9));
t30 = cos(pkin(9));
t12 = t29 * t19 + t30 * t20;
t10 = qJD(1) * pkin(7) + t12;
t32 = sin(qJ(5));
t33 = cos(qJ(5));
t7 = t32 * qJD(4) + t33 * t10;
t38 = cos(qJ(6));
t37 = qJD(1) * t32;
t36 = t33 * qJD(1);
t35 = qJD(1) * qJD(5);
t11 = t30 * t19 - t29 * t20;
t6 = t33 * qJD(4) - t32 * t10;
t31 = sin(qJ(6));
t24 = -pkin(1) * qJD(1) + qJD(2);
t21 = -qJD(6) + t36;
t15 = t31 * qJD(5) + t38 * t37;
t13 = -t38 * qJD(5) + t31 * t37;
t9 = -qJD(1) * pkin(4) - t11;
t5 = (-pkin(5) * t33 - pkin(8) * t32 - pkin(4)) * qJD(1) - t11;
t4 = qJD(5) * pkin(8) + t7;
t3 = -qJD(5) * pkin(5) - t6;
t2 = t31 * t5 + t38 * t4;
t1 = -t31 * t4 + t38 * t5;
t8 = [0, 0, 0, 0, 0, t27, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, t24 * qJD(1), t34 * qJ(2), qJ(2) ^ 2 * t27 + t24 ^ 2 / 0.2e1, 0, 0, 0, t27, 0, 0, t23 * qJD(1), 0, -t20 * qJD(1), t20 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t27, t11 * qJD(1), -t12 * qJD(1), 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t32 ^ 2 * t27, t32 * t34 * t33, t32 * t35, t33 ^ 2 * t27, t33 * t35, qJD(5) ^ 2 / 0.2e1, t6 * qJD(5) - t9 * t36, -t7 * qJD(5) + t9 * t37 (-t32 * t6 + t33 * t7) * qJD(1), t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, -t15 * t21, t13 ^ 2 / 0.2e1, t13 * t21, t21 ^ 2 / 0.2e1, -t1 * t21 + t3 * t13, t3 * t15 + t2 * t21, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
