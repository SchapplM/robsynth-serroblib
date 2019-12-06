% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:43
% EndTime: 2019-12-05 16:00:43
% DurationCPUTime: 0.11s
% Computational Cost: add. (124->29), mult. (260->81), div. (0->0), fcn. (125->6), ass. (0->31)
t24 = qJD(2) ^ 2;
t15 = t24 / 0.2e1;
t23 = cos(qJ(2));
t26 = -t23 * qJD(1) + qJD(3);
t10 = (-pkin(2) - pkin(6)) * qJD(2) + t26;
t33 = qJD(4) * t10;
t20 = sin(qJ(2));
t31 = t20 * qJD(1);
t12 = qJD(2) * qJ(3) + t31;
t32 = t12 * qJD(2);
t30 = t12 ^ 2 / 0.2e1;
t29 = qJD(1) * qJD(2);
t28 = qJD(2) * qJD(4);
t27 = -pkin(7) * qJD(2) + t10;
t25 = qJD(1) ^ 2;
t22 = cos(qJ(4));
t21 = cos(qJ(5));
t19 = sin(qJ(4));
t18 = sin(qJ(5));
t17 = t22 ^ 2;
t16 = t19 ^ 2;
t14 = qJD(4) + qJD(5);
t11 = -qJD(2) * pkin(2) + t26;
t8 = t31 + (pkin(4) * t19 + qJ(3)) * qJD(2);
t7 = (-t18 * t19 + t21 * t22) * qJD(2);
t5 = (-t18 * t22 - t19 * t21) * qJD(2);
t4 = t27 * t19;
t3 = qJD(4) * pkin(4) + t27 * t22;
t2 = t18 * t3 + t21 * t4;
t1 = -t18 * t4 + t21 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t25 / 0.2e1, 0, 0, 0, 0, 0, t15, t23 * t29, -t20 * t29, 0, (t20 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1) * t25, t15, 0, 0, 0, 0, 0, 0, t11 * qJD(2), t32, t30 + t11 ^ 2 / 0.2e1, t17 * t15, -t22 * t24 * t19, t22 * t28, t16 * t15, -t19 * t28, qJD(4) ^ 2 / 0.2e1, t19 * t32 + t22 * t33, -t19 * t33 + t22 * t32, (-t16 - t17) * t10 * qJD(2), t30 + (t16 / 0.2e1 + t17 / 0.2e1) * t10 ^ 2, t7 ^ 2 / 0.2e1, t7 * t5, t7 * t14, t5 ^ 2 / 0.2e1, t5 * t14, t14 ^ 2 / 0.2e1, t1 * t14 - t8 * t5, -t2 * t14 + t8 * t7, -t1 * t7 + t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg = t6;
