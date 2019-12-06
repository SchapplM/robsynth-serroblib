% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:14
% EndTime: 2019-12-05 17:38:14
% DurationCPUTime: 0.11s
% Computational Cost: add. (113->28), mult. (239->76), div. (0->0), fcn. (93->4), ass. (0->27)
t24 = qJD(1) ^ 2;
t16 = t24 / 0.2e1;
t30 = -pkin(1) - qJ(3);
t13 = qJD(1) * qJ(2) + qJD(3);
t10 = -qJD(1) * pkin(6) + t13;
t29 = qJD(4) * t10;
t11 = -t30 * qJD(1) - qJD(2);
t28 = t11 * qJD(1);
t27 = t11 ^ 2 / 0.2e1;
t26 = qJD(1) * qJD(4);
t25 = -pkin(7) * qJD(1) + t10;
t23 = cos(qJ(4));
t22 = cos(qJ(5));
t21 = sin(qJ(4));
t20 = sin(qJ(5));
t19 = t23 ^ 2;
t18 = t21 ^ 2;
t15 = qJD(4) + qJD(5);
t14 = -qJD(1) * pkin(1) + qJD(2);
t8 = -qJD(2) + (pkin(4) * t21 - t30) * qJD(1);
t7 = (-t20 * t21 + t22 * t23) * qJD(1);
t5 = (-t20 * t23 - t21 * t22) * qJD(1);
t4 = t25 * t21;
t3 = qJD(4) * pkin(4) + t25 * t23;
t2 = t20 * t3 + t22 * t4;
t1 = -t20 * t4 + t22 * t3;
t6 = [0, 0, 0, 0, 0, t16, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, t14 * qJD(1), t24 * qJ(2), qJ(2) ^ 2 * t16 + t14 ^ 2 / 0.2e1, t16, 0, 0, 0, 0, 0, 0, t13 * qJD(1), t28, t27 + t13 ^ 2 / 0.2e1, t19 * t16, -t23 * t24 * t21, t23 * t26, t18 * t16, -t21 * t26, qJD(4) ^ 2 / 0.2e1, t21 * t28 + t23 * t29, -t21 * t29 + t23 * t28, (-t18 - t19) * t10 * qJD(1), t27 + (t18 / 0.2e1 + t19 / 0.2e1) * t10 ^ 2, t7 ^ 2 / 0.2e1, t7 * t5, t7 * t15, t5 ^ 2 / 0.2e1, t5 * t15, t15 ^ 2 / 0.2e1, t1 * t15 - t8 * t5, -t2 * t15 + t8 * t7, -t1 * t7 + t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg = t6;
