% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:17
% EndTime: 2019-12-05 15:41:18
% DurationCPUTime: 0.10s
% Computational Cost: add. (79->26), mult. (178->64), div. (0->0), fcn. (62->4), ass. (0->30)
t20 = qJD(2) ^ 2;
t13 = t20 / 0.2e1;
t19 = cos(qJ(2));
t22 = -t19 * qJD(1) + qJD(3);
t5 = (-pkin(2) - pkin(6)) * qJD(2) + t22;
t31 = qJD(4) * t5;
t17 = sin(qJ(2));
t27 = t17 * qJD(1);
t7 = qJD(2) * qJ(3) + t27;
t30 = t7 ^ 2 / 0.2e1;
t16 = sin(qJ(4));
t29 = qJD(2) * t16;
t18 = cos(qJ(4));
t28 = qJD(2) * t18;
t26 = qJD(1) * qJD(2);
t25 = qJD(2) * qJD(4);
t24 = t18 * t20 * t16;
t23 = t16 * t25;
t21 = qJD(1) ^ 2;
t15 = t18 ^ 2;
t14 = t16 ^ 2;
t12 = qJD(4) ^ 2 / 0.2e1;
t11 = t18 * t25;
t10 = t15 * t13;
t9 = t14 * t13;
t6 = -qJD(2) * pkin(2) + t22;
t3 = qJD(4) * qJ(5) + t16 * t5;
t2 = t27 + (pkin(4) * t16 - qJ(5) * t18 + qJ(3)) * qJD(2);
t1 = -qJD(4) * pkin(4) - t18 * t5 + qJD(5);
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t21 / 0.2e1, 0, 0, 0, 0, 0, t13, t19 * t26, -t17 * t26, 0, (t17 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1) * t21, t13, 0, 0, 0, 0, 0, 0, t6 * qJD(2), t7 * qJD(2), t30 + t6 ^ 2 / 0.2e1, t10, -t24, t11, t9, -t23, t12, t18 * t31 + t7 * t29, -t16 * t31 + t7 * t28, (-t14 - t15) * t5 * qJD(2), t30 + (t14 / 0.2e1 + t15 / 0.2e1) * t5 ^ 2, t10, t11, t24, t12, t23, t9, -t1 * qJD(4) + t2 * t29, (t1 * t18 - t16 * t3) * qJD(2), t3 * qJD(4) - t2 * t28, t3 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t4;
