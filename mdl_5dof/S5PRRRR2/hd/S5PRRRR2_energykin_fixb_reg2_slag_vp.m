% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:50
% EndTime: 2019-12-05 17:04:50
% DurationCPUTime: 0.09s
% Computational Cost: add. (101->20), mult. (180->56), div. (0->0), fcn. (78->6), ass. (0->25)
t12 = qJD(2) + qJD(3);
t11 = qJD(4) + t12;
t10 = t11 ^ 2;
t29 = t10 / 0.2e1;
t15 = sin(qJ(4));
t18 = cos(qJ(3));
t26 = pkin(2) * qJD(2);
t22 = t18 * t26;
t21 = t12 * pkin(3) + t22;
t16 = sin(qJ(3));
t23 = t16 * t26;
t27 = cos(qJ(4));
t6 = t15 * t21 + t27 * t23;
t4 = t15 * t23 - t27 * t21;
t28 = t4 * t11;
t25 = t4 ^ 2 / 0.2e1;
t24 = qJD(5) * t11;
t19 = qJD(2) ^ 2;
t17 = cos(qJ(5));
t14 = sin(qJ(5));
t13 = qJD(1) ^ 2 / 0.2e1;
t3 = t11 * pkin(6) + t6;
t2 = t14 * qJD(1) + t17 * t3;
t1 = t17 * qJD(1) - t14 * t3;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, t19 / 0.2e1, 0, 0, 0, t13, 0, 0, 0, 0, 0, t12 ^ 2 / 0.2e1, t12 * t22, -t12 * t23, 0, t13 + (t16 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t19, 0, 0, 0, 0, 0, t29, -t28, -t6 * t11, 0, t6 ^ 2 / 0.2e1 + t25 + t13, t14 ^ 2 * t29, t14 * t10 * t17, t14 * t24, t17 ^ 2 * t29, t17 * t24, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t17 * t28, -t2 * qJD(5) + t14 * t28, (-t1 * t14 + t17 * t2) * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t25;];
T_reg = t5;
