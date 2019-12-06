% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:36
% EndTime: 2019-12-05 14:59:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (76->21), mult. (188->63), div. (0->0), fcn. (130->8), ass. (0->25)
t21 = qJD(4) ^ 2;
t27 = t21 / 0.2e1;
t13 = sin(pkin(9));
t15 = cos(pkin(9));
t14 = sin(pkin(8));
t24 = qJD(1) * t14;
t10 = t13 * qJD(2) + t15 * t24;
t16 = cos(pkin(8));
t12 = -t16 * qJD(1) + qJD(3);
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t6 = t20 * t10 + t18 * t12;
t5 = -t18 * t10 + t20 * t12;
t3 = -qJD(4) * pkin(4) - t5;
t26 = qJD(4) * t3;
t8 = -t15 * qJD(2) + t13 * t24;
t25 = t8 ^ 2 / 0.2e1;
t23 = qJD(4) * qJD(5);
t22 = qJD(1) ^ 2;
t19 = cos(qJ(5));
t17 = sin(qJ(5));
t4 = qJD(4) * pkin(6) + t6;
t2 = t17 * t8 + t19 * t4;
t1 = -t17 * t4 + t19 * t8;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t22 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t14 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1) * t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 / 0.2e1 + t25 + t12 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t27, t5 * qJD(4), -t6 * qJD(4), 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t25, t17 ^ 2 * t27, t17 * t21 * t19, t17 * t23, t19 ^ 2 * t27, t19 * t23, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t19 * t26, -t2 * qJD(5) + t17 * t26, (-t1 * t17 + t19 * t2) * qJD(4), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
