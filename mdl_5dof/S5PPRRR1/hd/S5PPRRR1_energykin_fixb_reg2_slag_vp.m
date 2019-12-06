% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:50
% EndTime: 2019-12-05 15:12:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (105->23), mult. (245->66), div. (0->0), fcn. (170->8), ass. (0->25)
t14 = qJD(3) + qJD(4);
t13 = t14 ^ 2;
t27 = t13 / 0.2e1;
t16 = sin(pkin(9));
t17 = cos(pkin(9));
t20 = sin(qJ(3));
t23 = cos(qJ(3));
t11 = (t16 * t23 + t17 * t20) * qJD(1);
t19 = sin(qJ(4));
t22 = cos(qJ(4));
t10 = (-t16 * t20 + t17 * t23) * qJD(1);
t9 = qJD(3) * pkin(3) + t10;
t6 = t22 * t11 + t19 * t9;
t5 = -t19 * t11 + t22 * t9;
t3 = -t14 * pkin(4) - t5;
t26 = t14 * t3;
t25 = qJD(5) * t14;
t24 = qJD(1) ^ 2;
t21 = cos(qJ(5));
t18 = sin(qJ(5));
t15 = qJD(2) ^ 2 / 0.2e1;
t4 = t14 * pkin(7) + t6;
t2 = t18 * qJD(2) + t21 * t4;
t1 = t21 * qJD(2) - t18 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t24 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 + (t16 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1) * t24, 0, 0, 0, 0, 0, qJD(3) ^ 2 / 0.2e1, t10 * qJD(3), -t11 * qJD(3), 0, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t15, 0, 0, 0, 0, 0, t27, t5 * t14, -t6 * t14, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15, t18 ^ 2 * t27, t18 * t13 * t21, t18 * t25, t21 ^ 2 * t27, t21 * t25, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t21 * t26, -t2 * qJD(5) + t18 * t26, (-t1 * t18 + t2 * t21) * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
