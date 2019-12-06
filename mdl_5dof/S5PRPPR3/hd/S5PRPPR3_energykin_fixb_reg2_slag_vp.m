% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:54
% EndTime: 2019-12-05 15:26:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (82->24), mult. (185->60), div. (0->0), fcn. (94->6), ass. (0->25)
t19 = qJD(2) ^ 2;
t12 = t19 / 0.2e1;
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t16 = sin(qJ(2));
t24 = qJD(1) * t16;
t18 = cos(qJ(2));
t9 = qJD(2) * pkin(2) + t18 * qJD(1);
t8 = t13 * t9 + t14 * t24;
t5 = qJD(2) * qJ(4) + t8;
t26 = t5 * qJD(2);
t25 = t5 ^ 2 / 0.2e1;
t23 = qJD(1) * qJD(2);
t22 = qJD(2) * qJD(5);
t7 = -t13 * t24 + t14 * t9;
t21 = qJD(4) - t7;
t20 = qJD(1) ^ 2;
t17 = cos(qJ(5));
t15 = sin(qJ(5));
t11 = qJD(3) ^ 2 / 0.2e1;
t4 = -qJD(2) * pkin(3) + t21;
t3 = (-pkin(3) - pkin(6)) * qJD(2) + t21;
t2 = t17 * qJD(3) + t15 * t3;
t1 = -t15 * qJD(3) + t17 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t20 / 0.2e1, 0, 0, 0, 0, 0, t12, t18 * t23, -t16 * t23, 0, (t16 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1) * t20, 0, 0, 0, 0, 0, t12, t7 * qJD(2), -t8 * qJD(2), 0, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t11, t12, 0, 0, 0, 0, 0, 0, t4 * qJD(2), t26, t11 + t25 + t4 ^ 2 / 0.2e1, t17 ^ 2 * t12, -t17 * t19 * t15, t17 * t22, t15 ^ 2 * t12, -t15 * t22, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t15 * t26, -t2 * qJD(5) + t17 * t26, (-t1 * t17 - t15 * t2) * qJD(2), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t25;];
T_reg = t6;
