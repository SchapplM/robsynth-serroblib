% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:11
% EndTime: 2019-12-05 16:40:11
% DurationCPUTime: 0.11s
% Computational Cost: add. (106->25), mult. (193->64), div. (0->0), fcn. (75->4), ass. (0->31)
t17 = qJD(2) + qJD(3);
t16 = t17 ^ 2;
t33 = t16 / 0.2e1;
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t21 = sin(qJ(3));
t30 = pkin(2) * qJD(2);
t27 = t21 * t30;
t7 = t17 * pkin(7) + t27;
t5 = t20 * qJD(1) + t22 * t7;
t32 = t17 * t20;
t31 = t17 * t22;
t29 = qJ(5) * t17;
t28 = qJD(4) * t17;
t23 = cos(qJ(3));
t26 = t23 * t30;
t24 = qJD(2) ^ 2;
t19 = qJD(1) ^ 2 / 0.2e1;
t18 = qJD(4) ^ 2 / 0.2e1;
t15 = t22 * qJD(1);
t13 = t22 * t28;
t12 = t20 * t28;
t11 = t22 ^ 2 * t33;
t10 = t20 ^ 2 * t33;
t9 = t20 * t16 * t22;
t8 = -t17 * pkin(3) - t26;
t4 = -t20 * t7 + t15;
t3 = -t26 + qJD(5) + (-pkin(4) * t22 - pkin(3)) * t17;
t2 = t22 * t29 + t5;
t1 = qJD(4) * pkin(4) + t15 + (-t7 - t29) * t20;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, t24 / 0.2e1, 0, 0, 0, t19, 0, 0, 0, 0, 0, t33, t17 * t26, -t17 * t27, 0, t19 + (t21 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t24, t10, t9, t12, t11, t13, t18, t4 * qJD(4) - t8 * t31, -t5 * qJD(4) + t8 * t32, (-t20 * t4 + t22 * t5) * t17, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t10, t9, t12, t11, t13, t18, t1 * qJD(4) - t3 * t31, -t2 * qJD(4) + t3 * t32, (-t1 * t20 + t2 * t22) * t17, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t6;
