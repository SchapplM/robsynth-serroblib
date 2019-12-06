% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:56
% EndTime: 2019-12-05 16:41:56
% DurationCPUTime: 0.11s
% Computational Cost: add. (105->26), mult. (190->63), div. (0->0), fcn. (72->4), ass. (0->29)
t14 = qJD(2) + qJD(3);
t13 = t14 ^ 2;
t31 = t13 / 0.2e1;
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t18 = sin(qJ(3));
t28 = pkin(2) * qJD(2);
t25 = t18 * t28;
t7 = t14 * pkin(7) + t25;
t5 = t17 * qJD(1) + t19 * t7;
t30 = t14 * t17;
t29 = t14 * t19;
t27 = qJD(4) * t14;
t26 = t17 * t13 * t19;
t20 = cos(qJ(3));
t24 = t20 * t28;
t23 = t19 * t27;
t4 = t19 * qJD(1) - t17 * t7;
t21 = qJD(2) ^ 2;
t16 = qJD(1) ^ 2 / 0.2e1;
t15 = qJD(4) ^ 2 / 0.2e1;
t11 = t17 * t27;
t10 = t19 ^ 2 * t31;
t9 = t17 ^ 2 * t31;
t8 = -t14 * pkin(3) - t24;
t3 = qJD(4) * qJ(5) + t5;
t2 = -qJD(4) * pkin(4) + qJD(5) - t4;
t1 = -t24 + (-pkin(4) * t19 - qJ(5) * t17 - pkin(3)) * t14;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, t21 / 0.2e1, 0, 0, 0, t16, 0, 0, 0, 0, 0, t31, t14 * t24, -t14 * t25, 0, t16 + (t18 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t21, t9, t26, t11, t10, t23, t15, t4 * qJD(4) - t8 * t29, -t5 * qJD(4) + t8 * t30, (-t17 * t4 + t19 * t5) * t14, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t9, t11, -t26, t15, -t23, t10, -t2 * qJD(4) - t1 * t29, (t17 * t2 + t19 * t3) * t14, t3 * qJD(4) - t1 * t30, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t6;
