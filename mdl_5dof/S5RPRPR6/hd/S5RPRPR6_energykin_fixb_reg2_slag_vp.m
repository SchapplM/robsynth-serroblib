% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:48
% EndTime: 2019-12-31 18:17:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (117->26), mult. (227->62), div. (0->0), fcn. (94->6), ass. (0->27)
t13 = qJD(1) + qJD(3);
t12 = t13 ^ 2;
t11 = t12 / 0.2e1;
t22 = qJD(1) ^ 2;
t29 = pkin(1) * t22;
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t16 = sin(pkin(8));
t25 = pkin(1) * qJD(1) * t16;
t17 = cos(pkin(8));
t9 = (pkin(1) * t17 + pkin(2)) * qJD(1);
t8 = t19 * t9 + t21 * t25;
t5 = t13 * qJ(4) + t8;
t28 = t5 * t13;
t27 = t5 ^ 2 / 0.2e1;
t26 = qJD(5) * t13;
t7 = -t19 * t25 + t21 * t9;
t24 = qJD(4) - t7;
t20 = cos(qJ(5));
t18 = sin(qJ(5));
t15 = t22 / 0.2e1;
t14 = qJD(2) ^ 2 / 0.2e1;
t4 = -t13 * pkin(3) + t24;
t3 = (-pkin(3) - pkin(7)) * t13 + t24;
t2 = t20 * qJD(2) + t18 * t3;
t1 = -t18 * qJD(2) + t20 * t3;
t6 = [0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t17 * t29, -t16 * t29, 0, t14 + (t16 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t22, 0, 0, 0, 0, 0, t11, t7 * t13, -t8 * t13, 0, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t14, t11, 0, 0, 0, 0, 0, 0, t4 * t13, t28, t14 + t27 + t4 ^ 2 / 0.2e1, t20 ^ 2 * t11, -t20 * t12 * t18, t20 * t26, t18 ^ 2 * t11, -t18 * t26, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t18 * t28, -t2 * qJD(5) + t20 * t28, (-t1 * t20 - t18 * t2) * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t27;];
T_reg = t6;
