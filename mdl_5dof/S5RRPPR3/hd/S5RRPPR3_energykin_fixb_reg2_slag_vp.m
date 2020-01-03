% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:35
% EndTime: 2019-12-31 19:26:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (136->25), mult. (225->62), div. (0->0), fcn. (94->6), ass. (0->27)
t13 = qJD(1) + qJD(2);
t12 = t13 ^ 2;
t11 = t12 / 0.2e1;
t15 = sin(pkin(8));
t16 = cos(pkin(8));
t18 = sin(qJ(2));
t28 = pkin(1) * qJD(1);
t25 = t18 * t28;
t20 = cos(qJ(2));
t24 = t20 * t28;
t9 = t13 * pkin(2) + t24;
t8 = t15 * t9 + t16 * t25;
t5 = t13 * qJ(4) + t8;
t29 = t5 * t13;
t27 = t5 ^ 2 / 0.2e1;
t26 = qJD(5) * t13;
t7 = -t15 * t25 + t16 * t9;
t23 = qJD(4) - t7;
t21 = qJD(1) ^ 2;
t19 = cos(qJ(5));
t17 = sin(qJ(5));
t14 = qJD(3) ^ 2 / 0.2e1;
t4 = -t13 * pkin(3) + t23;
t3 = (-pkin(3) - pkin(7)) * t13 + t23;
t2 = t19 * qJD(3) + t17 * t3;
t1 = -t17 * qJD(3) + t19 * t3;
t6 = [0, 0, 0, 0, 0, t21 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t13 * t24, -t13 * t25, 0, (t18 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t21, 0, 0, 0, 0, 0, t11, t7 * t13, -t8 * t13, 0, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t14, t11, 0, 0, 0, 0, 0, 0, t4 * t13, t29, t14 + t27 + t4 ^ 2 / 0.2e1, t19 ^ 2 * t11, -t19 * t12 * t17, t19 * t26, t17 ^ 2 * t11, -t17 * t26, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t17 * t29, -t2 * qJD(5) + t19 * t29, (-t1 * t19 - t17 * t2) * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t27;];
T_reg = t6;
