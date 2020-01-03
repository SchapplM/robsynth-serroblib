% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:36
% EndTime: 2019-12-31 17:38:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (100->25), mult. (187->66), div. (0->0), fcn. (82->6), ass. (0->23)
t19 = qJD(2) ^ 2;
t12 = t19 / 0.2e1;
t16 = sin(qJ(2));
t11 = qJD(2) * qJ(3) + t16 * qJD(1);
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t18 = cos(qJ(2));
t21 = -t18 * qJD(1) + qJD(3);
t8 = (-pkin(2) - pkin(3)) * qJD(2) + t21;
t6 = t14 * t11 + t13 * t8;
t5 = -t13 * t11 + t14 * t8;
t3 = qJD(2) * pkin(4) - t5;
t24 = qJD(2) * t3;
t23 = qJD(1) * qJD(2);
t22 = qJD(2) * qJD(5);
t20 = qJD(1) ^ 2;
t17 = cos(qJ(5));
t15 = sin(qJ(5));
t10 = -qJD(2) * pkin(2) + t21;
t4 = -qJD(2) * pkin(6) + t6;
t2 = t15 * qJD(4) + t17 * t4;
t1 = t17 * qJD(4) - t15 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t20 / 0.2e1, 0, 0, 0, 0, 0, t12, t18 * t23, -t16 * t23, 0, (t16 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1) * t20, 0, 0, 0, t12, 0, 0, -t10 * qJD(2), 0, t11 * qJD(2), t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t12, -t5 * qJD(2), t6 * qJD(2), 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t15 ^ 2 * t12, t15 * t19 * t17, -t15 * t22, t17 ^ 2 * t12, -t17 * t22, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t17 * t24, -t2 * qJD(5) - t15 * t24, (t1 * t15 - t17 * t2) * qJD(2), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
