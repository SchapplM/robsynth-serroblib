% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:21
% EndTime: 2019-12-05 15:05:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (91->23), mult. (224->64), div. (0->0), fcn. (148->8), ass. (0->26)
t24 = qJD(3) ^ 2;
t15 = t24 / 0.2e1;
t29 = qJD(1) ^ 2 / 0.2e1;
t21 = sin(qJ(3));
t23 = cos(qJ(3));
t17 = sin(pkin(8));
t27 = qJD(1) * t17;
t11 = t21 * qJD(2) + t23 * t27;
t16 = sin(pkin(9));
t18 = cos(pkin(9));
t10 = t23 * qJD(2) - t21 * t27;
t9 = qJD(3) * pkin(3) + t10;
t6 = t18 * t11 + t16 * t9;
t5 = -t16 * t11 + t18 * t9;
t3 = -qJD(3) * pkin(4) - t5;
t28 = qJD(3) * t3;
t26 = qJD(3) * qJD(5);
t22 = cos(qJ(5));
t20 = sin(qJ(5));
t19 = cos(pkin(8));
t13 = t19 ^ 2 * t29;
t12 = -t19 * qJD(1) + qJD(4);
t4 = qJD(3) * pkin(6) + t6;
t2 = t20 * t12 + t22 * t4;
t1 = t22 * t12 - t20 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 ^ 2 * t29 + t13 + qJD(2) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t15, t10 * qJD(3), -t11 * qJD(3), 0, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t13, 0, 0, 0, 0, 0, t15, t5 * qJD(3), -t6 * qJD(3), 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t20 ^ 2 * t15, t20 * t24 * t22, t20 * t26, t22 ^ 2 * t15, t22 * t26, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t22 * t28, -t2 * qJD(5) + t20 * t28, (-t1 * t20 + t2 * t22) * qJD(3), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
