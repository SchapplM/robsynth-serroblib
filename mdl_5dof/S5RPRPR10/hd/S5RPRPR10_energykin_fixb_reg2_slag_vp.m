% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:04
% EndTime: 2019-12-31 18:26:04
% DurationCPUTime: 0.10s
% Computational Cost: add. (182->25), mult. (276->67), div. (0->0), fcn. (108->6), ass. (0->26)
t17 = -qJD(1) + qJD(3);
t16 = t17 ^ 2;
t14 = t16 / 0.2e1;
t25 = qJD(1) ^ 2;
t18 = t25 / 0.2e1;
t13 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t26 = qJ(2) * qJD(1);
t11 = t22 * t13 + t24 * t26;
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t10 = t24 * t13 - t22 * t26;
t8 = t17 * pkin(3) + t10;
t6 = t20 * t11 + t19 * t8;
t5 = -t19 * t11 + t20 * t8;
t3 = -t17 * pkin(4) - t5;
t28 = t17 * t3;
t27 = qJD(5) * t17;
t23 = cos(qJ(5));
t21 = sin(qJ(5));
t15 = -pkin(1) * qJD(1) + qJD(2);
t4 = t17 * pkin(7) + t6;
t2 = t21 * qJD(4) + t23 * t4;
t1 = t23 * qJD(4) - t21 * t4;
t7 = [0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, -t15 * qJD(1), 0, t25 * qJ(2), qJ(2) ^ 2 * t18 + t15 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t14, t10 * t17, -t11 * t17, 0, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t14, t5 * t17, -t6 * t17, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t21 ^ 2 * t14, t21 * t16 * t23, t21 * t27, t23 ^ 2 * t14, t23 * t27, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t23 * t28, -t2 * qJD(5) + t21 * t28, (-t1 * t21 + t2 * t23) * t17, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
