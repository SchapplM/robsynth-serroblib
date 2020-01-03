% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:14
% EndTime: 2019-12-31 17:45:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (128->30), mult. (288->81), div. (0->0), fcn. (137->6), ass. (0->27)
t28 = qJD(1) ^ 2;
t20 = t28 / 0.2e1;
t33 = pkin(1) * t28;
t25 = cos(pkin(7));
t30 = -pkin(1) * t25 - pkin(2);
t13 = qJD(3) + (-qJ(4) + t30) * qJD(1);
t22 = sin(pkin(8));
t24 = cos(pkin(8));
t6 = t24 * qJD(2) + t22 * t13;
t23 = sin(pkin(7));
t16 = (-pkin(1) * t23 - qJ(3)) * qJD(1);
t32 = qJD(1) * t22;
t31 = qJD(1) * t24;
t14 = qJD(4) - t16;
t5 = -t22 * qJD(2) + t24 * t13;
t27 = cos(qJ(5));
t26 = sin(qJ(5));
t19 = qJD(2) ^ 2 / 0.2e1;
t15 = t30 * qJD(1) + qJD(3);
t12 = (-t22 * t26 + t24 * t27) * qJD(1);
t10 = (-t22 * t27 - t24 * t26) * qJD(1);
t9 = pkin(4) * t32 + t14;
t4 = -pkin(6) * t32 + t6;
t3 = -pkin(6) * t31 + t5;
t2 = t26 * t3 + t27 * t4;
t1 = -t26 * t4 + t27 * t3;
t7 = [0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t25 * t33, -t23 * t33, 0, t19 + (t23 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t28, t20, 0, 0, 0, 0, 0, 0, t15 * qJD(1), -t16 * qJD(1), t19 + t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t24 ^ 2 * t20, -t24 * t28 * t22, 0, t22 ^ 2 * t20, 0, 0, t14 * t32, t14 * t31, (-t22 * t6 - t24 * t5) * qJD(1), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t12 ^ 2 / 0.2e1, t12 * t10, t12 * qJD(5), t10 ^ 2 / 0.2e1, t10 * qJD(5), qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t9 * t10, -t2 * qJD(5) + t9 * t12, -t1 * t12 + t2 * t10, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg = t7;
