% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:36
% EndTime: 2019-12-31 17:00:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (64->31), mult. (208->66), div. (0->0), fcn. (67->2), ass. (0->27)
t17 = qJD(1) ^ 2;
t28 = t17 / 0.2e1;
t16 = cos(qJ(2));
t27 = t16 * t17;
t26 = -pkin(2) - qJ(4);
t15 = sin(qJ(2));
t25 = qJD(1) * t15;
t24 = qJD(1) * t16;
t23 = pkin(5) * t25 + qJD(3);
t22 = qJD(2) * qJ(3);
t21 = qJD(1) * qJD(2);
t10 = t15 * t21;
t20 = t16 * t21;
t19 = -qJ(3) * t15 - pkin(1);
t14 = t16 ^ 2;
t13 = t15 ^ 2;
t12 = qJD(2) ^ 2 / 0.2e1;
t9 = t14 * t28;
t8 = t13 * t28;
t7 = t15 * t27;
t6 = -pkin(5) * t24 - t22;
t5 = -qJD(2) * pkin(2) + t23;
t4 = (-pkin(2) * t16 + t19) * qJD(1);
t3 = t22 + qJD(4) + (pkin(3) + pkin(5)) * t24;
t2 = pkin(3) * t25 + t26 * qJD(2) + t23;
t1 = (t26 * t16 + t19) * qJD(1);
t11 = [0, 0, 0, 0, 0, t28, 0, 0, 0, 0, t8, t7, t10, t9, t20, t12, pkin(1) * t27 - pkin(5) * t10, -t17 * pkin(1) * t15 - pkin(5) * t20, (t13 + t14) * t17 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t14 / 0.2e1 + t13 / 0.2e1) * pkin(5) ^ 2) * t17, t12, -t10, -t20, t8, t7, t9, (t15 * t5 - t16 * t6) * qJD(1), t5 * qJD(2) + t4 * t24, -t6 * qJD(2) - t4 * t25, t4 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t12, -t20, t10, t9, -t7, t8, (t15 * t2 + t16 * t3) * qJD(1), t3 * qJD(2) - t1 * t25, -t2 * qJD(2) - t1 * t24, t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t11;
