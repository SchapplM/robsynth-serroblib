% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:41
% EndTime: 2019-12-31 18:16:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (103->35), mult. (237->70), div. (0->0), fcn. (67->2), ass. (0->30)
t23 = qJD(1) ^ 2;
t17 = t23 / 0.2e1;
t32 = -pkin(3) - pkin(4);
t21 = sin(qJ(3));
t9 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t6 = qJD(3) * qJ(4) + t21 * t9;
t31 = qJD(3) * t9;
t30 = t23 * qJ(2);
t29 = qJD(1) * t21;
t22 = cos(qJ(3));
t28 = qJD(1) * t22;
t27 = qJ(5) * qJD(1);
t26 = qJD(1) * qJD(3);
t25 = t21 * t26;
t24 = qJ(4) * t22 - qJ(2);
t20 = t22 ^ 2;
t19 = t21 ^ 2;
t16 = qJD(3) ^ 2 / 0.2e1;
t15 = qJ(2) ^ 2 * t17;
t14 = -qJD(1) * pkin(1) + qJD(2);
t13 = t22 * t26;
t12 = t20 * t17;
t11 = t19 * t17;
t10 = t22 * t23 * t21;
t5 = (pkin(3) * t21 - t24) * qJD(1);
t4 = -qJD(3) * pkin(3) - t22 * t9 + qJD(4);
t3 = t21 * t27 + t6;
t2 = qJD(5) + (t32 * t21 + t24) * qJD(1);
t1 = qJD(4) + (-t9 - t27) * t22 + t32 * qJD(3);
t7 = [0, 0, 0, 0, 0, t17, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, t14 * qJD(1), t30, t15 + t14 ^ 2 / 0.2e1, t12, -t10, t13, t11, -t25, t16, t21 * t30 + t22 * t31, -t21 * t31 + t22 * t30, (-t19 - t20) * t9 * qJD(1), t15 + (t19 / 0.2e1 + t20 / 0.2e1) * t9 ^ 2, t12, t13, t10, t16, t25, t11, -t4 * qJD(3) + t5 * t29, (-t21 * t6 + t22 * t4) * qJD(1), t6 * qJD(3) - t5 * t28, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t12, t10, -t13, t11, -t25, t16, -t1 * qJD(3) - t2 * t29, t3 * qJD(3) + t2 * t28, (-t1 * t22 + t21 * t3) * qJD(1), t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t7;
