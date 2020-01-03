% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:05
% EndTime: 2019-12-31 16:35:06
% DurationCPUTime: 0.10s
% Computational Cost: add. (78->25), mult. (226->78), div. (0->0), fcn. (121->6), ass. (0->29)
t21 = qJD(2) ^ 2;
t31 = t21 / 0.2e1;
t30 = cos(qJ(4));
t17 = sin(qJ(3));
t29 = qJD(2) * t17;
t19 = cos(qJ(3));
t28 = qJD(2) * t19;
t18 = sin(qJ(2));
t10 = qJD(2) * pkin(5) + t18 * qJD(1);
t27 = qJD(3) * t10;
t20 = cos(qJ(2));
t26 = t20 * qJD(1);
t25 = qJD(1) * qJD(2);
t24 = qJD(2) * qJD(3);
t23 = pkin(6) * qJD(2) + t10;
t22 = qJD(1) ^ 2;
t16 = sin(qJ(4));
t15 = t19 ^ 2;
t14 = t17 ^ 2;
t13 = qJD(3) + qJD(4);
t11 = -qJD(2) * pkin(2) - t26;
t8 = -t26 + (-pkin(3) * t19 - pkin(2)) * qJD(2);
t7 = (t16 * t19 + t30 * t17) * qJD(2);
t5 = t16 * t29 - t30 * t28;
t4 = t23 * t19;
t3 = qJD(3) * pkin(3) - t23 * t17;
t2 = t16 * t3 + t30 * t4;
t1 = -t16 * t4 + t30 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t22 / 0.2e1, 0, 0, 0, 0, 0, t31, t20 * t25, -t18 * t25, 0, (t18 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1) * t22, t14 * t31, t17 * t21 * t19, t17 * t24, t15 * t31, t19 * t24, qJD(3) ^ 2 / 0.2e1, -t11 * t28 - t17 * t27, t11 * t29 - t19 * t27, (t14 + t15) * t10 * qJD(2), t11 ^ 2 / 0.2e1 + (t15 / 0.2e1 + t14 / 0.2e1) * t10 ^ 2, t7 ^ 2 / 0.2e1, -t7 * t5, t7 * t13, t5 ^ 2 / 0.2e1, -t5 * t13, t13 ^ 2 / 0.2e1, t1 * t13 + t8 * t5, -t2 * t13 + t8 * t7, -t1 * t7 - t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg = t6;
