% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:43
% EndTime: 2019-12-31 16:40:43
% DurationCPUTime: 0.10s
% Computational Cost: add. (75->29), mult. (234->66), div. (0->0), fcn. (112->4), ass. (0->29)
t23 = qJD(1) ^ 2;
t31 = t23 / 0.2e1;
t30 = qJ(2) * t23;
t19 = sin(pkin(6));
t29 = qJD(1) * t19;
t20 = cos(pkin(6));
t28 = qJD(1) * t20;
t10 = qJ(2) * t29 + qJD(3);
t27 = t19 * t23 * t20;
t26 = qJ(2) ^ 2 * t31;
t25 = qJ(3) * t19 + pkin(1);
t22 = cos(qJ(4));
t21 = sin(qJ(4));
t18 = t20 ^ 2;
t17 = t19 ^ 2;
t16 = -qJD(1) * pkin(1) + qJD(2);
t15 = t18 * t30;
t13 = t18 * t31;
t12 = t17 * t31;
t11 = t18 * t26;
t9 = (-pkin(5) + qJ(2)) * t28;
t8 = -pkin(5) * t29 + t10;
t7 = (t19 * t22 - t20 * t21) * qJD(1);
t5 = (-t19 * t21 - t20 * t22) * qJD(1);
t4 = qJD(2) + (-pkin(2) * t20 - t25) * qJD(1);
t3 = -qJD(2) + ((pkin(2) + pkin(3)) * t20 + t25) * qJD(1);
t2 = t21 * t8 + t22 * t9;
t1 = -t21 * t9 + t22 * t8;
t6 = [0, 0, 0, 0, 0, t31, 0, 0, 0, 0, t12, t27, 0, t13, 0, 0, -t16 * t28, t16 * t29, t17 * t30 + t15, t11 + t17 * t26 + t16 ^ 2 / 0.2e1, t12, 0, -t27, 0, 0, t13, -t4 * t28, t10 * t29 + t15, -t4 * t29, t11 + t4 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t7 ^ 2 / 0.2e1, t7 * t5, t7 * qJD(4), t5 ^ 2 / 0.2e1, t5 * qJD(4), qJD(4) ^ 2 / 0.2e1, t1 * qJD(4) - t3 * t5, -t2 * qJD(4) + t3 * t7, -t1 * t7 + t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t6;
