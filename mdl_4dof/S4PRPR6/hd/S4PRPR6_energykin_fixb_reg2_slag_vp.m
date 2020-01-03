% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:40
% EndTime: 2019-12-31 16:24:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (70->23), mult. (206->69), div. (0->0), fcn. (115->6), ass. (0->26)
t20 = qJD(2) ^ 2;
t28 = t20 / 0.2e1;
t27 = cos(qJ(4));
t15 = sin(pkin(7));
t26 = qJD(2) * t15;
t16 = cos(pkin(7));
t25 = qJD(2) * t16;
t24 = qJD(1) * qJD(2);
t18 = sin(qJ(2));
t11 = qJD(2) * qJ(3) + t18 * qJD(1);
t23 = pkin(5) * qJD(2) + t11;
t19 = cos(qJ(2));
t22 = -t19 * qJD(1) + qJD(3);
t21 = qJD(1) ^ 2;
t17 = sin(qJ(4));
t14 = t16 ^ 2;
t13 = t15 ^ 2;
t9 = -qJD(2) * pkin(2) + t22;
t8 = (-pkin(3) * t16 - pkin(2)) * qJD(2) + t22;
t7 = (t27 * t15 + t16 * t17) * qJD(2);
t5 = t17 * t26 - t27 * t25;
t4 = t23 * t16;
t3 = t23 * t15;
t2 = -t17 * t3 + t27 * t4;
t1 = -t17 * t4 - t27 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t21 / 0.2e1, 0, 0, 0, 0, 0, t28, t19 * t24, -t18 * t24, 0, (t18 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1) * t21, t13 * t28, t15 * t20 * t16, 0, t14 * t28, 0, 0, -t9 * t25, t9 * t26, (t13 + t14) * t11 * qJD(2), t9 ^ 2 / 0.2e1 + (t14 / 0.2e1 + t13 / 0.2e1) * t11 ^ 2, t7 ^ 2 / 0.2e1, -t7 * t5, t7 * qJD(4), t5 ^ 2 / 0.2e1, -t5 * qJD(4), qJD(4) ^ 2 / 0.2e1, t1 * qJD(4) + t8 * t5, -t2 * qJD(4) + t8 * t7, -t1 * t7 - t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg = t6;
