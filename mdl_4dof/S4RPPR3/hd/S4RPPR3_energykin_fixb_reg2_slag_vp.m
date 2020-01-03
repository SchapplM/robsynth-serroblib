% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:55
% EndTime: 2019-12-31 16:37:56
% DurationCPUTime: 0.10s
% Computational Cost: add. (87->28), mult. (250->76), div. (0->0), fcn. (133->6), ass. (0->25)
t23 = qJD(1) ^ 2;
t17 = t23 / 0.2e1;
t29 = pkin(1) * t23;
t28 = cos(qJ(4));
t19 = sin(pkin(6));
t13 = (pkin(1) * t19 + qJ(3)) * qJD(1);
t18 = sin(pkin(7));
t20 = cos(pkin(7));
t6 = t18 * qJD(2) + t20 * t13;
t27 = qJD(1) * t18;
t26 = qJD(1) * t20;
t21 = cos(pkin(6));
t25 = -pkin(1) * t21 - pkin(2);
t22 = sin(qJ(4));
t16 = t20 * qJD(2);
t12 = t25 * qJD(1) + qJD(3);
t10 = (t28 * t18 + t20 * t22) * qJD(1);
t8 = t22 * t27 - t28 * t26;
t7 = qJD(3) + (-pkin(3) * t20 + t25) * qJD(1);
t5 = -t18 * t13 + t16;
t4 = pkin(5) * t26 + t6;
t3 = t16 + (-pkin(5) * qJD(1) - t13) * t18;
t2 = t22 * t3 + t28 * t4;
t1 = -t22 * t4 + t28 * t3;
t9 = [0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t21 * t29, -t19 * t29, 0, qJD(2) ^ 2 / 0.2e1 + (t19 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t23, t18 ^ 2 * t17, t18 * t23 * t20, 0, t20 ^ 2 * t17, 0, 0, -t12 * t26, t12 * t27, (-t18 * t5 + t20 * t6) * qJD(1), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * qJD(4), t8 ^ 2 / 0.2e1, -t8 * qJD(4), qJD(4) ^ 2 / 0.2e1, t1 * qJD(4) + t7 * t8, -t2 * qJD(4) + t7 * t10, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
