% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:40
% EndTime: 2019-12-31 16:41:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (81->20), mult. (199->64), div. (0->0), fcn. (89->4), ass. (0->21)
t21 = qJD(1) ^ 2;
t15 = t21 / 0.2e1;
t17 = sin(pkin(6));
t23 = qJD(1) * t17;
t11 = qJD(1) * qJ(2) + qJD(3);
t10 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t22 = -pkin(5) * qJD(1) + t10;
t20 = cos(qJ(4));
t19 = sin(qJ(4));
t18 = cos(pkin(6));
t14 = t18 ^ 2;
t13 = t17 ^ 2;
t12 = -qJD(1) * pkin(1) + qJD(2);
t8 = pkin(3) * t23 + t11;
t7 = (-t17 * t19 + t18 * t20) * qJD(1);
t5 = (-t17 * t20 - t18 * t19) * qJD(1);
t4 = t22 * t18;
t3 = t22 * t17;
t2 = t19 * t4 + t20 * t3;
t1 = -t19 * t3 + t20 * t4;
t6 = [0, 0, 0, 0, 0, t15, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, t12 * qJD(1), t21 * qJ(2), qJ(2) ^ 2 * t15 + t12 ^ 2 / 0.2e1, t14 * t15, -t18 * t21 * t17, 0, t13 * t15, 0, 0, t11 * t23, t11 * t18 * qJD(1), (-t13 - t14) * t10 * qJD(1), t11 ^ 2 / 0.2e1 + (t13 / 0.2e1 + t14 / 0.2e1) * t10 ^ 2, t7 ^ 2 / 0.2e1, t7 * t5, t7 * qJD(4), t5 ^ 2 / 0.2e1, t5 * qJD(4), qJD(4) ^ 2 / 0.2e1, t1 * qJD(4) - t8 * t5, -t2 * qJD(4) + t8 * t7, -t1 * t7 + t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg = t6;
