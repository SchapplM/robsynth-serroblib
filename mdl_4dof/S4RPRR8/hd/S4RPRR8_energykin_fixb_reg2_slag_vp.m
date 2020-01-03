% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:14
% EndTime: 2019-12-31 16:55:15
% DurationCPUTime: 0.09s
% Computational Cost: add. (93->24), mult. (220->68), div. (0->0), fcn. (93->4), ass. (0->24)
t21 = qJD(1) ^ 2;
t14 = t21 / 0.2e1;
t25 = t21 * qJ(2);
t10 = qJD(2) + (-pkin(1) - pkin(5)) * qJD(1);
t24 = qJD(3) * t10;
t23 = qJD(1) * qJD(3);
t22 = -pkin(6) * qJD(1) + t10;
t20 = cos(qJ(3));
t19 = cos(qJ(4));
t18 = sin(qJ(3));
t17 = sin(qJ(4));
t16 = t20 ^ 2;
t15 = t18 ^ 2;
t13 = qJD(3) + qJD(4);
t12 = qJ(2) ^ 2 * t14;
t11 = -qJD(1) * pkin(1) + qJD(2);
t8 = (pkin(3) * t18 + qJ(2)) * qJD(1);
t7 = (-t17 * t18 + t19 * t20) * qJD(1);
t5 = (-t17 * t20 - t18 * t19) * qJD(1);
t4 = t22 * t18;
t3 = qJD(3) * pkin(3) + t22 * t20;
t2 = t17 * t3 + t19 * t4;
t1 = -t17 * t4 + t19 * t3;
t6 = [0, 0, 0, 0, 0, t14, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, t11 * qJD(1), t25, t12 + t11 ^ 2 / 0.2e1, t16 * t14, -t20 * t21 * t18, t20 * t23, t15 * t14, -t18 * t23, qJD(3) ^ 2 / 0.2e1, t18 * t25 + t20 * t24, -t18 * t24 + t20 * t25, (-t15 - t16) * t10 * qJD(1), t12 + (t15 / 0.2e1 + t16 / 0.2e1) * t10 ^ 2, t7 ^ 2 / 0.2e1, t7 * t5, t7 * t13, t5 ^ 2 / 0.2e1, t5 * t13, t13 ^ 2 / 0.2e1, t1 * t13 - t8 * t5, -t2 * t13 + t8 * t7, -t1 * t7 + t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg = t6;
