% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRPR3
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:53
% EndTime: 2019-12-31 16:20:54
% DurationCPUTime: 0.09s
% Computational Cost: add. (62->23), mult. (188->62), div. (0->0), fcn. (107->4), ass. (0->22)
t20 = qJD(2) ^ 2;
t25 = t20 / 0.2e1;
t24 = cos(qJ(4));
t17 = sin(pkin(7));
t18 = cos(pkin(7));
t21 = qJ(3) * qJD(2);
t9 = t17 * qJD(1) + t18 * t21;
t23 = qJD(2) * t17;
t22 = qJD(2) * t18;
t19 = sin(qJ(4));
t16 = qJD(1) ^ 2 / 0.2e1;
t15 = t18 * qJD(1);
t13 = -qJD(2) * pkin(2) + qJD(3);
t10 = qJD(3) + (-pkin(3) * t18 - pkin(2)) * qJD(2);
t8 = -t17 * t21 + t15;
t7 = (t24 * t17 + t18 * t19) * qJD(2);
t5 = t19 * t23 - t24 * t22;
t4 = pkin(5) * t22 + t9;
t3 = t15 + (-pkin(5) - qJ(3)) * t23;
t2 = t19 * t3 + t24 * t4;
t1 = -t19 * t4 + t24 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, t25, 0, 0, 0, t16, t17 ^ 2 * t25, t17 * t20 * t18, 0, t18 ^ 2 * t25, 0, 0, -t13 * t22, t13 * t23, (-t17 * t8 + t18 * t9) * qJD(2), t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t7 ^ 2 / 0.2e1, -t7 * t5, t7 * qJD(4), t5 ^ 2 / 0.2e1, -t5 * qJD(4), qJD(4) ^ 2 / 0.2e1, qJD(4) * t1 + t10 * t5, -qJD(4) * t2 + t10 * t7, -t1 * t7 - t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t6;
