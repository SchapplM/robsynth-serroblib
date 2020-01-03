% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:47
% EndTime: 2019-12-31 16:42:47
% DurationCPUTime: 0.09s
% Computational Cost: add. (58->23), mult. (189->63), div. (0->0), fcn. (75->4), ass. (0->27)
t22 = qJD(1) ^ 2;
t17 = t22 / 0.2e1;
t29 = pkin(1) * t22;
t20 = sin(qJ(3));
t21 = cos(qJ(3));
t18 = sin(pkin(6));
t7 = (pkin(1) * t18 + pkin(5)) * qJD(1);
t4 = t20 * qJD(2) + t21 * t7;
t28 = qJD(1) * t20;
t27 = qJD(1) * t21;
t26 = qJ(4) * qJD(1);
t25 = qJD(1) * qJD(3);
t19 = cos(pkin(6));
t24 = -pkin(1) * t19 - pkin(2);
t16 = qJD(3) ^ 2 / 0.2e1;
t15 = t21 * qJD(2);
t13 = t21 * t25;
t12 = t20 * t25;
t11 = t21 ^ 2 * t17;
t10 = t20 ^ 2 * t17;
t9 = t20 * t22 * t21;
t8 = t24 * qJD(1);
t5 = qJD(4) + (-pkin(3) * t21 + t24) * qJD(1);
t3 = -t20 * t7 + t15;
t2 = t21 * t26 + t4;
t1 = qJD(3) * pkin(3) + t15 + (-t7 - t26) * t20;
t6 = [0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t19 * t29, -t18 * t29, 0, qJD(2) ^ 2 / 0.2e1 + (t18 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t22, t10, t9, t12, t11, t13, t16, t3 * qJD(3) - t8 * t27, -t4 * qJD(3) + t8 * t28, (-t20 * t3 + t21 * t4) * qJD(1), t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t10, t9, t12, t11, t13, t16, t1 * qJD(3) - t5 * t27, -t2 * qJD(3) + t5 * t28, (-t1 * t20 + t2 * t21) * qJD(1), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
