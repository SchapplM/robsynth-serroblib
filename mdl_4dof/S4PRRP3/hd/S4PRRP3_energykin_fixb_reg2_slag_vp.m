% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:55
% EndTime: 2019-12-31 16:26:55
% DurationCPUTime: 0.08s
% Computational Cost: add. (37->18), mult. (136->49), div. (0->0), fcn. (53->2), ass. (0->22)
t18 = qJD(2) ^ 2;
t23 = t18 / 0.2e1;
t17 = cos(qJ(3));
t22 = t17 * t18;
t16 = sin(qJ(3));
t20 = qJD(2) * t17;
t4 = pkin(5) * t20 + t16 * qJD(1);
t21 = qJD(2) * t16;
t19 = qJD(2) * qJD(3);
t15 = qJD(1) ^ 2 / 0.2e1;
t14 = qJD(3) ^ 2 / 0.2e1;
t13 = t17 * qJD(1);
t10 = t17 * t19;
t9 = t16 * t19;
t8 = t17 ^ 2 * t23;
t7 = t16 ^ 2 * t23;
t6 = t16 * t22;
t5 = qJD(4) + (-pkin(3) * t17 - pkin(2)) * qJD(2);
t3 = -pkin(5) * t21 + t13;
t2 = qJ(4) * t20 + t4;
t1 = qJD(3) * pkin(3) + t13 + (-pkin(5) - qJ(4)) * t21;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, t23, 0, 0, 0, t15, t7, t6, t9, t8, t10, t14, pkin(2) * t22 + qJD(3) * t3, -pkin(2) * t16 * t18 - qJD(3) * t4, (-t16 * t3 + t17 * t4) * qJD(2), t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + pkin(2) ^ 2 * t23, t7, t6, t9, t8, t10, t14, qJD(3) * t1 - t5 * t20, -qJD(3) * t2 + t5 * t21, (-t1 * t16 + t17 * t2) * qJD(2), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t11;
