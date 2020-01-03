% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRP4
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:57
% EndTime: 2019-12-31 16:27:57
% DurationCPUTime: 0.08s
% Computational Cost: add. (39->19), mult. (133->49), div. (0->0), fcn. (50->2), ass. (0->21)
t15 = qJD(2) ^ 2;
t22 = t15 / 0.2e1;
t13 = sin(qJ(3));
t14 = cos(qJ(3));
t19 = qJD(2) * t14;
t5 = pkin(5) * t19 + t13 * qJD(1);
t21 = t14 * t15;
t20 = qJD(2) * t13;
t18 = qJD(2) * qJD(3);
t17 = t13 * t21;
t16 = t14 * t18;
t4 = -pkin(5) * t20 + qJD(1) * t14;
t12 = qJD(1) ^ 2 / 0.2e1;
t11 = qJD(3) ^ 2 / 0.2e1;
t8 = t13 * t18;
t7 = t14 ^ 2 * t22;
t6 = t13 ^ 2 * t22;
t3 = qJD(3) * qJ(4) + t5;
t2 = (-pkin(3) * t14 - qJ(4) * t13 - pkin(2)) * qJD(2);
t1 = -qJD(3) * pkin(3) + qJD(4) - t4;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, t22, 0, 0, 0, t12, t6, t17, t8, t7, t16, t11, pkin(2) * t21 + qJD(3) * t4, -pkin(2) * t13 * t15 - qJD(3) * t5, (-t13 * t4 + t14 * t5) * qJD(2), t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + pkin(2) ^ 2 * t22, t6, t8, -t17, t11, -t16, t7, -qJD(3) * t1 - t2 * t19, (t1 * t13 + t14 * t3) * qJD(2), qJD(3) * t3 - t2 * t20, t3 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t9;
