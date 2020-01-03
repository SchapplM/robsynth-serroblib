% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRP5
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:21
% EndTime: 2019-12-31 16:29:21
% DurationCPUTime: 0.09s
% Computational Cost: add. (45->19), mult. (153->59), div. (0->0), fcn. (61->4), ass. (0->28)
t19 = qJD(2) ^ 2;
t28 = t19 / 0.2e1;
t16 = sin(qJ(2));
t5 = qJD(2) * pkin(5) + qJD(1) * t16;
t27 = qJD(3) * t5;
t18 = cos(qJ(2));
t26 = qJD(1) * t18;
t15 = sin(qJ(3));
t25 = qJD(2) * t15;
t17 = cos(qJ(3));
t24 = qJD(2) * t17;
t23 = qJD(1) * qJD(2);
t22 = qJD(2) * qJD(3);
t21 = qJ(4) * qJD(2) + t5;
t20 = qJD(1) ^ 2;
t14 = t17 ^ 2;
t13 = t15 ^ 2;
t12 = qJD(3) ^ 2 / 0.2e1;
t11 = t17 * t22;
t10 = t15 * t22;
t9 = t14 * t28;
t8 = t13 * t28;
t7 = t15 * t19 * t17;
t6 = -qJD(2) * pkin(2) - t26;
t3 = -t26 + qJD(4) + (-pkin(3) * t17 - pkin(2)) * qJD(2);
t2 = t21 * t17;
t1 = qJD(3) * pkin(3) - t21 * t15;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t20 / 0.2e1, 0, 0, 0, 0, 0, t28, t18 * t23, -t16 * t23, 0, (t16 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1) * t20, t8, t7, t10, t9, t11, t12, -t15 * t27 - t6 * t24, -t17 * t27 + t6 * t25, (t13 + t14) * t5 * qJD(2), t6 ^ 2 / 0.2e1 + (t14 / 0.2e1 + t13 / 0.2e1) * t5 ^ 2, t8, t7, t10, t9, t11, t12, qJD(3) * t1 - t3 * t24, -qJD(3) * t2 + t3 * t25, (-t1 * t15 + t17 * t2) * qJD(2), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t4;
