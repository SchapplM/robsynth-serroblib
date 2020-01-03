% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRP6
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:42
% EndTime: 2019-12-31 16:30:42
% DurationCPUTime: 0.09s
% Computational Cost: add. (47->22), mult. (150->60), div. (0->0), fcn. (58->4), ass. (0->27)
t17 = qJD(2) ^ 2;
t27 = t17 / 0.2e1;
t14 = sin(qJ(2));
t5 = qJD(2) * pkin(5) + t14 * qJD(1);
t26 = qJD(3) * t5;
t13 = sin(qJ(3));
t25 = qJD(2) * t13;
t15 = cos(qJ(3));
t24 = qJD(2) * t15;
t16 = cos(qJ(2));
t23 = t16 * qJD(1);
t22 = qJD(1) * qJD(2);
t21 = qJD(2) * qJD(3);
t20 = t13 * t17 * t15;
t19 = t15 * t21;
t18 = qJD(1) ^ 2;
t12 = t15 ^ 2;
t11 = t13 ^ 2;
t10 = qJD(3) ^ 2 / 0.2e1;
t9 = t13 * t21;
t8 = t12 * t27;
t7 = t11 * t27;
t6 = -qJD(2) * pkin(2) - t23;
t3 = qJD(3) * qJ(4) + t15 * t5;
t2 = -qJD(3) * pkin(3) + t13 * t5 + qJD(4);
t1 = -t23 + (-pkin(3) * t15 - qJ(4) * t13 - pkin(2)) * qJD(2);
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t18 / 0.2e1, 0, 0, 0, 0, 0, t27, t16 * t22, -t14 * t22, 0, (t14 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1) * t18, t7, t20, t9, t8, t19, t10, -t13 * t26 - t6 * t24, -t15 * t26 + t6 * t25, (t11 + t12) * t5 * qJD(2), t6 ^ 2 / 0.2e1 + (t12 / 0.2e1 + t11 / 0.2e1) * t5 ^ 2, t7, t9, -t20, t10, -t19, t8, -t2 * qJD(3) - t1 * t24, (t13 * t2 + t15 * t3) * qJD(2), t3 * qJD(3) - t1 * t25, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t4;
