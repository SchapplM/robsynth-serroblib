% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:40
% EndTime: 2019-12-31 17:15:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (126->32), mult. (375->77), div. (0->0), fcn. (208->4), ass. (0->33)
t28 = qJD(1) ^ 2;
t38 = t28 / 0.2e1;
t37 = -pkin(6) - pkin(5);
t26 = sin(qJ(2));
t34 = qJD(1) * t26;
t17 = qJD(2) * pkin(2) + t37 * t34;
t27 = cos(qJ(2));
t33 = qJD(1) * t27;
t18 = t37 * t33;
t25 = sin(qJ(3));
t36 = cos(qJ(3));
t4 = t25 * t17 - t36 * t18;
t35 = t27 * t28;
t32 = qJD(1) * qJD(2);
t31 = t26 * t32;
t30 = t27 * t32;
t3 = t36 * t17 + t25 * t18;
t19 = (-pkin(2) * t27 - pkin(1)) * qJD(1);
t24 = t27 ^ 2;
t23 = t26 ^ 2;
t22 = qJD(2) + qJD(3);
t21 = t22 ^ 2 / 0.2e1;
t15 = (t25 * t27 + t36 * t26) * qJD(1);
t13 = t25 * t34 - t36 * t33;
t12 = t15 ^ 2 / 0.2e1;
t11 = t13 ^ 2 / 0.2e1;
t8 = t15 * t22;
t7 = t13 * t22;
t6 = t13 * pkin(3) + qJD(4) + t19;
t5 = t15 * t13;
t2 = -t13 * qJ(4) + t4;
t1 = t22 * pkin(3) - t15 * qJ(4) + t3;
t9 = [0, 0, 0, 0, 0, t38, 0, 0, 0, 0, t23 * t38, t26 * t35, t31, t24 * t38, t30, qJD(2) ^ 2 / 0.2e1, pkin(1) * t35 - pkin(5) * t31, -t28 * pkin(1) * t26 - pkin(5) * t30, (t23 + t24) * t28 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t24 / 0.2e1 + t23 / 0.2e1) * pkin(5) ^ 2) * t28, t12, -t5, t8, t11, -t7, t21, t19 * t13 + t3 * t22, t19 * t15 - t4 * t22, -t4 * t13 - t3 * t15, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t12, -t5, t8, t11, -t7, t21, t1 * t22 + t6 * t13, t6 * t15 - t2 * t22, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg = t9;
