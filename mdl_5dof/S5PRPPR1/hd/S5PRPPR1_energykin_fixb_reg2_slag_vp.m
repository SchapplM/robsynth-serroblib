% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:08
% EndTime: 2019-12-05 15:22:08
% DurationCPUTime: 0.16s
% Computational Cost: add. (167->36), mult. (450->95), div. (0->0), fcn. (273->6), ass. (0->34)
t31 = qJD(2) ^ 2;
t40 = t31 / 0.2e1;
t26 = sin(pkin(8));
t39 = t26 ^ 2 * t31;
t27 = cos(pkin(9));
t38 = t26 * t27;
t28 = cos(pkin(8));
t13 = qJD(3) + (-pkin(3) * t28 - qJ(4) * t26 - pkin(2)) * qJD(2);
t35 = qJ(3) * qJD(2);
t17 = t26 * qJD(1) + t28 * t35;
t25 = sin(pkin(9));
t6 = t25 * t13 + t27 * t17;
t37 = qJD(2) * t26;
t36 = t28 * qJD(2);
t34 = t26 * t31 * t28;
t33 = t25 * t37;
t32 = t39 / 0.2e1;
t5 = t27 * t13 - t25 * t17;
t16 = t28 * qJD(1) - t26 * t35;
t15 = qJD(4) - t16;
t30 = cos(qJ(5));
t29 = sin(qJ(5));
t24 = qJD(1) ^ 2 / 0.2e1;
t22 = -qJD(2) * pkin(2) + qJD(3);
t20 = t28 ^ 2 * t40;
t18 = -qJD(5) + t36;
t10 = (-t25 * t29 + t27 * t30) * t37;
t8 = (-t25 * t30 - t27 * t29) * t37;
t7 = pkin(4) * t33 + t15;
t4 = -pkin(6) * t33 + t6;
t3 = (-pkin(4) * t28 - pkin(6) * t38) * qJD(2) + t5;
t2 = t29 * t3 + t30 * t4;
t1 = -t29 * t4 + t30 * t3;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, t40, 0, 0, 0, t24, t32, t34, 0, t20, 0, 0, -t22 * t36, t22 * t37, (-t16 * t26 + t17 * t28) * qJD(2), t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t27 ^ 2 * t32, -t27 * t25 * t39, -t27 * t34, t25 ^ 2 * t32, t25 * t34, t20, (t15 * t25 * t26 - t28 * t5) * qJD(2), (t15 * t38 + t28 * t6) * qJD(2), (-t25 * t6 - t27 * t5) * t37, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, t10 * t8, -t10 * t18, t8 ^ 2 / 0.2e1, -t8 * t18, t18 ^ 2 / 0.2e1, -t1 * t18 - t7 * t8, t7 * t10 + t2 * t18, -t1 * t10 + t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
