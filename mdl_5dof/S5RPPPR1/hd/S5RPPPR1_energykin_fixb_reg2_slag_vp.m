% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:29:04
% EndTime: 2019-12-05 17:29:04
% DurationCPUTime: 0.17s
% Computational Cost: add. (219->41), mult. (566->108), div. (0->0), fcn. (326->8), ass. (0->37)
t34 = qJD(1) ^ 2;
t25 = t34 / 0.2e1;
t44 = pkin(1) * t34;
t27 = sin(pkin(8));
t30 = cos(pkin(8));
t31 = cos(pkin(7));
t36 = -pkin(1) * t31 - pkin(2);
t10 = qJD(3) + (-pkin(3) * t30 - qJ(4) * t27 + t36) * qJD(1);
t28 = sin(pkin(7));
t20 = (pkin(1) * t28 + qJ(3)) * qJD(1);
t17 = t27 * qJD(2) + t30 * t20;
t26 = sin(pkin(9));
t29 = cos(pkin(9));
t6 = t26 * t10 + t29 * t17;
t43 = t27 ^ 2 * t34;
t42 = t27 * t29;
t41 = qJD(1) * t27;
t40 = t30 * qJD(1);
t39 = t27 * t34 * t30;
t38 = t26 * t41;
t37 = t43 / 0.2e1;
t5 = t29 * t10 - t26 * t17;
t16 = t30 * qJD(2) - t27 * t20;
t15 = qJD(4) - t16;
t33 = cos(qJ(5));
t32 = sin(qJ(5));
t23 = t30 ^ 2 * t25;
t21 = -qJD(5) + t40;
t19 = t36 * qJD(1) + qJD(3);
t14 = (-t26 * t32 + t29 * t33) * t41;
t12 = (-t26 * t33 - t29 * t32) * t41;
t7 = pkin(4) * t38 + t15;
t4 = -pkin(6) * t38 + t6;
t3 = (-pkin(4) * t30 - pkin(6) * t42) * qJD(1) + t5;
t2 = t32 * t3 + t33 * t4;
t1 = t33 * t3 - t32 * t4;
t8 = [0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t31 * t44, -t28 * t44, 0, qJD(2) ^ 2 / 0.2e1 + (t28 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t34, t37, t39, 0, t23, 0, 0, -t19 * t40, t19 * t41, (-t16 * t27 + t17 * t30) * qJD(1), t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t29 ^ 2 * t37, -t29 * t26 * t43, -t29 * t39, t26 ^ 2 * t37, t26 * t39, t23, (t15 * t26 * t27 - t30 * t5) * qJD(1), (t15 * t42 + t30 * t6) * qJD(1), (-t26 * t6 - t29 * t5) * t41, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, t14 * t12, -t14 * t21, t12 ^ 2 / 0.2e1, -t12 * t21, t21 ^ 2 / 0.2e1, -t1 * t21 - t7 * t12, t7 * t14 + t2 * t21, -t1 * t14 + t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
