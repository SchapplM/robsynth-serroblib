% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:58
% EndTime: 2020-01-03 11:55:58
% DurationCPUTime: 0.11s
% Computational Cost: add. (237->34), mult. (396->89), div. (0->0), fcn. (215->8), ass. (0->34)
t25 = qJD(1) + qJD(2);
t24 = t25 ^ 2;
t21 = t24 / 0.2e1;
t32 = cos(qJ(2));
t38 = pkin(1) * qJD(1);
t36 = t32 * t38;
t17 = t25 * pkin(2) + t36;
t27 = sin(pkin(8));
t29 = cos(pkin(8));
t31 = sin(qJ(2));
t37 = t31 * t38;
t12 = t27 * t17 + t29 * t37;
t10 = t25 * qJ(4) + t12;
t26 = sin(pkin(9));
t28 = cos(pkin(9));
t6 = t26 * qJD(3) + t28 * t10;
t41 = cos(qJ(5));
t40 = t25 * t26;
t39 = t25 * t28;
t11 = t29 * t17 - t27 * t37;
t35 = qJD(4) - t11;
t33 = qJD(1) ^ 2;
t30 = sin(qJ(5));
t23 = t28 * qJD(3);
t15 = (t41 * t26 + t28 * t30) * t25;
t13 = t30 * t40 - t41 * t39;
t9 = -t25 * pkin(3) + t35;
t7 = (-pkin(4) * t28 - pkin(3)) * t25 + t35;
t5 = -t26 * t10 + t23;
t4 = pkin(7) * t39 + t6;
t3 = t23 + (-pkin(7) * t25 - t10) * t26;
t2 = t30 * t3 + t41 * t4;
t1 = t41 * t3 - t30 * t4;
t8 = [0, 0, 0, 0, 0, t33 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t25 * t36, -t25 * t37, 0, (t31 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t33, 0, 0, 0, 0, 0, t21, t11 * t25, -t12 * t25, 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t26 ^ 2 * t21, t26 * t24 * t28, 0, t28 ^ 2 * t21, 0, 0, -t9 * t39, t9 * t40, (-t26 * t5 + t28 * t6) * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * qJD(5), t13 ^ 2 / 0.2e1, -t13 * qJD(5), qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t7 * t13, -t2 * qJD(5) + t7 * t15, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
