% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPPR2
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:40
% EndTime: 2019-12-05 15:24:40
% DurationCPUTime: 0.14s
% Computational Cost: add. (141->33), mult. (338->87), div. (0->0), fcn. (215->8), ass. (0->32)
t31 = qJD(2) ^ 2;
t23 = t31 / 0.2e1;
t30 = cos(qJ(2));
t17 = qJD(2) * pkin(2) + t30 * qJD(1);
t25 = sin(pkin(8));
t27 = cos(pkin(8));
t29 = sin(qJ(2));
t37 = qJD(1) * t29;
t12 = t25 * t17 + t27 * t37;
t10 = qJD(2) * qJ(4) + t12;
t24 = sin(pkin(9));
t26 = cos(pkin(9));
t6 = t24 * qJD(3) + t26 * t10;
t38 = cos(qJ(5));
t36 = qJD(2) * t24;
t35 = qJD(2) * t26;
t34 = qJD(1) * qJD(2);
t11 = t27 * t17 - t25 * t37;
t33 = qJD(4) - t11;
t32 = qJD(1) ^ 2;
t28 = sin(qJ(5));
t22 = t26 * qJD(3);
t15 = (t38 * t24 + t26 * t28) * qJD(2);
t13 = t28 * t36 - t38 * t35;
t9 = -qJD(2) * pkin(3) + t33;
t7 = (-pkin(4) * t26 - pkin(3)) * qJD(2) + t33;
t5 = -t24 * t10 + t22;
t4 = pkin(6) * t35 + t6;
t3 = t22 + (-pkin(6) * qJD(2) - t10) * t24;
t2 = t28 * t3 + t38 * t4;
t1 = -t28 * t4 + t38 * t3;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t32 / 0.2e1, 0, 0, 0, 0, 0, t23, t30 * t34, -t29 * t34, 0, (t29 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1) * t32, 0, 0, 0, 0, 0, t23, t11 * qJD(2), -t12 * qJD(2), 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t24 ^ 2 * t23, t24 * t31 * t26, 0, t26 ^ 2 * t23, 0, 0, -t9 * t35, t9 * t36, (-t24 * t5 + t26 * t6) * qJD(2), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * qJD(5), t13 ^ 2 / 0.2e1, -t13 * qJD(5), qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t7 * t13, -t2 * qJD(5) + t7 * t15, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
