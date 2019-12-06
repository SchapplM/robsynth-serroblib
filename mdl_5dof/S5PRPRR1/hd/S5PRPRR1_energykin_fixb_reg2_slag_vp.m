% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:03
% EndTime: 2019-12-05 15:43:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (210->39), mult. (555->95), div. (0->0), fcn. (378->6), ass. (0->32)
t34 = qJD(2) ^ 2;
t40 = t34 / 0.2e1;
t39 = cos(qJ(4));
t38 = cos(qJ(5));
t31 = cos(pkin(9));
t27 = t31 * qJD(1);
t30 = sin(pkin(9));
t37 = qJD(2) * t30;
t15 = t27 + (-pkin(6) - qJ(3)) * t37;
t35 = qJ(3) * qJD(2);
t21 = t30 * qJD(1) + t31 * t35;
t36 = qJD(2) * t31;
t16 = pkin(6) * t36 + t21;
t33 = sin(qJ(4));
t6 = t33 * t15 + t39 * t16;
t5 = t39 * t15 - t33 * t16;
t22 = qJD(3) + (-pkin(3) * t31 - pkin(2)) * qJD(2);
t32 = sin(qJ(5));
t29 = qJD(1) ^ 2 / 0.2e1;
t28 = qJD(4) + qJD(5);
t25 = -qJD(2) * pkin(2) + qJD(3);
t20 = -t30 * t35 + t27;
t19 = (t39 * t30 + t31 * t33) * qJD(2);
t17 = t33 * t37 - t39 * t36;
t10 = t17 * pkin(4) + t22;
t9 = -t32 * t17 + t38 * t19;
t7 = t38 * t17 + t32 * t19;
t4 = -t17 * pkin(7) + t6;
t3 = qJD(4) * pkin(4) - t19 * pkin(7) + t5;
t2 = t32 * t3 + t38 * t4;
t1 = t38 * t3 - t32 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, t40, 0, 0, 0, t29, t30 ^ 2 * t40, t30 * t34 * t31, 0, t31 ^ 2 * t40, 0, 0, -t25 * t36, t25 * t37, (-t20 * t30 + t21 * t31) * qJD(2), t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * qJD(4), t17 ^ 2 / 0.2e1, -t17 * qJD(4), qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) + t22 * t17, -t6 * qJD(4) + t22 * t19, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t28, t7 ^ 2 / 0.2e1, -t7 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t10 * t7, t10 * t9 - t2 * t28, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
