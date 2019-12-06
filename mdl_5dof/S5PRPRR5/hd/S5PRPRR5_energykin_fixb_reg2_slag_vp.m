% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRR5
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:40
% EndTime: 2019-12-05 15:54:40
% DurationCPUTime: 0.14s
% Computational Cost: add. (221->39), mult. (576->102), div. (0->0), fcn. (389->8), ass. (0->36)
t34 = qJD(2) ^ 2;
t43 = t34 / 0.2e1;
t42 = cos(qJ(4));
t41 = cos(qJ(5));
t28 = sin(pkin(9));
t32 = sin(qJ(2));
t23 = qJD(2) * qJ(3) + t32 * qJD(1);
t37 = pkin(6) * qJD(2) + t23;
t15 = t37 * t28;
t29 = cos(pkin(9));
t16 = t37 * t29;
t31 = sin(qJ(4));
t6 = -t31 * t15 + t42 * t16;
t40 = qJD(2) * t28;
t39 = qJD(2) * t29;
t38 = qJD(1) * qJD(2);
t5 = -t42 * t15 - t31 * t16;
t33 = cos(qJ(2));
t36 = -t33 * qJD(1) + qJD(3);
t20 = (-pkin(3) * t29 - pkin(2)) * qJD(2) + t36;
t35 = qJD(1) ^ 2;
t30 = sin(qJ(5));
t27 = qJD(4) + qJD(5);
t26 = t29 ^ 2;
t25 = t28 ^ 2;
t21 = -qJD(2) * pkin(2) + t36;
t19 = (t42 * t28 + t29 * t31) * qJD(2);
t17 = t31 * t40 - t42 * t39;
t10 = t17 * pkin(4) + t20;
t9 = -t30 * t17 + t41 * t19;
t7 = t41 * t17 + t30 * t19;
t4 = -t17 * pkin(7) + t6;
t3 = qJD(4) * pkin(4) - t19 * pkin(7) + t5;
t2 = t30 * t3 + t41 * t4;
t1 = t41 * t3 - t30 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t35 / 0.2e1, 0, 0, 0, 0, 0, t43, t33 * t38, -t32 * t38, 0, (t32 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1) * t35, t25 * t43, t28 * t34 * t29, 0, t26 * t43, 0, 0, -t21 * t39, t21 * t40, (t25 + t26) * t23 * qJD(2), t21 ^ 2 / 0.2e1 + (t26 / 0.2e1 + t25 / 0.2e1) * t23 ^ 2, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * qJD(4), t17 ^ 2 / 0.2e1, -t17 * qJD(4), qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) + t20 * t17, -t6 * qJD(4) + t20 * t19, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t27, t7 ^ 2 / 0.2e1, -t7 * t27, t27 ^ 2 / 0.2e1, t1 * t27 + t10 * t7, t10 * t9 - t2 * t27, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
