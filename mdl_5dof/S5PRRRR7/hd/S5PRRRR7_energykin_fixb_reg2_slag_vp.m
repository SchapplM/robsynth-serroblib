% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:57
% EndTime: 2019-12-05 17:12:58
% DurationCPUTime: 0.16s
% Computational Cost: add. (249->41), mult. (608->111), div. (0->0), fcn. (395->8), ass. (0->39)
t35 = qJD(2) ^ 2;
t46 = t35 / 0.2e1;
t45 = cos(qJ(4));
t44 = cos(qJ(5));
t31 = sin(qJ(3));
t32 = sin(qJ(2));
t22 = qJD(2) * pkin(6) + t32 * qJD(1);
t37 = pkin(7) * qJD(2) + t22;
t14 = qJD(3) * pkin(3) - t31 * t37;
t33 = cos(qJ(3));
t16 = t37 * t33;
t30 = sin(qJ(4));
t6 = t30 * t14 + t45 * t16;
t43 = qJD(2) * t31;
t42 = qJD(2) * t33;
t41 = qJD(3) * t22;
t34 = cos(qJ(2));
t40 = t34 * qJD(1);
t39 = qJD(1) * qJD(2);
t38 = qJD(2) * qJD(3);
t26 = qJD(3) + qJD(4);
t5 = t45 * t14 - t30 * t16;
t20 = -t40 + (-pkin(3) * t33 - pkin(2)) * qJD(2);
t36 = qJD(1) ^ 2;
t29 = sin(qJ(5));
t28 = t33 ^ 2;
t27 = t31 ^ 2;
t25 = qJD(5) + t26;
t23 = -qJD(2) * pkin(2) - t40;
t19 = (t30 * t33 + t45 * t31) * qJD(2);
t17 = t30 * t43 - t45 * t42;
t10 = t17 * pkin(4) + t20;
t9 = -t29 * t17 + t44 * t19;
t7 = t44 * t17 + t29 * t19;
t4 = -t17 * pkin(8) + t6;
t3 = t26 * pkin(4) - t19 * pkin(8) + t5;
t2 = t29 * t3 + t44 * t4;
t1 = -t29 * t4 + t44 * t3;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t36 / 0.2e1, 0, 0, 0, 0, 0, t46, t34 * t39, -t32 * t39, 0, (t32 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1) * t36, t27 * t46, t31 * t35 * t33, t31 * t38, t28 * t46, t33 * t38, qJD(3) ^ 2 / 0.2e1, -t23 * t42 - t31 * t41, t23 * t43 - t33 * t41, (t27 + t28) * t22 * qJD(2), t23 ^ 2 / 0.2e1 + (t28 / 0.2e1 + t27 / 0.2e1) * t22 ^ 2, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * t26, t17 ^ 2 / 0.2e1, -t17 * t26, t26 ^ 2 / 0.2e1, t20 * t17 + t5 * t26, t20 * t19 - t6 * t26, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t25, t7 ^ 2 / 0.2e1, -t7 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t10 * t7, t10 * t9 - t2 * t25, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
