% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:20
% EndTime: 2019-12-05 16:23:21
% DurationCPUTime: 0.18s
% Computational Cost: add. (239->41), mult. (608->108), div. (0->0), fcn. (395->8), ass. (0->39)
t35 = qJD(2) ^ 2;
t46 = t35 / 0.2e1;
t45 = cos(qJ(5));
t31 = sin(qJ(3));
t32 = sin(qJ(2));
t22 = qJD(2) * pkin(6) + t32 * qJD(1);
t37 = qJ(4) * qJD(2) + t22;
t14 = qJD(3) * pkin(3) - t37 * t31;
t33 = cos(qJ(3));
t16 = t37 * t33;
t29 = sin(pkin(9));
t44 = cos(pkin(9));
t6 = t29 * t14 + t44 * t16;
t43 = qJD(2) * t31;
t42 = qJD(2) * t33;
t41 = qJD(3) * t22;
t34 = cos(qJ(2));
t40 = t34 * qJD(1);
t39 = qJD(1) * qJD(2);
t38 = qJD(2) * qJD(3);
t5 = t44 * t14 - t29 * t16;
t20 = -t40 + qJD(4) + (-pkin(3) * t33 - pkin(2)) * qJD(2);
t36 = qJD(1) ^ 2;
t30 = sin(qJ(5));
t28 = t33 ^ 2;
t27 = t31 ^ 2;
t26 = qJD(3) ^ 2 / 0.2e1;
t25 = qJD(3) + qJD(5);
t23 = -qJD(2) * pkin(2) - t40;
t19 = (t29 * t33 + t44 * t31) * qJD(2);
t17 = t29 * t43 - t44 * t42;
t10 = t17 * pkin(4) + t20;
t9 = -t30 * t17 + t45 * t19;
t7 = t45 * t17 + t30 * t19;
t4 = -t17 * pkin(7) + t6;
t3 = qJD(3) * pkin(4) - t19 * pkin(7) + t5;
t2 = t30 * t3 + t45 * t4;
t1 = t45 * t3 - t30 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t36 / 0.2e1, 0, 0, 0, 0, 0, t46, t34 * t39, -t32 * t39, 0, (t32 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1) * t36, t27 * t46, t31 * t35 * t33, t31 * t38, t28 * t46, t33 * t38, t26, -t23 * t42 - t31 * t41, t23 * t43 - t33 * t41, (t27 + t28) * t22 * qJD(2), t23 ^ 2 / 0.2e1 + (t28 / 0.2e1 + t27 / 0.2e1) * t22 ^ 2, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * qJD(3), t17 ^ 2 / 0.2e1, -t17 * qJD(3), t26, t5 * qJD(3) + t20 * t17, -t6 * qJD(3) + t20 * t19, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t25, t7 ^ 2 / 0.2e1, -t7 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t10 * t7, t10 * t9 - t2 * t25, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
