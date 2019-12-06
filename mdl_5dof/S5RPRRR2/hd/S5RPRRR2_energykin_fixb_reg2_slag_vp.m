% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:12
% EndTime: 2019-12-05 18:12:12
% DurationCPUTime: 0.17s
% Computational Cost: add. (489->51), mult. (1359->122), div. (0->0), fcn. (997->8), ass. (0->40)
t44 = qJD(1) ^ 2;
t52 = t44 / 0.2e1;
t51 = cos(qJ(3));
t50 = cos(qJ(4));
t49 = cos(qJ(5));
t48 = pkin(6) + qJ(2);
t39 = sin(pkin(9));
t47 = qJD(1) * t39;
t30 = t48 * t47;
t40 = cos(pkin(9));
t46 = qJD(1) * t40;
t31 = t48 * t46;
t43 = sin(qJ(3));
t20 = -t51 * t30 - t43 * t31;
t29 = (t51 * t39 + t40 * t43) * qJD(1);
t14 = qJD(3) * pkin(3) - t29 * pkin(7) + t20;
t21 = -t43 * t30 + t51 * t31;
t27 = t43 * t47 - t51 * t46;
t15 = -t27 * pkin(7) + t21;
t42 = sin(qJ(4));
t6 = t42 * t14 + t50 * t15;
t38 = qJD(3) + qJD(4);
t5 = t50 * t14 - t42 * t15;
t32 = qJD(2) + (-pkin(2) * t40 - pkin(1)) * qJD(1);
t22 = t27 * pkin(3) + t32;
t41 = sin(qJ(5));
t37 = t40 ^ 2;
t36 = t39 ^ 2;
t35 = qJD(5) + t38;
t34 = -qJD(1) * pkin(1) + qJD(2);
t19 = -t42 * t27 + t50 * t29;
t17 = t50 * t27 + t42 * t29;
t10 = t17 * pkin(4) + t22;
t9 = -t41 * t17 + t49 * t19;
t7 = t49 * t17 + t41 * t19;
t4 = -t17 * pkin(8) + t6;
t3 = t38 * pkin(4) - t19 * pkin(8) + t5;
t2 = t41 * t3 + t49 * t4;
t1 = t49 * t3 - t41 * t4;
t8 = [0, 0, 0, 0, 0, t52, 0, 0, 0, 0, t36 * t52, t39 * t44 * t40, 0, t37 * t52, 0, 0, -t34 * t46, t34 * t47, (t36 + t37) * t44 * qJ(2), t34 ^ 2 / 0.2e1 + (t37 / 0.2e1 + t36 / 0.2e1) * qJ(2) ^ 2 * t44, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * qJD(3), t27 ^ 2 / 0.2e1, -t27 * qJD(3), qJD(3) ^ 2 / 0.2e1, t20 * qJD(3) + t32 * t27, -t21 * qJD(3) + t32 * t29, -t20 * t29 - t21 * t27, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * t38, t17 ^ 2 / 0.2e1, -t17 * t38, t38 ^ 2 / 0.2e1, t22 * t17 + t5 * t38, t22 * t19 - t6 * t38, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t35, t7 ^ 2 / 0.2e1, -t7 * t35, t35 ^ 2 / 0.2e1, t1 * t35 + t10 * t7, t10 * t9 - t2 * t35, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
