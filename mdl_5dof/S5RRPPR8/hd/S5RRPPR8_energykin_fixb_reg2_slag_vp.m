% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:14
% EndTime: 2019-12-31 19:39:14
% DurationCPUTime: 0.17s
% Computational Cost: add. (276->50), mult. (672->115), div. (0->0), fcn. (374->6), ass. (0->41)
t44 = qJD(1) ^ 2;
t55 = t44 / 0.2e1;
t54 = cos(qJ(5));
t43 = cos(qJ(2));
t53 = t43 * t44;
t42 = sin(qJ(2));
t49 = qJ(4) * qJD(1);
t52 = qJD(1) * t42;
t50 = pkin(6) * t52 + qJD(3);
t14 = -t42 * t49 + (-pkin(2) - pkin(3)) * qJD(2) + t50;
t51 = qJD(1) * t43;
t23 = pkin(6) * t51 + qJD(2) * qJ(3);
t20 = -t43 * t49 + t23;
t38 = sin(pkin(8));
t39 = cos(pkin(8));
t6 = t38 * t14 + t39 * t20;
t48 = qJD(1) * qJD(2);
t47 = t42 * t53;
t21 = -qJD(1) * pkin(1) - pkin(2) * t51 - qJ(3) * t52;
t26 = t42 * t48;
t46 = t43 * t48;
t5 = t39 * t14 - t38 * t20;
t13 = pkin(3) * t51 + qJD(4) - t21;
t41 = sin(qJ(5));
t37 = t43 ^ 2;
t36 = t42 ^ 2;
t34 = qJD(2) ^ 2 / 0.2e1;
t32 = qJD(2) - qJD(5);
t25 = t37 * t55;
t24 = t36 * t55;
t22 = -qJD(2) * pkin(2) + t50;
t19 = (-t38 * t43 + t39 * t42) * qJD(1);
t17 = (-t38 * t42 - t39 * t43) * qJD(1);
t10 = t41 * t17 + t54 * t19;
t8 = -t54 * t17 + t41 * t19;
t7 = -t17 * pkin(4) + t13;
t4 = t17 * pkin(7) + t6;
t3 = -qJD(2) * pkin(4) - t19 * pkin(7) + t5;
t2 = t41 * t3 + t54 * t4;
t1 = t54 * t3 - t41 * t4;
t9 = [0, 0, 0, 0, 0, t55, 0, 0, 0, 0, t24, t47, t26, t25, t46, t34, pkin(1) * t53 - pkin(6) * t26, -t44 * pkin(1) * t42 - pkin(6) * t46, (t36 + t37) * t44 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t37 / 0.2e1 + t36 / 0.2e1) * pkin(6) ^ 2) * t44, t24, t26, -t47, t34, -t46, t25, -t22 * qJD(2) - t21 * t51, (t22 * t42 + t23 * t43) * qJD(1), t23 * qJD(2) - t21 * t52, t23 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, t19 * t17, -t19 * qJD(2), t17 ^ 2 / 0.2e1, -t17 * qJD(2), t34, -t5 * qJD(2) - t13 * t17, t6 * qJD(2) + t13 * t19, t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, -t10 * t32, t8 ^ 2 / 0.2e1, t8 * t32, t32 ^ 2 / 0.2e1, -t1 * t32 + t7 * t8, t7 * t10 + t2 * t32, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
