% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:46
% EndTime: 2019-12-31 22:25:46
% DurationCPUTime: 0.17s
% Computational Cost: add. (523->52), mult. (1178->130), div. (0->0), fcn. (801->8), ass. (0->44)
t45 = qJD(1) ^ 2;
t56 = t45 / 0.2e1;
t55 = -pkin(7) - pkin(6);
t54 = cos(qJ(4));
t53 = cos(qJ(5));
t44 = cos(qJ(2));
t52 = t44 * t45;
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t50 = qJD(1) * t44;
t42 = sin(qJ(2));
t51 = qJD(1) * t42;
t26 = t41 * t51 - t43 * t50;
t28 = (t41 * t44 + t42 * t43) * qJD(1);
t33 = (-pkin(2) * t44 - pkin(1)) * qJD(1);
t13 = t26 * pkin(3) - t28 * pkin(8) + t33;
t31 = qJD(2) * pkin(2) + t55 * t51;
t32 = t55 * t50;
t18 = t41 * t31 - t43 * t32;
t36 = qJD(2) + qJD(3);
t16 = t36 * pkin(8) + t18;
t40 = sin(qJ(4));
t6 = t40 * t13 + t54 * t16;
t49 = qJD(1) * qJD(2);
t48 = t42 * t49;
t47 = t44 * t49;
t5 = t54 * t13 - t40 * t16;
t17 = t43 * t31 + t41 * t32;
t15 = -t36 * pkin(3) - t17;
t25 = qJD(4) + t26;
t39 = sin(qJ(5));
t38 = t44 ^ 2;
t37 = t42 ^ 2;
t23 = qJD(5) + t25;
t22 = t54 * t28 + t40 * t36;
t20 = t40 * t28 - t54 * t36;
t10 = -t39 * t20 + t53 * t22;
t8 = t53 * t20 + t39 * t22;
t7 = t20 * pkin(4) + t15;
t4 = -t20 * pkin(9) + t6;
t3 = t25 * pkin(4) - t22 * pkin(9) + t5;
t2 = t39 * t3 + t53 * t4;
t1 = t53 * t3 - t39 * t4;
t9 = [0, 0, 0, 0, 0, t56, 0, 0, 0, 0, t37 * t56, t42 * t52, t48, t38 * t56, t47, qJD(2) ^ 2 / 0.2e1, pkin(1) * t52 - pkin(6) * t48, -t45 * pkin(1) * t42 - pkin(6) * t47, (t37 + t38) * t45 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t38 / 0.2e1 + t37 / 0.2e1) * pkin(6) ^ 2) * t45, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * t36, t26 ^ 2 / 0.2e1, -t26 * t36, t36 ^ 2 / 0.2e1, t17 * t36 + t33 * t26, -t18 * t36 + t33 * t28, -t17 * t28 - t18 * t26, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t25, t20 ^ 2 / 0.2e1, -t20 * t25, t25 ^ 2 / 0.2e1, t15 * t20 + t5 * t25, t15 * t22 - t6 * t25, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t23, t8 ^ 2 / 0.2e1, -t8 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t7 * t8, t7 * t10 - t2 * t23, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
