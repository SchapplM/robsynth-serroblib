% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:57
% EndTime: 2019-12-31 20:21:57
% DurationCPUTime: 0.17s
% Computational Cost: add. (470->52), mult. (1178->127), div. (0->0), fcn. (801->8), ass. (0->44)
t45 = qJD(1) ^ 2;
t56 = t45 / 0.2e1;
t55 = cos(qJ(4));
t54 = cos(qJ(5));
t44 = cos(qJ(2));
t53 = t44 * t45;
t52 = pkin(6) + qJ(3);
t39 = sin(pkin(9));
t40 = cos(pkin(9));
t50 = qJD(1) * t44;
t43 = sin(qJ(2));
t51 = qJD(1) * t43;
t26 = t39 * t51 - t40 * t50;
t28 = (t39 * t44 + t40 * t43) * qJD(1);
t33 = qJD(3) + (-pkin(2) * t44 - pkin(1)) * qJD(1);
t13 = t26 * pkin(3) - t28 * pkin(7) + t33;
t31 = qJD(2) * pkin(2) - t52 * t51;
t32 = t52 * t50;
t18 = t39 * t31 + t40 * t32;
t16 = qJD(2) * pkin(7) + t18;
t42 = sin(qJ(4));
t6 = t42 * t13 + t55 * t16;
t49 = qJD(1) * qJD(2);
t48 = t43 * t49;
t47 = t44 * t49;
t5 = t55 * t13 - t42 * t16;
t17 = t40 * t31 - t39 * t32;
t25 = qJD(4) + t26;
t15 = -qJD(2) * pkin(3) - t17;
t41 = sin(qJ(5));
t38 = t44 ^ 2;
t37 = t43 ^ 2;
t36 = qJD(2) ^ 2 / 0.2e1;
t23 = qJD(5) + t25;
t22 = t42 * qJD(2) + t55 * t28;
t20 = -t55 * qJD(2) + t42 * t28;
t10 = -t41 * t20 + t54 * t22;
t8 = t54 * t20 + t41 * t22;
t7 = t20 * pkin(4) + t15;
t4 = -t20 * pkin(8) + t6;
t3 = t25 * pkin(4) - t22 * pkin(8) + t5;
t2 = t41 * t3 + t54 * t4;
t1 = t54 * t3 - t41 * t4;
t9 = [0, 0, 0, 0, 0, t56, 0, 0, 0, 0, t37 * t56, t43 * t53, t48, t38 * t56, t47, t36, pkin(1) * t53 - pkin(6) * t48, -t45 * pkin(1) * t43 - pkin(6) * t47, (t37 + t38) * t45 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t38 / 0.2e1 + t37 / 0.2e1) * pkin(6) ^ 2) * t45, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * qJD(2), t26 ^ 2 / 0.2e1, -t26 * qJD(2), t36, t17 * qJD(2) + t33 * t26, -t18 * qJD(2) + t33 * t28, -t17 * t28 - t18 * t26, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t25, t20 ^ 2 / 0.2e1, -t20 * t25, t25 ^ 2 / 0.2e1, t15 * t20 + t5 * t25, t15 * t22 - t6 * t25, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t23, t8 ^ 2 / 0.2e1, -t8 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t7 * t8, t7 * t10 - t2 * t23, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
