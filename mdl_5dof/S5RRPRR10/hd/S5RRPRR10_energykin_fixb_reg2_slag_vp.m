% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:38
% EndTime: 2019-12-31 20:26:38
% DurationCPUTime: 0.18s
% Computational Cost: add. (616->57), mult. (1792->135), div. (0->0), fcn. (1361->10), ass. (0->47)
t50 = cos(qJ(2));
t57 = cos(pkin(5)) * qJD(1);
t55 = pkin(1) * t57;
t40 = t50 * t55;
t41 = qJD(2) + t57;
t49 = sin(qJ(2));
t44 = sin(pkin(5));
t58 = qJD(1) * t44;
t54 = t49 * t58;
t23 = t41 * pkin(2) + t40 + (-pkin(7) - qJ(3)) * t54;
t53 = t50 * t58;
t33 = pkin(7) * t53 + t49 * t55;
t26 = qJ(3) * t53 + t33;
t43 = sin(pkin(10));
t45 = cos(pkin(10));
t12 = t43 * t23 + t45 * t26;
t10 = t41 * pkin(8) + t12;
t29 = t43 * t54 - t45 * t53;
t31 = (t43 * t50 + t45 * t49) * t58;
t34 = qJD(3) + (-pkin(2) * t50 - pkin(1)) * t58;
t17 = t29 * pkin(3) - t31 * pkin(8) + t34;
t48 = sin(qJ(4));
t61 = cos(qJ(4));
t7 = t61 * t10 + t48 * t17;
t60 = cos(qJ(5));
t51 = qJD(1) ^ 2;
t59 = t44 ^ 2 * t51;
t56 = t50 * t59;
t52 = t59 / 0.2e1;
t11 = t45 * t23 - t43 * t26;
t20 = t48 * t31 - t61 * t41;
t9 = -t41 * pkin(3) - t11;
t6 = -t48 * t10 + t61 * t17;
t47 = sin(qJ(5));
t37 = t41 ^ 2 / 0.2e1;
t32 = -pkin(7) * t54 + t40;
t28 = qJD(4) + t29;
t22 = t61 * t31 + t48 * t41;
t19 = qJD(5) + t20;
t15 = t60 * t22 + t47 * t28;
t13 = t47 * t22 - t60 * t28;
t5 = t20 * pkin(4) - t22 * pkin(9) + t9;
t4 = t28 * pkin(9) + t7;
t3 = -t28 * pkin(4) - t6;
t2 = t60 * t4 + t47 * t5;
t1 = -t47 * t4 + t60 * t5;
t8 = [0, 0, 0, 0, 0, t51 / 0.2e1, 0, 0, 0, 0, t49 ^ 2 * t52, t49 * t56, t41 * t54, t50 ^ 2 * t52, t41 * t53, t37, pkin(1) * t56 + t32 * t41, -pkin(1) * t49 * t59 - t33 * t41, (-t32 * t49 + t33 * t50) * t58, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t52, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t41, t29 ^ 2 / 0.2e1, -t29 * t41, t37, t11 * t41 + t34 * t29, -t12 * t41 + t34 * t31, -t11 * t31 - t12 * t29, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t28, t20 ^ 2 / 0.2e1, -t20 * t28, t28 ^ 2 / 0.2e1, t9 * t20 + t6 * t28, t9 * t22 - t7 * t28, -t7 * t20 - t6 * t22, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t19, t13 ^ 2 / 0.2e1, -t13 * t19, t19 ^ 2 / 0.2e1, t1 * t19 + t3 * t13, t3 * t15 - t2 * t19, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
