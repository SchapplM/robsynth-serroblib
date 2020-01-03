% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR14
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR14_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:36
% EndTime: 2019-12-31 20:38:36
% DurationCPUTime: 0.20s
% Computational Cost: add. (715->56), mult. (1819->134), div. (0->0), fcn. (1380->10), ass. (0->46)
t49 = sin(qJ(2));
t50 = cos(qJ(2));
t45 = sin(pkin(5));
t58 = qJD(1) * t45;
t53 = t50 * t58;
t57 = cos(pkin(5)) * qJD(1);
t55 = pkin(1) * t57;
t33 = pkin(7) * t53 + t49 * t55;
t42 = qJD(2) + t57;
t26 = t42 * qJ(3) + t33;
t28 = (-pkin(2) * t50 - qJ(3) * t49 - pkin(1)) * t58;
t44 = sin(pkin(10));
t59 = cos(pkin(10));
t16 = t59 * t26 + t44 * t28;
t54 = t49 * t58;
t29 = -t59 * t42 + t44 * t54;
t11 = -t29 * pkin(8) + t16;
t48 = sin(qJ(4));
t62 = cos(qJ(4));
t15 = -t44 * t26 + t59 * t28;
t31 = t44 * t42 + t59 * t54;
t9 = -pkin(3) * t53 - t31 * pkin(8) + t15;
t6 = t62 * t11 + t48 * t9;
t61 = cos(qJ(5));
t51 = qJD(1) ^ 2;
t60 = t45 ^ 2 * t51;
t56 = t50 * t60;
t52 = t60 / 0.2e1;
t19 = t62 * t29 + t48 * t31;
t32 = -pkin(7) * t54 + t50 * t55;
t5 = -t48 * t11 + t62 * t9;
t23 = -t42 * pkin(2) + qJD(3) - t32;
t18 = t29 * pkin(3) + t23;
t47 = sin(qJ(5));
t38 = t50 ^ 2 * t52;
t36 = -qJD(4) + t53;
t21 = -t48 * t29 + t62 * t31;
t17 = qJD(5) + t19;
t14 = t61 * t21 - t47 * t36;
t12 = t47 * t21 + t61 * t36;
t7 = t19 * pkin(4) - t21 * pkin(9) + t18;
t4 = -t36 * pkin(9) + t6;
t3 = t36 * pkin(4) - t5;
t2 = t61 * t4 + t47 * t7;
t1 = -t47 * t4 + t61 * t7;
t8 = [0, 0, 0, 0, 0, t51 / 0.2e1, 0, 0, 0, 0, t49 ^ 2 * t52, t49 * t56, t42 * t54, t38, t42 * t53, t42 ^ 2 / 0.2e1, pkin(1) * t56 + t32 * t42, -pkin(1) * t49 * t60 - t33 * t42, (-t32 * t49 + t33 * t50) * t58, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t52, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t53, t29 ^ 2 / 0.2e1, t29 * t53, t38, -t15 * t53 + t23 * t29, t16 * t53 + t23 * t31, -t15 * t31 - t16 * t29, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, -t21 * t36, t19 ^ 2 / 0.2e1, t19 * t36, t36 ^ 2 / 0.2e1, t18 * t19 - t5 * t36, t18 * t21 + t6 * t36, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t17, t12 ^ 2 / 0.2e1, -t12 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t12, t3 * t14 - t2 * t17, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
