% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:33
% EndTime: 2019-12-31 21:40:34
% DurationCPUTime: 0.20s
% Computational Cost: add. (733->56), mult. (1789->134), div. (0->0), fcn. (1355->10), ass. (0->46)
t63 = cos(qJ(3));
t62 = cos(qJ(5));
t45 = sin(pkin(5));
t51 = qJD(1) ^ 2;
t61 = t45 ^ 2 * t51;
t50 = cos(qJ(2));
t49 = sin(qJ(2));
t59 = qJD(1) * t45;
t54 = t49 * t59;
t58 = cos(pkin(5)) * qJD(1);
t55 = pkin(1) * t58;
t33 = -pkin(7) * t54 + t50 * t55;
t42 = qJD(2) + t58;
t26 = -t42 * pkin(2) - t33;
t48 = sin(qJ(3));
t30 = -t63 * t42 + t48 * t54;
t32 = t48 * t42 + t63 * t54;
t13 = t30 * pkin(3) - t32 * qJ(4) + t26;
t53 = t50 * t59;
t34 = pkin(7) * t53 + t49 * t55;
t27 = t42 * pkin(8) + t34;
t28 = (-pkin(2) * t50 - pkin(8) * t49 - pkin(1)) * t59;
t18 = t63 * t27 + t48 * t28;
t37 = -qJD(3) + t53;
t16 = -t37 * qJ(4) + t18;
t44 = sin(pkin(10));
t60 = cos(pkin(10));
t6 = t44 * t13 + t60 * t16;
t57 = t30 ^ 2 / 0.2e1;
t56 = t50 * t61;
t52 = t61 / 0.2e1;
t5 = t60 * t13 - t44 * t16;
t17 = -t48 * t27 + t63 * t28;
t15 = t37 * pkin(3) + qJD(4) - t17;
t47 = sin(qJ(5));
t29 = qJD(5) + t30;
t22 = t60 * t32 - t44 * t37;
t20 = t44 * t32 + t60 * t37;
t10 = -t47 * t20 + t62 * t22;
t8 = t62 * t20 + t47 * t22;
t7 = t20 * pkin(4) + t15;
t4 = -t20 * pkin(9) + t6;
t3 = t30 * pkin(4) - t22 * pkin(9) + t5;
t2 = t47 * t3 + t62 * t4;
t1 = t62 * t3 - t47 * t4;
t9 = [0, 0, 0, 0, 0, t51 / 0.2e1, 0, 0, 0, 0, t49 ^ 2 * t52, t49 * t56, t42 * t54, t50 ^ 2 * t52, t42 * t53, t42 ^ 2 / 0.2e1, pkin(1) * t56 + t33 * t42, -pkin(1) * t49 * t61 - t34 * t42, (-t33 * t49 + t34 * t50) * t59, t34 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t52, t32 ^ 2 / 0.2e1, -t32 * t30, -t32 * t37, t57, t30 * t37, t37 ^ 2 / 0.2e1, -t17 * t37 + t26 * t30, t18 * t37 + t26 * t32, -t17 * t32 - t18 * t30, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t30, t20 ^ 2 / 0.2e1, -t20 * t30, t57, t15 * t20 + t5 * t30, t15 * t22 - t6 * t30, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t29, t8 ^ 2 / 0.2e1, -t8 * t29, t29 ^ 2 / 0.2e1, t1 * t29 + t7 * t8, t7 * t10 - t2 * t29, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
