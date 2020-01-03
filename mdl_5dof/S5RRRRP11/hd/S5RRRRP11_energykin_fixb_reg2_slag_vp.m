% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:30
% EndTime: 2019-12-31 22:18:30
% DurationCPUTime: 0.20s
% Computational Cost: add. (532->50), mult. (1319->117), div. (0->0), fcn. (966->8), ass. (0->45)
t44 = sin(qJ(2));
t45 = cos(qJ(2));
t40 = sin(pkin(5));
t54 = qJD(1) * t40;
t48 = t45 * t54;
t53 = cos(pkin(5)) * qJD(1);
t50 = pkin(1) * t53;
t30 = pkin(7) * t48 + t44 * t50;
t38 = qJD(2) + t53;
t22 = t38 * pkin(8) + t30;
t24 = (-pkin(2) * t45 - pkin(8) * t44 - pkin(1)) * t54;
t43 = sin(qJ(3));
t59 = cos(qJ(3));
t13 = t59 * t22 + t43 * t24;
t33 = -qJD(3) + t48;
t10 = -t33 * pkin(9) + t13;
t42 = sin(qJ(4));
t58 = cos(qJ(4));
t49 = t44 * t54;
t29 = -pkin(7) * t49 + t45 * t50;
t21 = -t38 * pkin(2) - t29;
t26 = -t59 * t38 + t43 * t49;
t28 = t43 * t38 + t59 * t49;
t7 = t26 * pkin(3) - t28 * pkin(9) + t21;
t5 = t58 * t10 + t42 * t7;
t15 = t42 * t28 + t58 * t33;
t17 = t58 * t28 - t42 * t33;
t57 = t17 * t15;
t25 = qJD(4) + t26;
t56 = t25 * t15;
t46 = qJD(1) ^ 2;
t55 = t40 ^ 2 * t46;
t52 = t15 ^ 2 / 0.2e1;
t51 = t45 * t55;
t47 = t55 / 0.2e1;
t12 = -t43 * t22 + t59 * t24;
t4 = -t42 * t10 + t58 * t7;
t9 = t33 * pkin(3) - t12;
t23 = t25 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t11 = t17 * t25;
t3 = t15 * pkin(4) - t17 * qJ(5) + t9;
t2 = t25 * qJ(5) + t5;
t1 = -t25 * pkin(4) + qJD(5) - t4;
t6 = [0, 0, 0, 0, 0, t46 / 0.2e1, 0, 0, 0, 0, t44 ^ 2 * t47, t44 * t51, t38 * t49, t45 ^ 2 * t47, t38 * t48, t38 ^ 2 / 0.2e1, pkin(1) * t51 + t29 * t38, -pkin(1) * t44 * t55 - t30 * t38, (-t29 * t44 + t30 * t45) * t54, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t47, t28 ^ 2 / 0.2e1, -t28 * t26, -t28 * t33, t26 ^ 2 / 0.2e1, t26 * t33, t33 ^ 2 / 0.2e1, -t12 * t33 + t21 * t26, t13 * t33 + t21 * t28, -t12 * t28 - t13 * t26, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14, -t57, t11, t52, -t56, t23, t9 * t15 + t4 * t25, t9 * t17 - t5 * t25, -t5 * t15 - t4 * t17, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t14, t11, t57, t23, t56, t52, -t1 * t25 + t3 * t15, t1 * t17 - t2 * t15, -t3 * t17 + t2 * t25, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
