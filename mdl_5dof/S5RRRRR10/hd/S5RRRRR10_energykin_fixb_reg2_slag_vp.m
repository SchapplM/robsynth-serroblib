% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:35:56
% EndTime: 2019-12-31 22:35:57
% DurationCPUTime: 0.20s
% Computational Cost: add. (745->56), mult. (1817->137), div. (0->0), fcn. (1380->10), ass. (0->46)
t50 = sin(qJ(2));
t51 = cos(qJ(2));
t45 = sin(pkin(5));
t59 = qJD(1) * t45;
t54 = t51 * t59;
t58 = cos(pkin(5)) * qJD(1);
t56 = pkin(1) * t58;
t33 = pkin(7) * t54 + t50 * t56;
t43 = qJD(2) + t58;
t26 = t43 * pkin(8) + t33;
t28 = (-pkin(2) * t51 - pkin(8) * t50 - pkin(1)) * t59;
t49 = sin(qJ(3));
t63 = cos(qJ(3));
t16 = t63 * t26 + t49 * t28;
t55 = t50 * t59;
t29 = -t63 * t43 + t49 * t55;
t11 = -t29 * pkin(9) + t16;
t48 = sin(qJ(4));
t62 = cos(qJ(4));
t15 = -t49 * t26 + t63 * t28;
t31 = t49 * t43 + t63 * t55;
t38 = -qJD(3) + t54;
t9 = -t38 * pkin(3) - t31 * pkin(9) + t15;
t6 = t62 * t11 + t48 * t9;
t61 = cos(qJ(5));
t52 = qJD(1) ^ 2;
t60 = t45 ^ 2 * t52;
t57 = t51 * t60;
t53 = t60 / 0.2e1;
t18 = t62 * t29 + t48 * t31;
t32 = -pkin(7) * t55 + t51 * t56;
t5 = -t48 * t11 + t62 * t9;
t25 = -t43 * pkin(2) - t32;
t21 = t29 * pkin(3) + t25;
t47 = sin(qJ(5));
t35 = -qJD(4) + t38;
t20 = -t48 * t29 + t62 * t31;
t17 = qJD(5) + t18;
t14 = t61 * t20 - t47 * t35;
t12 = t47 * t20 + t61 * t35;
t7 = t18 * pkin(4) - t20 * pkin(10) + t21;
t4 = -t35 * pkin(10) + t6;
t3 = t35 * pkin(4) - t5;
t2 = t61 * t4 + t47 * t7;
t1 = -t47 * t4 + t61 * t7;
t8 = [0, 0, 0, 0, 0, t52 / 0.2e1, 0, 0, 0, 0, t50 ^ 2 * t53, t50 * t57, t43 * t55, t51 ^ 2 * t53, t43 * t54, t43 ^ 2 / 0.2e1, pkin(1) * t57 + t32 * t43, -pkin(1) * t50 * t60 - t33 * t43, (-t32 * t50 + t33 * t51) * t59, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t53, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t38, t29 ^ 2 / 0.2e1, t29 * t38, t38 ^ 2 / 0.2e1, -t15 * t38 + t25 * t29, t16 * t38 + t25 * t31, -t15 * t31 - t16 * t29, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, -t20 * t35, t18 ^ 2 / 0.2e1, t18 * t35, t35 ^ 2 / 0.2e1, t21 * t18 - t5 * t35, t21 * t20 + t6 * t35, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t17, t12 ^ 2 / 0.2e1, -t12 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t12, t3 * t14 - t2 * t17, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
