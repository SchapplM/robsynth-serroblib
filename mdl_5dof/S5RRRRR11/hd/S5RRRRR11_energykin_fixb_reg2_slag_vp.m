% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR11
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:40
% EndTime: 2019-12-31 22:43:41
% DurationCPUTime: 0.20s
% Computational Cost: add. (743->56), mult. (1789->137), div. (0->0), fcn. (1355->10), ass. (0->46)
t63 = cos(qJ(4));
t62 = cos(qJ(5));
t45 = sin(pkin(5));
t53 = qJD(1) ^ 2;
t61 = t45 ^ 2 * t53;
t52 = cos(qJ(2));
t50 = sin(qJ(2));
t60 = qJD(1) * t45;
t56 = t50 * t60;
t59 = cos(pkin(5)) * qJD(1);
t57 = pkin(1) * t59;
t34 = -pkin(7) * t56 + t52 * t57;
t43 = qJD(2) + t59;
t26 = -t43 * pkin(2) - t34;
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t31 = -t51 * t43 + t49 * t56;
t33 = t49 * t43 + t51 * t56;
t13 = t31 * pkin(3) - t33 * pkin(9) + t26;
t55 = t52 * t60;
t35 = pkin(7) * t55 + t50 * t57;
t27 = t43 * pkin(8) + t35;
t29 = (-pkin(2) * t52 - pkin(8) * t50 - pkin(1)) * t60;
t18 = t51 * t27 + t49 * t29;
t38 = -qJD(3) + t55;
t16 = -t38 * pkin(9) + t18;
t48 = sin(qJ(4));
t6 = t48 * t13 + t63 * t16;
t58 = t52 * t61;
t54 = t61 / 0.2e1;
t5 = t63 * t13 - t48 * t16;
t17 = -t49 * t27 + t51 * t29;
t15 = t38 * pkin(3) - t17;
t30 = qJD(4) + t31;
t47 = sin(qJ(5));
t28 = qJD(5) + t30;
t22 = t63 * t33 - t48 * t38;
t20 = t48 * t33 + t63 * t38;
t10 = -t47 * t20 + t62 * t22;
t8 = t62 * t20 + t47 * t22;
t7 = t20 * pkin(4) + t15;
t4 = -t20 * pkin(10) + t6;
t3 = t30 * pkin(4) - t22 * pkin(10) + t5;
t2 = t47 * t3 + t62 * t4;
t1 = t62 * t3 - t47 * t4;
t9 = [0, 0, 0, 0, 0, t53 / 0.2e1, 0, 0, 0, 0, t50 ^ 2 * t54, t50 * t58, t43 * t56, t52 ^ 2 * t54, t43 * t55, t43 ^ 2 / 0.2e1, pkin(1) * t58 + t34 * t43, -pkin(1) * t50 * t61 - t35 * t43, (-t34 * t50 + t35 * t52) * t60, t35 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t54, t33 ^ 2 / 0.2e1, -t33 * t31, -t33 * t38, t31 ^ 2 / 0.2e1, t31 * t38, t38 ^ 2 / 0.2e1, -t17 * t38 + t26 * t31, t18 * t38 + t26 * t33, -t17 * t33 - t18 * t31, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t30, t20 ^ 2 / 0.2e1, -t20 * t30, t30 ^ 2 / 0.2e1, t15 * t20 + t5 * t30, t15 * t22 - t6 * t30, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t28, t8 ^ 2 / 0.2e1, -t8 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t8, t7 * t10 - t2 * t28, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
