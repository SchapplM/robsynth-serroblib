% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP14_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:09:49
% EndTime: 2019-03-09 13:09:49
% DurationCPUTime: 0.20s
% Computational Cost: add. (706->65), mult. (1678->136), div. (0->0), fcn. (1185->8), ass. (0->58)
t70 = -pkin(2) - pkin(9);
t51 = sin(qJ(2));
t53 = cos(qJ(2));
t47 = sin(pkin(6));
t65 = qJD(1) * t47;
t60 = t53 * t65;
t48 = cos(pkin(6));
t64 = t48 * qJD(1);
t61 = pkin(1) * t64;
t33 = pkin(8) * t60 + t51 * t61;
t45 = qJD(2) + t64;
t25 = -t45 * qJ(3) - t33;
t22 = pkin(3) * t60 - t25;
t50 = sin(qJ(4));
t52 = cos(qJ(4));
t29 = t50 * t45 + t52 * t60;
t31 = t52 * t45 - t50 * t60;
t12 = t29 * pkin(4) - t31 * pkin(10) + t22;
t49 = sin(qJ(5));
t69 = cos(qJ(5));
t59 = t51 * t65;
t41 = pkin(8) * t59;
t15 = qJD(3) + t41 + t70 * t45 + (-pkin(1) * t48 * t53 + pkin(3) * t47 * t51) * qJD(1);
t57 = -qJ(3) * t51 - pkin(1);
t23 = (t70 * t53 + t57) * t65;
t10 = t50 * t15 + t52 * t23;
t36 = qJD(4) + t59;
t8 = t36 * pkin(10) + t10;
t4 = t49 * t12 + t69 * t8;
t17 = t49 * t31 - t69 * t36;
t19 = t69 * t31 + t49 * t36;
t68 = t19 * t17;
t28 = qJD(5) + t29;
t67 = t28 * t17;
t54 = qJD(1) ^ 2;
t66 = t47 ^ 2 * t54;
t63 = t17 ^ 2 / 0.2e1;
t62 = t53 * t66;
t58 = t66 / 0.2e1;
t9 = t52 * t15 - t50 * t23;
t56 = t45 * t59;
t55 = t45 * t60;
t32 = t53 * t61 - t41;
t7 = -t36 * pkin(4) - t9;
t3 = t69 * t12 - t49 * t8;
t39 = t53 ^ 2 * t58;
t38 = t51 ^ 2 * t58;
t37 = t45 ^ 2 / 0.2e1;
t35 = t51 * t62;
t27 = (-pkin(2) * t53 + t57) * t65;
t26 = t28 ^ 2 / 0.2e1;
t24 = -t45 * pkin(2) + qJD(3) - t32;
t16 = t19 ^ 2 / 0.2e1;
t13 = t19 * t28;
t5 = t17 * pkin(5) - t19 * qJ(6) + t7;
t2 = t28 * qJ(6) + t4;
t1 = -t28 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t54 / 0.2e1, 0, 0, 0, 0, t38, t35, t56, t39, t55, t37, pkin(1) * t62 + t32 * t45, -pkin(1) * t51 * t66 - t33 * t45 (-t32 * t51 + t33 * t53) * t65, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t58, t37, -t56, -t55, t38, t35, t39 (t24 * t51 - t25 * t53) * t65, t24 * t45 + t27 * t60, -t25 * t45 - t27 * t59, t27 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t36, t29 ^ 2 / 0.2e1, -t29 * t36, t36 ^ 2 / 0.2e1, t22 * t29 + t9 * t36, -t10 * t36 + t22 * t31, -t10 * t29 - t9 * t31, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t16, -t68, t13, t63, -t67, t26, t7 * t17 + t3 * t28, t7 * t19 - t4 * t28, -t4 * t17 - t3 * t19, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t16, t13, t68, t26, t67, t63, -t1 * t28 + t5 * t17, t1 * t19 - t2 * t17, -t5 * t19 + t2 * t28, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
