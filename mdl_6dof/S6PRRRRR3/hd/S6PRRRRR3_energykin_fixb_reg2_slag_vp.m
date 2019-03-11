% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:51:31
% EndTime: 2019-03-09 00:51:31
% DurationCPUTime: 0.18s
% Computational Cost: add. (721->61), mult. (1606->152), div. (0->0), fcn. (1185->12), ass. (0->51)
t57 = qJD(2) ^ 2;
t69 = t57 / 0.2e1;
t54 = sin(qJ(2));
t48 = sin(pkin(6));
t65 = qJD(1) * t48;
t37 = qJD(2) * pkin(8) + t54 * t65;
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t49 = cos(pkin(6));
t64 = qJD(1) * t49;
t29 = t55 * t37 + t53 * t64;
t25 = qJD(3) * pkin(9) + t29;
t56 = cos(qJ(2));
t60 = t56 * t65;
t30 = -t60 + (-pkin(3) * t55 - pkin(9) * t53 - pkin(2)) * qJD(2);
t52 = sin(qJ(4));
t68 = cos(qJ(4));
t17 = t68 * t25 + t52 * t30;
t63 = qJD(2) * t53;
t34 = -t68 * qJD(3) + t52 * t63;
t15 = -t34 * pkin(10) + t17;
t51 = sin(qJ(5));
t67 = cos(qJ(5));
t16 = -t52 * t25 + t68 * t30;
t36 = t52 * qJD(3) + t68 * t63;
t62 = t55 * qJD(2);
t44 = -qJD(4) + t62;
t9 = -t44 * pkin(4) - t36 * pkin(10) + t16;
t6 = t67 * t15 + t51 * t9;
t66 = cos(qJ(6));
t61 = qJD(2) * qJD(3);
t5 = -t51 * t15 + t67 * t9;
t59 = qJD(2) * t65;
t28 = -t53 * t37 + t55 * t64;
t41 = -qJD(5) + t44;
t24 = -qJD(3) * pkin(3) - t28;
t18 = t34 * pkin(4) + t24;
t58 = qJD(1) ^ 2;
t50 = sin(qJ(6));
t39 = -qJD(6) + t41;
t38 = -qJD(2) * pkin(2) - t60;
t22 = -t51 * t34 + t67 * t36;
t20 = t67 * t34 + t51 * t36;
t14 = -t50 * t20 + t66 * t22;
t12 = t66 * t20 + t50 * t22;
t10 = t20 * pkin(5) + t18;
t4 = -t20 * pkin(11) + t6;
t3 = -t41 * pkin(5) - t22 * pkin(11) + t5;
t2 = t50 * t3 + t66 * t4;
t1 = t66 * t3 - t50 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t58 / 0.2e1, 0, 0, 0, 0, 0, t69, t56 * t59, -t54 * t59, 0 (t49 ^ 2 / 0.2e1 + (t54 ^ 2 / 0.2e1 + t56 ^ 2 / 0.2e1) * t48 ^ 2) * t58, t53 ^ 2 * t69, t53 * t57 * t55, t53 * t61, t55 ^ 2 * t69, t55 * t61, qJD(3) ^ 2 / 0.2e1, t28 * qJD(3) - t38 * t62, -t29 * qJD(3) + t38 * t63 (-t28 * t53 + t29 * t55) * qJD(2), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t36 ^ 2 / 0.2e1, -t36 * t34, -t36 * t44, t34 ^ 2 / 0.2e1, t34 * t44, t44 ^ 2 / 0.2e1, -t16 * t44 + t24 * t34, t17 * t44 + t24 * t36, -t16 * t36 - t17 * t34, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, -t22 * t41, t20 ^ 2 / 0.2e1, t20 * t41, t41 ^ 2 / 0.2e1, t18 * t20 - t5 * t41, t18 * t22 + t6 * t41, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, -t14 * t39, t12 ^ 2 / 0.2e1, t12 * t39, t39 ^ 2 / 0.2e1, -t1 * t39 + t10 * t12, t10 * t14 + t2 * t39, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t7;
