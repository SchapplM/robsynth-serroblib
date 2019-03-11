% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:31:25
% EndTime: 2019-03-10 01:31:25
% DurationCPUTime: 0.21s
% Computational Cost: add. (997->65), mult. (2165->144), div. (0->0), fcn. (1533->8), ass. (0->52)
t56 = qJD(1) ^ 2;
t70 = t56 / 0.2e1;
t51 = sin(qJ(5));
t67 = cos(qJ(5));
t54 = sin(qJ(2));
t55 = cos(qJ(2));
t33 = (-pkin(2) * t55 - pkin(8) * t54 - pkin(1)) * qJD(1);
t62 = t55 * qJD(1);
t42 = pkin(7) * t62 + qJD(2) * pkin(8);
t53 = sin(qJ(3));
t69 = cos(qJ(3));
t27 = t69 * t33 - t53 * t42;
t63 = qJD(1) * t54;
t36 = t53 * qJD(2) + t69 * t63;
t45 = -qJD(3) + t62;
t20 = -t45 * pkin(3) - t36 * pkin(9) + t27;
t28 = t53 * t33 + t69 * t42;
t34 = -t69 * qJD(2) + t53 * t63;
t22 = -t34 * pkin(9) + t28;
t52 = sin(qJ(4));
t68 = cos(qJ(4));
t10 = t68 * t20 - t52 * t22;
t26 = -t52 * t34 + t68 * t36;
t43 = -qJD(4) + t45;
t7 = -t43 * pkin(4) - t26 * pkin(10) + t10;
t11 = t52 * t20 + t68 * t22;
t24 = t68 * t34 + t52 * t36;
t9 = -t24 * pkin(10) + t11;
t4 = t51 * t7 + t67 * t9;
t14 = t67 * t24 + t51 * t26;
t16 = -t51 * t24 + t67 * t26;
t66 = t16 * t14;
t39 = -qJD(5) + t43;
t65 = t39 * t14;
t64 = t55 * t56;
t61 = t14 ^ 2 / 0.2e1;
t60 = qJD(1) * qJD(2);
t59 = t54 * t60;
t58 = t55 * t60;
t41 = -qJD(2) * pkin(2) + pkin(7) * t63;
t29 = t34 * pkin(3) + t41;
t3 = -t51 * t9 + t67 * t7;
t17 = t24 * pkin(4) + t29;
t50 = t55 ^ 2;
t49 = t54 ^ 2;
t38 = t39 ^ 2 / 0.2e1;
t13 = t16 ^ 2 / 0.2e1;
t12 = t16 * t39;
t5 = t14 * pkin(5) - t16 * qJ(6) + t17;
t2 = -t39 * qJ(6) + t4;
t1 = t39 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t70, 0, 0, 0, 0, t49 * t70, t54 * t64, t59, t50 * t70, t58, qJD(2) ^ 2 / 0.2e1, pkin(1) * t64 - pkin(7) * t59, -t56 * pkin(1) * t54 - pkin(7) * t58 (t49 + t50) * t56 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t56, t36 ^ 2 / 0.2e1, -t36 * t34, -t36 * t45, t34 ^ 2 / 0.2e1, t34 * t45, t45 ^ 2 / 0.2e1, -t27 * t45 + t41 * t34, t28 * t45 + t41 * t36, -t27 * t36 - t28 * t34, t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, -t26 * t43, t24 ^ 2 / 0.2e1, t24 * t43, t43 ^ 2 / 0.2e1, -t10 * t43 + t29 * t24, t11 * t43 + t29 * t26, -t10 * t26 - t11 * t24, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t13, -t66, -t12, t61, t65, t38, t17 * t14 - t3 * t39, t17 * t16 + t4 * t39, -t4 * t14 - t3 * t16, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t13, -t12, t66, t38, -t65, t61, t1 * t39 + t5 * t14, t1 * t16 - t2 * t14, -t5 * t16 - t2 * t39, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
