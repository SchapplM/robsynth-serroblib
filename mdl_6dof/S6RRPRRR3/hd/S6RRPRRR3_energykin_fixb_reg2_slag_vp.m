% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:17
% EndTime: 2019-03-09 13:25:17
% DurationCPUTime: 0.24s
% Computational Cost: add. (1157->68), mult. (2743->160), div. (0->0), fcn. (2027->10), ass. (0->54)
t59 = qJD(1) ^ 2;
t71 = t59 / 0.2e1;
t52 = sin(pkin(11));
t53 = cos(pkin(11));
t58 = cos(qJ(2));
t64 = qJD(1) * t58;
t57 = sin(qJ(2));
t65 = qJD(1) * t57;
t39 = t52 * t65 - t53 * t64;
t41 = (t52 * t58 + t53 * t57) * qJD(1);
t46 = qJD(3) + (-pkin(2) * t58 - pkin(1)) * qJD(1);
t25 = t39 * pkin(3) - t41 * pkin(8) + t46;
t66 = pkin(7) + qJ(3);
t44 = qJD(2) * pkin(2) - t66 * t65;
t45 = t66 * t64;
t30 = t52 * t44 + t53 * t45;
t28 = qJD(2) * pkin(8) + t30;
t56 = sin(qJ(4));
t70 = cos(qJ(4));
t17 = t56 * t25 + t70 * t28;
t32 = -t70 * qJD(2) + t56 * t41;
t15 = -t32 * pkin(9) + t17;
t55 = sin(qJ(5));
t69 = cos(qJ(5));
t16 = t70 * t25 - t56 * t28;
t34 = t56 * qJD(2) + t70 * t41;
t38 = qJD(4) + t39;
t9 = t38 * pkin(4) - t34 * pkin(9) + t16;
t6 = t69 * t15 + t55 * t9;
t68 = cos(qJ(6));
t67 = t58 * t59;
t63 = qJD(1) * qJD(2);
t5 = -t55 * t15 + t69 * t9;
t62 = t57 * t63;
t61 = t58 * t63;
t29 = t53 * t44 - t52 * t45;
t27 = -qJD(2) * pkin(3) - t29;
t36 = qJD(5) + t38;
t18 = t32 * pkin(4) + t27;
t54 = sin(qJ(6));
t51 = t58 ^ 2;
t50 = t57 ^ 2;
t49 = qJD(2) ^ 2 / 0.2e1;
t35 = qJD(6) + t36;
t22 = -t55 * t32 + t69 * t34;
t20 = t69 * t32 + t55 * t34;
t14 = -t54 * t20 + t68 * t22;
t12 = t68 * t20 + t54 * t22;
t10 = t20 * pkin(5) + t18;
t4 = -t20 * pkin(10) + t6;
t3 = t36 * pkin(5) - t22 * pkin(10) + t5;
t2 = t54 * t3 + t68 * t4;
t1 = t68 * t3 - t54 * t4;
t7 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t50 * t71, t57 * t67, t62, t51 * t71, t61, t49, pkin(1) * t67 - pkin(7) * t62, -t59 * pkin(1) * t57 - pkin(7) * t61 (t50 + t51) * t59 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t51 / 0.2e1 + t50 / 0.2e1) * pkin(7) ^ 2) * t59, t41 ^ 2 / 0.2e1, -t41 * t39, t41 * qJD(2), t39 ^ 2 / 0.2e1, -t39 * qJD(2), t49, t29 * qJD(2) + t46 * t39, -t30 * qJD(2) + t46 * t41, -t29 * t41 - t30 * t39, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t38, t32 ^ 2 / 0.2e1, -t32 * t38, t38 ^ 2 / 0.2e1, t16 * t38 + t27 * t32, -t17 * t38 + t27 * t34, -t16 * t34 - t17 * t32, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t36, t20 ^ 2 / 0.2e1, -t20 * t36, t36 ^ 2 / 0.2e1, t18 * t20 + t5 * t36, t18 * t22 - t6 * t36, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t35, t12 ^ 2 / 0.2e1, -t12 * t35, t35 ^ 2 / 0.2e1, t1 * t35 + t10 * t12, t10 * t14 - t2 * t35, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t7;
