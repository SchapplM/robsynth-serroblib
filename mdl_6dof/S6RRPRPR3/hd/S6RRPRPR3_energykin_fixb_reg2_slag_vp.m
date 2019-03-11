% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:06
% EndTime: 2019-03-09 10:19:06
% DurationCPUTime: 0.24s
% Computational Cost: add. (1147->68), mult. (2743->157), div. (0->0), fcn. (2027->10), ass. (0->54)
t59 = qJD(1) ^ 2;
t71 = t59 / 0.2e1;
t53 = sin(pkin(10));
t54 = cos(pkin(10));
t58 = cos(qJ(2));
t64 = qJD(1) * t58;
t57 = sin(qJ(2));
t65 = qJD(1) * t57;
t39 = t53 * t65 - t54 * t64;
t41 = (t53 * t58 + t54 * t57) * qJD(1);
t46 = qJD(3) + (-pkin(2) * t58 - pkin(1)) * qJD(1);
t25 = t39 * pkin(3) - t41 * pkin(8) + t46;
t67 = pkin(7) + qJ(3);
t44 = qJD(2) * pkin(2) - t67 * t65;
t45 = t67 * t64;
t30 = t53 * t44 + t54 * t45;
t28 = qJD(2) * pkin(8) + t30;
t56 = sin(qJ(4));
t70 = cos(qJ(4));
t17 = t56 * t25 + t70 * t28;
t32 = -t70 * qJD(2) + t56 * t41;
t15 = -t32 * qJ(5) + t17;
t52 = sin(pkin(11));
t66 = cos(pkin(11));
t16 = t70 * t25 - t56 * t28;
t34 = t56 * qJD(2) + t70 * t41;
t38 = qJD(4) + t39;
t9 = t38 * pkin(4) - t34 * qJ(5) + t16;
t6 = t66 * t15 + t52 * t9;
t69 = cos(qJ(6));
t68 = t58 * t59;
t63 = qJD(1) * qJD(2);
t5 = -t52 * t15 + t66 * t9;
t62 = t57 * t63;
t61 = t58 * t63;
t29 = t54 * t44 - t53 * t45;
t27 = -qJD(2) * pkin(3) - t29;
t18 = t32 * pkin(4) + qJD(5) + t27;
t55 = sin(qJ(6));
t51 = t58 ^ 2;
t50 = t57 ^ 2;
t49 = qJD(2) ^ 2 / 0.2e1;
t36 = qJD(6) + t38;
t35 = t38 ^ 2 / 0.2e1;
t22 = -t52 * t32 + t66 * t34;
t20 = t66 * t32 + t52 * t34;
t13 = -t55 * t20 + t69 * t22;
t11 = t69 * t20 + t55 * t22;
t10 = t20 * pkin(5) + t18;
t4 = -t20 * pkin(9) + t6;
t3 = t38 * pkin(5) - t22 * pkin(9) + t5;
t2 = t55 * t3 + t69 * t4;
t1 = t69 * t3 - t55 * t4;
t7 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t50 * t71, t57 * t68, t62, t51 * t71, t61, t49, pkin(1) * t68 - pkin(7) * t62, -t59 * pkin(1) * t57 - pkin(7) * t61 (t50 + t51) * t59 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t51 / 0.2e1 + t50 / 0.2e1) * pkin(7) ^ 2) * t59, t41 ^ 2 / 0.2e1, -t41 * t39, t41 * qJD(2), t39 ^ 2 / 0.2e1, -t39 * qJD(2), t49, t29 * qJD(2) + t46 * t39, -t30 * qJD(2) + t46 * t41, -t29 * t41 - t30 * t39, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t38, t32 ^ 2 / 0.2e1, -t32 * t38, t35, t16 * t38 + t27 * t32, -t17 * t38 + t27 * t34, -t16 * t34 - t17 * t32, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t38, t20 ^ 2 / 0.2e1, -t20 * t38, t35, t18 * t20 + t5 * t38, t18 * t22 - t6 * t38, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t36, t11 ^ 2 / 0.2e1, -t11 * t36, t36 ^ 2 / 0.2e1, t1 * t36 + t10 * t11, t10 * t13 - t2 * t36, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t7;
