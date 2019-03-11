% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:18
% EndTime: 2019-03-09 22:09:18
% DurationCPUTime: 0.24s
% Computational Cost: add. (1300->68), mult. (2743->160), div. (0->0), fcn. (2027->10), ass. (0->54)
t59 = qJD(1) ^ 2;
t71 = t59 / 0.2e1;
t70 = -pkin(8) - pkin(7);
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t58 = cos(qJ(2));
t64 = qJD(1) * t58;
t56 = sin(qJ(2));
t65 = qJD(1) * t56;
t39 = t55 * t65 - t57 * t64;
t41 = (t55 * t58 + t56 * t57) * qJD(1);
t46 = (-pkin(2) * t58 - pkin(1)) * qJD(1);
t25 = t39 * pkin(3) - t41 * pkin(9) + t46;
t44 = qJD(2) * pkin(2) + t70 * t65;
t45 = t70 * t64;
t30 = t55 * t44 - t57 * t45;
t49 = qJD(2) + qJD(3);
t28 = t49 * pkin(9) + t30;
t54 = sin(qJ(4));
t69 = cos(qJ(4));
t17 = t54 * t25 + t69 * t28;
t32 = t54 * t41 - t69 * t49;
t15 = -t32 * qJ(5) + t17;
t52 = sin(pkin(11));
t66 = cos(pkin(11));
t16 = t69 * t25 - t54 * t28;
t34 = t69 * t41 + t54 * t49;
t38 = qJD(4) + t39;
t9 = t38 * pkin(4) - t34 * qJ(5) + t16;
t6 = t66 * t15 + t52 * t9;
t68 = cos(qJ(6));
t67 = t58 * t59;
t63 = qJD(1) * qJD(2);
t5 = -t52 * t15 + t66 * t9;
t62 = t56 * t63;
t61 = t58 * t63;
t29 = t57 * t44 + t55 * t45;
t27 = -t49 * pkin(3) - t29;
t18 = t32 * pkin(4) + qJD(5) + t27;
t53 = sin(qJ(6));
t51 = t58 ^ 2;
t50 = t56 ^ 2;
t36 = qJD(6) + t38;
t35 = t38 ^ 2 / 0.2e1;
t22 = -t52 * t32 + t66 * t34;
t20 = t66 * t32 + t52 * t34;
t13 = -t53 * t20 + t68 * t22;
t11 = t68 * t20 + t53 * t22;
t10 = t20 * pkin(5) + t18;
t4 = -t20 * pkin(10) + t6;
t3 = t38 * pkin(5) - t22 * pkin(10) + t5;
t2 = t53 * t3 + t68 * t4;
t1 = t68 * t3 - t53 * t4;
t7 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t50 * t71, t56 * t67, t62, t51 * t71, t61, qJD(2) ^ 2 / 0.2e1, pkin(1) * t67 - pkin(7) * t62, -t59 * pkin(1) * t56 - pkin(7) * t61 (t50 + t51) * t59 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t51 / 0.2e1 + t50 / 0.2e1) * pkin(7) ^ 2) * t59, t41 ^ 2 / 0.2e1, -t41 * t39, t41 * t49, t39 ^ 2 / 0.2e1, -t39 * t49, t49 ^ 2 / 0.2e1, t29 * t49 + t46 * t39, -t30 * t49 + t46 * t41, -t29 * t41 - t30 * t39, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t38, t32 ^ 2 / 0.2e1, -t32 * t38, t35, t16 * t38 + t27 * t32, -t17 * t38 + t27 * t34, -t16 * t34 - t17 * t32, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t38, t20 ^ 2 / 0.2e1, -t20 * t38, t35, t18 * t20 + t5 * t38, t18 * t22 - t6 * t38, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t36, t11 ^ 2 / 0.2e1, -t11 * t36, t36 ^ 2 / 0.2e1, t1 * t36 + t10 * t11, t10 * t13 - t2 * t36, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t7;
