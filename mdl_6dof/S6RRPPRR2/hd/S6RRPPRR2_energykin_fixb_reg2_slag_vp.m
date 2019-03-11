% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:34
% EndTime: 2019-03-09 08:52:35
% DurationCPUTime: 0.22s
% Computational Cost: add. (1127->68), mult. (2743->157), div. (0->0), fcn. (2027->10), ass. (0->54)
t58 = qJD(1) ^ 2;
t71 = t58 / 0.2e1;
t52 = sin(pkin(10));
t53 = cos(pkin(10));
t57 = cos(qJ(2));
t64 = qJD(1) * t57;
t56 = sin(qJ(2));
t65 = qJD(1) * t56;
t38 = t52 * t65 - t53 * t64;
t40 = (t52 * t57 + t53 * t56) * qJD(1);
t45 = qJD(3) + (-pkin(2) * t57 - pkin(1)) * qJD(1);
t25 = t38 * pkin(3) - t40 * qJ(4) + t45;
t67 = pkin(7) + qJ(3);
t43 = qJD(2) * pkin(2) - t67 * t65;
t44 = t67 * t64;
t30 = t52 * t43 + t53 * t44;
t28 = qJD(2) * qJ(4) + t30;
t51 = sin(pkin(11));
t66 = cos(pkin(11));
t17 = t51 * t25 + t66 * t28;
t32 = -t66 * qJD(2) + t51 * t40;
t15 = -t32 * pkin(8) + t17;
t55 = sin(qJ(5));
t70 = cos(qJ(5));
t16 = t66 * t25 - t51 * t28;
t34 = t51 * qJD(2) + t66 * t40;
t9 = t38 * pkin(4) - t34 * pkin(8) + t16;
t6 = t70 * t15 + t55 * t9;
t69 = cos(qJ(6));
t68 = t57 * t58;
t63 = t38 ^ 2 / 0.2e1;
t62 = qJD(1) * qJD(2);
t5 = -t55 * t15 + t70 * t9;
t61 = t56 * t62;
t60 = t57 * t62;
t29 = t53 * t43 - t52 * t44;
t37 = qJD(5) + t38;
t27 = -qJD(2) * pkin(3) + qJD(4) - t29;
t18 = t32 * pkin(4) + t27;
t54 = sin(qJ(6));
t50 = t57 ^ 2;
t49 = t56 ^ 2;
t48 = qJD(2) ^ 2 / 0.2e1;
t35 = qJD(6) + t37;
t22 = -t55 * t32 + t70 * t34;
t20 = t70 * t32 + t55 * t34;
t13 = -t54 * t20 + t69 * t22;
t11 = t69 * t20 + t54 * t22;
t10 = t20 * pkin(5) + t18;
t4 = -t20 * pkin(9) + t6;
t3 = t37 * pkin(5) - t22 * pkin(9) + t5;
t2 = t54 * t3 + t69 * t4;
t1 = t69 * t3 - t54 * t4;
t7 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t49 * t71, t56 * t68, t61, t50 * t71, t60, t48, pkin(1) * t68 - pkin(7) * t61, -t58 * pkin(1) * t56 - pkin(7) * t60 (t49 + t50) * t58 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t58, t40 ^ 2 / 0.2e1, -t40 * t38, t40 * qJD(2), t63, -t38 * qJD(2), t48, t29 * qJD(2) + t45 * t38, -t30 * qJD(2) + t45 * t40, -t29 * t40 - t30 * t38, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t38, t32 ^ 2 / 0.2e1, -t32 * t38, t63, t16 * t38 + t27 * t32, -t17 * t38 + t27 * t34, -t16 * t34 - t17 * t32, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t37, t20 ^ 2 / 0.2e1, -t20 * t37, t37 ^ 2 / 0.2e1, t18 * t20 + t5 * t37, t18 * t22 - t6 * t37, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t35, t11 ^ 2 / 0.2e1, -t11 * t35, t35 ^ 2 / 0.2e1, t1 * t35 + t10 * t11, t10 * t13 - t2 * t35, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t7;
