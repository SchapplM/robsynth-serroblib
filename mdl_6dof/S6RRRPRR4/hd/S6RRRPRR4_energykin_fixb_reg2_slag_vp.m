% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:17:44
% EndTime: 2019-03-09 18:17:45
% DurationCPUTime: 0.22s
% Computational Cost: add. (1280->68), mult. (2743->160), div. (0->0), fcn. (2027->10), ass. (0->54)
t58 = qJD(1) ^ 2;
t71 = t58 / 0.2e1;
t70 = -pkin(8) - pkin(7);
t54 = sin(qJ(3));
t56 = cos(qJ(3));
t57 = cos(qJ(2));
t64 = qJD(1) * t57;
t55 = sin(qJ(2));
t65 = qJD(1) * t55;
t38 = t54 * t65 - t56 * t64;
t40 = (t54 * t57 + t55 * t56) * qJD(1);
t45 = (-pkin(2) * t57 - pkin(1)) * qJD(1);
t25 = t38 * pkin(3) - t40 * qJ(4) + t45;
t43 = qJD(2) * pkin(2) + t70 * t65;
t44 = t70 * t64;
t30 = t54 * t43 - t56 * t44;
t48 = qJD(2) + qJD(3);
t28 = t48 * qJ(4) + t30;
t51 = sin(pkin(11));
t66 = cos(pkin(11));
t17 = t51 * t25 + t66 * t28;
t32 = t51 * t40 - t66 * t48;
t15 = -t32 * pkin(9) + t17;
t53 = sin(qJ(5));
t69 = cos(qJ(5));
t16 = t66 * t25 - t51 * t28;
t34 = t66 * t40 + t51 * t48;
t9 = t38 * pkin(4) - t34 * pkin(9) + t16;
t6 = t69 * t15 + t53 * t9;
t68 = cos(qJ(6));
t67 = t57 * t58;
t63 = t38 ^ 2 / 0.2e1;
t62 = qJD(1) * qJD(2);
t5 = -t53 * t15 + t69 * t9;
t61 = t55 * t62;
t60 = t57 * t62;
t29 = t56 * t43 + t54 * t44;
t37 = qJD(5) + t38;
t27 = -t48 * pkin(3) + qJD(4) - t29;
t18 = t32 * pkin(4) + t27;
t52 = sin(qJ(6));
t50 = t57 ^ 2;
t49 = t55 ^ 2;
t35 = qJD(6) + t37;
t22 = -t53 * t32 + t69 * t34;
t20 = t69 * t32 + t53 * t34;
t13 = -t52 * t20 + t68 * t22;
t11 = t68 * t20 + t52 * t22;
t10 = t20 * pkin(5) + t18;
t4 = -t20 * pkin(10) + t6;
t3 = t37 * pkin(5) - t22 * pkin(10) + t5;
t2 = t52 * t3 + t68 * t4;
t1 = t68 * t3 - t52 * t4;
t7 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t49 * t71, t55 * t67, t61, t50 * t71, t60, qJD(2) ^ 2 / 0.2e1, pkin(1) * t67 - pkin(7) * t61, -t58 * pkin(1) * t55 - pkin(7) * t60 (t49 + t50) * t58 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t58, t40 ^ 2 / 0.2e1, -t40 * t38, t40 * t48, t63, -t38 * t48, t48 ^ 2 / 0.2e1, t29 * t48 + t45 * t38, -t30 * t48 + t45 * t40, -t29 * t40 - t30 * t38, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t38, t32 ^ 2 / 0.2e1, -t32 * t38, t63, t16 * t38 + t27 * t32, -t17 * t38 + t27 * t34, -t16 * t34 - t17 * t32, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t37, t20 ^ 2 / 0.2e1, -t20 * t37, t37 ^ 2 / 0.2e1, t18 * t20 + t5 * t37, t18 * t22 - t6 * t37, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t35, t11 ^ 2 / 0.2e1, -t11 * t35, t35 ^ 2 / 0.2e1, t1 * t35 + t10 * t11, t10 * t13 - t2 * t35, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t7;
