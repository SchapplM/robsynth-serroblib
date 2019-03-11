% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:05
% EndTime: 2019-03-09 15:23:05
% DurationCPUTime: 0.22s
% Computational Cost: add. (1167->68), mult. (2781->157), div. (0->0), fcn. (2086->10), ass. (0->54)
t57 = qJD(1) ^ 2;
t71 = t57 / 0.2e1;
t70 = -pkin(8) - pkin(7);
t55 = sin(qJ(2));
t64 = qJD(1) * t55;
t42 = qJD(2) * pkin(2) + t70 * t64;
t56 = cos(qJ(2));
t63 = qJD(1) * t56;
t43 = t70 * t63;
t54 = sin(qJ(3));
t69 = cos(qJ(3));
t32 = t69 * t42 + t54 * t43;
t40 = (t54 * t56 + t69 * t55) * qJD(1);
t48 = qJD(2) + qJD(3);
t20 = t48 * pkin(3) - t40 * qJ(4) + t32;
t33 = t54 * t42 - t69 * t43;
t38 = t54 * t64 - t69 * t63;
t24 = -t38 * qJ(4) + t33;
t52 = sin(pkin(10));
t66 = cos(pkin(10));
t12 = t52 * t20 + t66 * t24;
t10 = t48 * qJ(5) + t12;
t29 = t66 * t38 + t52 * t40;
t31 = -t52 * t38 + t66 * t40;
t44 = (-pkin(2) * t56 - pkin(1)) * qJD(1);
t34 = t38 * pkin(3) + qJD(4) + t44;
t15 = t29 * pkin(4) - t31 * qJ(5) + t34;
t51 = sin(pkin(11));
t65 = cos(pkin(11));
t6 = t65 * t10 + t51 * t15;
t68 = cos(qJ(6));
t67 = t56 * t57;
t62 = t29 ^ 2 / 0.2e1;
t61 = qJD(1) * qJD(2);
t60 = t55 * t61;
t59 = t56 * t61;
t5 = -t51 * t10 + t65 * t15;
t11 = t66 * t20 - t52 * t24;
t9 = -t48 * pkin(4) + qJD(5) - t11;
t53 = sin(qJ(6));
t50 = t56 ^ 2;
t49 = t55 ^ 2;
t47 = t48 ^ 2 / 0.2e1;
t28 = qJD(6) + t29;
t27 = t65 * t31 + t51 * t48;
t25 = t51 * t31 - t65 * t48;
t18 = -t53 * t25 + t68 * t27;
t16 = t68 * t25 + t53 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(9) + t6;
t3 = t29 * pkin(5) - t27 * pkin(9) + t5;
t2 = t53 * t3 + t68 * t4;
t1 = t68 * t3 - t53 * t4;
t8 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t49 * t71, t55 * t67, t60, t50 * t71, t59, qJD(2) ^ 2 / 0.2e1, pkin(1) * t67 - pkin(7) * t60, -t57 * pkin(1) * t55 - pkin(7) * t59 (t49 + t50) * t57 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t57, t40 ^ 2 / 0.2e1, -t40 * t38, t40 * t48, t38 ^ 2 / 0.2e1, -t38 * t48, t47, t32 * t48 + t44 * t38, -t33 * t48 + t44 * t40, -t32 * t40 - t33 * t38, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t48, t62, -t29 * t48, t47, t11 * t48 + t34 * t29, -t12 * t48 + t34 * t31, -t11 * t31 - t12 * t29, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t29, t25 ^ 2 / 0.2e1, -t25 * t29, t62, t9 * t25 + t5 * t29, t9 * t27 - t6 * t29, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t28, t16 ^ 2 / 0.2e1, -t16 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t16, t7 * t18 - t2 * t28, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
