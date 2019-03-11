% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR1
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
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:09:55
% EndTime: 2019-03-09 10:09:56
% DurationCPUTime: 0.24s
% Computational Cost: add. (1135->68), mult. (2781->157), div. (0->0), fcn. (2086->10), ass. (0->54)
t57 = qJD(1) ^ 2;
t71 = t57 / 0.2e1;
t55 = sin(qJ(2));
t64 = qJD(1) * t55;
t67 = pkin(7) + qJ(3);
t42 = qJD(2) * pkin(2) - t67 * t64;
t56 = cos(qJ(2));
t63 = qJD(1) * t56;
t43 = t67 * t63;
t52 = sin(pkin(10));
t66 = cos(pkin(10));
t32 = t66 * t42 - t52 * t43;
t40 = (t52 * t56 + t66 * t55) * qJD(1);
t22 = qJD(2) * pkin(3) - t40 * pkin(8) + t32;
t33 = t52 * t42 + t66 * t43;
t38 = t52 * t64 - t66 * t63;
t24 = -t38 * pkin(8) + t33;
t54 = sin(qJ(4));
t70 = cos(qJ(4));
t12 = t54 * t22 + t70 * t24;
t47 = qJD(2) + qJD(4);
t10 = t47 * qJ(5) + t12;
t29 = t70 * t38 + t54 * t40;
t31 = -t54 * t38 + t70 * t40;
t44 = qJD(3) + (-pkin(2) * t56 - pkin(1)) * qJD(1);
t34 = t38 * pkin(3) + t44;
t15 = t29 * pkin(4) - t31 * qJ(5) + t34;
t51 = sin(pkin(11));
t65 = cos(pkin(11));
t6 = t65 * t10 + t51 * t15;
t69 = cos(qJ(6));
t68 = t56 * t57;
t62 = t29 ^ 2 / 0.2e1;
t61 = qJD(1) * qJD(2);
t60 = t55 * t61;
t59 = t56 * t61;
t5 = -t51 * t10 + t65 * t15;
t11 = t70 * t22 - t54 * t24;
t9 = -t47 * pkin(4) + qJD(5) - t11;
t53 = sin(qJ(6));
t50 = t56 ^ 2;
t49 = t55 ^ 2;
t48 = qJD(2) ^ 2 / 0.2e1;
t28 = qJD(6) + t29;
t27 = t65 * t31 + t51 * t47;
t25 = t51 * t31 - t65 * t47;
t18 = -t53 * t25 + t69 * t27;
t16 = t69 * t25 + t53 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(9) + t6;
t3 = t29 * pkin(5) - t27 * pkin(9) + t5;
t2 = t53 * t3 + t69 * t4;
t1 = t69 * t3 - t53 * t4;
t8 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t49 * t71, t55 * t68, t60, t50 * t71, t59, t48, pkin(1) * t68 - pkin(7) * t60, -t57 * pkin(1) * t55 - pkin(7) * t59 (t49 + t50) * t57 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t57, t40 ^ 2 / 0.2e1, -t40 * t38, t40 * qJD(2), t38 ^ 2 / 0.2e1, -t38 * qJD(2), t48, t32 * qJD(2) + t44 * t38, -t33 * qJD(2) + t44 * t40, -t32 * t40 - t33 * t38, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t47, t62, -t29 * t47, t47 ^ 2 / 0.2e1, t11 * t47 + t34 * t29, -t12 * t47 + t34 * t31, -t11 * t31 - t12 * t29, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t29, t25 ^ 2 / 0.2e1, -t25 * t29, t62, t9 * t25 + t5 * t29, t9 * t27 - t6 * t29, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t28, t16 ^ 2 / 0.2e1, -t16 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t16, t7 * t18 - t2 * t28, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
