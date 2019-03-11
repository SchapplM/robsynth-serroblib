% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRR2
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:08
% EndTime: 2019-03-09 13:19:08
% DurationCPUTime: 0.25s
% Computational Cost: add. (1145->68), mult. (2781->160), div. (0->0), fcn. (2086->10), ass. (0->54)
t59 = qJD(1) ^ 2;
t71 = t59 / 0.2e1;
t56 = sin(qJ(2));
t65 = qJD(1) * t56;
t67 = pkin(7) + qJ(3);
t43 = qJD(2) * pkin(2) - t67 * t65;
t58 = cos(qJ(2));
t64 = qJD(1) * t58;
t44 = t67 * t64;
t52 = sin(pkin(11));
t66 = cos(pkin(11));
t33 = t66 * t43 - t52 * t44;
t41 = (t52 * t58 + t66 * t56) * qJD(1);
t22 = qJD(2) * pkin(3) - t41 * pkin(8) + t33;
t34 = t52 * t43 + t66 * t44;
t39 = t52 * t65 - t66 * t64;
t24 = -t39 * pkin(8) + t34;
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t12 = t55 * t22 + t57 * t24;
t48 = qJD(2) + qJD(4);
t10 = t48 * pkin(9) + t12;
t30 = t57 * t39 + t55 * t41;
t32 = -t55 * t39 + t57 * t41;
t45 = qJD(3) + (-pkin(2) * t58 - pkin(1)) * qJD(1);
t35 = t39 * pkin(3) + t45;
t15 = t30 * pkin(4) - t32 * pkin(9) + t35;
t54 = sin(qJ(5));
t70 = cos(qJ(5));
t6 = t70 * t10 + t54 * t15;
t69 = cos(qJ(6));
t68 = t58 * t59;
t63 = qJD(1) * qJD(2);
t62 = t56 * t63;
t61 = t58 * t63;
t5 = -t54 * t10 + t70 * t15;
t11 = t57 * t22 - t55 * t24;
t29 = qJD(5) + t30;
t9 = -t48 * pkin(4) - t11;
t53 = sin(qJ(6));
t51 = t58 ^ 2;
t50 = t56 ^ 2;
t49 = qJD(2) ^ 2 / 0.2e1;
t28 = qJD(6) + t29;
t27 = t70 * t32 + t54 * t48;
t25 = t54 * t32 - t70 * t48;
t18 = -t53 * t25 + t69 * t27;
t16 = t69 * t25 + t53 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(10) + t6;
t3 = t29 * pkin(5) - t27 * pkin(10) + t5;
t2 = t53 * t3 + t69 * t4;
t1 = t69 * t3 - t53 * t4;
t8 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t50 * t71, t56 * t68, t62, t51 * t71, t61, t49, pkin(1) * t68 - pkin(7) * t62, -t59 * pkin(1) * t56 - pkin(7) * t61 (t50 + t51) * t59 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t51 / 0.2e1 + t50 / 0.2e1) * pkin(7) ^ 2) * t59, t41 ^ 2 / 0.2e1, -t41 * t39, t41 * qJD(2), t39 ^ 2 / 0.2e1, -t39 * qJD(2), t49, t33 * qJD(2) + t45 * t39, -t34 * qJD(2) + t45 * t41, -t33 * t41 - t34 * t39, t34 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t30, t32 * t48, t30 ^ 2 / 0.2e1, -t30 * t48, t48 ^ 2 / 0.2e1, t11 * t48 + t35 * t30, -t12 * t48 + t35 * t32, -t11 * t32 - t12 * t30, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t29, t25 ^ 2 / 0.2e1, -t25 * t29, t29 ^ 2 / 0.2e1, t9 * t25 + t5 * t29, t9 * t27 - t6 * t29, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t28, t16 ^ 2 / 0.2e1, -t16 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t16, t7 * t18 - t2 * t28, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
