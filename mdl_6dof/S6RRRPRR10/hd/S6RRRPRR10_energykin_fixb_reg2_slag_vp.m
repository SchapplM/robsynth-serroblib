% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:21:55
% EndTime: 2019-03-09 19:21:55
% DurationCPUTime: 0.19s
% Computational Cost: add. (747->66), mult. (1550->144), div. (0->0), fcn. (1020->8), ass. (0->53)
t56 = qJD(1) ^ 2;
t71 = t56 / 0.2e1;
t70 = cos(qJ(3));
t69 = cos(qJ(5));
t68 = cos(qJ(6));
t53 = sin(qJ(3));
t54 = sin(qJ(2));
t64 = qJD(1) * t54;
t33 = -t70 * qJD(2) + t53 * t64;
t35 = t53 * qJD(2) + t70 * t64;
t67 = t35 * t33;
t55 = cos(qJ(2));
t63 = t55 * qJD(1);
t44 = -qJD(3) + t63;
t66 = t44 * t33;
t65 = t55 * t56;
t29 = (-pkin(2) * t55 - pkin(8) * t54 - pkin(1)) * qJD(1);
t40 = pkin(7) * t63 + qJD(2) * pkin(8);
t24 = t70 * t29 - t53 * t40;
t58 = qJD(4) - t24;
t13 = -t35 * pkin(9) + (pkin(3) + pkin(4)) * t44 + t58;
t25 = t53 * t29 + t70 * t40;
t18 = -t44 * qJ(4) + t25;
t15 = t33 * pkin(9) + t18;
t52 = sin(qJ(5));
t6 = t52 * t13 + t69 * t15;
t39 = -qJD(2) * pkin(2) + pkin(7) * t64;
t62 = t33 ^ 2 / 0.2e1;
t61 = qJD(1) * qJD(2);
t60 = t54 * t61;
t59 = t55 * t61;
t5 = t69 * t13 - t52 * t15;
t19 = t33 * pkin(3) - t35 * qJ(4) + t39;
t43 = qJD(5) + t44;
t16 = -t33 * pkin(4) - t19;
t51 = sin(qJ(6));
t49 = t55 ^ 2;
t48 = t54 ^ 2;
t41 = t44 ^ 2 / 0.2e1;
t38 = qJD(6) + t43;
t30 = t35 ^ 2 / 0.2e1;
t26 = t35 * t44;
t23 = t52 * t33 + t69 * t35;
t21 = -t69 * t33 + t52 * t35;
t17 = t44 * pkin(3) + t58;
t10 = -t51 * t21 + t68 * t23;
t8 = t68 * t21 + t51 * t23;
t7 = t21 * pkin(5) + t16;
t4 = -t21 * pkin(10) + t6;
t3 = t43 * pkin(5) - t23 * pkin(10) + t5;
t2 = t51 * t3 + t68 * t4;
t1 = t68 * t3 - t51 * t4;
t9 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t48 * t71, t54 * t65, t60, t49 * t71, t59, qJD(2) ^ 2 / 0.2e1, pkin(1) * t65 - pkin(7) * t60, -t56 * pkin(1) * t54 - pkin(7) * t59 (t48 + t49) * t56 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t49 / 0.2e1 + t48 / 0.2e1) * pkin(7) ^ 2) * t56, t30, -t67, -t26, t62, t66, t41, -t24 * t44 + t39 * t33, t25 * t44 + t39 * t35, -t24 * t35 - t25 * t33, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, t30, -t26, t67, t41, -t66, t62, t17 * t44 + t19 * t33, t17 * t35 - t18 * t33, -t18 * t44 - t19 * t35, t18 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * t43, t21 ^ 2 / 0.2e1, -t21 * t43, t43 ^ 2 / 0.2e1, t16 * t21 + t5 * t43, t16 * t23 - t6 * t43, -t6 * t21 - t5 * t23, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t38, t8 ^ 2 / 0.2e1, -t8 * t38, t38 ^ 2 / 0.2e1, t1 * t38 + t7 * t8, t7 * t10 - t2 * t38, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t9;
