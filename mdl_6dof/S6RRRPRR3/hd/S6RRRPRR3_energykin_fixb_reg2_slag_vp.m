% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR3
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:21
% EndTime: 2019-03-09 18:13:21
% DurationCPUTime: 0.20s
% Computational Cost: add. (700->63), mult. (1545->143), div. (0->0), fcn. (1045->8), ass. (0->54)
t55 = qJD(1) ^ 2;
t71 = t55 / 0.2e1;
t70 = -pkin(8) - pkin(7);
t53 = sin(qJ(2));
t63 = qJD(1) * t53;
t36 = qJD(2) * pkin(2) + t70 * t63;
t54 = cos(qJ(2));
t62 = qJD(1) * t54;
t37 = t70 * t62;
t52 = sin(qJ(3));
t69 = cos(qJ(3));
t24 = t52 * t36 - t69 * t37;
t46 = qJD(2) + qJD(3);
t19 = t46 * qJ(4) + t24;
t31 = t52 * t63 - t69 * t62;
t12 = t31 * pkin(9) + t19;
t51 = sin(qJ(5));
t68 = cos(qJ(5));
t33 = (t52 * t54 + t69 * t53) * qJD(1);
t23 = t69 * t36 + t52 * t37;
t57 = qJD(4) - t23;
t9 = -t33 * pkin(9) + (-pkin(3) - pkin(4)) * t46 + t57;
t7 = t68 * t12 + t51 * t9;
t67 = cos(qJ(6));
t66 = t33 * t31;
t65 = t46 * t31;
t64 = t54 * t55;
t38 = -qJD(1) * pkin(1) - pkin(2) * t62;
t61 = t31 ^ 2 / 0.2e1;
t60 = qJD(1) * qJD(2);
t59 = t53 * t60;
t58 = t54 * t60;
t20 = -t68 * t31 + t51 * t33;
t16 = t31 * pkin(3) - t33 * qJ(4) + t38;
t10 = -t31 * pkin(4) - t16;
t6 = -t51 * t12 + t68 * t9;
t50 = sin(qJ(6));
t48 = t54 ^ 2;
t47 = t53 ^ 2;
t44 = -qJD(5) + t46;
t41 = t46 ^ 2 / 0.2e1;
t28 = t33 ^ 2 / 0.2e1;
t25 = t33 * t46;
t22 = t51 * t31 + t68 * t33;
t18 = qJD(6) + t20;
t17 = -t46 * pkin(3) + t57;
t15 = t67 * t22 - t50 * t44;
t13 = t50 * t22 + t67 * t44;
t5 = -t44 * pkin(10) + t7;
t4 = t44 * pkin(5) - t6;
t3 = t20 * pkin(5) - t22 * pkin(10) + t10;
t2 = t50 * t3 + t67 * t5;
t1 = t67 * t3 - t50 * t5;
t8 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t47 * t71, t53 * t64, t59, t48 * t71, t58, qJD(2) ^ 2 / 0.2e1, pkin(1) * t64 - pkin(7) * t59, -t55 * pkin(1) * t53 - pkin(7) * t58 (t47 + t48) * t55 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t48 / 0.2e1 + t47 / 0.2e1) * pkin(7) ^ 2) * t55, t28, -t66, t25, t61, -t65, t41, t23 * t46 + t38 * t31, -t24 * t46 + t38 * t33, -t23 * t33 - t24 * t31, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t28, t25, t66, t41, t65, t61, t16 * t31 - t17 * t46, t17 * t33 - t19 * t31, -t16 * t33 + t19 * t46, t19 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, -t22 * t44, t20 ^ 2 / 0.2e1, t20 * t44, t44 ^ 2 / 0.2e1, t10 * t20 - t6 * t44, t10 * t22 + t7 * t44, -t7 * t20 - t6 * t22, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t18, t13 ^ 2 / 0.2e1, -t13 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t4 * t13, t4 * t15 - t2 * t18, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t8;
