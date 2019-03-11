% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:09:29
% EndTime: 2019-03-09 16:09:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (658->67), mult. (1570->133), div. (0->0), fcn. (1133->8), ass. (0->52)
t70 = -pkin(4) - pkin(10);
t69 = cos(qJ(6));
t64 = cos(pkin(6)) * qJD(1);
t46 = qJD(2) + t64;
t51 = sin(qJ(3));
t53 = cos(qJ(3));
t52 = sin(qJ(2));
t48 = sin(pkin(6));
t65 = qJD(1) * t48;
t61 = t52 * t65;
t29 = -t53 * t46 + t51 * t61;
t54 = cos(qJ(2));
t60 = t54 * t65;
t39 = -qJD(3) + t60;
t68 = t29 * t39;
t31 = t51 * t46 + t53 * t61;
t67 = t31 * t29;
t17 = t31 * t39;
t55 = qJD(1) ^ 2;
t66 = t48 ^ 2 * t55;
t62 = pkin(1) * t64;
t34 = pkin(8) * t60 + t52 * t62;
t22 = t46 * pkin(9) + t34;
t23 = (-pkin(2) * t54 - pkin(9) * t52 - pkin(1)) * t65;
t13 = t53 * t22 + t51 * t23;
t33 = -pkin(8) * t61 + t54 * t62;
t24 = t29 ^ 2 / 0.2e1;
t25 = t31 ^ 2 / 0.2e1;
t36 = t39 ^ 2 / 0.2e1;
t63 = t54 * t66;
t11 = -t39 * qJ(4) + t13;
t21 = -t46 * pkin(2) - t33;
t59 = t66 / 0.2e1;
t12 = -t51 * t22 + t53 * t23;
t9 = t29 * pkin(3) - t31 * qJ(4) + t21;
t58 = qJD(4) - t12;
t8 = -t29 * qJ(5) - t11;
t57 = qJD(5) - t9;
t56 = -t31 * qJ(5) + t58;
t50 = sin(qJ(6));
t28 = qJD(6) + t31;
t16 = t69 * t29 + t50 * t39;
t14 = t50 * t29 - t69 * t39;
t10 = t39 * pkin(3) + t58;
t7 = -t29 * pkin(4) + t57;
t6 = -t39 * pkin(5) - t8;
t5 = (pkin(3) + pkin(4)) * t39 + t56;
t4 = (pkin(3) - t70) * t39 + t56;
t3 = t31 * pkin(5) + t70 * t29 + t57;
t2 = t50 * t3 + t69 * t4;
t1 = t69 * t3 - t50 * t4;
t15 = [0, 0, 0, 0, 0, t55 / 0.2e1, 0, 0, 0, 0, t52 ^ 2 * t59, t52 * t63, t46 * t61, t54 ^ 2 * t59, t46 * t60, t46 ^ 2 / 0.2e1, pkin(1) * t63 + t33 * t46, -pkin(1) * t52 * t66 - t34 * t46 (-t33 * t52 + t34 * t54) * t65, t34 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t59, t25, -t67, -t17, t24, t68, t36, -t12 * t39 + t21 * t29, t13 * t39 + t21 * t31, -t12 * t31 - t13 * t29, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t25, -t17, t67, t36, -t68, t24, t10 * t39 + t9 * t29, t10 * t31 - t11 * t29, -t11 * t39 - t9 * t31, t11 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t24, -t67, t68, t25, -t17, t36, t7 * t31 + t8 * t39, t7 * t29 - t5 * t39, -t8 * t29 - t5 * t31, t5 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t28, t14 ^ 2 / 0.2e1, -t14 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t6 * t14, t6 * t16 - t2 * t28, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg  = t15;
