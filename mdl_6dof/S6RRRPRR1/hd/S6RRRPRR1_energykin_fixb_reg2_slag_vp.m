% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR1
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
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:04:05
% EndTime: 2019-03-09 18:04:05
% DurationCPUTime: 0.23s
% Computational Cost: add. (1240->68), mult. (3051->160), div. (0->0), fcn. (2311->10), ass. (0->54)
t57 = qJD(1) ^ 2;
t70 = t57 / 0.2e1;
t69 = -pkin(8) - pkin(7);
t55 = sin(qJ(2));
t63 = qJD(1) * t55;
t41 = qJD(2) * pkin(2) + t63 * t69;
t56 = cos(qJ(2));
t62 = qJD(1) * t56;
t42 = t69 * t62;
t54 = sin(qJ(3));
t68 = cos(qJ(3));
t31 = t41 * t68 + t42 * t54;
t39 = (t54 * t56 + t55 * t68) * qJD(1);
t48 = qJD(2) + qJD(3);
t24 = pkin(3) * t48 - qJ(4) * t39 + t31;
t32 = t41 * t54 - t42 * t68;
t37 = t54 * t63 - t62 * t68;
t26 = -qJ(4) * t37 + t32;
t51 = sin(pkin(11));
t64 = cos(pkin(11));
t13 = t24 * t51 + t26 * t64;
t28 = t37 * t64 + t39 * t51;
t11 = -pkin(9) * t28 + t13;
t53 = sin(qJ(5));
t67 = cos(qJ(5));
t12 = t24 * t64 - t26 * t51;
t30 = -t37 * t51 + t39 * t64;
t9 = pkin(4) * t48 - pkin(9) * t30 + t12;
t6 = t11 * t67 + t53 * t9;
t66 = cos(qJ(6));
t65 = t56 * t57;
t61 = qJD(1) * qJD(2);
t60 = t55 * t61;
t59 = t56 * t61;
t18 = t28 * t67 + t30 * t53;
t43 = (-pkin(2) * t56 - pkin(1)) * qJD(1);
t5 = -t11 * t53 + t67 * t9;
t33 = t37 * pkin(3) + qJD(4) + t43;
t21 = t28 * pkin(4) + t33;
t52 = sin(qJ(6));
t50 = t56 ^ 2;
t49 = t55 ^ 2;
t47 = qJD(5) + t48;
t46 = t48 ^ 2 / 0.2e1;
t20 = -t28 * t53 + t30 * t67;
t17 = qJD(6) + t18;
t16 = t20 * t66 + t47 * t52;
t14 = t20 * t52 - t47 * t66;
t7 = t18 * pkin(5) - t20 * pkin(10) + t21;
t4 = pkin(10) * t47 + t6;
t3 = -pkin(5) * t47 - t5;
t2 = t4 * t66 + t52 * t7;
t1 = -t4 * t52 + t66 * t7;
t8 = [0, 0, 0, 0, 0, t70, 0, 0, 0, 0, t49 * t70, t55 * t65, t60, t50 * t70, t59, qJD(2) ^ 2 / 0.2e1, pkin(1) * t65 - pkin(7) * t60, -pkin(1) * t55 * t57 - pkin(7) * t59 (t49 + t50) * t57 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t57, t39 ^ 2 / 0.2e1, -t39 * t37, t39 * t48, t37 ^ 2 / 0.2e1, -t37 * t48, t46, t31 * t48 + t37 * t43, -t32 * t48 + t39 * t43, -t31 * t39 - t32 * t37, t32 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * t48, t28 ^ 2 / 0.2e1, -t28 * t48, t46, t12 * t48 + t28 * t33, -t13 * t48 + t30 * t33, -t12 * t30 - t13 * t28, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t47, t18 ^ 2 / 0.2e1, -t18 * t47, t47 ^ 2 / 0.2e1, t18 * t21 + t47 * t5, t20 * t21 - t47 * t6, -t18 * t6 - t20 * t5, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t17, t14 ^ 2 / 0.2e1, -t14 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t14 * t3, t16 * t3 - t17 * t2, -t1 * t16 - t14 * t2, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
