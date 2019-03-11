% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:24
% EndTime: 2019-03-09 13:53:24
% DurationCPUTime: 0.18s
% Computational Cost: add. (669->66), mult. (1439->149), div. (0->0), fcn. (916->8), ass. (0->51)
t59 = qJD(1) ^ 2;
t70 = t59 / 0.2e1;
t56 = sin(qJ(2));
t66 = qJD(1) * t56;
t64 = pkin(7) * t66 + qJD(3);
t25 = -pkin(8) * t66 + (-pkin(2) - pkin(3)) * qJD(2) + t64;
t58 = cos(qJ(2));
t65 = qJD(1) * t58;
t34 = pkin(7) * t65 + qJD(2) * qJ(3);
t31 = -pkin(8) * t65 + t34;
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t16 = t55 * t25 + t57 * t31;
t28 = (-t55 * t56 - t57 * t58) * qJD(1);
t11 = t28 * pkin(9) + t16;
t54 = sin(qJ(5));
t69 = cos(qJ(5));
t15 = t57 * t25 - t55 * t31;
t30 = (-t55 * t58 + t56 * t57) * qJD(1);
t46 = qJD(2) - qJD(4);
t9 = -t46 * pkin(4) - t30 * pkin(9) + t15;
t6 = t69 * t11 + t54 * t9;
t68 = cos(qJ(6));
t67 = t58 * t59;
t63 = qJD(1) * qJD(2);
t62 = t56 * t67;
t32 = -qJD(1) * pkin(1) - pkin(2) * t65 - qJ(3) * t66;
t38 = t56 * t63;
t61 = t58 * t63;
t18 = -t69 * t28 + t54 * t30;
t24 = pkin(3) * t65 - t32;
t21 = -t28 * pkin(4) + t24;
t5 = -t54 * t11 + t69 * t9;
t53 = sin(qJ(6));
t51 = t58 ^ 2;
t50 = t56 ^ 2;
t48 = qJD(2) ^ 2 / 0.2e1;
t44 = -qJD(5) + t46;
t37 = t51 * t70;
t36 = t50 * t70;
t33 = -qJD(2) * pkin(2) + t64;
t20 = t54 * t28 + t69 * t30;
t17 = qJD(6) + t18;
t14 = t68 * t20 - t53 * t44;
t12 = t53 * t20 + t68 * t44;
t7 = t18 * pkin(5) - t20 * pkin(10) + t21;
t4 = -t44 * pkin(10) + t6;
t3 = t44 * pkin(5) - t5;
t2 = t68 * t4 + t53 * t7;
t1 = -t53 * t4 + t68 * t7;
t8 = [0, 0, 0, 0, 0, t70, 0, 0, 0, 0, t36, t62, t38, t37, t61, t48, pkin(1) * t67 - pkin(7) * t38, -t59 * pkin(1) * t56 - pkin(7) * t61 (t50 + t51) * t59 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t51 / 0.2e1 + t50 / 0.2e1) * pkin(7) ^ 2) * t59, t36, t38, -t62, t48, -t61, t37, -t33 * qJD(2) - t32 * t65 (t33 * t56 + t34 * t58) * qJD(1), t34 * qJD(2) - t32 * t66, t34 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, t30 * t28, -t30 * t46, t28 ^ 2 / 0.2e1, -t28 * t46, t46 ^ 2 / 0.2e1, -t15 * t46 - t24 * t28, t16 * t46 + t24 * t30, -t15 * t30 + t16 * t28, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, -t20 * t44, t18 ^ 2 / 0.2e1, t18 * t44, t44 ^ 2 / 0.2e1, t21 * t18 - t5 * t44, t21 * t20 + t6 * t44, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t17, t12 ^ 2 / 0.2e1, -t12 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t12, t3 * t14 - t2 * t17, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
