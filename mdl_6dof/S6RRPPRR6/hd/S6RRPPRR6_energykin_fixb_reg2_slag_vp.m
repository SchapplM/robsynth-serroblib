% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:14
% EndTime: 2019-03-09 09:15:14
% DurationCPUTime: 0.20s
% Computational Cost: add. (639->66), mult. (1439->148), div. (0->0), fcn. (916->8), ass. (0->51)
t57 = qJD(1) ^ 2;
t69 = t57 / 0.2e1;
t55 = sin(qJ(2));
t62 = qJ(4) * qJD(1);
t65 = qJD(1) * t55;
t63 = pkin(7) * t65 + qJD(3);
t25 = -t55 * t62 + (-pkin(2) - pkin(3)) * qJD(2) + t63;
t56 = cos(qJ(2));
t64 = qJD(1) * t56;
t34 = pkin(7) * t64 + qJD(2) * qJ(3);
t31 = -t56 * t62 + t34;
t50 = sin(pkin(10));
t51 = cos(pkin(10));
t15 = t51 * t25 - t50 * t31;
t30 = (-t50 * t56 + t51 * t55) * qJD(1);
t10 = -qJD(2) * pkin(4) - t30 * pkin(8) + t15;
t16 = t50 * t25 + t51 * t31;
t28 = (-t50 * t55 - t51 * t56) * qJD(1);
t11 = t28 * pkin(8) + t16;
t54 = sin(qJ(5));
t68 = cos(qJ(5));
t6 = t54 * t10 + t68 * t11;
t67 = cos(qJ(6));
t66 = t56 * t57;
t61 = qJD(1) * qJD(2);
t60 = t55 * t66;
t32 = -qJD(1) * pkin(1) - pkin(2) * t64 - qJ(3) * t65;
t38 = t55 * t61;
t59 = t56 * t61;
t19 = -t68 * t28 + t54 * t30;
t24 = pkin(3) * t64 + qJD(4) - t32;
t5 = t68 * t10 - t54 * t11;
t18 = -t28 * pkin(4) + t24;
t53 = sin(qJ(6));
t49 = t56 ^ 2;
t48 = t55 ^ 2;
t46 = qJD(2) ^ 2 / 0.2e1;
t44 = qJD(2) - qJD(5);
t37 = t49 * t69;
t36 = t48 * t69;
t33 = -qJD(2) * pkin(2) + t63;
t21 = t54 * t28 + t68 * t30;
t17 = qJD(6) + t19;
t14 = t67 * t21 - t53 * t44;
t12 = t53 * t21 + t67 * t44;
t7 = t19 * pkin(5) - t21 * pkin(9) + t18;
t4 = -t44 * pkin(9) + t6;
t3 = t44 * pkin(5) - t5;
t2 = t67 * t4 + t53 * t7;
t1 = -t53 * t4 + t67 * t7;
t8 = [0, 0, 0, 0, 0, t69, 0, 0, 0, 0, t36, t60, t38, t37, t59, t46, pkin(1) * t66 - pkin(7) * t38, -t57 * pkin(1) * t55 - pkin(7) * t59 (t48 + t49) * t57 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t49 / 0.2e1 + t48 / 0.2e1) * pkin(7) ^ 2) * t57, t36, t38, -t60, t46, -t59, t37, -t33 * qJD(2) - t32 * t64 (t33 * t55 + t34 * t56) * qJD(1), t34 * qJD(2) - t32 * t65, t34 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, t30 * t28, -t30 * qJD(2), t28 ^ 2 / 0.2e1, -t28 * qJD(2), t46, -t15 * qJD(2) - t24 * t28, t16 * qJD(2) + t24 * t30, -t15 * t30 + t16 * t28, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, -t21 * t44, t19 ^ 2 / 0.2e1, t19 * t44, t44 ^ 2 / 0.2e1, t18 * t19 - t5 * t44, t18 * t21 + t6 * t44, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t17, t12 ^ 2 / 0.2e1, -t12 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t12, t3 * t14 - t2 * t17, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
