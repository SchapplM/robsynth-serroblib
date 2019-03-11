% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:47:42
% EndTime: 2019-03-09 10:47:42
% DurationCPUTime: 0.18s
% Computational Cost: add. (654->66), mult. (1439->146), div. (0->0), fcn. (916->8), ass. (0->51)
t58 = qJD(1) ^ 2;
t69 = t58 / 0.2e1;
t55 = sin(qJ(2));
t65 = qJD(1) * t55;
t63 = pkin(7) * t65 + qJD(3);
t25 = -pkin(8) * t65 + (-pkin(2) - pkin(3)) * qJD(2) + t63;
t57 = cos(qJ(2));
t64 = qJD(1) * t57;
t34 = pkin(7) * t64 + qJD(2) * qJ(3);
t31 = -pkin(8) * t64 + t34;
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t16 = t54 * t25 + t56 * t31;
t28 = (-t54 * t55 - t56 * t57) * qJD(1);
t11 = t28 * qJ(5) + t16;
t51 = sin(pkin(10));
t66 = cos(pkin(10));
t15 = t56 * t25 - t54 * t31;
t30 = (-t54 * t57 + t55 * t56) * qJD(1);
t45 = qJD(2) - qJD(4);
t9 = -t45 * pkin(4) - t30 * qJ(5) + t15;
t6 = t66 * t11 + t51 * t9;
t68 = cos(qJ(6));
t67 = t57 * t58;
t62 = qJD(1) * qJD(2);
t61 = t55 * t67;
t32 = -qJD(1) * pkin(1) - pkin(2) * t64 - qJ(3) * t65;
t38 = t55 * t62;
t60 = t57 * t62;
t19 = -t66 * t28 + t51 * t30;
t24 = pkin(3) * t64 - t32;
t5 = -t51 * t11 + t66 * t9;
t18 = -t28 * pkin(4) + qJD(5) + t24;
t53 = sin(qJ(6));
t50 = t57 ^ 2;
t49 = t55 ^ 2;
t47 = qJD(2) ^ 2 / 0.2e1;
t40 = t45 ^ 2 / 0.2e1;
t37 = t50 * t69;
t36 = t49 * t69;
t33 = -qJD(2) * pkin(2) + t63;
t21 = t51 * t28 + t66 * t30;
t17 = qJD(6) + t19;
t14 = t68 * t21 - t53 * t45;
t12 = t53 * t21 + t68 * t45;
t7 = t19 * pkin(5) - t21 * pkin(9) + t18;
t4 = -t45 * pkin(9) + t6;
t3 = t45 * pkin(5) - t5;
t2 = t68 * t4 + t53 * t7;
t1 = -t53 * t4 + t68 * t7;
t8 = [0, 0, 0, 0, 0, t69, 0, 0, 0, 0, t36, t61, t38, t37, t60, t47, pkin(1) * t67 - pkin(7) * t38, -t58 * pkin(1) * t55 - pkin(7) * t60 (t49 + t50) * t58 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t58, t36, t38, -t61, t47, -t60, t37, -t33 * qJD(2) - t32 * t64 (t33 * t55 + t34 * t57) * qJD(1), t34 * qJD(2) - t32 * t65, t34 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, t30 * t28, -t30 * t45, t28 ^ 2 / 0.2e1, -t28 * t45, t40, -t15 * t45 - t24 * t28, t16 * t45 + t24 * t30, -t15 * t30 + t16 * t28, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, -t21 * t45, t19 ^ 2 / 0.2e1, t19 * t45, t40, t18 * t19 - t5 * t45, t18 * t21 + t6 * t45, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t17, t12 ^ 2 / 0.2e1, -t12 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t12, t3 * t14 - t2 * t17, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
