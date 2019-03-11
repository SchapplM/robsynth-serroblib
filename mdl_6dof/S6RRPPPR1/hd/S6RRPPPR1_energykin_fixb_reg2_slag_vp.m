% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:16
% EndTime: 2019-03-09 08:08:17
% DurationCPUTime: 0.22s
% Computational Cost: add. (661->63), mult. (1651->138), div. (0->0), fcn. (1145->8), ass. (0->55)
t52 = qJD(1) ^ 2;
t71 = t52 / 0.2e1;
t70 = -pkin(4) - pkin(5);
t69 = cos(qJ(6));
t48 = sin(pkin(9));
t50 = sin(qJ(2));
t51 = cos(qJ(2));
t64 = cos(pkin(9));
t36 = (t48 * t51 + t64 * t50) * qJD(1);
t47 = sin(pkin(10));
t63 = cos(pkin(10));
t25 = -t63 * qJD(2) + t47 * t36;
t27 = t47 * qJD(2) + t63 * t36;
t68 = t27 * t25;
t61 = qJD(1) * t51;
t62 = qJD(1) * t50;
t34 = t48 * t62 - t64 * t61;
t67 = t34 * t25;
t66 = t51 * t52;
t65 = pkin(7) + qJ(3);
t41 = qJD(3) + (-pkin(2) * t51 - pkin(1)) * qJD(1);
t15 = t34 * pkin(3) - t36 * qJ(4) + t41;
t39 = qJD(2) * pkin(2) - t65 * t62;
t40 = t65 * t61;
t22 = t48 * t39 + t64 * t40;
t20 = qJD(2) * qJ(4) + t22;
t10 = t47 * t15 + t63 * t20;
t21 = t64 * t39 - t48 * t40;
t60 = t25 ^ 2 / 0.2e1;
t30 = t34 ^ 2 / 0.2e1;
t59 = qJD(1) * qJD(2);
t7 = t34 * qJ(5) + t10;
t58 = t50 * t59;
t57 = t51 * t59;
t9 = t63 * t15 - t47 * t20;
t56 = qJD(2) * pkin(3) - qJD(4) + t21;
t55 = qJD(5) - t9;
t54 = t27 * qJ(5) + t56;
t49 = sin(qJ(6));
t46 = t51 ^ 2;
t45 = t50 ^ 2;
t44 = qJD(2) ^ 2 / 0.2e1;
t32 = -qJD(6) + t34;
t24 = t27 ^ 2 / 0.2e1;
t18 = t27 * t34;
t13 = t49 * t25 + t69 * t27;
t11 = -t69 * t25 + t49 * t27;
t8 = t25 * pkin(4) - t54;
t6 = -t34 * pkin(4) + t55;
t5 = t70 * t25 + t54;
t4 = t25 * pkin(8) + t7;
t3 = -t27 * pkin(8) + t70 * t34 + t55;
t2 = t49 * t3 + t69 * t4;
t1 = t69 * t3 - t49 * t4;
t12 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t45 * t71, t50 * t66, t58, t46 * t71, t57, t44, pkin(1) * t66 - pkin(7) * t58, -t52 * pkin(1) * t50 - pkin(7) * t57 (t45 + t46) * t52 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t46 / 0.2e1 + t45 / 0.2e1) * pkin(7) ^ 2) * t52, t36 ^ 2 / 0.2e1, -t36 * t34, t36 * qJD(2), t30, -t34 * qJD(2), t44, t21 * qJD(2) + t41 * t34, -t22 * qJD(2) + t41 * t36, -t21 * t36 - t22 * t34, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1, t24, -t68, t18, t60, -t67, t30, -t25 * t56 + t9 * t34, -t10 * t34 - t27 * t56, -t10 * t25 - t9 * t27, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t56 ^ 2 / 0.2e1, t24, t18, t68, t30, t67, t60, t8 * t25 - t6 * t34, -t7 * t25 + t6 * t27, -t8 * t27 + t7 * t34, t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, -t13 * t32, t11 ^ 2 / 0.2e1, t11 * t32, t32 ^ 2 / 0.2e1, -t1 * t32 + t5 * t11, t5 * t13 + t2 * t32, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
