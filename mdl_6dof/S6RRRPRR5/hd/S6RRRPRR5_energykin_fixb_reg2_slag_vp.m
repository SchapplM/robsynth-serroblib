% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR5
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:08
% EndTime: 2019-03-09 18:23:09
% DurationCPUTime: 0.22s
% Computational Cost: add. (735->65), mult. (1578->143), div. (0->0), fcn. (1072->8), ass. (0->56)
t47 = sin(qJ(3));
t48 = sin(qJ(2));
t49 = cos(qJ(3));
t50 = cos(qJ(2));
t32 = (t47 * t50 + t48 * t49) * qJD(1);
t51 = qJD(1) ^ 2;
t71 = t51 / 0.2e1;
t70 = pkin(3) + pkin(9);
t69 = -pkin(8) - pkin(7);
t42 = qJD(2) + qJD(3);
t62 = qJD(1) * t48;
t36 = qJD(2) * pkin(2) + t69 * t62;
t61 = qJD(1) * t50;
t37 = t69 * t61;
t20 = t49 * t36 + t47 * t37;
t55 = qJD(4) - t20;
t11 = t32 * pkin(4) - t70 * t42 + t55;
t30 = t47 * t62 - t49 * t61;
t38 = (-pkin(2) * t50 - pkin(1)) * qJD(1);
t53 = -t32 * qJ(4) + t38;
t12 = t70 * t30 + t53;
t46 = sin(qJ(5));
t68 = cos(qJ(5));
t6 = t46 * t11 + t68 * t12;
t67 = cos(qJ(6));
t66 = t32 * t30;
t65 = t32 * t42;
t64 = t42 * t30;
t63 = t50 * t51;
t21 = t47 * t36 - t49 * t37;
t60 = t30 ^ 2 / 0.2e1;
t59 = t32 ^ 2 / 0.2e1;
t58 = qJD(1) * qJD(2);
t19 = -t42 * qJ(4) - t21;
t5 = t68 * t11 - t46 * t12;
t57 = t48 * t58;
t56 = t50 * t58;
t16 = -t30 * pkin(4) - t19;
t29 = qJD(5) + t32;
t45 = sin(qJ(6));
t44 = t50 ^ 2;
t43 = t48 ^ 2;
t40 = t42 ^ 2 / 0.2e1;
t27 = qJD(6) + t29;
t25 = t46 * t30 + t68 * t42;
t23 = -t68 * t30 + t46 * t42;
t18 = -t42 * pkin(3) + t55;
t17 = t30 * pkin(3) + t53;
t15 = -t45 * t23 + t67 * t25;
t13 = t67 * t23 + t45 * t25;
t7 = t23 * pkin(5) + t16;
t4 = -t23 * pkin(10) + t6;
t3 = t29 * pkin(5) - t25 * pkin(10) + t5;
t2 = t45 * t3 + t67 * t4;
t1 = t67 * t3 - t45 * t4;
t8 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t43 * t71, t48 * t63, t57, t44 * t71, t56, qJD(2) ^ 2 / 0.2e1, pkin(1) * t63 - pkin(7) * t57, -t51 * pkin(1) * t48 - pkin(7) * t56 (t43 + t44) * t51 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t44 / 0.2e1 + t43 / 0.2e1) * pkin(7) ^ 2) * t51, t59, -t66, t65, t60, -t64, t40, t20 * t42 + t38 * t30, -t21 * t42 + t38 * t32, -t20 * t32 - t21 * t30, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t40, -t65, t64, t59, -t66, t60, t18 * t32 + t19 * t30, -t17 * t30 + t18 * t42, -t17 * t32 - t19 * t42, t17 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t29, t23 ^ 2 / 0.2e1, -t23 * t29, t29 ^ 2 / 0.2e1, t16 * t23 + t5 * t29, t16 * t25 - t6 * t29, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t27, t13 ^ 2 / 0.2e1, -t13 * t27, t27 ^ 2 / 0.2e1, t1 * t27 + t7 * t13, t7 * t15 - t2 * t27, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
