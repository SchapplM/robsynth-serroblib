% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:26
% EndTime: 2019-03-08 22:14:27
% DurationCPUTime: 0.17s
% Computational Cost: add. (374->57), mult. (839->133), div. (0->0), fcn. (542->10), ass. (0->50)
t47 = sin(qJ(5));
t48 = sin(qJ(3));
t50 = cos(qJ(5));
t51 = cos(qJ(3));
t22 = (t47 * t48 + t50 * t51) * qJD(2);
t53 = qJD(2) ^ 2;
t66 = t53 / 0.2e1;
t49 = sin(qJ(2));
t43 = sin(pkin(6));
t64 = qJD(1) * t43;
t27 = qJD(2) * pkin(8) + t49 * t64;
t44 = cos(pkin(6));
t63 = qJD(1) * t44;
t19 = t51 * t27 + t48 * t63;
t14 = qJD(3) * qJ(4) + t19;
t61 = qJD(2) * t51;
t11 = -pkin(9) * t61 + t14;
t18 = -t48 * t27 + t51 * t63;
t55 = qJD(4) - t18;
t62 = qJD(2) * t48;
t9 = -pkin(9) * t62 + (-pkin(3) - pkin(4)) * qJD(3) + t55;
t6 = t50 * t11 + t47 * t9;
t65 = cos(qJ(6));
t52 = cos(qJ(2));
t28 = -qJD(2) * pkin(2) - t52 * t64;
t60 = qJD(2) * qJD(3);
t59 = t48 * t53 * t51;
t58 = qJD(2) * t64;
t57 = t51 * t60;
t20 = -pkin(3) * t61 - qJ(4) * t62 + t28;
t12 = pkin(4) * t61 - t20;
t5 = -t47 * t11 + t50 * t9;
t54 = qJD(1) ^ 2;
t46 = sin(qJ(6));
t41 = qJD(3) ^ 2 / 0.2e1;
t39 = qJD(3) - qJD(5);
t34 = t48 * t60;
t33 = t51 ^ 2 * t66;
t32 = t48 ^ 2 * t66;
t24 = (-t47 * t51 + t48 * t50) * qJD(2);
t21 = qJD(6) + t22;
t17 = t65 * t24 - t46 * t39;
t15 = t46 * t24 + t65 * t39;
t13 = -qJD(3) * pkin(3) + t55;
t7 = t22 * pkin(5) - t24 * pkin(10) + t12;
t4 = -t39 * pkin(10) + t6;
t3 = t39 * pkin(5) - t5;
t2 = t65 * t4 + t46 * t7;
t1 = -t46 * t4 + t65 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t54 / 0.2e1, 0, 0, 0, 0, 0, t66, t52 * t58, -t49 * t58, 0 (t44 ^ 2 / 0.2e1 + (t49 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1) * t43 ^ 2) * t54, t32, t59, t34, t33, t57, t41, t18 * qJD(3) - t28 * t61, -t19 * qJD(3) + t28 * t62 (-t18 * t48 + t19 * t51) * qJD(2), t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t32, t34, -t59, t41, -t57, t33, -t13 * qJD(3) - t20 * t61 (t13 * t48 + t14 * t51) * qJD(2), t14 * qJD(3) - t20 * t62, t14 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, -t24 * t39, t22 ^ 2 / 0.2e1, t22 * t39, t39 ^ 2 / 0.2e1, t12 * t22 - t5 * t39, t12 * t24 + t6 * t39, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t21, t15 ^ 2 / 0.2e1, -t15 * t21, t21 ^ 2 / 0.2e1, t1 * t21 + t3 * t15, t3 * t17 - t2 * t21, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
