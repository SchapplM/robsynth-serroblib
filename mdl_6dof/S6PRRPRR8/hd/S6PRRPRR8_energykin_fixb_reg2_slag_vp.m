% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:36
% EndTime: 2019-03-08 22:40:36
% DurationCPUTime: 0.24s
% Computational Cost: add. (559->66), mult. (1393->146), div. (0->0), fcn. (1048->12), ass. (0->63)
t44 = sin(pkin(7));
t47 = cos(pkin(6));
t68 = qJD(1) * t47;
t54 = cos(qJ(2));
t45 = sin(pkin(6));
t69 = qJD(1) * t45;
t31 = qJD(2) * pkin(2) + t54 * t69;
t46 = cos(pkin(7));
t72 = t31 * t46;
t76 = t44 * t68 + t72;
t53 = cos(qJ(3));
t75 = t76 * t53;
t74 = -pkin(3) - pkin(10);
t40 = t46 * t68;
t50 = sin(qJ(3));
t70 = qJ(4) * t50;
t14 = t40 + (-t31 + (t74 * t53 - t70) * qJD(2)) * t44;
t49 = sin(qJ(5));
t52 = cos(qJ(5));
t41 = t46 * qJD(2) + qJD(3);
t51 = sin(qJ(2));
t67 = qJD(2) * t44;
t29 = pkin(9) * t67 + t51 * t69;
t27 = t50 * t29;
t65 = qJD(4) + t27;
t66 = qJD(2) * t50;
t9 = -t53 * t72 + (pkin(4) * t66 - t53 * t68) * t44 + t74 * t41 + t65;
t6 = t52 * t14 + t49 * t9;
t73 = cos(qJ(6));
t55 = qJD(2) ^ 2;
t71 = t44 ^ 2 * t55;
t16 = t53 * t29 + t76 * t50;
t63 = t53 * t67;
t62 = t44 * t66;
t61 = t71 / 0.2e1;
t60 = qJD(2) * t69;
t12 = -t41 * qJ(4) - t16;
t59 = t41 * t62;
t58 = t41 * t63;
t10 = pkin(4) * t63 - t12;
t5 = -t49 * t14 + t52 * t9;
t23 = t49 * t41 + t52 * t63;
t56 = qJD(1) ^ 2;
t48 = sin(qJ(6));
t37 = t53 ^ 2 * t61;
t36 = t50 ^ 2 * t61;
t35 = t41 ^ 2 / 0.2e1;
t34 = qJD(5) + t62;
t33 = t50 * t53 * t71;
t25 = t52 * t41 - t49 * t63;
t22 = qJD(6) + t23;
t21 = -t44 * t31 + t40;
t20 = t73 * t25 + t48 * t34;
t18 = t48 * t25 - t73 * t34;
t17 = t40 + (-t31 + (-pkin(3) * t53 - t70) * qJD(2)) * t44;
t15 = -t27 + t75;
t11 = -t41 * pkin(3) + t65 - t75;
t7 = t23 * pkin(5) - t25 * pkin(11) + t10;
t4 = t34 * pkin(11) + t6;
t3 = -t34 * pkin(5) - t5;
t2 = t73 * t4 + t48 * t7;
t1 = -t48 * t4 + t73 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t56 / 0.2e1, 0, 0, 0, 0, 0, t55 / 0.2e1, t54 * t60, -t51 * t60, 0 (t47 ^ 2 / 0.2e1 + (t51 ^ 2 / 0.2e1 + t54 ^ 2 / 0.2e1) * t45 ^ 2) * t56, t36, t33, t59, t37, t58, t35, t15 * t41 - t21 * t63, -t16 * t41 + t21 * t62 (-t15 * t50 + t16 * t53) * t67, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t35, -t59, -t58, t36, t33, t37 (t11 * t50 - t12 * t53) * t67, t11 * t41 + t17 * t63, -t12 * t41 - t17 * t62, t17 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t34, t23 ^ 2 / 0.2e1, -t23 * t34, t34 ^ 2 / 0.2e1, t10 * t23 + t5 * t34, t10 * t25 - t6 * t34, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t22, t18 ^ 2 / 0.2e1, -t18 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t3 * t18, -t2 * t22 + t3 * t20, -t1 * t20 - t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
