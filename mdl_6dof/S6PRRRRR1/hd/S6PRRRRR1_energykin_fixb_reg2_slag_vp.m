% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:39:50
% EndTime: 2019-03-09 00:39:50
% DurationCPUTime: 0.17s
% Computational Cost: add. (689->62), mult. (1614->154), div. (0->0), fcn. (1223->12), ass. (0->52)
t53 = qJD(2) ^ 2;
t65 = t53 / 0.2e1;
t50 = sin(qJ(2));
t44 = sin(pkin(6));
t61 = qJD(1) * t44;
t35 = qJD(2) * pkin(8) + t50 * t61;
t51 = cos(qJ(3));
t45 = cos(pkin(6));
t60 = qJD(1) * t45;
t39 = t51 * t60;
t49 = sin(qJ(3));
t25 = qJD(3) * pkin(3) + t39 + (-pkin(9) * qJD(2) - t35) * t49;
t28 = t51 * t35 + t49 * t60;
t58 = qJD(2) * t51;
t26 = pkin(9) * t58 + t28;
t48 = sin(qJ(4));
t64 = cos(qJ(4));
t13 = t48 * t25 + t64 * t26;
t59 = qJD(2) * t49;
t31 = t48 * t59 - t64 * t58;
t11 = -t31 * pkin(10) + t13;
t47 = sin(qJ(5));
t63 = cos(qJ(5));
t12 = t64 * t25 - t48 * t26;
t33 = (t48 * t51 + t64 * t49) * qJD(2);
t43 = qJD(3) + qJD(4);
t9 = t43 * pkin(4) - t33 * pkin(10) + t12;
t6 = t63 * t11 + t47 * t9;
t62 = cos(qJ(6));
t57 = qJD(2) * qJD(3);
t52 = cos(qJ(2));
t56 = t52 * t61;
t55 = qJD(2) * t61;
t18 = t63 * t31 + t47 * t33;
t5 = -t47 * t11 + t63 * t9;
t30 = -t56 + (-pkin(3) * t51 - pkin(2)) * qJD(2);
t23 = t31 * pkin(4) + t30;
t54 = qJD(1) ^ 2;
t46 = sin(qJ(6));
t41 = qJD(5) + t43;
t36 = -qJD(2) * pkin(2) - t56;
t27 = -t49 * t35 + t39;
t20 = -t47 * t31 + t63 * t33;
t17 = qJD(6) + t18;
t16 = t62 * t20 + t46 * t41;
t14 = t46 * t20 - t62 * t41;
t7 = t18 * pkin(5) - t20 * pkin(11) + t23;
t4 = t41 * pkin(11) + t6;
t3 = -t41 * pkin(5) - t5;
t2 = t62 * t4 + t46 * t7;
t1 = -t46 * t4 + t62 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t54 / 0.2e1, 0, 0, 0, 0, 0, t65, t52 * t55, -t50 * t55, 0 (t45 ^ 2 / 0.2e1 + (t50 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1) * t44 ^ 2) * t54, t49 ^ 2 * t65, t49 * t53 * t51, t49 * t57, t51 ^ 2 * t65, t51 * t57, qJD(3) ^ 2 / 0.2e1, t27 * qJD(3) - t36 * t58, -t28 * qJD(3) + t36 * t59 (-t27 * t49 + t28 * t51) * qJD(2), t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * t43, t31 ^ 2 / 0.2e1, -t31 * t43, t43 ^ 2 / 0.2e1, t12 * t43 + t30 * t31, -t13 * t43 + t30 * t33, -t12 * t33 - t13 * t31, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t41, t18 ^ 2 / 0.2e1, -t18 * t41, t41 ^ 2 / 0.2e1, t23 * t18 + t5 * t41, t23 * t20 - t6 * t41, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t17, t14 ^ 2 / 0.2e1, -t14 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t14, t3 * t16 - t2 * t17, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
