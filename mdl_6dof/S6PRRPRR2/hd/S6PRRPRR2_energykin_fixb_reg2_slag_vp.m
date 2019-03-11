% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:15
% EndTime: 2019-03-08 22:00:15
% DurationCPUTime: 0.17s
% Computational Cost: add. (629->62), mult. (1509->151), div. (0->0), fcn. (1124->12), ass. (0->53)
t55 = qJD(2) ^ 2;
t67 = t55 / 0.2e1;
t52 = sin(qJ(2));
t46 = sin(pkin(6));
t64 = qJD(1) * t46;
t37 = qJD(2) * pkin(8) + t52 * t64;
t53 = cos(qJ(3));
t48 = cos(pkin(6));
t63 = qJD(1) * t48;
t41 = t53 * t63;
t51 = sin(qJ(3));
t60 = qJ(4) * qJD(2);
t22 = qJD(3) * pkin(3) + t41 + (-t37 - t60) * t51;
t29 = t53 * t37 + t51 * t63;
t23 = t53 * t60 + t29;
t45 = sin(pkin(12));
t47 = cos(pkin(12));
t12 = t45 * t22 + t47 * t23;
t10 = qJD(3) * pkin(9) + t12;
t54 = cos(qJ(2));
t58 = t54 * t64;
t31 = -t58 + qJD(4) + (-pkin(3) * t53 - pkin(2)) * qJD(2);
t61 = qJD(2) * t53;
t62 = qJD(2) * t51;
t33 = t45 * t62 - t47 * t61;
t35 = (t45 * t53 + t47 * t51) * qJD(2);
t18 = t33 * pkin(4) - t35 * pkin(9) + t31;
t50 = sin(qJ(5));
t66 = cos(qJ(5));
t6 = t66 * t10 + t50 * t18;
t65 = cos(qJ(6));
t59 = qJD(2) * qJD(3);
t57 = qJD(2) * t64;
t5 = -t50 * t10 + t66 * t18;
t11 = t47 * t22 - t45 * t23;
t32 = qJD(5) + t33;
t9 = -qJD(3) * pkin(4) - t11;
t56 = qJD(1) ^ 2;
t49 = sin(qJ(6));
t44 = qJD(3) ^ 2 / 0.2e1;
t38 = -qJD(2) * pkin(2) - t58;
t30 = qJD(6) + t32;
t28 = -t51 * t37 + t41;
t27 = t50 * qJD(3) + t66 * t35;
t25 = -t66 * qJD(3) + t50 * t35;
t15 = -t49 * t25 + t65 * t27;
t13 = t65 * t25 + t49 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(10) + t6;
t3 = t32 * pkin(5) - t27 * pkin(10) + t5;
t2 = t49 * t3 + t65 * t4;
t1 = t65 * t3 - t49 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t56 / 0.2e1, 0, 0, 0, 0, 0, t67, t54 * t57, -t52 * t57, 0 (t48 ^ 2 / 0.2e1 + (t52 ^ 2 / 0.2e1 + t54 ^ 2 / 0.2e1) * t46 ^ 2) * t56, t51 ^ 2 * t67, t51 * t55 * t53, t51 * t59, t53 ^ 2 * t67, t53 * t59, t44, t28 * qJD(3) - t38 * t61, -t29 * qJD(3) + t38 * t62 (-t28 * t51 + t29 * t53) * qJD(2), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * qJD(3), t33 ^ 2 / 0.2e1, -t33 * qJD(3), t44, t11 * qJD(3) + t31 * t33, -t12 * qJD(3) + t31 * t35, -t11 * t35 - t12 * t33, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t32, t25 ^ 2 / 0.2e1, -t25 * t32, t32 ^ 2 / 0.2e1, t9 * t25 + t5 * t32, t9 * t27 - t6 * t32, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t30, t13 ^ 2 / 0.2e1, -t13 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t13, t7 * t15 - t2 * t30, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
