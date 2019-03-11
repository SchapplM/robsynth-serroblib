% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:02
% EndTime: 2019-03-09 08:48:02
% DurationCPUTime: 0.18s
% Computational Cost: add. (630->63), mult. (1545->141), div. (0->0), fcn. (1045->8), ass. (0->53)
t54 = qJD(1) ^ 2;
t70 = t54 / 0.2e1;
t48 = sin(pkin(10));
t52 = sin(qJ(2));
t53 = cos(qJ(2));
t64 = cos(pkin(10));
t33 = (t48 * t53 + t64 * t52) * qJD(1);
t63 = qJD(1) * t52;
t65 = pkin(7) + qJ(3);
t36 = qJD(2) * pkin(2) - t65 * t63;
t62 = qJD(1) * t53;
t37 = t65 * t62;
t23 = t64 * t36 - t48 * t37;
t56 = qJD(4) - t23;
t10 = -t33 * pkin(8) + (-pkin(3) - pkin(4)) * qJD(2) + t56;
t24 = t48 * t36 + t64 * t37;
t22 = qJD(2) * qJ(4) + t24;
t31 = t48 * t63 - t64 * t62;
t12 = t31 * pkin(8) + t22;
t51 = sin(qJ(5));
t69 = cos(qJ(5));
t7 = t51 * t10 + t69 * t12;
t68 = cos(qJ(6));
t67 = t33 * t31;
t66 = t53 * t54;
t61 = qJD(2) * t31;
t60 = t31 ^ 2 / 0.2e1;
t59 = qJD(1) * qJD(2);
t38 = -qJD(1) * pkin(1) - pkin(2) * t62 + qJD(3);
t58 = t52 * t59;
t57 = t53 * t59;
t19 = -t69 * t31 + t51 * t33;
t16 = t31 * pkin(3) - t33 * qJ(4) + t38;
t6 = t69 * t10 - t51 * t12;
t8 = -t31 * pkin(4) - t16;
t50 = sin(qJ(6));
t47 = t53 ^ 2;
t46 = t52 ^ 2;
t44 = qJD(2) ^ 2 / 0.2e1;
t42 = qJD(2) - qJD(5);
t28 = t33 * qJD(2);
t27 = t33 ^ 2 / 0.2e1;
t21 = t51 * t31 + t69 * t33;
t18 = -qJD(2) * pkin(3) + t56;
t17 = qJD(6) + t19;
t15 = t68 * t21 - t50 * t42;
t13 = t50 * t21 + t68 * t42;
t5 = -t42 * pkin(9) + t7;
t4 = t42 * pkin(5) - t6;
t3 = t19 * pkin(5) - t21 * pkin(9) + t8;
t2 = t50 * t3 + t68 * t5;
t1 = t68 * t3 - t50 * t5;
t9 = [0, 0, 0, 0, 0, t70, 0, 0, 0, 0, t46 * t70, t52 * t66, t58, t47 * t70, t57, t44, pkin(1) * t66 - pkin(7) * t58, -t54 * pkin(1) * t52 - pkin(7) * t57 (t46 + t47) * t54 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t47 / 0.2e1 + t46 / 0.2e1) * pkin(7) ^ 2) * t54, t27, -t67, t28, t60, -t61, t44, t23 * qJD(2) + t38 * t31, -t24 * qJD(2) + t38 * t33, -t23 * t33 - t24 * t31, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t27, t28, t67, t44, t61, t60, -t18 * qJD(2) + t16 * t31, t18 * t33 - t22 * t31, t22 * qJD(2) - t16 * t33, t22 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, -t21 * t42, t19 ^ 2 / 0.2e1, t19 * t42, t42 ^ 2 / 0.2e1, t8 * t19 - t6 * t42, t8 * t21 + t7 * t42, -t7 * t19 - t6 * t21, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t17, t13 ^ 2 / 0.2e1, -t13 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t4 * t13, t4 * t15 - t2 * t17, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t9;
