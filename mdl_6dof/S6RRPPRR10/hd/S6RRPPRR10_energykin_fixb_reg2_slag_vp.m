% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR10
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:36:59
% EndTime: 2019-03-09 09:37:00
% DurationCPUTime: 0.21s
% Computational Cost: add. (728->66), mult. (1564->146), div. (0->0), fcn. (971->8), ass. (0->52)
t54 = qJD(1) ^ 2;
t67 = t54 / 0.2e1;
t66 = cos(qJ(5));
t65 = cos(qJ(6));
t53 = cos(qJ(2));
t64 = t53 * t54;
t63 = -pkin(2) - qJ(4);
t52 = sin(qJ(2));
t56 = -qJ(3) * t52 - pkin(1);
t26 = (t63 * t53 + t56) * qJD(1);
t61 = t52 * qJD(1);
t60 = pkin(7) * t61 + qJD(3);
t27 = pkin(3) * t61 + t63 * qJD(2) + t60;
t48 = sin(pkin(10));
t49 = cos(pkin(10));
t16 = -t48 * t26 + t49 * t27;
t62 = qJD(1) * t53;
t33 = t49 * qJD(2) - t48 * t62;
t12 = pkin(4) * t61 - t33 * pkin(8) + t16;
t17 = t49 * t26 + t48 * t27;
t31 = t48 * qJD(2) + t49 * t62;
t15 = -t31 * pkin(8) + t17;
t51 = sin(qJ(5));
t6 = t51 * t12 + t66 * t15;
t35 = -pkin(7) * t62 - qJD(2) * qJ(3);
t59 = qJD(1) * qJD(2);
t58 = t52 * t59;
t57 = t53 * t59;
t5 = t66 * t12 - t51 * t15;
t29 = pkin(3) * t62 + qJD(4) - t35;
t38 = qJD(5) + t61;
t22 = t31 * pkin(4) + t29;
t50 = sin(qJ(6));
t47 = t53 ^ 2;
t46 = t52 ^ 2;
t44 = qJD(2) ^ 2 / 0.2e1;
t40 = t47 * t67;
t39 = t46 * t67;
t37 = t52 * t64;
t36 = qJD(6) + t38;
t34 = -qJD(2) * pkin(2) + t60;
t30 = (-pkin(2) * t53 + t56) * qJD(1);
t21 = -t51 * t31 + t66 * t33;
t19 = t66 * t31 + t51 * t33;
t13 = t19 * pkin(5) + t22;
t9 = -t50 * t19 + t65 * t21;
t7 = t65 * t19 + t50 * t21;
t4 = -t19 * pkin(9) + t6;
t3 = t38 * pkin(5) - t21 * pkin(9) + t5;
t2 = t50 * t3 + t65 * t4;
t1 = t65 * t3 - t50 * t4;
t8 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t39, t37, t58, t40, t57, t44, pkin(1) * t64 - pkin(7) * t58, -t54 * pkin(1) * t52 - pkin(7) * t57 (t46 + t47) * t54 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t47 / 0.2e1 + t46 / 0.2e1) * pkin(7) ^ 2) * t54, t44, -t58, -t57, t39, t37, t40 (t34 * t52 - t35 * t53) * qJD(1), t34 * qJD(2) + t30 * t62, -t35 * qJD(2) - t30 * t61, t30 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * t61, t31 ^ 2 / 0.2e1, -t31 * t61, t39, t16 * t61 + t29 * t31, -t17 * t61 + t29 * t33, -t16 * t33 - t17 * t31, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * t38, t19 ^ 2 / 0.2e1, -t19 * t38, t38 ^ 2 / 0.2e1, t22 * t19 + t5 * t38, t22 * t21 - t6 * t38, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t36, t7 ^ 2 / 0.2e1, -t7 * t36, t36 ^ 2 / 0.2e1, t1 * t36 + t13 * t7, t13 * t9 - t2 * t36, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg  = t8;
