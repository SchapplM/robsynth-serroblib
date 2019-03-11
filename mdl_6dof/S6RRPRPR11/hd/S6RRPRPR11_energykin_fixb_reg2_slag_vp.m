% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR11
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:02
% EndTime: 2019-03-09 11:15:02
% DurationCPUTime: 0.19s
% Computational Cost: add. (748->66), mult. (1563->145), div. (0->0), fcn. (971->8), ass. (0->53)
t55 = qJD(1) ^ 2;
t68 = t55 / 0.2e1;
t67 = -pkin(2) - pkin(8);
t66 = cos(qJ(6));
t54 = cos(qJ(2));
t65 = t54 * t55;
t52 = sin(qJ(2));
t57 = -qJ(3) * t52 - pkin(1);
t26 = (t67 * t54 + t57) * qJD(1);
t62 = t52 * qJD(1);
t61 = pkin(7) * t62 + qJD(3);
t27 = pkin(3) * t62 + t67 * qJD(2) + t61;
t51 = sin(qJ(4));
t53 = cos(qJ(4));
t16 = -t51 * t26 + t53 * t27;
t63 = qJD(1) * t54;
t33 = t53 * qJD(2) - t51 * t63;
t39 = qJD(4) + t62;
t12 = t39 * pkin(4) - t33 * qJ(5) + t16;
t17 = t53 * t26 + t51 * t27;
t31 = t51 * qJD(2) + t53 * t63;
t15 = -t31 * qJ(5) + t17;
t49 = sin(pkin(10));
t64 = cos(pkin(10));
t6 = t49 * t12 + t64 * t15;
t35 = -pkin(7) * t63 - qJD(2) * qJ(3);
t60 = qJD(1) * qJD(2);
t29 = pkin(3) * t63 - t35;
t59 = t52 * t60;
t58 = t54 * t60;
t5 = t64 * t12 - t49 * t15;
t22 = t31 * pkin(4) + qJD(5) + t29;
t50 = sin(qJ(6));
t48 = t54 ^ 2;
t47 = t52 ^ 2;
t45 = qJD(2) ^ 2 / 0.2e1;
t41 = t48 * t68;
t40 = t47 * t68;
t38 = t52 * t65;
t37 = qJD(6) + t39;
t36 = t39 ^ 2 / 0.2e1;
t34 = -qJD(2) * pkin(2) + t61;
t30 = (-pkin(2) * t54 + t57) * qJD(1);
t21 = -t49 * t31 + t64 * t33;
t19 = t64 * t31 + t49 * t33;
t13 = t19 * pkin(5) + t22;
t9 = -t50 * t19 + t66 * t21;
t7 = t66 * t19 + t50 * t21;
t4 = -t19 * pkin(9) + t6;
t3 = t39 * pkin(5) - t21 * pkin(9) + t5;
t2 = t50 * t3 + t66 * t4;
t1 = t66 * t3 - t50 * t4;
t8 = [0, 0, 0, 0, 0, t68, 0, 0, 0, 0, t40, t38, t59, t41, t58, t45, pkin(1) * t65 - pkin(7) * t59, -t55 * pkin(1) * t52 - pkin(7) * t58 (t47 + t48) * t55 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t48 / 0.2e1 + t47 / 0.2e1) * pkin(7) ^ 2) * t55, t45, -t59, -t58, t40, t38, t41 (t34 * t52 - t35 * t54) * qJD(1), t34 * qJD(2) + t30 * t63, -t35 * qJD(2) - t30 * t62, t30 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * t39, t31 ^ 2 / 0.2e1, -t31 * t39, t36, t16 * t39 + t29 * t31, -t17 * t39 + t29 * t33, -t16 * t33 - t17 * t31, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * t39, t19 ^ 2 / 0.2e1, -t19 * t39, t36, t22 * t19 + t5 * t39, t22 * t21 - t6 * t39, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t37, t7 ^ 2 / 0.2e1, -t7 * t37, t37 ^ 2 / 0.2e1, t1 * t37 + t13 * t7, t13 * t9 - t2 * t37, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg  = t8;
