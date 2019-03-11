% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:27:36
% EndTime: 2019-03-08 22:27:37
% DurationCPUTime: 0.26s
% Computational Cost: add. (958->64), mult. (2375->157), div. (0->0), fcn. (1925->14), ass. (0->57)
t60 = cos(qJ(2));
t52 = sin(pkin(6));
t71 = qJD(1) * t52;
t40 = qJD(2) * pkin(2) + t60 * t71;
t51 = sin(pkin(7));
t53 = cos(pkin(7));
t54 = cos(pkin(6));
t70 = qJD(1) * t54;
t77 = t40 * t53 + t51 * t70;
t58 = sin(qJ(2));
t69 = qJD(2) * t51;
t38 = pkin(9) * t69 + t58 * t71;
t57 = sin(qJ(3));
t59 = cos(qJ(3));
t25 = -t57 * t38 + t77 * t59;
t26 = t59 * t38 + t77 * t57;
t47 = t53 * qJD(2) + qJD(3);
t20 = t47 * qJ(4) + t26;
t46 = t53 * t70;
t29 = t46 + (-t40 + (-pkin(3) * t59 - qJ(4) * t57) * qJD(2)) * t51;
t50 = sin(pkin(13));
t72 = cos(pkin(13));
t13 = t72 * t20 + t50 * t29;
t67 = t57 * t69;
t32 = -t72 * t47 + t50 * t67;
t11 = -t32 * pkin(10) + t13;
t56 = sin(qJ(5));
t76 = cos(qJ(5));
t12 = -t50 * t20 + t72 * t29;
t34 = t50 * t47 + t72 * t67;
t66 = t59 * t69;
t9 = -pkin(4) * t66 - t34 * pkin(10) + t12;
t6 = t76 * t11 + t56 * t9;
t75 = cos(qJ(6));
t61 = qJD(2) ^ 2;
t73 = t51 ^ 2 * t61;
t65 = t73 / 0.2e1;
t64 = qJD(2) * t71;
t22 = t76 * t32 + t56 * t34;
t5 = -t56 * t11 + t76 * t9;
t19 = -t47 * pkin(3) + qJD(4) - t25;
t14 = t32 * pkin(4) + t19;
t62 = qJD(1) ^ 2;
t55 = sin(qJ(6));
t45 = t59 ^ 2 * t65;
t43 = -qJD(5) + t66;
t31 = -t51 * t40 + t46;
t24 = -t56 * t32 + t76 * t34;
t21 = qJD(6) + t22;
t17 = t75 * t24 - t55 * t43;
t15 = t55 * t24 + t75 * t43;
t7 = t22 * pkin(5) - t24 * pkin(11) + t14;
t4 = -t43 * pkin(11) + t6;
t3 = t43 * pkin(5) - t5;
t2 = t75 * t4 + t55 * t7;
t1 = -t55 * t4 + t75 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t62 / 0.2e1, 0, 0, 0, 0, 0, t61 / 0.2e1, t60 * t64, -t58 * t64, 0 (t54 ^ 2 / 0.2e1 + (t58 ^ 2 / 0.2e1 + t60 ^ 2 / 0.2e1) * t52 ^ 2) * t62, t57 ^ 2 * t65, t57 * t59 * t73, t47 * t67, t45, t47 * t66, t47 ^ 2 / 0.2e1, t25 * t47 - t31 * t66, -t26 * t47 + t31 * t67 (-t25 * t57 + t26 * t59) * t69, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, -t34 * t66, t32 ^ 2 / 0.2e1, t32 * t66, t45, -t12 * t66 + t19 * t32, t13 * t66 + t19 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, -t24 * t43, t22 ^ 2 / 0.2e1, t22 * t43, t43 ^ 2 / 0.2e1, t14 * t22 - t5 * t43, t14 * t24 + t6 * t43, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t21, t15 ^ 2 / 0.2e1, -t15 * t21, t21 ^ 2 / 0.2e1, t1 * t21 + t3 * t15, t3 * t17 - t2 * t21, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
