% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR8
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:06
% EndTime: 2019-03-09 09:26:07
% DurationCPUTime: 0.18s
% Computational Cost: add. (677->65), mult. (1552->143), div. (0->0), fcn. (1020->8), ass. (0->51)
t53 = qJD(1) ^ 2;
t68 = t53 / 0.2e1;
t67 = cos(qJ(5));
t66 = cos(qJ(6));
t47 = sin(pkin(10));
t51 = sin(qJ(2));
t62 = qJD(1) * t51;
t63 = cos(pkin(10));
t32 = -t63 * qJD(2) + t47 * t62;
t34 = t47 * qJD(2) + t63 * t62;
t65 = t34 * t32;
t52 = cos(qJ(2));
t64 = t52 * t53;
t30 = (-pkin(2) * t52 - qJ(3) * t51 - pkin(1)) * qJD(1);
t61 = t52 * qJD(1);
t38 = pkin(7) * t61 + qJD(2) * qJ(3);
t24 = t63 * t30 - t47 * t38;
t18 = pkin(3) * t61 + qJD(4) - t24;
t13 = pkin(4) * t61 - t34 * pkin(8) + t18;
t25 = t47 * t30 + t63 * t38;
t19 = -qJ(4) * t61 + t25;
t16 = t32 * pkin(8) + t19;
t50 = sin(qJ(5));
t6 = t50 * t13 + t67 * t16;
t60 = t32 ^ 2 / 0.2e1;
t59 = qJD(1) * qJD(2);
t58 = t32 * t61;
t57 = t34 * t61;
t37 = -qJD(2) * pkin(2) + pkin(7) * t62 + qJD(3);
t56 = t51 * t59;
t55 = t52 * t59;
t5 = t67 * t13 - t50 * t16;
t40 = qJD(5) + t61;
t17 = t32 * pkin(3) - t34 * qJ(4) + t37;
t15 = -t32 * pkin(4) - t17;
t49 = sin(qJ(6));
t46 = t52 ^ 2;
t45 = t51 ^ 2;
t41 = t46 * t68;
t39 = qJD(6) + t40;
t28 = t34 ^ 2 / 0.2e1;
t23 = t50 * t32 + t67 * t34;
t21 = -t67 * t32 + t50 * t34;
t10 = -t49 * t21 + t66 * t23;
t8 = t66 * t21 + t49 * t23;
t7 = t21 * pkin(5) + t15;
t4 = -t21 * pkin(9) + t6;
t3 = t40 * pkin(5) - t23 * pkin(9) + t5;
t2 = t49 * t3 + t66 * t4;
t1 = t66 * t3 - t49 * t4;
t9 = [0, 0, 0, 0, 0, t68, 0, 0, 0, 0, t45 * t68, t51 * t64, t56, t41, t55, qJD(2) ^ 2 / 0.2e1, pkin(1) * t64 - pkin(7) * t56, -t53 * pkin(1) * t51 - pkin(7) * t55 (t45 + t46) * t53 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t46 / 0.2e1 + t45 / 0.2e1) * pkin(7) ^ 2) * t53, t28, -t65, -t57, t60, t58, t41, -t24 * t61 + t37 * t32, t25 * t61 + t37 * t34, -t24 * t34 - t25 * t32, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t28, -t57, t65, t41, -t58, t60, t17 * t32 + t18 * t61, t18 * t34 - t19 * t32, -t17 * t34 - t19 * t61, t19 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * t40, t21 ^ 2 / 0.2e1, -t21 * t40, t40 ^ 2 / 0.2e1, t15 * t21 + t5 * t40, t15 * t23 - t6 * t40, -t6 * t21 - t5 * t23, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t39, t8 ^ 2 / 0.2e1, -t8 * t39, t39 ^ 2 / 0.2e1, t1 * t39 + t7 * t8, t7 * t10 - t2 * t39, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t9;
