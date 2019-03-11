% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:44:36
% EndTime: 2019-03-08 23:44:36
% DurationCPUTime: 0.25s
% Computational Cost: add. (973->64), mult. (2342->157), div. (0->0), fcn. (1897->14), ass. (0->57)
t60 = cos(qJ(2));
t52 = sin(pkin(6));
t72 = qJD(1) * t52;
t41 = qJD(2) * pkin(2) + t60 * t72;
t51 = sin(pkin(7));
t53 = cos(pkin(7));
t54 = cos(pkin(6));
t71 = qJD(1) * t54;
t78 = t41 * t53 + t51 * t71;
t58 = sin(qJ(2));
t70 = qJD(2) * t51;
t39 = pkin(9) * t70 + t58 * t72;
t57 = sin(qJ(3));
t59 = cos(qJ(3));
t23 = -t57 * t39 + t78 * t59;
t24 = t59 * t39 + t78 * t57;
t47 = t53 * qJD(2) + qJD(3);
t22 = t47 * pkin(10) + t24;
t46 = t53 * t71;
t27 = t46 + (-t41 + (-pkin(3) * t59 - pkin(10) * t57) * qJD(2)) * t51;
t56 = sin(qJ(4));
t77 = cos(qJ(4));
t14 = t77 * t22 + t56 * t27;
t66 = t59 * t70;
t44 = -qJD(4) + t66;
t10 = -t44 * qJ(5) + t14;
t21 = -t47 * pkin(3) - t23;
t67 = t57 * t70;
t33 = -t77 * t47 + t56 * t67;
t35 = t56 * t47 + t77 * t67;
t15 = t33 * pkin(4) - t35 * qJ(5) + t21;
t50 = sin(pkin(13));
t73 = cos(pkin(13));
t6 = t73 * t10 + t50 * t15;
t76 = cos(qJ(6));
t61 = qJD(2) ^ 2;
t74 = t51 ^ 2 * t61;
t69 = t33 ^ 2 / 0.2e1;
t65 = t74 / 0.2e1;
t64 = qJD(2) * t72;
t5 = -t50 * t10 + t73 * t15;
t13 = -t56 * t22 + t77 * t27;
t9 = t44 * pkin(4) + qJD(5) - t13;
t62 = qJD(1) ^ 2;
t55 = sin(qJ(6));
t32 = qJD(6) + t33;
t31 = -t51 * t41 + t46;
t30 = t73 * t35 - t50 * t44;
t28 = t50 * t35 + t73 * t44;
t18 = -t55 * t28 + t76 * t30;
t16 = t76 * t28 + t55 * t30;
t7 = t28 * pkin(5) + t9;
t4 = -t28 * pkin(11) + t6;
t3 = t33 * pkin(5) - t30 * pkin(11) + t5;
t2 = t55 * t3 + t76 * t4;
t1 = t76 * t3 - t55 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t62 / 0.2e1, 0, 0, 0, 0, 0, t61 / 0.2e1, t60 * t64, -t58 * t64, 0 (t54 ^ 2 / 0.2e1 + (t58 ^ 2 / 0.2e1 + t60 ^ 2 / 0.2e1) * t52 ^ 2) * t62, t57 ^ 2 * t65, t57 * t59 * t74, t47 * t67, t59 ^ 2 * t65, t47 * t66, t47 ^ 2 / 0.2e1, t23 * t47 - t31 * t66, -t24 * t47 + t31 * t67 (-t23 * t57 + t24 * t59) * t70, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, -t35 * t44, t69, t33 * t44, t44 ^ 2 / 0.2e1, -t13 * t44 + t21 * t33, t14 * t44 + t21 * t35, -t13 * t35 - t14 * t33, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * t33, t28 ^ 2 / 0.2e1, -t28 * t33, t69, t9 * t28 + t5 * t33, t9 * t30 - t6 * t33, -t6 * t28 - t5 * t30, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t32, t16 ^ 2 / 0.2e1, -t16 * t32, t32 ^ 2 / 0.2e1, t1 * t32 + t7 * t16, t7 * t18 - t2 * t32, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
