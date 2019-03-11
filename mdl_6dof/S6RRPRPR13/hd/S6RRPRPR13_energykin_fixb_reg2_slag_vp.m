% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR13_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:29:32
% EndTime: 2019-03-09 11:29:32
% DurationCPUTime: 0.22s
% Computational Cost: add. (943->71), mult. (2214->153), div. (0->0), fcn. (1616->10), ass. (0->59)
t74 = -pkin(2) - pkin(9);
t56 = sin(qJ(2));
t52 = sin(pkin(6));
t70 = qJD(1) * t52;
t64 = t56 * t70;
t45 = pkin(8) * t64;
t53 = cos(pkin(6));
t69 = t53 * qJD(1);
t49 = qJD(2) + t69;
t58 = cos(qJ(2));
t20 = qJD(3) + t45 + t74 * t49 + (-pkin(1) * t53 * t58 + pkin(3) * t52 * t56) * qJD(1);
t62 = -qJ(3) * t56 - pkin(1);
t28 = (t74 * t58 + t62) * t70;
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t17 = t55 * t20 + t57 * t28;
t40 = qJD(4) + t64;
t10 = t40 * qJ(5) + t17;
t65 = t58 * t70;
t66 = pkin(1) * t69;
t37 = pkin(8) * t65 + t56 * t66;
t30 = -t49 * qJ(3) - t37;
t27 = pkin(3) * t65 - t30;
t33 = t55 * t49 + t57 * t65;
t35 = t57 * t49 - t55 * t65;
t18 = t33 * pkin(4) - t35 * qJ(5) + t27;
t51 = sin(pkin(11));
t71 = cos(pkin(11));
t6 = t71 * t10 + t51 * t18;
t73 = cos(qJ(6));
t59 = qJD(1) ^ 2;
t72 = t52 ^ 2 * t59;
t68 = t33 ^ 2 / 0.2e1;
t67 = t58 * t72;
t63 = t72 / 0.2e1;
t5 = -t51 * t10 + t71 * t18;
t16 = t57 * t20 - t55 * t28;
t61 = t49 * t64;
t60 = t49 * t65;
t36 = t58 * t66 - t45;
t9 = -t40 * pkin(4) + qJD(5) - t16;
t54 = sin(qJ(6));
t43 = t58 ^ 2 * t63;
t42 = t56 ^ 2 * t63;
t41 = t49 ^ 2 / 0.2e1;
t39 = t56 * t67;
t32 = qJD(6) + t33;
t31 = (-pkin(2) * t58 + t62) * t70;
t29 = -t49 * pkin(2) + qJD(3) - t36;
t24 = t71 * t35 + t51 * t40;
t22 = t51 * t35 - t71 * t40;
t13 = -t54 * t22 + t73 * t24;
t11 = t73 * t22 + t54 * t24;
t7 = t22 * pkin(5) + t9;
t4 = -t22 * pkin(10) + t6;
t3 = t33 * pkin(5) - t24 * pkin(10) + t5;
t2 = t54 * t3 + t73 * t4;
t1 = t73 * t3 - t54 * t4;
t8 = [0, 0, 0, 0, 0, t59 / 0.2e1, 0, 0, 0, 0, t42, t39, t61, t43, t60, t41, pkin(1) * t67 + t36 * t49, -pkin(1) * t56 * t72 - t37 * t49 (-t36 * t56 + t37 * t58) * t70, t37 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t63, t41, -t61, -t60, t42, t39, t43 (t29 * t56 - t30 * t58) * t70, t29 * t49 + t31 * t65, -t30 * t49 - t31 * t64, t31 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * t40, t68, -t33 * t40, t40 ^ 2 / 0.2e1, t16 * t40 + t27 * t33, -t17 * t40 + t27 * t35, -t16 * t35 - t17 * t33, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, t24 * t33, t22 ^ 2 / 0.2e1, -t22 * t33, t68, t9 * t22 + t5 * t33, t9 * t24 - t6 * t33, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t32, t11 ^ 2 / 0.2e1, -t11 * t32, t32 ^ 2 / 0.2e1, t1 * t32 + t7 * t11, t7 * t13 - t2 * t32, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
