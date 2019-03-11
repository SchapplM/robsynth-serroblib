% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR13_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:53:18
% EndTime: 2019-03-09 14:53:18
% DurationCPUTime: 0.24s
% Computational Cost: add. (953->71), mult. (2214->156), div. (0->0), fcn. (1616->10), ass. (0->59)
t74 = -pkin(2) - pkin(9);
t57 = sin(qJ(2));
t52 = sin(pkin(6));
t70 = qJD(1) * t52;
t65 = t57 * t70;
t46 = pkin(8) * t65;
t53 = cos(pkin(6));
t69 = t53 * qJD(1);
t50 = qJD(2) + t69;
t59 = cos(qJ(2));
t20 = qJD(3) + t46 + t74 * t50 + (-pkin(1) * t53 * t59 + pkin(3) * t52 * t57) * qJD(1);
t63 = -qJ(3) * t57 - pkin(1);
t28 = (t74 * t59 + t63) * t70;
t56 = sin(qJ(4));
t58 = cos(qJ(4));
t12 = t56 * t20 + t58 * t28;
t41 = qJD(4) + t65;
t10 = t41 * pkin(10) + t12;
t66 = t59 * t70;
t67 = pkin(1) * t69;
t38 = pkin(8) * t66 + t57 * t67;
t30 = -t50 * qJ(3) - t38;
t27 = pkin(3) * t66 - t30;
t34 = t56 * t50 + t58 * t66;
t36 = t58 * t50 - t56 * t66;
t18 = t34 * pkin(4) - t36 * pkin(10) + t27;
t55 = sin(qJ(5));
t73 = cos(qJ(5));
t6 = t73 * t10 + t55 * t18;
t72 = cos(qJ(6));
t60 = qJD(1) ^ 2;
t71 = t52 ^ 2 * t60;
t68 = t59 * t71;
t64 = t71 / 0.2e1;
t5 = -t55 * t10 + t73 * t18;
t11 = t58 * t20 - t56 * t28;
t62 = t50 * t65;
t61 = t50 * t66;
t37 = t59 * t67 - t46;
t9 = -t41 * pkin(4) - t11;
t33 = qJD(5) + t34;
t54 = sin(qJ(6));
t44 = t59 ^ 2 * t64;
t43 = t57 ^ 2 * t64;
t42 = t50 ^ 2 / 0.2e1;
t40 = t57 * t68;
t32 = (-pkin(2) * t59 + t63) * t70;
t31 = qJD(6) + t33;
t29 = -t50 * pkin(2) + qJD(3) - t37;
t24 = t73 * t36 + t55 * t41;
t22 = t55 * t36 - t73 * t41;
t17 = -t54 * t22 + t72 * t24;
t15 = t72 * t22 + t54 * t24;
t7 = t22 * pkin(5) + t9;
t4 = -t22 * pkin(11) + t6;
t3 = t33 * pkin(5) - t24 * pkin(11) + t5;
t2 = t54 * t3 + t72 * t4;
t1 = t72 * t3 - t54 * t4;
t8 = [0, 0, 0, 0, 0, t60 / 0.2e1, 0, 0, 0, 0, t43, t40, t62, t44, t61, t42, pkin(1) * t68 + t37 * t50, -pkin(1) * t57 * t71 - t38 * t50 (-t37 * t57 + t38 * t59) * t70, t38 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t64, t42, -t62, -t61, t43, t40, t44 (t29 * t57 - t30 * t59) * t70, t29 * t50 + t32 * t66, -t30 * t50 - t32 * t65, t32 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t36 ^ 2 / 0.2e1, -t36 * t34, t36 * t41, t34 ^ 2 / 0.2e1, -t34 * t41, t41 ^ 2 / 0.2e1, t11 * t41 + t27 * t34, -t12 * t41 + t27 * t36, -t11 * t36 - t12 * t34, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, t24 * t33, t22 ^ 2 / 0.2e1, -t22 * t33, t33 ^ 2 / 0.2e1, t9 * t22 + t5 * t33, t9 * t24 - t6 * t33, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t31, t15 ^ 2 / 0.2e1, -t15 * t31, t31 ^ 2 / 0.2e1, t1 * t31 + t7 * t15, t7 * t17 - t2 * t31, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
