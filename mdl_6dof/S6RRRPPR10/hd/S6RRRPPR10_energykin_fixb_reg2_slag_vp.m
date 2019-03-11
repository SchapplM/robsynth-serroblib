% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:27:52
% EndTime: 2019-03-09 16:27:52
% DurationCPUTime: 0.22s
% Computational Cost: add. (1021->69), mult. (2362->148), div. (0->0), fcn. (1790->10), ass. (0->57)
t68 = cos(pkin(6)) * qJD(1);
t49 = qJD(2) + t68;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t56 = sin(qJ(2));
t52 = sin(pkin(6));
t69 = qJD(1) * t52;
t64 = t56 * t69;
t37 = t55 * t49 + t57 * t64;
t58 = cos(qJ(2));
t63 = t58 * t69;
t44 = -qJD(3) + t63;
t65 = pkin(1) * t68;
t40 = pkin(8) * t63 + t56 * t65;
t30 = t49 * pkin(9) + t40;
t32 = (-pkin(2) * t58 - pkin(9) * t56 - pkin(1)) * t69;
t20 = -t55 * t30 + t57 * t32;
t61 = qJD(4) - t20;
t71 = pkin(3) + qJ(5);
t10 = t37 * pkin(4) + t71 * t44 + t61;
t35 = -t57 * t49 + t55 * t64;
t39 = -pkin(8) * t64 + t58 * t65;
t29 = -t49 * pkin(2) - t39;
t60 = -t37 * qJ(4) + t29;
t12 = t71 * t35 + t60;
t51 = sin(pkin(11));
t70 = cos(pkin(11));
t6 = t51 * t10 + t70 * t12;
t76 = cos(qJ(6));
t75 = t37 * t35;
t74 = t37 * t44;
t73 = t44 * t35;
t59 = qJD(1) ^ 2;
t72 = t52 ^ 2 * t59;
t21 = t57 * t30 + t55 * t32;
t67 = t35 ^ 2 / 0.2e1;
t33 = t37 ^ 2 / 0.2e1;
t66 = t58 * t72;
t19 = t44 * qJ(4) - t21;
t62 = t72 / 0.2e1;
t5 = t70 * t10 - t51 * t12;
t13 = -t35 * pkin(4) + qJD(5) - t19;
t54 = sin(qJ(6));
t41 = t44 ^ 2 / 0.2e1;
t34 = qJD(6) + t37;
t25 = t51 * t35 - t70 * t44;
t23 = -t70 * t35 - t51 * t44;
t18 = t44 * pkin(3) + t61;
t17 = t35 * pkin(3) + t60;
t16 = -t54 * t23 + t76 * t25;
t14 = t76 * t23 + t54 * t25;
t7 = t23 * pkin(5) + t13;
t4 = -t23 * pkin(10) + t6;
t3 = t37 * pkin(5) - t25 * pkin(10) + t5;
t2 = t54 * t3 + t76 * t4;
t1 = t76 * t3 - t54 * t4;
t8 = [0, 0, 0, 0, 0, t59 / 0.2e1, 0, 0, 0, 0, t56 ^ 2 * t62, t56 * t66, t49 * t64, t58 ^ 2 * t62, t49 * t63, t49 ^ 2 / 0.2e1, pkin(1) * t66 + t39 * t49, -pkin(1) * t56 * t72 - t40 * t49 (-t39 * t56 + t40 * t58) * t69, t40 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t62, t33, -t75, -t74, t67, t73, t41, -t20 * t44 + t29 * t35, t21 * t44 + t29 * t37, -t20 * t37 - t21 * t35, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t41, t74, -t73, t33, -t75, t67, t18 * t37 + t19 * t35, -t17 * t35 - t18 * t44, -t17 * t37 + t19 * t44, t17 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t37, t23 ^ 2 / 0.2e1, -t23 * t37, t33, t13 * t23 + t5 * t37, t13 * t25 - t6 * t37, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t34, t14 ^ 2 / 0.2e1, -t14 * t34, t34 ^ 2 / 0.2e1, t1 * t34 + t7 * t14, t7 * t16 - t2 * t34, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
