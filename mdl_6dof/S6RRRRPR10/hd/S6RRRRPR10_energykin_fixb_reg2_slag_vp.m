% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:08:30
% EndTime: 2019-03-09 23:08:30
% DurationCPUTime: 0.25s
% Computational Cost: add. (1125->69), mult. (2637->150), div. (0->0), fcn. (2044->10), ass. (0->58)
t77 = pkin(4) + pkin(11);
t76 = cos(qJ(3));
t75 = cos(qJ(6));
t69 = cos(pkin(6)) * qJD(1);
t49 = qJD(2) + t69;
t55 = sin(qJ(3));
t56 = sin(qJ(2));
t51 = sin(pkin(6));
t70 = qJD(1) * t51;
t64 = t56 * t70;
t35 = -t76 * t49 + t55 * t64;
t37 = t55 * t49 + t76 * t64;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t23 = t57 * t35 + t54 * t37;
t25 = -t54 * t35 + t57 * t37;
t74 = t25 * t23;
t58 = cos(qJ(2));
t63 = t58 * t70;
t44 = -qJD(3) + t63;
t41 = -qJD(4) + t44;
t73 = t25 * t41;
t72 = t41 * t23;
t59 = qJD(1) ^ 2;
t71 = t51 ^ 2 * t59;
t65 = pkin(1) * t69;
t39 = pkin(8) * t63 + t56 * t65;
t32 = t49 * pkin(9) + t39;
t34 = (-pkin(2) * t58 - pkin(9) * t56 - pkin(1)) * t70;
t19 = -t55 * t32 + t76 * t34;
t12 = -t44 * pkin(3) - t37 * pkin(10) + t19;
t20 = t76 * t32 + t55 * t34;
t15 = -t35 * pkin(10) + t20;
t9 = t54 * t12 + t57 * t15;
t68 = t23 ^ 2 / 0.2e1;
t67 = t25 ^ 2 / 0.2e1;
t66 = t58 * t71;
t62 = t71 / 0.2e1;
t8 = t57 * t12 - t54 * t15;
t7 = t41 * qJ(5) - t9;
t61 = qJD(5) - t8;
t38 = -pkin(8) * t64 + t58 * t65;
t31 = -t49 * pkin(2) - t38;
t27 = t35 * pkin(3) + t31;
t60 = -t25 * qJ(5) + t27;
t53 = sin(qJ(6));
t40 = t41 ^ 2 / 0.2e1;
t22 = qJD(6) + t25;
t18 = t53 * t23 - t75 * t41;
t16 = -t75 * t23 - t53 * t41;
t10 = t23 * pkin(4) + t60;
t6 = t41 * pkin(4) + t61;
t5 = t77 * t23 + t60;
t4 = -t23 * pkin(5) - t7;
t3 = t25 * pkin(5) + t77 * t41 + t61;
t2 = t53 * t3 + t75 * t5;
t1 = t75 * t3 - t53 * t5;
t11 = [0, 0, 0, 0, 0, t59 / 0.2e1, 0, 0, 0, 0, t56 ^ 2 * t62, t56 * t66, t49 * t64, t58 ^ 2 * t62, t49 * t63, t49 ^ 2 / 0.2e1, pkin(1) * t66 + t38 * t49, -pkin(1) * t56 * t71 - t39 * t49 (-t38 * t56 + t39 * t58) * t70, t39 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t62, t37 ^ 2 / 0.2e1, -t37 * t35, -t37 * t44, t35 ^ 2 / 0.2e1, t35 * t44, t44 ^ 2 / 0.2e1, -t19 * t44 + t31 * t35, t20 * t44 + t31 * t37, -t19 * t37 - t20 * t35, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t67, -t74, -t73, t68, t72, t40, t27 * t23 - t8 * t41, t27 * t25 + t9 * t41, -t9 * t23 - t8 * t25, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t40, t73, -t72, t67, -t74, t68, t7 * t23 + t6 * t25, -t10 * t23 - t6 * t41, -t10 * t25 + t7 * t41, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t22, t16 ^ 2 / 0.2e1, -t16 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t4 * t16, t4 * t18 - t2 * t22, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
