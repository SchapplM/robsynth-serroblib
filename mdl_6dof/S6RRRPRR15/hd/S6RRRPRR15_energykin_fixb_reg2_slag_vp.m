% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR15_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:38:38
% EndTime: 2019-03-09 20:38:38
% DurationCPUTime: 0.31s
% Computational Cost: add. (1615->78), mult. (4364->168), div. (0->0), fcn. (3549->12), ass. (0->67)
t79 = cos(pkin(6)) * qJD(1);
t53 = qJD(2) + t79;
t60 = sin(qJ(3));
t61 = sin(qJ(2));
t62 = cos(qJ(2));
t81 = cos(pkin(7));
t68 = t62 * t81;
t56 = sin(pkin(6));
t80 = qJD(1) * t56;
t55 = sin(pkin(7));
t82 = t55 * t60;
t89 = cos(qJ(3));
t36 = (t60 * t68 + t89 * t61) * t80 + t53 * t82;
t90 = pkin(3) + pkin(11);
t67 = t68 * t80;
t70 = t55 * t89;
t73 = t61 * t80;
t34 = -t53 * t70 + t60 * t73 - t89 * t67;
t74 = pkin(1) * t79;
t52 = t62 * t74;
t33 = t53 * pkin(2) + t52 + (-t81 * pkin(10) - pkin(9)) * t73;
t40 = (-pkin(10) * t55 * t61 - pkin(2) * t62 - pkin(1)) * t80;
t21 = -t55 * t33 + t81 * t40;
t66 = -t36 * qJ(4) + t21;
t11 = t90 * t34 + t66;
t59 = sin(qJ(5));
t88 = cos(qJ(5));
t72 = t62 * t80;
t43 = -t81 * t53 + t55 * t72 - qJD(3);
t46 = pkin(9) * t72 + t61 * t74;
t30 = (t53 * t55 + t67) * pkin(10) + t46;
t69 = t33 * t81;
t16 = -t60 * t30 + t40 * t70 + t89 * t69;
t64 = qJD(4) - t16;
t9 = t36 * pkin(4) + t90 * t43 + t64;
t6 = t88 * t11 + t59 * t9;
t87 = cos(qJ(6));
t86 = t36 * t34;
t85 = t36 * t43;
t84 = t43 * t34;
t63 = qJD(1) ^ 2;
t83 = t56 ^ 2 * t63;
t78 = t34 ^ 2 / 0.2e1;
t77 = t36 ^ 2 / 0.2e1;
t76 = t62 * t83;
t17 = t89 * t30 + t40 * t82 + t60 * t69;
t71 = t83 / 0.2e1;
t23 = -t88 * t34 - t59 * t43;
t15 = t43 * qJ(4) - t17;
t12 = -t34 * pkin(4) - t15;
t5 = -t59 * t11 + t88 * t9;
t58 = sin(qJ(6));
t45 = -pkin(9) * t73 + t52;
t41 = t43 ^ 2 / 0.2e1;
t32 = qJD(5) + t36;
t25 = t59 * t34 - t88 * t43;
t22 = qJD(6) + t23;
t20 = t87 * t25 + t58 * t32;
t18 = t58 * t25 - t87 * t32;
t14 = t43 * pkin(3) + t64;
t13 = t34 * pkin(3) + t66;
t7 = t23 * pkin(5) - t25 * pkin(12) + t12;
t4 = t32 * pkin(12) + t6;
t3 = -t32 * pkin(5) - t5;
t2 = t87 * t4 + t58 * t7;
t1 = -t58 * t4 + t87 * t7;
t8 = [0, 0, 0, 0, 0, t63 / 0.2e1, 0, 0, 0, 0, t61 ^ 2 * t71, t61 * t76, t53 * t73, t62 ^ 2 * t71, t53 * t72, t53 ^ 2 / 0.2e1, pkin(1) * t76 + t45 * t53, -pkin(1) * t61 * t83 - t46 * t53 (-t45 * t61 + t46 * t62) * t80, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t71, t77, -t86, -t85, t78, t84, t41, -t16 * t43 + t21 * t34, t17 * t43 + t21 * t36, -t16 * t36 - t17 * t34, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t41, t85, -t84, t77, -t86, t78, t14 * t36 + t15 * t34, -t13 * t34 - t14 * t43, -t13 * t36 + t15 * t43, t13 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t32, t23 ^ 2 / 0.2e1, -t23 * t32, t32 ^ 2 / 0.2e1, t12 * t23 + t5 * t32, t12 * t25 - t6 * t32, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t22, t18 ^ 2 / 0.2e1, -t18 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t3 * t18, -t2 * t22 + t3 * t20, -t1 * t20 - t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
