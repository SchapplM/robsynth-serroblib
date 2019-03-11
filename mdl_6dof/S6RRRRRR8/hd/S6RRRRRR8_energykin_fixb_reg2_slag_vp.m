% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:06:42
% EndTime: 2019-03-10 05:06:43
% DurationCPUTime: 0.42s
% Computational Cost: add. (2727->80), mult. (7202->186), div. (0->0), fcn. (6009->14), ass. (0->63)
t81 = cos(pkin(6)) * qJD(1);
t60 = qJD(2) + t81;
t62 = sin(pkin(7));
t64 = cos(pkin(7));
t72 = cos(qJ(2));
t63 = sin(pkin(6));
t82 = qJD(1) * t63;
t77 = t72 * t82;
t90 = t60 * t62 + t64 * t77;
t70 = sin(qJ(2));
t79 = pkin(1) * t81;
t53 = pkin(9) * t77 + t70 * t79;
t39 = t90 * pkin(10) + t53;
t59 = t72 * t79;
t78 = t70 * t82;
t42 = t60 * pkin(2) + t59 + (-pkin(10) * t64 - pkin(9)) * t78;
t49 = (-pkin(10) * t62 * t70 - pkin(2) * t72 - pkin(1)) * t82;
t69 = sin(qJ(3));
t71 = cos(qJ(3));
t28 = -t69 * t39 + (t42 * t64 + t49 * t62) * t71;
t30 = -t62 * t42 + t64 * t49;
t43 = t69 * t78 - t90 * t71;
t83 = t64 * t69;
t84 = t62 * t69;
t45 = t60 * t84 + (t70 * t71 + t72 * t83) * t82;
t20 = t43 * pkin(3) - t45 * pkin(11) + t30;
t29 = t71 * t39 + t42 * t83 + t49 * t84;
t50 = -t64 * t60 + t62 * t77 - qJD(3);
t23 = -t50 * pkin(11) + t29;
t68 = sin(qJ(4));
t89 = cos(qJ(4));
t13 = t68 * t20 + t89 * t23;
t32 = t68 * t45 + t89 * t50;
t11 = -t32 * pkin(12) + t13;
t67 = sin(qJ(5));
t88 = cos(qJ(5));
t12 = t89 * t20 - t68 * t23;
t34 = t89 * t45 - t68 * t50;
t41 = qJD(4) + t43;
t9 = t41 * pkin(4) - t34 * pkin(12) + t12;
t6 = t88 * t11 + t67 * t9;
t87 = cos(qJ(6));
t73 = qJD(1) ^ 2;
t85 = t63 ^ 2 * t73;
t80 = t72 * t85;
t76 = t85 / 0.2e1;
t25 = t88 * t32 + t67 * t34;
t5 = -t67 * t11 + t88 * t9;
t22 = t50 * pkin(3) - t28;
t14 = t32 * pkin(4) + t22;
t66 = sin(qJ(6));
t52 = -pkin(9) * t78 + t59;
t40 = qJD(5) + t41;
t27 = -t67 * t32 + t88 * t34;
t24 = qJD(6) + t25;
t17 = t87 * t27 + t66 * t40;
t15 = t66 * t27 - t87 * t40;
t7 = t25 * pkin(5) - t27 * pkin(13) + t14;
t4 = t40 * pkin(13) + t6;
t3 = -t40 * pkin(5) - t5;
t2 = t87 * t4 + t66 * t7;
t1 = -t66 * t4 + t87 * t7;
t8 = [0, 0, 0, 0, 0, t73 / 0.2e1, 0, 0, 0, 0, t70 ^ 2 * t76, t70 * t80, t60 * t78, t72 ^ 2 * t76, t60 * t77, t60 ^ 2 / 0.2e1, pkin(1) * t80 + t52 * t60, -pkin(1) * t70 * t85 - t53 * t60 (-t52 * t70 + t53 * t72) * t82, t53 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t76, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t50, t43 ^ 2 / 0.2e1, t43 * t50, t50 ^ 2 / 0.2e1, -t28 * t50 + t30 * t43, t29 * t50 + t30 * t45, -t28 * t45 - t29 * t43, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t41, t32 ^ 2 / 0.2e1, -t32 * t41, t41 ^ 2 / 0.2e1, t12 * t41 + t22 * t32, -t13 * t41 + t22 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t40, t25 ^ 2 / 0.2e1, -t25 * t40, t40 ^ 2 / 0.2e1, t14 * t25 + t5 * t40, t14 * t27 - t6 * t40, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t24, t15 ^ 2 / 0.2e1, -t15 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t3 * t15, t3 * t17 - t2 * t24, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
