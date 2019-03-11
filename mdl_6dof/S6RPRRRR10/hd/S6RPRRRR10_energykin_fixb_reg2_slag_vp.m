% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:34:07
% EndTime: 2019-03-09 07:34:08
% DurationCPUTime: 0.37s
% Computational Cost: add. (2246->82), mult. (7201->197), div. (0->0), fcn. (6011->14), ass. (0->66)
t62 = sin(pkin(13));
t65 = cos(pkin(13));
t64 = sin(pkin(6));
t81 = qJD(1) * t64;
t75 = qJ(2) * t81;
t67 = cos(pkin(6));
t80 = qJD(1) * t67;
t78 = pkin(1) * t80;
t53 = t62 * t78 + t65 * t75;
t63 = sin(pkin(7));
t66 = cos(pkin(7));
t83 = t64 * t65;
t39 = (t63 * t67 + t66 * t83) * qJD(1) * pkin(9) + t53;
t58 = t65 * t78;
t85 = t62 * t64;
t42 = t58 + (pkin(2) * t67 + (-pkin(9) * t66 - qJ(2)) * t85) * qJD(1);
t48 = qJD(2) + (-pkin(9) * t62 * t63 - pkin(2) * t65 - pkin(1)) * t81;
t71 = sin(qJ(3));
t72 = cos(qJ(3));
t28 = -t71 * t39 + (t42 * t66 + t48 * t63) * t72;
t73 = qJD(1) ^ 2;
t90 = t73 / 0.2e1;
t30 = -t63 * t42 + t66 * t48;
t77 = t65 * t81;
t43 = t71 * t62 * t81 + (-t63 * t80 - t66 * t77) * t72;
t82 = t66 * t71;
t84 = t63 * t71;
t45 = (t67 * t84 + (t62 * t72 + t65 * t82) * t64) * qJD(1);
t20 = t43 * pkin(3) - t45 * pkin(10) + t30;
t29 = t72 * t39 + t42 * t82 + t48 * t84;
t50 = t63 * t77 - t66 * t80 - qJD(3);
t23 = -t50 * pkin(10) + t29;
t70 = sin(qJ(4));
t89 = cos(qJ(4));
t13 = t70 * t20 + t89 * t23;
t32 = t70 * t45 + t89 * t50;
t11 = -t32 * pkin(11) + t13;
t69 = sin(qJ(5));
t88 = cos(qJ(5));
t12 = t89 * t20 - t70 * t23;
t34 = t89 * t45 - t70 * t50;
t41 = qJD(4) + t43;
t9 = t41 * pkin(4) - t34 * pkin(11) + t12;
t6 = t88 * t11 + t69 * t9;
t87 = cos(qJ(6));
t86 = t64 ^ 2 * t73;
t79 = t64 * t67 * t73;
t76 = t86 / 0.2e1;
t25 = t88 * t32 + t69 * t34;
t5 = -t69 * t11 + t88 * t9;
t22 = t50 * pkin(3) - t28;
t14 = t32 * pkin(4) + t22;
t68 = sin(qJ(6));
t59 = -pkin(1) * t81 + qJD(2);
t52 = -t62 * t75 + t58;
t40 = qJD(5) + t41;
t27 = -t69 * t32 + t88 * t34;
t24 = qJD(6) + t25;
t17 = t87 * t27 + t68 * t40;
t15 = t68 * t27 - t87 * t40;
t7 = t25 * pkin(5) - t27 * pkin(12) + t14;
t4 = t40 * pkin(12) + t6;
t3 = -t40 * pkin(5) - t5;
t2 = t87 * t4 + t68 * t7;
t1 = -t68 * t4 + t87 * t7;
t8 = [0, 0, 0, 0, 0, t90, 0, 0, 0, 0, t62 ^ 2 * t76, t62 * t65 * t86, t62 * t79, t65 ^ 2 * t76, t65 * t79, t67 ^ 2 * t90 (t52 * t67 - t59 * t83) * qJD(1) (-t53 * t67 + t59 * t85) * qJD(1) (-t52 * t62 + t53 * t65) * t81, t53 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t50, t43 ^ 2 / 0.2e1, t43 * t50, t50 ^ 2 / 0.2e1, -t28 * t50 + t30 * t43, t29 * t50 + t30 * t45, -t28 * t45 - t29 * t43, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t41, t32 ^ 2 / 0.2e1, -t32 * t41, t41 ^ 2 / 0.2e1, t12 * t41 + t22 * t32, -t13 * t41 + t22 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t40, t25 ^ 2 / 0.2e1, -t25 * t40, t40 ^ 2 / 0.2e1, t14 * t25 + t5 * t40, t14 * t27 - t6 * t40, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t24, t15 ^ 2 / 0.2e1, -t15 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t3 * t15, t3 * t17 - t2 * t24, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
