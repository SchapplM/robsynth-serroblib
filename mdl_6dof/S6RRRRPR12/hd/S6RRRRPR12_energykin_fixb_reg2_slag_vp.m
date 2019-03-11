% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:44:38
% EndTime: 2019-03-09 23:44:39
% DurationCPUTime: 0.42s
% Computational Cost: add. (2712->80), mult. (7202->183), div. (0->0), fcn. (6009->14), ass. (0->64)
t82 = cos(pkin(6)) * qJD(1);
t60 = qJD(2) + t82;
t63 = sin(pkin(7));
t71 = cos(qJ(2));
t85 = cos(pkin(7));
t76 = t71 * t85;
t64 = sin(pkin(6));
t83 = qJD(1) * t64;
t91 = t60 * t63 + t76 * t83;
t69 = sin(qJ(2));
t78 = t71 * t83;
t80 = pkin(1) * t82;
t53 = pkin(9) * t78 + t69 * t80;
t40 = t91 * pkin(10) + t53;
t49 = (-pkin(10) * t63 * t69 - pkin(2) * t71 - pkin(1)) * t83;
t68 = sin(qJ(3));
t70 = cos(qJ(3));
t59 = t71 * t80;
t79 = t69 * t83;
t42 = t60 * pkin(2) + t59 + (-t85 * pkin(10) - pkin(9)) * t79;
t75 = t85 * t42;
t28 = -t68 * t40 + t70 * (t49 * t63 + t75);
t30 = -t63 * t42 + t85 * t49;
t43 = t68 * t79 - t91 * t70;
t86 = t63 * t68;
t45 = t60 * t86 + (t68 * t76 + t69 * t70) * t83;
t20 = t43 * pkin(3) - t45 * pkin(11) + t30;
t29 = t70 * t40 + t49 * t86 + t68 * t75;
t50 = -t85 * t60 + t63 * t78 - qJD(3);
t23 = -t50 * pkin(11) + t29;
t67 = sin(qJ(4));
t90 = cos(qJ(4));
t13 = t67 * t20 + t90 * t23;
t32 = t67 * t45 + t90 * t50;
t11 = -t32 * qJ(5) + t13;
t62 = sin(pkin(13));
t84 = cos(pkin(13));
t12 = t90 * t20 - t67 * t23;
t34 = t90 * t45 - t67 * t50;
t41 = qJD(4) + t43;
t9 = t41 * pkin(4) - t34 * qJ(5) + t12;
t6 = t84 * t11 + t62 * t9;
t89 = cos(qJ(6));
t72 = qJD(1) ^ 2;
t87 = t64 ^ 2 * t72;
t81 = t71 * t87;
t77 = t87 / 0.2e1;
t25 = t84 * t32 + t62 * t34;
t5 = -t62 * t11 + t84 * t9;
t22 = t50 * pkin(3) - t28;
t14 = t32 * pkin(4) + qJD(5) + t22;
t66 = sin(qJ(6));
t52 = -pkin(9) * t79 + t59;
t39 = t41 ^ 2 / 0.2e1;
t27 = -t62 * t32 + t84 * t34;
t24 = qJD(6) + t25;
t17 = t89 * t27 + t66 * t41;
t15 = t66 * t27 - t89 * t41;
t7 = t25 * pkin(5) - t27 * pkin(12) + t14;
t4 = t41 * pkin(12) + t6;
t3 = -t41 * pkin(5) - t5;
t2 = t89 * t4 + t66 * t7;
t1 = -t66 * t4 + t89 * t7;
t8 = [0, 0, 0, 0, 0, t72 / 0.2e1, 0, 0, 0, 0, t69 ^ 2 * t77, t69 * t81, t60 * t79, t71 ^ 2 * t77, t60 * t78, t60 ^ 2 / 0.2e1, pkin(1) * t81 + t52 * t60, -pkin(1) * t69 * t87 - t53 * t60 (-t52 * t69 + t53 * t71) * t83, t53 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t77, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t50, t43 ^ 2 / 0.2e1, t43 * t50, t50 ^ 2 / 0.2e1, -t28 * t50 + t30 * t43, t29 * t50 + t30 * t45, -t28 * t45 - t29 * t43, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t41, t32 ^ 2 / 0.2e1, -t32 * t41, t39, t12 * t41 + t22 * t32, -t13 * t41 + t22 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t41, t25 ^ 2 / 0.2e1, -t25 * t41, t39, t14 * t25 + t5 * t41, t14 * t27 - t6 * t41, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t24, t15 ^ 2 / 0.2e1, -t15 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t3 * t15, t3 * t17 - t2 * t24, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
