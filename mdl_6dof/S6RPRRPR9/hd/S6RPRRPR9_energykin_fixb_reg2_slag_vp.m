% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:31:32
% EndTime: 2019-03-09 05:31:32
% DurationCPUTime: 0.44s
% Computational Cost: add. (2231->81), mult. (7201->192), div. (0->0), fcn. (6011->14), ass. (0->66)
t64 = sin(pkin(7));
t65 = sin(pkin(6));
t67 = cos(pkin(6));
t66 = cos(pkin(12));
t85 = cos(pkin(7));
t76 = t66 * t85;
t94 = qJD(1) * (t64 * t67 + t65 * t76);
t63 = sin(pkin(12));
t83 = qJD(1) * t65;
t77 = qJ(2) * t83;
t82 = qJD(1) * t67;
t79 = pkin(1) * t82;
t53 = t63 * t79 + t66 * t77;
t39 = pkin(9) * t94 + t53;
t48 = qJD(2) + (-pkin(9) * t63 * t64 - pkin(2) * t66 - pkin(1)) * t83;
t70 = sin(qJ(3));
t71 = cos(qJ(3));
t58 = t66 * t79;
t88 = t63 * t65;
t42 = t58 + (pkin(2) * t67 + (-t85 * pkin(9) - qJ(2)) * t88) * qJD(1);
t75 = t85 * t42;
t28 = -t70 * t39 + (t48 * t64 + t75) * t71;
t72 = qJD(1) ^ 2;
t92 = t72 / 0.2e1;
t30 = -t64 * t42 + t85 * t48;
t43 = t70 * t63 * t83 - t71 * t94;
t86 = t64 * t70;
t45 = (t67 * t86 + (t63 * t71 + t70 * t76) * t65) * qJD(1);
t20 = t43 * pkin(3) - t45 * pkin(10) + t30;
t29 = t71 * t39 + t48 * t86 + t70 * t75;
t50 = t64 * t66 * t83 - t85 * t82 - qJD(3);
t23 = -t50 * pkin(10) + t29;
t69 = sin(qJ(4));
t91 = cos(qJ(4));
t13 = t69 * t20 + t91 * t23;
t32 = t69 * t45 + t91 * t50;
t11 = -t32 * qJ(5) + t13;
t62 = sin(pkin(13));
t84 = cos(pkin(13));
t12 = t91 * t20 - t69 * t23;
t34 = t91 * t45 - t69 * t50;
t41 = qJD(4) + t43;
t9 = t41 * pkin(4) - t34 * qJ(5) + t12;
t6 = t84 * t11 + t62 * t9;
t90 = cos(qJ(6));
t89 = t65 ^ 2 * t72;
t80 = t65 * t67 * t72;
t78 = t89 / 0.2e1;
t25 = t84 * t32 + t62 * t34;
t5 = -t62 * t11 + t84 * t9;
t22 = t50 * pkin(3) - t28;
t14 = t32 * pkin(4) + qJD(5) + t22;
t68 = sin(qJ(6));
t59 = -pkin(1) * t83 + qJD(2);
t52 = -t63 * t77 + t58;
t40 = t41 ^ 2 / 0.2e1;
t27 = -t62 * t32 + t84 * t34;
t24 = qJD(6) + t25;
t17 = t90 * t27 + t68 * t41;
t15 = t68 * t27 - t90 * t41;
t7 = t25 * pkin(5) - t27 * pkin(11) + t14;
t4 = t41 * pkin(11) + t6;
t3 = -t41 * pkin(5) - t5;
t2 = t90 * t4 + t68 * t7;
t1 = -t68 * t4 + t90 * t7;
t8 = [0, 0, 0, 0, 0, t92, 0, 0, 0, 0, t63 ^ 2 * t78, t63 * t66 * t89, t63 * t80, t66 ^ 2 * t78, t66 * t80, t67 ^ 2 * t92 (-t59 * t65 * t66 + t52 * t67) * qJD(1) (-t53 * t67 + t59 * t88) * qJD(1) (-t52 * t63 + t53 * t66) * t83, t53 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t50, t43 ^ 2 / 0.2e1, t43 * t50, t50 ^ 2 / 0.2e1, -t28 * t50 + t30 * t43, t29 * t50 + t30 * t45, -t28 * t45 - t29 * t43, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t41, t32 ^ 2 / 0.2e1, -t32 * t41, t40, t12 * t41 + t22 * t32, -t13 * t41 + t22 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t41, t25 ^ 2 / 0.2e1, -t25 * t41, t40, t14 * t25 + t5 * t41, t14 * t27 - t6 * t41, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t24, t15 ^ 2 / 0.2e1, -t15 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t3 * t15, t3 * t17 - t2 * t24, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
