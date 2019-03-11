% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPR11
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:44:30
% EndTime: 2019-03-09 05:44:30
% DurationCPUTime: 0.37s
% Computational Cost: add. (2213->82), mult. (7110->194), div. (0->0), fcn. (5932->14), ass. (0->66)
t63 = sin(pkin(12));
t66 = cos(pkin(12));
t65 = sin(pkin(6));
t82 = qJD(1) * t65;
t75 = qJ(2) * t82;
t68 = cos(pkin(6));
t81 = qJD(1) * t68;
t78 = pkin(1) * t81;
t53 = t63 * t78 + t66 * t75;
t64 = sin(pkin(7));
t67 = cos(pkin(7));
t85 = t65 * t66;
t40 = (t64 * t68 + t67 * t85) * qJD(1) * pkin(9) + t53;
t58 = t66 * t78;
t87 = t63 * t65;
t42 = t58 + (pkin(2) * t68 + (-pkin(9) * t67 - qJ(2)) * t87) * qJD(1);
t48 = qJD(2) + (-pkin(9) * t63 * t64 - pkin(2) * t66 - pkin(1)) * t82;
t71 = sin(qJ(3));
t72 = cos(qJ(3));
t25 = -t71 * t40 + (t42 * t67 + t48 * t64) * t72;
t73 = qJD(1) ^ 2;
t91 = t73 / 0.2e1;
t31 = -t64 * t42 + t67 * t48;
t77 = t66 * t82;
t43 = t71 * t63 * t82 + (-t64 * t81 - t67 * t77) * t72;
t84 = t67 * t71;
t86 = t64 * t71;
t45 = (t68 * t86 + (t63 * t72 + t66 * t84) * t65) * qJD(1);
t20 = t43 * pkin(3) - t45 * pkin(10) + t31;
t26 = t72 * t40 + t42 * t84 + t48 * t86;
t50 = t64 * t77 - t67 * t81 - qJD(3);
t24 = -t50 * pkin(10) + t26;
t70 = sin(qJ(4));
t90 = cos(qJ(4));
t12 = t70 * t20 + t90 * t24;
t41 = qJD(4) + t43;
t10 = t41 * qJ(5) + t12;
t23 = t50 * pkin(3) - t25;
t33 = t70 * t45 + t90 * t50;
t35 = t90 * t45 - t70 * t50;
t15 = t33 * pkin(4) - t35 * qJ(5) + t23;
t62 = sin(pkin(13));
t83 = cos(pkin(13));
t6 = t83 * t10 + t62 * t15;
t89 = cos(qJ(6));
t88 = t65 ^ 2 * t73;
t80 = t33 ^ 2 / 0.2e1;
t79 = t65 * t68 * t73;
t76 = t88 / 0.2e1;
t5 = -t62 * t10 + t83 * t15;
t11 = t90 * t20 - t70 * t24;
t9 = -t41 * pkin(4) + qJD(5) - t11;
t69 = sin(qJ(6));
t59 = -pkin(1) * t82 + qJD(2);
t52 = -t63 * t75 + t58;
t32 = qJD(6) + t33;
t30 = t83 * t35 + t62 * t41;
t28 = t62 * t35 - t83 * t41;
t18 = -t69 * t28 + t89 * t30;
t16 = t89 * t28 + t69 * t30;
t7 = t28 * pkin(5) + t9;
t4 = -t28 * pkin(11) + t6;
t3 = t33 * pkin(5) - t30 * pkin(11) + t5;
t2 = t69 * t3 + t89 * t4;
t1 = t89 * t3 - t69 * t4;
t8 = [0, 0, 0, 0, 0, t91, 0, 0, 0, 0, t63 ^ 2 * t76, t63 * t66 * t88, t63 * t79, t66 ^ 2 * t76, t66 * t79, t68 ^ 2 * t91 (t52 * t68 - t59 * t85) * qJD(1) (-t53 * t68 + t59 * t87) * qJD(1) (-t52 * t63 + t53 * t66) * t82, t53 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t50, t43 ^ 2 / 0.2e1, t43 * t50, t50 ^ 2 / 0.2e1, -t25 * t50 + t31 * t43, t26 * t50 + t31 * t45, -t25 * t45 - t26 * t43, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * t41, t80, -t33 * t41, t41 ^ 2 / 0.2e1, t11 * t41 + t23 * t33, -t12 * t41 + t23 * t35, -t11 * t35 - t12 * t33, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * t33, t28 ^ 2 / 0.2e1, -t28 * t33, t80, t9 * t28 + t5 * t33, t9 * t30 - t6 * t33, -t6 * t28 - t5 * t30, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t32, t16 ^ 2 / 0.2e1, -t16 * t32, t32 ^ 2 / 0.2e1, t1 * t32 + t7 * t16, t7 * t18 - t2 * t32, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
