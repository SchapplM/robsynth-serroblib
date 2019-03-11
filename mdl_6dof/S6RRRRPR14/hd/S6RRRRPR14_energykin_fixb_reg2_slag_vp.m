% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR14
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR14_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:25:28
% EndTime: 2019-03-10 00:25:28
% DurationCPUTime: 0.42s
% Computational Cost: add. (2686->80), mult. (7111->183), div. (0->0), fcn. (5930->14), ass. (0->63)
t82 = cos(pkin(6)) * qJD(1);
t60 = qJD(2) + t82;
t63 = sin(pkin(7));
t65 = cos(pkin(7));
t72 = cos(qJ(2));
t64 = sin(pkin(6));
t83 = qJD(1) * t64;
t77 = t72 * t83;
t91 = t60 * t63 + t65 * t77;
t70 = sin(qJ(2));
t79 = pkin(1) * t82;
t53 = pkin(9) * t77 + t70 * t79;
t40 = t91 * pkin(10) + t53;
t59 = t72 * t79;
t78 = t70 * t83;
t42 = t60 * pkin(2) + t59 + (-pkin(10) * t65 - pkin(9)) * t78;
t49 = (-pkin(10) * t63 * t70 - pkin(2) * t72 - pkin(1)) * t83;
t69 = sin(qJ(3));
t71 = cos(qJ(3));
t25 = -t69 * t40 + (t42 * t65 + t49 * t63) * t71;
t31 = -t63 * t42 + t65 * t49;
t43 = t69 * t78 - t91 * t71;
t85 = t65 * t69;
t86 = t63 * t69;
t45 = t60 * t86 + (t70 * t71 + t72 * t85) * t83;
t20 = t43 * pkin(3) - t45 * pkin(11) + t31;
t26 = t71 * t40 + t42 * t85 + t49 * t86;
t50 = -t65 * t60 + t63 * t77 - qJD(3);
t24 = -t50 * pkin(11) + t26;
t68 = sin(qJ(4));
t90 = cos(qJ(4));
t12 = t68 * t20 + t90 * t24;
t41 = qJD(4) + t43;
t10 = t41 * qJ(5) + t12;
t23 = t50 * pkin(3) - t25;
t33 = t68 * t45 + t90 * t50;
t35 = t90 * t45 - t68 * t50;
t15 = t33 * pkin(4) - t35 * qJ(5) + t23;
t62 = sin(pkin(13));
t84 = cos(pkin(13));
t6 = t84 * t10 + t62 * t15;
t89 = cos(qJ(6));
t73 = qJD(1) ^ 2;
t87 = t64 ^ 2 * t73;
t81 = t33 ^ 2 / 0.2e1;
t80 = t72 * t87;
t76 = t87 / 0.2e1;
t5 = -t62 * t10 + t84 * t15;
t11 = t90 * t20 - t68 * t24;
t9 = -t41 * pkin(4) + qJD(5) - t11;
t67 = sin(qJ(6));
t52 = -pkin(9) * t78 + t59;
t32 = qJD(6) + t33;
t30 = t84 * t35 + t62 * t41;
t28 = t62 * t35 - t84 * t41;
t18 = -t67 * t28 + t89 * t30;
t16 = t89 * t28 + t67 * t30;
t7 = t28 * pkin(5) + t9;
t4 = -t28 * pkin(12) + t6;
t3 = t33 * pkin(5) - t30 * pkin(12) + t5;
t2 = t67 * t3 + t89 * t4;
t1 = t89 * t3 - t67 * t4;
t8 = [0, 0, 0, 0, 0, t73 / 0.2e1, 0, 0, 0, 0, t70 ^ 2 * t76, t70 * t80, t60 * t78, t72 ^ 2 * t76, t60 * t77, t60 ^ 2 / 0.2e1, pkin(1) * t80 + t52 * t60, -pkin(1) * t70 * t87 - t53 * t60 (-t52 * t70 + t53 * t72) * t83, t53 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t76, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t50, t43 ^ 2 / 0.2e1, t43 * t50, t50 ^ 2 / 0.2e1, -t25 * t50 + t31 * t43, t26 * t50 + t31 * t45, -t25 * t45 - t26 * t43, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * t41, t81, -t33 * t41, t41 ^ 2 / 0.2e1, t11 * t41 + t23 * t33, -t12 * t41 + t23 * t35, -t11 * t35 - t12 * t33, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t33 * t30, t28 ^ 2 / 0.2e1, -t33 * t28, t81, t9 * t28 + t5 * t33, t9 * t30 - t6 * t33, -t6 * t28 - t5 * t30, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t32, t16 ^ 2 / 0.2e1, -t16 * t32, t32 ^ 2 / 0.2e1, t1 * t32 + t7 * t16, t7 * t18 - t2 * t32, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
