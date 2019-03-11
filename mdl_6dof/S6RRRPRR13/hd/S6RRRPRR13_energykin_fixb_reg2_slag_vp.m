% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR13_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:05:07
% EndTime: 2019-03-09 20:05:07
% DurationCPUTime: 0.42s
% Computational Cost: add. (2695->80), mult. (7202->183), div. (0->0), fcn. (6009->14), ass. (0->64)
t82 = cos(pkin(6)) * qJD(1);
t59 = qJD(2) + t82;
t62 = sin(pkin(7));
t70 = cos(qJ(2));
t85 = cos(pkin(7));
t75 = t70 * t85;
t63 = sin(pkin(6));
t83 = qJD(1) * t63;
t91 = t59 * t62 + t75 * t83;
t68 = sin(qJ(2));
t77 = t70 * t83;
t79 = pkin(1) * t82;
t52 = pkin(9) * t77 + t68 * t79;
t39 = t91 * pkin(10) + t52;
t48 = (-pkin(10) * t62 * t68 - pkin(2) * t70 - pkin(1)) * t83;
t67 = sin(qJ(3));
t69 = cos(qJ(3));
t58 = t70 * t79;
t78 = t68 * t83;
t41 = t59 * pkin(2) + t58 + (-t85 * pkin(10) - pkin(9)) * t78;
t74 = t85 * t41;
t28 = -t67 * t39 + t69 * (t48 * t62 + t74);
t30 = -t62 * t41 + t85 * t48;
t42 = t67 * t78 - t91 * t69;
t86 = t62 * t67;
t44 = t59 * t86 + (t67 * t75 + t68 * t69) * t83;
t20 = t42 * pkin(3) - t44 * qJ(4) + t30;
t29 = t69 * t39 + t48 * t86 + t67 * t74;
t49 = -t85 * t59 + t62 * t77 - qJD(3);
t23 = -t49 * qJ(4) + t29;
t61 = sin(pkin(13));
t84 = cos(pkin(13));
t13 = t61 * t20 + t84 * t23;
t32 = t61 * t44 + t84 * t49;
t11 = -t32 * pkin(11) + t13;
t66 = sin(qJ(5));
t12 = t84 * t20 - t61 * t23;
t34 = t84 * t44 - t61 * t49;
t9 = t42 * pkin(4) - t34 * pkin(11) + t12;
t90 = cos(qJ(5));
t6 = t90 * t11 + t66 * t9;
t89 = cos(qJ(6));
t71 = qJD(1) ^ 2;
t87 = t63 ^ 2 * t71;
t81 = t42 ^ 2 / 0.2e1;
t80 = t70 * t87;
t76 = t87 / 0.2e1;
t25 = t90 * t32 + t66 * t34;
t5 = -t66 * t11 + t90 * t9;
t22 = t49 * pkin(3) + qJD(4) - t28;
t14 = t32 * pkin(4) + t22;
t65 = sin(qJ(6));
t51 = -pkin(9) * t78 + t58;
t40 = qJD(5) + t42;
t27 = -t66 * t32 + t90 * t34;
t24 = qJD(6) + t25;
t17 = t89 * t27 + t65 * t40;
t15 = t65 * t27 - t89 * t40;
t7 = t25 * pkin(5) - t27 * pkin(12) + t14;
t4 = t40 * pkin(12) + t6;
t3 = -t40 * pkin(5) - t5;
t2 = t89 * t4 + t65 * t7;
t1 = -t65 * t4 + t89 * t7;
t8 = [0, 0, 0, 0, 0, t71 / 0.2e1, 0, 0, 0, 0, t68 ^ 2 * t76, t68 * t80, t59 * t78, t70 ^ 2 * t76, t59 * t77, t59 ^ 2 / 0.2e1, pkin(1) * t80 + t51 * t59, -pkin(1) * t68 * t87 - t52 * t59 (-t51 * t68 + t52 * t70) * t83, t52 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t76, t44 ^ 2 / 0.2e1, -t44 * t42, -t44 * t49, t81, t42 * t49, t49 ^ 2 / 0.2e1, -t28 * t49 + t30 * t42, t29 * t49 + t30 * t44, -t28 * t44 - t29 * t42, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t42, t32 ^ 2 / 0.2e1, -t32 * t42, t81, t12 * t42 + t22 * t32, -t13 * t42 + t22 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t40, t25 ^ 2 / 0.2e1, -t25 * t40, t40 ^ 2 / 0.2e1, t14 * t25 + t5 * t40, t14 * t27 - t6 * t40, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t24, t15 ^ 2 / 0.2e1, -t15 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t3 * t15, t3 * t17 - t2 * t24, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
