% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:38:26
% EndTime: 2019-05-06 16:38:30
% DurationCPUTime: 0.90s
% Computational Cost: add. (510->145), mult. (1183->217), div. (0->0), fcn. (1429->12), ass. (0->78)
t50 = sin(qJ(2));
t51 = sin(qJ(1));
t53 = cos(qJ(2));
t54 = cos(qJ(1));
t83 = cos(pkin(6));
t76 = t54 * t83;
t26 = t50 * t76 + t51 * t53;
t77 = t51 * t83;
t28 = -t50 * t77 + t54 * t53;
t107 = -g(1) * t28 - g(2) * t26;
t44 = pkin(11) + qJ(6);
t41 = sin(t44);
t42 = cos(t44);
t25 = t51 * t50 - t53 * t76;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t46 = sin(pkin(6));
t88 = t46 * t54;
t66 = -t25 * t49 + t52 * t88;
t106 = t26 * t42 + t41 * t66;
t105 = -t26 * t41 + t42 * t66;
t102 = g(3) * t46;
t101 = t25 * pkin(9);
t27 = t54 * t50 + t53 * t77;
t100 = t27 * pkin(9);
t45 = sin(pkin(11));
t99 = t25 * t45;
t96 = t27 * t45;
t95 = t41 * t49;
t94 = t42 * t49;
t93 = t45 * t49;
t92 = t45 * t53;
t91 = t46 * t50;
t90 = t46 * t51;
t89 = t46 * t53;
t47 = cos(pkin(11));
t87 = t47 * t49;
t86 = t49 * t50;
t85 = pkin(2) * t89 + qJ(3) * t91;
t84 = t54 * pkin(1) + pkin(8) * t90;
t82 = pkin(9) * t89 + t85;
t81 = pkin(5) * t45 + pkin(9);
t80 = -t51 * pkin(1) + pkin(8) * t88;
t19 = t25 * pkin(2);
t79 = -t19 - t101;
t21 = t27 * pkin(2);
t78 = -t21 - t100;
t75 = t26 * qJ(3) - t19;
t74 = t28 * qJ(3) - t21;
t73 = g(3) * t82;
t10 = -t27 * t52 + t49 * t90;
t65 = t25 * t52 + t49 * t88;
t72 = g(1) * t65 + g(2) * t10;
t71 = g(1) * t25 - g(2) * t27;
t9 = g(1) * t26 - g(2) * t28;
t70 = g(1) * t54 + g(2) * t51;
t69 = pkin(4) * t49 - qJ(5) * t52;
t40 = t47 * pkin(5) + pkin(4);
t48 = -pkin(10) - qJ(5);
t68 = t40 * t49 + t48 * t52;
t67 = t28 * pkin(2) + t27 * qJ(3) + t84;
t62 = pkin(3) * t90 + t67;
t61 = -t26 * pkin(2) - t25 * qJ(3) + t80;
t23 = t83 * t49 + t52 * t89;
t60 = g(1) * t10 - g(2) * t65 + g(3) * t23;
t11 = t27 * t49 + t52 * t90;
t24 = -t49 * t89 + t83 * t52;
t59 = g(1) * t11 - g(2) * t66 + g(3) * t24;
t58 = pkin(3) * t88 + t61;
t7 = -g(1) * t27 - g(2) * t25 + g(3) * t89;
t57 = g(3) * t91 - t107;
t56 = t28 * pkin(9) + t62;
t55 = -t26 * pkin(9) + t58;
t29 = t70 * t46;
t6 = t57 * t52;
t5 = t11 * t42 + t28 * t41;
t4 = -t11 * t41 + t28 * t42;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t51 - g(2) * t54, t70, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t71, -t29, -g(1) * t80 - g(2) * t84, 0, 0, 0, 0, 0, 0, -t29, -t9, t71, -g(1) * t61 - g(2) * t67, 0, 0, 0, 0, 0, 0, -g(1) * t66 - g(2) * t11, t72, t9, -g(1) * t55 - g(2) * t56, 0, 0, 0, 0, 0, 0, -g(1) * (-t26 * t45 + t47 * t66) - g(2) * (t11 * t47 + t28 * t45) -g(1) * (-t26 * t47 - t45 * t66) - g(2) * (-t11 * t45 + t28 * t47) -t72, -g(1) * (pkin(4) * t66 + qJ(5) * t65 + t55) - g(2) * (t11 * pkin(4) + t10 * qJ(5) + t56) 0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t5, g(1) * t106 - g(2) * t4, -t72, -g(1) * (-t81 * t26 + t40 * t66 - t48 * t65 + t58) - g(2) * (-t10 * t48 + t11 * t40 + t81 * t28 + t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t57, -g(1) * t74 - g(2) * t75 - g(3) * t85, 0, 0, 0, 0, 0, 0, -t57 * t49, -t6, -t7, -g(1) * (t74 - t100) - g(2) * (t75 - t101) - t73, 0, 0, 0, 0, 0, 0, -g(1) * (t28 * t87 - t96) - g(2) * (t26 * t87 - t99) - (t47 * t86 + t92) * t102, -g(1) * (-t27 * t47 - t28 * t93) - g(2) * (-t25 * t47 - t26 * t93) - (-t45 * t86 + t47 * t53) * t102, t6, -g(1) * t78 - g(2) * t79 - g(3) * (t69 * t91 + t82) + t107 * (qJ(3) + t69) 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t41 + t28 * t94) - g(2) * (-t25 * t41 + t26 * t94) - (t41 * t53 + t42 * t86) * t102, -g(1) * (-t27 * t42 - t28 * t95) - g(2) * (-t25 * t42 - t26 * t95) - (-t41 * t86 + t42 * t53) * t102, t6, -g(1) * (-pkin(5) * t96 + t78) - g(2) * (-pkin(5) * t99 + t79) - t73 - (pkin(5) * t92 + t68 * t50) * t102 + t107 * (qJ(3) + t68); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t59, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t47, -t60 * t45, -t59, -g(1) * (-t10 * pkin(4) + t11 * qJ(5)) - g(2) * (pkin(4) * t65 - qJ(5) * t66) - g(3) * (-t23 * pkin(4) + t24 * qJ(5)) 0, 0, 0, 0, 0, 0, t60 * t42, -t60 * t41, -t59, -g(1) * (-t10 * t40 - t11 * t48) - g(2) * (t40 * t65 + t48 * t66) - g(3) * (-t23 * t40 - t24 * t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t106 - g(3) * (-t24 * t41 + t42 * t91) g(1) * t5 - g(2) * t105 - g(3) * (-t24 * t42 - t41 * t91) 0, 0;];
taug_reg  = t1;
