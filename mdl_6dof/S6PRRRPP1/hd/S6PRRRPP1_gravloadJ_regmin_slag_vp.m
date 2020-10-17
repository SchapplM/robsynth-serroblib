% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:37:38
% EndTime: 2019-05-05 06:37:40
% DurationCPUTime: 0.51s
% Computational Cost: add. (447->117), mult. (967->186), div. (0->0), fcn. (1187->12), ass. (0->75)
t56 = -qJ(5) - pkin(9);
t58 = sin(qJ(3));
t101 = -t56 * t58 + pkin(2);
t60 = cos(qJ(4));
t51 = t60 * pkin(4) + pkin(3);
t61 = cos(qJ(3));
t100 = -t51 * t61 - t101;
t55 = sin(pkin(6));
t99 = g(3) * t55;
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t78 = cos(pkin(10));
t79 = cos(pkin(6));
t72 = t79 * t78;
t77 = sin(pkin(10));
t37 = t59 * t72 + t77 * t62;
t75 = t55 * t78;
t19 = t37 * t61 - t58 * t75;
t57 = sin(qJ(4));
t98 = t19 * t57;
t71 = t79 * t77;
t39 = -t59 * t71 + t78 * t62;
t74 = t55 * t77;
t21 = t39 * t61 + t58 * t74;
t97 = t21 * t57;
t36 = t77 * t59 - t62 * t72;
t96 = t36 * t60;
t95 = t37 * t57;
t38 = t78 * t59 + t62 * t71;
t94 = t38 * t60;
t93 = t39 * t57;
t54 = qJ(4) + pkin(11);
t52 = sin(t54);
t91 = t52 * t61;
t53 = cos(t54);
t90 = t53 * t61;
t89 = t55 * t59;
t88 = t55 * t62;
t86 = t57 * t59;
t85 = t57 * t61;
t84 = t60 * t61;
t83 = t61 * t62;
t18 = t37 * t58 + t61 * t75;
t82 = -t18 * t51 - t19 * t56;
t20 = t39 * t58 - t61 * t74;
t81 = -t20 * t51 - t21 * t56;
t40 = t58 * t89 - t79 * t61;
t41 = t79 * t58 + t61 * t89;
t80 = -t40 * t51 - t41 * t56;
t76 = t60 * t88;
t73 = -pkin(5) * t53 - qJ(6) * t52;
t70 = -t41 * t57 - t76;
t14 = t41 * t52 + t53 * t88;
t3 = t19 * t52 - t36 * t53;
t5 = t21 * t52 - t38 * t53;
t69 = g(1) * t5 + g(2) * t3 + g(3) * t14;
t68 = g(1) * t20 + g(2) * t18 + g(3) * t40;
t67 = g(1) * t21 + g(2) * t19 + g(3) * t41;
t66 = -g(1) * t38 - g(2) * t36 + g(3) * t88;
t65 = pkin(8) * t89 + t101 * t88 + (pkin(4) * t86 + t51 * t83) * t55;
t64 = pkin(4) * t95 + t37 * pkin(8) + t100 * t36;
t63 = pkin(4) * t93 + t39 * pkin(8) + t100 * t38;
t31 = pkin(4) * t94;
t29 = pkin(4) * t96;
t23 = (t52 * t59 + t53 * t83) * t55;
t22 = (t52 * t83 - t53 * t59) * t55;
t15 = t41 * t53 - t52 * t88;
t11 = -t38 * t90 + t39 * t52;
t10 = -t38 * t91 - t39 * t53;
t9 = -t36 * t90 + t37 * t52;
t8 = -t36 * t91 - t37 * t53;
t7 = t66 * t58;
t6 = t21 * t53 + t38 * t52;
t4 = t19 * t53 + t36 * t52;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, -t66, g(1) * t39 + g(2) * t37 + g(3) * t89, 0, 0, 0, 0, 0, -t66 * t61, t7, 0, 0, 0, 0, 0, -g(1) * (-t38 * t84 + t93) - g(2) * (-t36 * t84 + t95) - (t60 * t83 + t86) * t99, -g(1) * (t38 * t85 + t39 * t60) - g(2) * (t36 * t85 + t37 * t60) - (-t57 * t83 + t59 * t60) * t99, -t7, -g(1) * t63 - g(2) * t64 - g(3) * t65, -g(1) * t11 - g(2) * t9 - g(3) * t23, -t7, -g(1) * t10 - g(2) * t8 - g(3) * t22, -g(1) * (t11 * pkin(5) + t10 * qJ(6) + t63) - g(2) * (t9 * pkin(5) + t8 * qJ(6) + t64) - g(3) * (t23 * pkin(5) + t22 * qJ(6) + t65); 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t67, 0, 0, 0, 0, 0, t68 * t60, -t68 * t57, -t67, -g(1) * t81 - g(2) * t82 - g(3) * t80, t68 * t53, -t67, t68 * t52, -g(1) * (t73 * t20 + t81) - g(2) * (t73 * t18 + t82) - g(3) * (t73 * t40 + t80); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t94 - t97) - g(2) * (t96 - t98) - g(3) * t70, -g(1) * (-t21 * t60 - t38 * t57) - g(2) * (-t19 * t60 - t36 * t57) - g(3) * (-t41 * t60 + t57 * t88) 0, -g(1) * t31 - g(2) * t29 + (g(3) * t76 + t67 * t57) * pkin(4), t69, 0, -g(1) * t6 - g(2) * t4 - g(3) * t15, -g(1) * (-pkin(4) * t97 - t5 * pkin(5) + t6 * qJ(6) + t31) - g(2) * (-pkin(4) * t98 - t3 * pkin(5) + t4 * qJ(6) + t29) - g(3) * (t70 * pkin(4) - t14 * pkin(5) + t15 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, 0, 0, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69;];
taug_reg  = t1;
