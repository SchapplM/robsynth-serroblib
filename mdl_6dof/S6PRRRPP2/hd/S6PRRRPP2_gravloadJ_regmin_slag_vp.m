% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t69 = sin(qJ(4));
t110 = qJ(5) * t69 + pkin(3);
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t109 = -pkin(3) * t73 - pkin(9) * t70 - pkin(2);
t67 = sin(pkin(10));
t71 = sin(qJ(2));
t74 = cos(qJ(2));
t93 = cos(pkin(10));
t94 = cos(pkin(6));
t82 = t94 * t93;
t53 = t67 * t74 + t71 * t82;
t68 = sin(pkin(6));
t86 = t68 * t93;
t30 = -t53 * t70 - t73 * t86;
t72 = cos(qJ(4));
t106 = t30 * t72;
t102 = t68 * t73;
t87 = t67 * t94;
t55 = -t71 * t87 + t93 * t74;
t32 = t67 * t102 - t55 * t70;
t105 = t32 * t72;
t103 = t68 * t71;
t56 = -t70 * t103 + t94 * t73;
t104 = t56 * t72;
t101 = t68 * t74;
t100 = t69 * t73;
t99 = t72 * t73;
t98 = t73 * t74;
t97 = pkin(9) - qJ(6);
t95 = qJ(6) * t70;
t92 = t70 * t101;
t91 = t69 * t101;
t90 = pkin(4) * t106 + t110 * t30;
t89 = pkin(4) * t105 + t110 * t32;
t88 = pkin(4) * t104 + t110 * t56;
t31 = t53 * t73 - t70 * t86;
t52 = t67 * t71 - t74 * t82;
t12 = t31 * t69 - t52 * t72;
t13 = t31 * t72 + t52 * t69;
t85 = -t12 * pkin(4) + t13 * qJ(5);
t33 = t67 * t68 * t70 + t55 * t73;
t54 = t93 * t71 + t74 * t87;
t14 = t33 * t69 - t54 * t72;
t15 = t33 * t72 + t54 * t69;
t84 = -t14 * pkin(4) + t15 * qJ(5);
t57 = t71 * t102 + t94 * t70;
t34 = t72 * t101 + t57 * t69;
t35 = t57 * t72 - t91;
t83 = -t34 * pkin(4) + t35 * qJ(5);
t2 = g(1) * t14 + g(2) * t12 + g(3) * t34;
t81 = g(1) * t15 + g(2) * t13 + g(3) * t35;
t19 = -t52 * t100 - t53 * t72;
t21 = -t54 * t100 - t55 * t72;
t37 = -t72 * t103 + t73 * t91;
t80 = g(1) * t21 + g(2) * t19 + g(3) * t37;
t79 = g(1) * t32 + g(2) * t30 + g(3) * t56;
t9 = g(1) * t33 + g(2) * t31 + g(3) * t57;
t38 = (t69 * t71 + t72 * t98) * t68;
t78 = t68 * pkin(3) * t98 + pkin(2) * t101 + t38 * pkin(4) + pkin(8) * t103 + pkin(9) * t92 + t37 * qJ(5);
t77 = -g(1) * t54 - g(2) * t52 + g(3) * t101;
t20 = -t52 * t99 + t53 * t69;
t76 = t20 * pkin(4) + t53 * pkin(8) + t19 * qJ(5) + t109 * t52;
t22 = -t54 * t99 + t55 * t69;
t75 = t22 * pkin(4) + t55 * pkin(8) + t21 * qJ(5) + t109 * t54;
t16 = t77 * t70;
t7 = t79 * t72;
t6 = t79 * t69;
t5 = -g(1) * t22 - g(2) * t20 - g(3) * t38;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, -t77, g(1) * t55 + g(2) * t53 + g(3) * t103, 0, 0, 0, 0, 0, -t77 * t73, t16, 0, 0, 0, 0, 0, t5, t80, t5, -t16, -t80, -g(1) * t75 - g(2) * t76 - g(3) * t78, t5, -t80, t16, -g(1) * (t22 * pkin(5) + t54 * t95 + t75) - g(2) * (t20 * pkin(5) + t52 * t95 + t76) - g(3) * (t38 * pkin(5) - qJ(6) * t92 + t78); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t9, 0, 0, 0, 0, 0, -t7, t6, -t7, -t9, -t6, -g(1) * (t33 * pkin(9) + t89) - g(2) * (t31 * pkin(9) + t90) - g(3) * (t57 * pkin(9) + t88) -t7, -t6, t9, -g(1) * (pkin(5) * t105 + t97 * t33 + t89) - g(2) * (pkin(5) * t106 + t97 * t31 + t90) - g(3) * (pkin(5) * t104 + t97 * t57 + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t81, t2, 0, -t81, -g(1) * t84 - g(2) * t85 - g(3) * t83, t2, -t81, 0, -g(1) * (-t14 * pkin(5) + t84) - g(2) * (-t12 * pkin(5) + t85) - g(3) * (-t34 * pkin(5) + t83); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79;];
taug_reg  = t1;
