% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t105 = cos(qJ(3));
t107 = cos(qJ(1));
t103 = sin(qJ(2));
t104 = sin(qJ(1));
t106 = cos(qJ(2));
t98 = cos(pkin(6));
t87 = t98 * t106;
t67 = t104 * t103 - t107 * t87;
t95 = sin(pkin(7));
t96 = sin(pkin(6));
t80 = t96 * t95;
t97 = cos(pkin(7));
t117 = t107 * t80 + t67 * t97;
t86 = t98 * t103;
t41 = t104 * t106 + t107 * t86;
t56 = sin(qJ(3));
t16 = t117 * t105 + t41 * t56;
t52 = pkin(13) + qJ(6);
t50 = sin(t52);
t51 = cos(t52);
t81 = t97 * t96;
t111 = -t107 * t81 + t67 * t95;
t19 = -t41 * t105 + t117 * t56;
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t7 = -t111 * t55 + t19 * t57;
t121 = t16 * t51 + t7 * t50;
t120 = -t16 * t50 + t7 * t51;
t6 = t111 * t57 + t19 * t55;
t63 = t107 * t103 + t104 * t87;
t113 = -t104 * t80 + t63 * t97;
t112 = t106 * t81 + t98 * t95;
t42 = -t104 * t86 + t107 * t106;
t20 = t113 * t105 + t42 * t56;
t83 = t96 * t103;
t31 = -t112 * t105 + t56 * t83;
t75 = g(1) * t20 + g(2) * t16 + g(3) * t31;
t102 = t50 * t57;
t101 = t51 * t57;
t53 = sin(pkin(13));
t100 = t53 * t57;
t54 = cos(pkin(13));
t99 = t54 * t57;
t94 = pkin(9) * t96;
t93 = t95 * pkin(10);
t92 = t55 * t95;
t91 = t56 * t97;
t90 = t57 * t95;
t21 = t42 * t105 - t113 * t56;
t58 = -t104 * t81 - t63 * t95;
t8 = t21 * t55 + t58 * t57;
t89 = g(1) * t6 + g(2) * t8;
t85 = t97 * t105;
t84 = t106 * t96;
t32 = t105 * t83 + t112 * t56;
t62 = -t106 * t80 + t98 * t97;
t14 = t32 * t55 - t62 * t57;
t78 = g(1) * t8 - g(2) * t6 + g(3) * t14;
t15 = t32 * t57 + t62 * t55;
t9 = t21 * t57 - t58 * t55;
t77 = g(1) * t9 - g(2) * t7 + g(3) * t15;
t23 = -t67 * t105 - t41 * t91;
t10 = t23 * t55 - t41 * t90;
t25 = -t63 * t105 - t42 * t91;
t12 = t25 * t55 - t42 * t90;
t73 = t103 * t81;
t39 = t105 * t84 - t56 * t73;
t69 = t103 * t80;
t26 = t39 * t55 - t57 * t69;
t76 = g(1) * t12 + g(2) * t10 + g(3) * t26;
t74 = g(1) * t21 - g(2) * t19 + g(3) * t32;
t38 = t105 * t73 + t56 * t84;
t27 = t39 * t57 + t55 * t69;
t24 = t42 * t85 - t63 * t56;
t22 = t41 * t85 - t67 * t56;
t13 = t25 * t57 + t42 * t92;
t11 = t23 * t57 + t41 * t92;
t3 = t75 * t55;
t2 = t20 * t50 + t9 * t51;
t1 = t20 * t51 - t9 * t50;
t4 = [0, g(1) * t104 - g(2) * t107, g(1) * t107 + g(2) * t104, 0, 0, 0, 0, 0, g(1) * t41 - g(2) * t42, -g(1) * t67 + g(2) * t63, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t21, -g(1) * t16 + g(2) * t20, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, t89, -g(1) * (-t16 * t53 + t7 * t54) - g(2) * (t20 * t53 + t9 * t54) -g(1) * (-t16 * t54 - t7 * t53) - g(2) * (t20 * t54 - t9 * t53) -t89, -g(1) * (-t104 * pkin(1) - t41 * pkin(2) + t19 * pkin(3) + t7 * pkin(4) - pkin(11) * t16 + t6 * qJ(5) + t107 * t94) - g(2) * (t107 * pkin(1) + t42 * pkin(2) + t21 * pkin(3) + t9 * pkin(4) + t20 * pkin(11) + t8 * qJ(5) + t104 * t94) + (g(1) * t111 + g(2) * t58) * pkin(10), 0, 0, 0, 0, 0, -g(1) * t120 - g(2) * t2, g(1) * t121 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t63 + g(2) * t67 - g(3) * t84, g(1) * t42 + g(2) * t41 + g(3) * t83, 0, 0, 0, 0, 0, -g(1) * t25 - g(2) * t23 - g(3) * t39, g(1) * t24 + g(2) * t22 + g(3) * t38, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t27, t76, -g(1) * (t13 * t54 + t24 * t53) - g(2) * (t11 * t54 + t22 * t53) - g(3) * (t27 * t54 + t38 * t53) -g(1) * (-t13 * t53 + t24 * t54) - g(2) * (-t11 * t53 + t22 * t54) - g(3) * (-t27 * t53 + t38 * t54) -t76, -g(1) * (-t63 * pkin(2) + t25 * pkin(3) + t13 * pkin(4) + t24 * pkin(11) + t12 * qJ(5) + t42 * t93) - g(2) * (-t67 * pkin(2) + t23 * pkin(3) + t11 * pkin(4) + t22 * pkin(11) + t10 * qJ(5) + t41 * t93) - g(3) * (pkin(2) * t84 + t39 * pkin(3) + t27 * pkin(4) + pkin(10) * t69 + t38 * pkin(11) + t26 * qJ(5)) 0, 0, 0, 0, 0, -g(1) * (t13 * t51 + t24 * t50) - g(2) * (t11 * t51 + t22 * t50) - g(3) * (t27 * t51 + t38 * t50) -g(1) * (-t13 * t50 + t24 * t51) - g(2) * (-t11 * t50 + t22 * t51) - g(3) * (-t27 * t50 + t38 * t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, 0, 0, 0, 0, 0, t75 * t57, -t3, -g(1) * (-t20 * t99 + t21 * t53) - g(2) * (-t16 * t99 - t19 * t53) - g(3) * (-t31 * t99 + t32 * t53) -g(1) * (t20 * t100 + t21 * t54) - g(2) * (t16 * t100 - t19 * t54) - g(3) * (t31 * t100 + t32 * t54) t3, -t74 * pkin(11) + t75 * (pkin(4) * t57 + qJ(5) * t55 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (-t20 * t101 + t21 * t50) - g(2) * (-t16 * t101 - t19 * t50) - g(3) * (-t31 * t101 + t32 * t50) -g(1) * (t20 * t102 + t21 * t51) - g(2) * (t16 * t102 - t19 * t51) - g(3) * (t31 * t102 + t32 * t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t77, t78 * t54, -t78 * t53, -t77, -g(1) * (-t8 * pkin(4) + t9 * qJ(5)) - g(2) * (pkin(4) * t6 - qJ(5) * t7) - g(3) * (-t14 * pkin(4) + t15 * qJ(5)) 0, 0, 0, 0, 0, t78 * t51, -t78 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t121 - g(3) * (-t15 * t50 + t31 * t51) g(1) * t2 - g(2) * t120 - g(3) * (-t15 * t51 - t31 * t50);];
taug_reg  = t4;
