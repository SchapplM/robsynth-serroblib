% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR15_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t57 = sin(qJ(2));
t60 = cos(qJ(1));
t89 = cos(pkin(6));
t81 = t60 * t89;
t97 = sin(qJ(1));
t99 = cos(qJ(2));
t43 = t97 * t57 - t99 * t81;
t44 = t57 * t81 + t97 * t99;
t56 = sin(qJ(3));
t88 = cos(pkin(7));
t98 = cos(qJ(3));
t74 = t88 * t98;
t52 = sin(pkin(7));
t53 = sin(pkin(6));
t93 = t53 * t60;
t86 = t52 * t93;
t18 = t43 * t74 + t44 * t56 + t98 * t86;
t54 = sin(qJ(6));
t58 = cos(qJ(6));
t82 = t56 * t88;
t19 = -t43 * t82 + t44 * t98 - t56 * t86;
t55 = sin(qJ(4));
t59 = cos(qJ(4));
t83 = t53 * t88;
t65 = t43 * t52 - t60 * t83;
t6 = t19 * t55 - t65 * t59;
t110 = t18 * t58 + t6 * t54;
t109 = -t18 * t54 + t6 * t58;
t7 = t19 * t59 + t65 * t55;
t75 = t89 * t97;
t64 = t60 * t57 + t99 * t75;
t84 = t53 * t97;
t104 = -t52 * t84 + t64 * t88;
t45 = -t57 * t75 + t60 * t99;
t22 = t104 * t98 + t45 * t56;
t80 = t89 * t52;
t85 = t53 * t99;
t94 = t53 * t57;
t33 = t56 * t94 - t74 * t85 - t98 * t80;
t68 = g(1) * t22 + g(2) * t18 + g(3) * t33;
t103 = pkin(10) * t52;
t96 = t52 * t55;
t95 = t52 * t59;
t92 = t54 * t55;
t91 = t55 * t58;
t90 = t57 * t52;
t87 = t53 * t90;
t23 = -t104 * t56 + t45 * t98;
t61 = -t64 * t52 - t97 * t83;
t10 = t23 * t55 + t61 * t59;
t78 = -g(1) * t6 + g(2) * t10;
t11 = t23 * t59 - t61 * t55;
t77 = -g(1) * t7 + g(2) * t11;
t76 = -g(1) * t18 + g(2) * t22;
t34 = t56 * t80 + (t98 * t57 + t99 * t82) * t53;
t63 = -t52 * t85 + t89 * t88;
t16 = t34 * t55 - t63 * t59;
t72 = g(1) * t10 + g(2) * t6 + g(3) * t16;
t17 = t34 * t59 + t63 * t55;
t71 = g(1) * t11 + g(2) * t7 + g(3) * t17;
t25 = -t43 * t98 - t44 * t82;
t12 = t25 * t55 - t44 * t95;
t27 = -t45 * t82 - t64 * t98;
t14 = t27 * t55 - t45 * t95;
t41 = (-t57 * t82 + t98 * t99) * t53;
t28 = t41 * t55 - t59 * t87;
t70 = g(1) * t14 + g(2) * t12 + g(3) * t28;
t13 = t25 * t59 + t44 * t96;
t15 = t27 * t59 + t45 * t96;
t29 = t41 * t59 + t55 * t87;
t69 = g(1) * t15 + g(2) * t13 + g(3) * t29;
t67 = g(1) * t23 + g(2) * t19 + g(3) * t34;
t24 = -t43 * t56 + t44 * t74;
t26 = t45 * t74 - t64 * t56;
t40 = (t99 * t56 + t57 * t74) * t53;
t66 = g(1) * t26 + g(2) * t24 + g(3) * t40;
t5 = t10 * t54 + t22 * t58;
t4 = t10 * t58 - t22 * t54;
t3 = t68 * t59;
t2 = t68 * t55;
t1 = [0, g(1) * t97 - g(2) * t60, g(1) * t60 + g(2) * t97, 0, 0, 0, 0, 0, g(1) * t44 - g(2) * t45, -g(1) * t43 + g(2) * t64, 0, 0, 0, 0, 0, g(1) * t19 - g(2) * t23, t76, 0, 0, 0, 0, 0, -t77, t78, -t76, t77, -t78, -g(1) * (-t97 * pkin(1) - t44 * pkin(2) - pkin(3) * t19 - pkin(4) * t7 + pkin(9) * t93 - pkin(11) * t18 - qJ(5) * t6) - g(2) * (t60 * pkin(1) + t45 * pkin(2) + t23 * pkin(3) + t11 * pkin(4) + pkin(9) * t84 + t22 * pkin(11) + t10 * qJ(5)) + (g(1) * t65 + g(2) * t61) * pkin(10), 0, 0, 0, 0, 0, g(1) * t110 - g(2) * t5, g(1) * t109 - g(2) * t4; 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t64 + g(2) * t43 - g(3) * t85, g(1) * t45 + g(2) * t44 + g(3) * t94, 0, 0, 0, 0, 0, -g(1) * t27 - g(2) * t25 - g(3) * t41, t66, 0, 0, 0, 0, 0, -t69, t70, -t66, t69, -t70, -g(1) * (-t64 * pkin(2) + t27 * pkin(3) + t15 * pkin(4) + t26 * pkin(11) + t14 * qJ(5) + t45 * t103) - g(2) * (-t43 * pkin(2) + t25 * pkin(3) + t13 * pkin(4) + t24 * pkin(11) + t12 * qJ(5) + t44 * t103) - g(3) * (t41 * pkin(3) + t29 * pkin(4) + t40 * pkin(11) + t28 * qJ(5) + (t99 * pkin(2) + pkin(10) * t90) * t53) 0, 0, 0, 0, 0, -g(1) * (t14 * t54 + t26 * t58) - g(2) * (t12 * t54 + t24 * t58) - g(3) * (t28 * t54 + t40 * t58) -g(1) * (t14 * t58 - t26 * t54) - g(2) * (t12 * t58 - t24 * t54) - g(3) * (t28 * t58 - t40 * t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t67, 0, 0, 0, 0, 0, t3, -t2, -t67, -t3, t2, -t67 * pkin(11) + t68 * (pkin(4) * t59 + qJ(5) * t55 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (-t22 * t92 + t23 * t58) - g(2) * (-t18 * t92 + t19 * t58) - g(3) * (-t33 * t92 + t34 * t58) -g(1) * (-t22 * t91 - t23 * t54) - g(2) * (-t18 * t91 - t19 * t54) - g(3) * (-t33 * t91 - t34 * t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t71, 0, -t72, -t71, -g(1) * (-t10 * pkin(4) + t11 * qJ(5)) - g(2) * (-t6 * pkin(4) + t7 * qJ(5)) - g(3) * (-t16 * pkin(4) + t17 * qJ(5)) 0, 0, 0, 0, 0, -t71 * t54, -t71 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t109 - g(3) * (t16 * t58 - t33 * t54) g(1) * t5 + g(2) * t110 - g(3) * (-t16 * t54 - t33 * t58);];
taug_reg  = t1;
