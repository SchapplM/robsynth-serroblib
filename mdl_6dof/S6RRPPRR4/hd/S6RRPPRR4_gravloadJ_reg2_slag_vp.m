% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t52 = cos(pkin(6));
t55 = sin(qJ(2));
t59 = cos(qJ(2));
t89 = sin(pkin(11));
t90 = cos(pkin(11));
t67 = -t55 * t89 + t59 * t90;
t34 = t67 * t52;
t42 = -t55 * t90 - t59 * t89;
t56 = sin(qJ(1));
t60 = cos(qJ(1));
t24 = -t56 * t34 + t60 * t42;
t93 = t60 * t55;
t95 = t56 * t59;
t38 = -t52 * t95 - t93;
t117 = t38 * pkin(2) + t24 * pkin(3);
t115 = t42 * t52;
t25 = t115 * t56 + t60 * t67;
t20 = t115 * t60 - t56 * t67;
t51 = sin(pkin(6));
t101 = t51 * t59;
t114 = -g(1) * t38 - g(3) * t101;
t33 = t42 * t51;
t70 = -g(1) * t25 + g(2) * t20 + g(3) * t33;
t53 = sin(qJ(6));
t57 = cos(qJ(6));
t100 = t51 * t60;
t21 = t60 * t34 + t56 * t42;
t54 = sin(qJ(5));
t58 = cos(qJ(5));
t74 = t58 * t100 + t21 * t54;
t113 = -t20 * t57 + t53 * t74;
t112 = t20 * t53 + t57 * t74;
t107 = t21 * pkin(9);
t106 = t24 * pkin(9);
t32 = t67 * t51;
t105 = t32 * pkin(9);
t102 = t51 * t56;
t99 = t53 * t54;
t98 = t54 * t57;
t96 = t56 * t55;
t92 = t60 * t59;
t91 = pkin(2) * t101 + t32 * pkin(3);
t86 = t52 * t92;
t35 = t52 * t55 * pkin(2) + (-pkin(8) - qJ(3)) * t51;
t50 = t59 * pkin(2) + pkin(1);
t85 = -t56 * t35 + t60 * t50;
t7 = t54 * t102 + t24 * t58;
t73 = t54 * t100 - t21 * t58;
t82 = g(1) * t73 + g(2) * t7;
t81 = -t33 * qJ(4) + t91;
t80 = g(1) * t21 - g(2) * t24;
t6 = -g(1) * t20 - g(2) * t25;
t79 = g(1) * t60 + g(2) * t56;
t78 = g(1) * t56 - g(2) * t60;
t77 = -t60 * t35 - t56 * t50;
t43 = pkin(2) * t86;
t76 = -pkin(2) * t96 + t21 * pkin(3) + t43;
t26 = -t32 * t58 - t52 * t54;
t72 = g(1) * t7 - g(2) * t73 - g(3) * t26;
t27 = -t32 * t54 + t52 * t58;
t8 = t58 * t102 - t24 * t54;
t71 = g(1) * t8 - g(2) * t74 + g(3) * t27;
t3 = g(1) * t24 + g(2) * t21 + g(3) * t32;
t69 = t25 * pkin(3) - t24 * qJ(4) + t85;
t68 = -t20 * qJ(4) + t76;
t66 = pkin(3) * t20 + t21 * qJ(4) + t77;
t63 = pkin(4) * t102 + t25 * pkin(9) + t69;
t62 = qJ(4) * t25 + t117;
t61 = pkin(4) * t100 + pkin(9) * t20 + t66;
t40 = t79 * t51;
t39 = -t52 * t96 + t92;
t37 = -t52 * t93 - t95;
t36 = -t86 + t96;
t31 = -g(3) * t52 - t78 * t51;
t5 = t25 * t53 + t8 * t57;
t4 = t25 * t57 - t8 * t53;
t1 = t70 * t58;
t2 = [0, 0, 0, 0, 0, 0, t78, t79, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t37 - g(2) * t39, -g(1) * t36 - g(2) * t38, -t40, -g(1) * (-t56 * pkin(1) + pkin(8) * t100) - g(2) * (t60 * pkin(1) + pkin(8) * t102) 0, 0, 0, 0, 0, 0, t6, t80, -t40, -g(1) * t77 - g(2) * t85, 0, 0, 0, 0, 0, 0, -t40, -t6, -t80, -g(1) * t66 - g(2) * t69, 0, 0, 0, 0, 0, 0, -g(1) * t74 - g(2) * t8, t82, t6, -g(1) * t61 - g(2) * t63, 0, 0, 0, 0, 0, 0, -g(1) * t112 - g(2) * t5, g(1) * t113 - g(2) * t4, -t82, -g(1) * (pkin(5) * t74 + pkin(10) * t73 + t61) - g(2) * (t8 * pkin(5) + t7 * pkin(10) + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t36 + t114, g(3) * t51 * t55 + g(1) * t39 - g(2) * t37, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t70, 0, -g(2) * t43 + (g(2) * t96 + t114) * pkin(2), 0, 0, 0, 0, 0, 0, 0, t3, t70, -g(1) * t62 - g(2) * t68 - g(3) * t81, 0, 0, 0, 0, 0, 0, t70 * t54, t1, -t3, -g(1) * (t62 + t106) - g(2) * (t68 + t107) - g(3) * (t81 + t105) 0, 0, 0, 0, 0, 0, -g(1) * (t24 * t53 + t25 * t98) - g(2) * (-t20 * t98 + t21 * t53) - g(3) * (t32 * t53 - t33 * t98) -g(1) * (t24 * t57 - t25 * t99) - g(2) * (t20 * t99 + t21 * t57) - g(3) * (t32 * t57 + t33 * t99) -t1, -g(1) * (t106 + t117) - g(2) * (t76 + t107) - g(3) * (t91 + t105) + t70 * (pkin(5) * t54 - pkin(10) * t58 + qJ(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t71, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t57, -t72 * t53, -t71, -g(1) * (-t7 * pkin(5) + t8 * pkin(10)) - g(2) * (pkin(5) * t73 - pkin(10) * t74) - g(3) * (t26 * pkin(5) + t27 * pkin(10)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t113 - g(3) * (-t27 * t53 - t33 * t57) g(1) * t5 - g(2) * t112 - g(3) * (-t27 * t57 + t33 * t53) 0, 0;];
taug_reg  = t2;
