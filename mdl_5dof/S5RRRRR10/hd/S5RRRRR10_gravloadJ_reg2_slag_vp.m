% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:35:58
% EndTime: 2019-12-31 22:36:00
% DurationCPUTime: 0.67s
% Computational Cost: add. (519->126), mult. (953->202), div. (0->0), fcn. (1143->12), ass. (0->77)
t57 = sin(qJ(2));
t58 = sin(qJ(1));
t61 = cos(qJ(2));
t62 = cos(qJ(1));
t85 = cos(pkin(5));
t73 = t62 * t85;
t33 = t57 * t73 + t58 * t61;
t53 = qJ(3) + qJ(4);
t50 = sin(t53);
t51 = cos(t53);
t54 = sin(pkin(5));
t92 = t54 * t62;
t14 = t33 * t51 - t50 * t92;
t32 = t58 * t57 - t61 * t73;
t55 = sin(qJ(5));
t59 = cos(qJ(5));
t104 = t14 * t55 - t32 * t59;
t103 = t14 * t59 + t32 * t55;
t102 = g(3) * t54;
t74 = t58 * t85;
t35 = -t57 * t74 + t62 * t61;
t56 = sin(qJ(3));
t99 = t35 * t56;
t98 = t51 * t55;
t97 = t51 * t59;
t96 = t54 * t57;
t95 = t54 * t58;
t60 = cos(qJ(3));
t94 = t54 * t60;
t93 = t54 * t61;
t91 = t55 * t61;
t63 = -pkin(9) - pkin(8);
t90 = t57 * t63;
t89 = t59 * t61;
t49 = t60 * pkin(3) + pkin(2);
t88 = -t32 * t49 - t33 * t63;
t34 = t62 * t57 + t61 * t74;
t87 = -t34 * t49 - t35 * t63;
t86 = t62 * pkin(1) + pkin(7) * t95;
t84 = t56 * t96;
t83 = t56 * t95;
t82 = t58 * t94;
t44 = t56 * t92;
t81 = t60 * t92;
t80 = -t58 * pkin(1) + pkin(7) * t92;
t76 = -t33 * t50 - t51 * t92;
t79 = pkin(4) * t76 + t14 * pkin(10);
t17 = t35 * t50 - t51 * t95;
t18 = t35 * t51 + t50 * t95;
t78 = -t17 * pkin(4) + t18 * pkin(10);
t26 = -t50 * t96 + t85 * t51;
t27 = t85 * t50 + t51 * t96;
t77 = t26 * pkin(4) + t27 * pkin(10);
t75 = t33 * t60 - t44;
t72 = t85 * t60;
t71 = pkin(3) * t83 - t34 * t63 + t35 * t49 + t86;
t70 = pkin(4) * t51 + pkin(10) * t50;
t69 = g(1) * t76 + g(2) * t17;
t10 = g(1) * t32 - g(2) * t34;
t68 = g(1) * t62 + g(2) * t58;
t67 = t33 * t56 + t81;
t66 = pkin(3) * t44 + t32 * t63 - t33 * t49 + t80;
t3 = g(1) * t17 - g(2) * t76 - g(3) * t26;
t5 = g(1) * t18 + g(2) * t14 + g(3) * t27;
t65 = -g(1) * t34 - g(2) * t32 + g(3) * t93;
t64 = g(1) * t35 + g(2) * t33 + g(3) * t96;
t48 = pkin(3) * t72;
t41 = pkin(3) * t82;
t36 = t49 * t93;
t20 = t35 * t60 + t83;
t19 = t82 - t99;
t8 = t18 * t59 + t34 * t55;
t7 = -t18 * t55 + t34 * t59;
t6 = t65 * t50;
t2 = t3 * t59;
t1 = t3 * t55;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t58 - g(2) * t62, t68, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t33 - g(2) * t35, -t10, -t68 * t54, -g(1) * t80 - g(2) * t86, 0, 0, 0, 0, 0, 0, g(1) * t75 - g(2) * t20, -g(1) * t67 - g(2) * t19, t10, -g(1) * (-t33 * pkin(2) - t32 * pkin(8) + t80) - g(2) * (t35 * pkin(2) + t34 * pkin(8) + t86), 0, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t18, t69, t10, -g(1) * t66 - g(2) * t71, 0, 0, 0, 0, 0, 0, g(1) * t103 - g(2) * t8, -g(1) * t104 - g(2) * t7, -t69, -g(1) * (-pkin(4) * t14 + pkin(10) * t76 + t66) - g(2) * (t18 * pkin(4) + t17 * pkin(10) + t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t64, 0, 0, 0, 0, 0, 0, 0, 0, -t65 * t60, t65 * t56, -t64, -g(1) * (-t34 * pkin(2) + t35 * pkin(8)) - g(2) * (-t32 * pkin(2) + t33 * pkin(8)) - (pkin(2) * t61 + pkin(8) * t57) * t102, 0, 0, 0, 0, 0, 0, -t65 * t51, t6, -t64, -g(1) * t87 - g(2) * t88 - g(3) * (-t54 * t90 + t36), 0, 0, 0, 0, 0, 0, -g(1) * (-t34 * t97 + t35 * t55) - g(2) * (-t32 * t97 + t33 * t55) - (t51 * t89 + t55 * t57) * t102, -g(1) * (t34 * t98 + t35 * t59) - g(2) * (t32 * t98 + t33 * t59) - (-t51 * t91 + t57 * t59) * t102, -t6, -g(1) * (-t34 * t70 + t87) - g(2) * (-t32 * t70 + t88) - g(3) * t36 - (t61 * t70 - t90) * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t19 + g(2) * t67 - g(3) * (t72 - t84), g(1) * t20 + g(2) * t75 - g(3) * (-t85 * t56 - t57 * t94), 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t41 - g(3) * t48 + (g(2) * t81 + t56 * t64) * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (-pkin(3) * t99 + t41 + t78) - g(2) * (-pkin(3) * t67 + t79) - g(3) * (-pkin(3) * t84 + t48 + t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * t78 - g(2) * t79 - g(3) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t104 - g(3) * (-t27 * t55 - t54 * t89), g(1) * t8 + g(2) * t103 - g(3) * (-t27 * t59 + t54 * t91), 0, 0;];
taug_reg = t4;
