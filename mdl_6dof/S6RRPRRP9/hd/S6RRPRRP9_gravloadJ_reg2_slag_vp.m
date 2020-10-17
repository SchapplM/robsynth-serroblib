% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:34:30
% EndTime: 2019-05-06 18:34:33
% DurationCPUTime: 0.78s
% Computational Cost: add. (706->137), mult. (1263->209), div. (0->0), fcn. (1524->12), ass. (0->76)
t59 = sin(qJ(2));
t60 = sin(qJ(1));
t62 = cos(qJ(2));
t79 = cos(pkin(6));
t96 = cos(qJ(1));
t70 = t79 * t96;
t34 = t59 * t70 + t60 * t62;
t52 = pkin(11) + qJ(4);
t49 = sin(t52);
t50 = cos(t52);
t54 = sin(pkin(6));
t77 = t54 * t96;
t18 = t34 * t50 - t49 * t77;
t33 = t60 * t59 - t62 * t70;
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t101 = t18 * t58 - t33 * t61;
t95 = t33 * t58;
t100 = t18 * t61 + t95;
t75 = t60 * t79;
t36 = -t59 * t75 + t96 * t62;
t87 = t54 * t60;
t22 = t36 * t50 + t49 * t87;
t35 = t96 * t59 + t62 * t75;
t12 = -t22 * t58 + t35 * t61;
t88 = t54 * t59;
t28 = t79 * t49 + t50 * t88;
t83 = t61 * t62;
t1 = g(2) * t101 - g(3) * (-t28 * t58 - t54 * t83) - g(1) * t12;
t55 = cos(pkin(11));
t47 = t55 * pkin(3) + pkin(2);
t86 = t54 * t62;
t37 = t47 * t86;
t98 = g(3) * t37;
t97 = g(3) * t54;
t93 = t34 * t58;
t92 = t35 * t58;
t91 = t36 * t58;
t90 = t50 * t58;
t89 = t50 * t61;
t57 = -pkin(9) - qJ(3);
t85 = t57 * t59;
t84 = t58 * t62;
t82 = -t33 * t47 - t34 * t57;
t81 = -t35 * t47 - t36 * t57;
t80 = t96 * pkin(1) + pkin(8) * t87;
t53 = sin(pkin(11));
t78 = t53 * t87;
t76 = -t60 * pkin(1) + pkin(8) * t77;
t74 = t53 * t77;
t73 = pkin(3) * t78 - t35 * t57 + t36 * t47 + t80;
t72 = pkin(4) * t50 + pkin(10) * t49;
t17 = t34 * t49 + t50 * t77;
t21 = t36 * t49 - t50 * t87;
t71 = -g(1) * t17 + g(2) * t21;
t16 = g(1) * t33 - g(2) * t35;
t48 = t61 * pkin(5) + pkin(4);
t56 = -qJ(6) - pkin(10);
t69 = t48 * t50 - t49 * t56;
t68 = g(1) * t96 + g(2) * t60;
t67 = pkin(3) * t74 + t33 * t57 - t34 * t47 + t76;
t27 = t49 * t88 - t79 * t50;
t66 = g(1) * t21 + g(2) * t17 + g(3) * t27;
t65 = g(1) * t22 + g(2) * t18 + g(3) * t28;
t14 = -g(1) * t35 - g(2) * t33 + g(3) * t86;
t64 = g(1) * t36 + g(2) * t34 + g(3) * t88;
t13 = t22 * t61 + t92;
t11 = t14 * t49;
t8 = t66 * t61;
t7 = t66 * t58;
t6 = -g(1) * (-t35 * t89 + t91) - g(2) * (-t33 * t89 + t93) - (t50 * t83 + t58 * t59) * t97;
t5 = -g(1) * (t35 * t90 + t36 * t61) - g(2) * (t33 * t90 + t34 * t61) - (-t50 * t84 + t59 * t61) * t97;
t4 = g(1) * t100 - g(2) * t13;
t3 = -g(1) * t101 - g(2) * t12;
t2 = g(1) * t13 + g(2) * t100 - g(3) * (-t28 * t61 + t54 * t84);
t9 = [0, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t96, t68, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t34 - g(2) * t36, -t16, -t68 * t54, -g(1) * t76 - g(2) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (-t34 * t55 + t74) - g(2) * (t36 * t55 + t78) -g(1) * (t34 * t53 + t55 * t77) - g(2) * (-t36 * t53 + t55 * t87) t16, -g(1) * (-t34 * pkin(2) - t33 * qJ(3) + t76) - g(2) * (t36 * pkin(2) + t35 * qJ(3) + t80) 0, 0, 0, 0, 0, 0, g(1) * t18 - g(2) * t22, t71, t16, -g(1) * t67 - g(2) * t73, 0, 0, 0, 0, 0, 0, t4, t3, -t71, -g(1) * (-pkin(4) * t18 - pkin(10) * t17 + t67) - g(2) * (t22 * pkin(4) + t21 * pkin(10) + t73) 0, 0, 0, 0, 0, 0, t4, t3, -t71, -g(1) * (-pkin(5) * t95 + t17 * t56 - t18 * t48 + t67) - g(2) * (pkin(5) * t92 - t21 * t56 + t22 * t48 + t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t64, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t55, t14 * t53, -t64, -g(1) * (-t35 * pkin(2) + t36 * qJ(3)) - g(2) * (-t33 * pkin(2) + t34 * qJ(3)) - (pkin(2) * t62 + qJ(3) * t59) * t97, 0, 0, 0, 0, 0, 0, -t14 * t50, t11, -t64, -g(1) * t81 - g(2) * t82 - g(3) * (-t54 * t85 + t37) 0, 0, 0, 0, 0, 0, t6, t5, -t11, -g(1) * (-t72 * t35 + t81) - g(2) * (-t72 * t33 + t82) - t98 - (t72 * t62 - t85) * t97, 0, 0, 0, 0, 0, 0, t6, t5, -t11, -g(1) * (pkin(5) * t91 - t69 * t35 + t81) - g(2) * (pkin(5) * t93 - t69 * t33 + t82) - t98 - (t69 * t62 + (pkin(5) * t58 - t57) * t59) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t65, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t65, -g(1) * (-t21 * pkin(4) + t22 * pkin(10)) - g(2) * (-t17 * pkin(4) + t18 * pkin(10)) - g(3) * (-t27 * pkin(4) + t28 * pkin(10)) 0, 0, 0, 0, 0, 0, t8, -t7, -t65, -g(1) * (-t21 * t48 - t22 * t56) - g(2) * (-t17 * t48 - t18 * t56) - g(3) * (-t27 * t48 - t28 * t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66;];
taug_reg  = t9;
