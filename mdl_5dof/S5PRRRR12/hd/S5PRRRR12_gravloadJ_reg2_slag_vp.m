% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:20
% EndTime: 2024-09-28 18:08:21
% DurationCPUTime: 1.03s
% Computational Cost: add. (656->140), mult. (881->237), div. (0->0), fcn. (925->26), ass. (0->101)
t67 = sin(pkin(11));
t72 = cos(pkin(5));
t107 = t67 * t72;
t70 = cos(pkin(11));
t68 = sin(pkin(6));
t117 = pkin(10) * t68;
t94 = t67 * t117;
t35 = t70 * pkin(4) + t72 * t94;
t93 = t70 * t117;
t37 = pkin(4) * t107 - t93;
t74 = sin(qJ(4));
t78 = cos(qJ(4));
t118 = t35 * t74 + t37 * t78;
t10 = -pkin(3) * t107 - t118;
t16 = -t35 * t78 + t37 * t74;
t12 = t70 * pkin(3) - t16;
t75 = sin(qJ(3));
t79 = cos(qJ(3));
t122 = t10 * t79 - t12 * t75;
t102 = t70 * t72;
t36 = -t67 * pkin(4) + t72 * t93;
t38 = pkin(4) * t102 + t94;
t87 = t36 * t74 + t38 * t78;
t11 = pkin(3) * t102 + t87;
t17 = t36 * t78 - t38 * t74;
t15 = pkin(3) * t67 - t17;
t121 = t11 * t79 - t15 * t75;
t120 = -t11 * t75 - t15 * t79;
t119 = -t10 * t75 - t12 * t79;
t69 = sin(pkin(5));
t116 = g(3) * t69;
t44 = -t74 * pkin(4) + t78 * t117;
t109 = t44 * t75;
t66 = qJ(2) + qJ(3);
t65 = qJ(4) + t66;
t58 = cos(t65);
t108 = t67 * t58;
t106 = t68 * t69;
t76 = sin(qJ(2));
t105 = t69 * t76;
t80 = cos(qJ(2));
t104 = t69 * t80;
t103 = t70 * t58;
t71 = cos(pkin(6));
t73 = sin(qJ(5));
t101 = t71 * t73;
t77 = cos(qJ(5));
t100 = t71 * t77;
t99 = t72 * t73;
t98 = t72 * t76;
t97 = t72 * t77;
t96 = t72 * t80;
t95 = t75 * t76;
t92 = t73 * t106;
t91 = t77 * t106;
t90 = t71 * t99;
t89 = t71 * t97;
t61 = pkin(5) - t66;
t60 = pkin(5) + t66;
t43 = -pkin(4) * t78 - t74 * t117;
t42 = pkin(3) - t43;
t84 = -t42 * t79 - t109;
t82 = -pkin(3) * t95 + (t79 * pkin(3) + pkin(2)) * t80;
t81 = -g(1) * (-t67 * t96 - t70 * t76) - g(2) * (-t67 * t76 + t70 * t96) - g(3) * t104;
t63 = cos(t66);
t62 = sin(t66);
t57 = sin(t65);
t56 = -qJ(4) + t61;
t55 = qJ(4) + t60;
t54 = cos(t61);
t53 = sin(t60);
t52 = cos(t56);
t51 = sin(t55);
t49 = cos(t60) / 0.2e1;
t48 = sin(t61) / 0.2e1;
t47 = cos(t55) / 0.2e1;
t46 = sin(t56) / 0.2e1;
t45 = -t76 * pkin(2) - pkin(3) * t62;
t41 = t54 / 0.2e1 + t49;
t40 = t48 - t53 / 0.2e1;
t39 = t44 * t79;
t34 = t52 / 0.2e1 + t47;
t33 = t46 - t51 / 0.2e1;
t28 = -t70 * t100 + t67 * t99;
t27 = t70 * t101 + t67 * t97;
t26 = -t67 * t100 - t70 * t99;
t25 = t67 * t101 - t70 * t97;
t24 = -t67 * t73 + t70 * t89;
t23 = t67 * t77 + t70 * t90;
t22 = -t67 * t89 - t70 * t73;
t21 = t67 * t90 - t70 * t77;
t19 = t82 * t72;
t18 = (-t75 * t42 + t39) * t105;
t7 = -g(1) * (-t67 * t40 - t70 * t63) - g(2) * (t70 * t40 - t67 * t63) - g(3) * (t49 - t54 / 0.2e1);
t6 = -g(1) * (-t67 * t41 - t70 * t62) - g(2) * (t70 * t41 - t67 * t62) - g(3) * (t53 / 0.2e1 + t48);
t5 = (-g(1) * (-t57 * t107 + t103) - g(2) * (t57 * t102 + t108) - t57 * t116) * t68;
t4 = -g(1) * (-t67 * t33 - t103) - g(2) * (t70 * t33 - t108) - g(3) * (t47 - t52 / 0.2e1);
t3 = -g(1) * (-t67 * t34 - t70 * t57) - g(2) * (t70 * t34 - t67 * t57) - g(3) * (t51 / 0.2e1 + t46);
t2 = -g(1) * (-t22 * t57 + t28 * t58) - g(2) * (-t24 * t57 + t26 * t58) - (-t57 * t100 - t58 * t73) * t116;
t1 = -g(1) * (t21 * t57 - t27 * t58) - g(2) * (-t23 * t57 - t25 * t58) - (-t57 * t101 + t58 * t77) * t116;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -g(1) * (t67 * t98 - t70 * t80) - g(2) * (-t67 * t80 - t70 * t98) + g(3) * t105, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, t81 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (-t67 * t19 + t70 * t45) - g(2) * (t70 * t19 + t67 * t45) - t82 * t116, 0, 0, 0, 0, 0, 0, t1, t2, t5, -g(1) * (-(t70 * pkin(2) - t119) * t76 + (-pkin(2) * t107 + t122) * t80) - g(2) * (-(t67 * pkin(2) - t120) * t76 + (pkin(2) * t102 + t121) * t80) - g(3) * ((pkin(2) - t84) * t104 + t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, ((g(1) * t70 + g(2) * t67) * t62 + ((g(1) * t67 - g(2) * t70) * t72 - t116) * (t79 * t80 - t95)) * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, t5, -g(1) * (t119 * t76 + t122 * t80) - g(2) * (t120 * t76 + t121 * t80) - g(3) * (-t84 * t104 + t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t5, -g(1) * ((-t118 * t79 + t16 * t75) * t80 + (t118 * t75 + t16 * t79) * t76) - g(2) * ((t17 * t75 + t79 * t87) * t80 + (t17 * t79 - t87 * t75) * t76) - ((t43 * t75 + t39) * t76 - (t43 * t79 - t109) * t80) * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t22 * t58 + t28 * t57 + t67 * t91) - g(2) * (t24 * t58 + t26 * t57 - t70 * t91) - g(3) * (t68 * t97 + (t58 * t100 - t57 * t73) * t69), -g(1) * (t21 * t58 + t27 * t57 - t67 * t92) - g(2) * (-t23 * t58 + t25 * t57 + t70 * t92) - g(3) * (-t68 * t99 + (-t58 * t101 - t57 * t77) * t69), 0, 0;];
taug_reg = t8;
