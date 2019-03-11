% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t106 = sin(qJ(1));
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t107 = cos(qJ(1));
t87 = cos(pkin(6));
t74 = t87 * t107;
t27 = t106 * t59 + t56 * t74;
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t53 = sin(pkin(6));
t83 = t53 * t107;
t15 = t27 * t58 - t55 * t83;
t26 = t106 * t56 - t59 * t74;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t119 = t15 * t54 - t26 * t57;
t118 = t15 * t57 + t26 * t54;
t52 = qJ(4) + qJ(5);
t46 = sin(t52);
t47 = cos(t52);
t117 = t15 * t46 - t26 * t47;
t116 = t15 * t47 + t26 * t46;
t48 = qJ(6) + t52;
t43 = sin(t48);
t44 = cos(t48);
t115 = t15 * t43 - t26 * t44;
t114 = t15 * t44 + t26 * t43;
t73 = t87 * t106;
t29 = t107 * t59 - t56 * t73;
t82 = t53 * t106;
t19 = t29 * t58 + t55 * t82;
t28 = t107 * t56 + t59 * t73;
t11 = -t19 * t54 + t28 * t57;
t95 = t53 * t56;
t25 = t87 * t55 + t58 * t95;
t94 = t53 * t59;
t113 = g(2) * t119 - g(3) * (-t25 * t54 - t57 * t94) - g(1) * t11;
t9 = -t19 * t46 + t28 * t47;
t3 = g(2) * t117 - g(3) * (-t25 * t46 - t47 * t94) - g(1) * t9;
t60 = -pkin(11) - pkin(10);
t110 = g(3) * t53;
t109 = t54 * pkin(4);
t31 = pkin(5) * t46 + t109;
t108 = pkin(9) + t31;
t99 = t43 * t58;
t98 = t44 * t58;
t97 = t46 * t58;
t96 = t47 * t58;
t93 = t54 * t56;
t92 = t54 * t58;
t91 = t57 * t58;
t90 = t58 * t59;
t89 = pkin(2) * t94 + pkin(9) * t95;
t49 = t57 * pkin(4);
t32 = pkin(5) * t47 + t49;
t88 = t107 * pkin(1) + pkin(8) * t82;
t86 = t29 * pkin(2) + t88;
t85 = pkin(9) + t109;
t84 = g(3) * t89;
t20 = t26 * pkin(2);
t81 = t27 * pkin(9) - t20;
t22 = t28 * pkin(2);
t80 = t29 * pkin(9) - t22;
t79 = -t27 * t55 - t58 * t83;
t78 = -t106 * pkin(1) + pkin(8) * t83;
t77 = pkin(3) * t58 + pkin(10) * t55;
t18 = t29 * t55 - t58 * t82;
t76 = g(1) * t79 + g(2) * t18;
t75 = g(1) * t26 - g(2) * t28;
t30 = pkin(3) + t32;
t51 = -pkin(12) + t60;
t72 = t30 * t58 - t51 * t55;
t45 = t49 + pkin(3);
t71 = t45 * t58 - t55 * t60;
t70 = t28 * pkin(9) + t86;
t69 = -t27 * pkin(2) + t78;
t24 = -t55 * t95 + t87 * t58;
t68 = g(1) * t18 - g(2) * t79 - g(3) * t24;
t67 = g(1) * t19 + g(2) * t15 + g(3) * t25;
t66 = g(1) * t107 + g(2) * t106;
t65 = -t26 * pkin(9) + t69;
t64 = -g(1) * t28 - g(2) * t26 + g(3) * t94;
t63 = g(1) * t29 + g(2) * t27 + g(3) * t95;
t13 = t64 * t55;
t12 = t19 * t57 + t28 * t54;
t10 = t19 * t47 + t28 * t46;
t8 = t19 * t44 + t28 * t43;
t7 = -t19 * t43 + t28 * t44;
t4 = g(1) * t10 + g(2) * t116 - g(3) * (-t25 * t47 + t46 * t94);
t2 = g(1) * t8 + g(2) * t114 - g(3) * (-t25 * t44 + t43 * t94);
t1 = -g(1) * t7 + g(2) * t115 - g(3) * (-t25 * t43 - t44 * t94);
t5 = [0, 0, 0, 0, 0, 0, g(1) * t106 - g(2) * t107, t66, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t27 - g(2) * t29, -t75, -t66 * t53, -g(1) * t78 - g(2) * t88, 0, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t19, t76, t75, -g(1) * t65 - g(2) * t70, 0, 0, 0, 0, 0, 0, g(1) * t118 - g(2) * t12, -g(1) * t119 - g(2) * t11, -t76, -g(1) * (-pkin(3) * t15 + pkin(10) * t79 + t65) - g(2) * (t19 * pkin(3) + t18 * pkin(10) + t70) 0, 0, 0, 0, 0, 0, g(1) * t116 - g(2) * t10, -g(1) * t117 - g(2) * t9, -t76, -g(1) * (-t15 * t45 - t85 * t26 - t60 * t79 + t69) - g(2) * (-t18 * t60 + t19 * t45 + t28 * t85 + t86) 0, 0, 0, 0, 0, 0, g(1) * t114 - g(2) * t8, -g(1) * t115 - g(2) * t7, -t76, -g(1) * (-t108 * t26 - t15 * t30 - t51 * t79 + t69) - g(2) * (t108 * t28 - t18 * t51 + t19 * t30 + t86); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t63, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t58, t13, -t63, -g(1) * t80 - g(2) * t81 - t84, 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t91 + t29 * t54) - g(2) * (-t26 * t91 + t27 * t54) - (t57 * t90 + t93) * t110, -g(1) * (t28 * t92 + t29 * t57) - g(2) * (t26 * t92 + t27 * t57) - (-t54 * t90 + t56 * t57) * t110, -t13, -g(1) * (-t28 * t77 + t80) - g(2) * (-t77 * t26 + t81) - g(3) * (t77 * t94 + t89) 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t96 + t29 * t46) - g(2) * (-t26 * t96 + t27 * t46) - (t46 * t56 + t47 * t90) * t110, -g(1) * (t28 * t97 + t29 * t47) - g(2) * (t26 * t97 + t27 * t47) - (-t46 * t90 + t47 * t56) * t110, -t13, -g(1) * (-t28 * t71 + t29 * t85 - t22) - g(2) * (-t71 * t26 + t27 * t85 - t20) - t84 - (pkin(4) * t93 + t59 * t71) * t110, 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t98 + t29 * t43) - g(2) * (-t26 * t98 + t27 * t43) - (t43 * t56 + t44 * t90) * t110, -g(1) * (t28 * t99 + t29 * t44) - g(2) * (t26 * t99 + t27 * t44) - (-t43 * t90 + t44 * t56) * t110, -t13, -g(1) * (t108 * t29 - t72 * t28 - t22) - g(2) * (t108 * t27 - t72 * t26 - t20) - t84 - (t31 * t56 + t59 * t72) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t67, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t57, -t68 * t54, -t67, -g(1) * (-t18 * pkin(3) + t19 * pkin(10)) - g(2) * (pkin(3) * t79 + t15 * pkin(10)) - g(3) * (t24 * pkin(3) + t25 * pkin(10)) 0, 0, 0, 0, 0, 0, t68 * t47, -t68 * t46, -t67, -g(1) * (-t18 * t45 - t19 * t60) - g(2) * (-t15 * t60 + t45 * t79) - g(3) * (t24 * t45 - t25 * t60) 0, 0, 0, 0, 0, 0, t68 * t44, -t68 * t43, -t67, -g(1) * (-t18 * t30 - t19 * t51) - g(2) * (-t15 * t51 + t30 * t79) - g(3) * (t24 * t30 - t25 * t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, g(1) * t12 + g(2) * t118 - g(3) * (-t25 * t57 + t54 * t94) 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t113 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t19 * t31 + t28 * t32) - g(2) * (-t15 * t31 + t26 * t32) - g(3) * (-t25 * t31 - t32 * t94); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t5;
