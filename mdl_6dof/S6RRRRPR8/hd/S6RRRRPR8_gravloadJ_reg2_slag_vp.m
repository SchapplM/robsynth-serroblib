% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t54 = sin(qJ(1));
t108 = g(2) * t54;
t58 = cos(qJ(1));
t32 = g(1) * t58 + t108;
t53 = sin(qJ(2));
t57 = cos(qJ(2));
t123 = -g(3) * t57 + t32 * t53;
t104 = g(3) * t53;
t50 = qJ(3) + qJ(4);
t45 = sin(t50);
t46 = cos(t50);
t51 = sin(qJ(6));
t55 = cos(qJ(6));
t69 = t45 * t55 - t46 * t51;
t94 = t58 * t45;
t23 = -t54 * t46 + t57 * t94;
t93 = t58 * t46;
t24 = t54 * t45 + t57 * t93;
t7 = t23 * t55 - t24 * t51;
t95 = t54 * t57;
t21 = t45 * t95 + t93;
t22 = t46 * t95 - t94;
t81 = -t21 * t55 + t22 * t51;
t122 = -g(1) * t7 + g(2) * t81 - t69 * t104;
t119 = -t22 * pkin(4) - t21 * qJ(5);
t110 = g(1) * t54;
t72 = -g(2) * t58 + t110;
t68 = t45 * t51 + t46 * t55;
t71 = t21 * t51 + t22 * t55;
t8 = t23 * t51 + t24 * t55;
t117 = g(1) * t8 + g(2) * t71 + t68 * t104;
t113 = -pkin(4) - pkin(5);
t52 = sin(qJ(3));
t112 = pkin(3) * t52;
t111 = pkin(4) * t46;
t59 = -pkin(9) - pkin(8);
t102 = pkin(10) + t59;
t100 = t45 * t53;
t99 = t46 * t53;
t98 = t53 * t59;
t97 = t54 * t52;
t56 = cos(qJ(3));
t96 = t54 * t56;
t44 = t56 * pkin(3) + pkin(2);
t35 = t57 * t44;
t92 = t58 * t52;
t91 = t58 * t56;
t90 = t58 * pkin(1) + t54 * pkin(7);
t89 = qJ(5) * t45;
t87 = t58 * t98;
t86 = t57 * t92;
t48 = t58 * pkin(7);
t85 = pkin(3) * t92 + t54 * t98 + t48;
t84 = t113 * t45;
t83 = t102 * t58;
t82 = -pkin(1) - t35;
t80 = -t21 * pkin(4) + t22 * qJ(5);
t79 = -t23 * pkin(4) + t24 * qJ(5);
t78 = -t44 - t89;
t77 = pkin(3) * t97 + t58 * t35 + t90;
t76 = g(3) * (t35 + (t111 + t89) * t57);
t75 = t57 * pkin(2) + t53 * pkin(8);
t73 = g(1) * t21 - g(2) * t23;
t26 = t52 * t95 + t91;
t66 = t32 * t57;
t64 = t24 * pkin(4) + t23 * qJ(5) + t77;
t4 = g(1) * t23 + g(2) * t21 + g(3) * t100;
t6 = g(1) * t24 + g(2) * t22 + g(3) * t99;
t63 = t82 * t54 + t85;
t41 = pkin(3) * t96;
t62 = -pkin(3) * t86 + t41 + t79;
t25 = t66 + t104;
t60 = -t26 * pkin(3) + t80;
t33 = qJ(5) * t99;
t30 = t72 * t53;
t29 = t57 * t91 + t97;
t28 = -t86 + t96;
t27 = -t56 * t95 + t92;
t18 = t23 * pkin(5);
t15 = t21 * pkin(5);
t11 = t123 * t46;
t10 = t123 * t45;
t9 = g(1) * t22 - g(2) * t24;
t1 = [0, 0, 0, 0, 0, 0, t72, t32, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t57, -t30, -t32, -g(1) * (-t54 * pkin(1) + t48) - g(2) * t90, 0, 0, 0, 0, 0, 0, -g(1) * t27 - g(2) * t29, -g(1) * t26 - g(2) * t28, t30, -g(1) * t48 - g(2) * (t75 * t58 + t90) - (-pkin(1) - t75) * t110, 0, 0, 0, 0, 0, 0, t9, -t73, t30, -g(1) * t63 - g(2) * (t77 - t87) 0, 0, 0, 0, 0, 0, t9, t30, t73, -g(1) * (t63 + t119) - g(2) * (t64 - t87) 0, 0, 0, 0, 0, 0, g(1) * t71 - g(2) * t8, -g(1) * t81 - g(2) * t7, -t30, -g(1) * (-t22 * pkin(5) + t119 + t85) - g(2) * (t24 * pkin(5) - t53 * t83 + t64) - (t53 * pkin(10) + t82) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t25, 0, 0, 0, 0, 0, 0, 0, 0, t123 * t56, -t123 * t52, -t25, -g(3) * t75 + t32 * (pkin(2) * t53 - pkin(8) * t57) 0, 0, 0, 0, 0, 0, t11, -t10, -t25, -g(3) * (t35 - t98) + t32 * (t44 * t53 + t57 * t59) 0, 0, 0, 0, 0, 0, t11, -t25, t10, -t76 + t59 * t66 + (g(3) * t59 + t32 * (-t78 + t111)) * t53, 0, 0, 0, 0, 0, 0, t123 * t68, t123 * t69, t25, -t76 + (-g(3) * pkin(5) * t46 + g(1) * t83 + t102 * t108) * t57 + (g(3) * t102 + t32 * (-t113 * t46 - t78)) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t28 + g(2) * t26 + t52 * t104, g(1) * t29 - g(2) * t27 + t56 * t104, 0, 0, 0, 0, 0, 0, 0, 0, t4, t6, 0, -g(1) * t41 + (g(2) * t91 + t25 * t52) * pkin(3), 0, 0, 0, 0, 0, 0, t4, 0, -t6, -g(1) * t62 - g(2) * t60 - g(3) * (t33 + (-pkin(4) * t45 - t112) * t53) 0, 0, 0, 0, 0, 0, -t122, -t117, 0, -g(1) * (-t18 + t62) - g(2) * (-t15 + t60) - g(3) * t33 - (t84 - t112) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, -t6, -g(1) * t79 - g(2) * t80 - g(3) * (-pkin(4) * t100 + t33) 0, 0, 0, 0, 0, 0, -t122, -t117, 0, -g(1) * (-t18 + t79) - g(2) * (-t15 + t80) - g(3) * (t53 * t84 + t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t117, 0, 0;];
taug_reg  = t1;
