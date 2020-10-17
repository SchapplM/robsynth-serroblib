% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:11:28
% EndTime: 2019-05-08 06:11:32
% DurationCPUTime: 1.05s
% Computational Cost: add. (813->171), mult. (1708->258), div. (0->0), fcn. (2087->12), ass. (0->87)
t100 = sin(qJ(1));
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t101 = cos(qJ(1));
t85 = cos(pkin(6));
t73 = t85 * t101;
t29 = t100 * t58 + t55 * t73;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t52 = sin(pkin(6));
t81 = t52 * t101;
t17 = t29 * t57 - t54 * t81;
t28 = t100 * t55 - t58 * t73;
t51 = qJ(4) + qJ(5);
t46 = sin(t51);
t47 = cos(t51);
t111 = t17 * t46 - t28 * t47;
t110 = t17 * t47 + t28 * t46;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t109 = t17 * t53 - t28 * t56;
t108 = t17 * t56 + t28 * t53;
t72 = t85 * t100;
t31 = t101 * t58 - t55 * t72;
t80 = t52 * t100;
t21 = t31 * t57 + t54 * t80;
t30 = t101 * t55 + t58 * t72;
t13 = -t21 * t53 + t30 * t56;
t93 = t52 * t55;
t27 = t85 * t54 + t57 * t93;
t92 = t52 * t58;
t107 = g(2) * t109 - g(3) * (-t27 * t53 - t56 * t92) - g(1) * t13;
t11 = -t21 * t46 + t30 * t47;
t1 = g(2) * t111 - g(3) * (-t27 * t46 - t47 * t92) - g(1) * t11;
t59 = -pkin(11) - pkin(10);
t104 = g(3) * t52;
t103 = t53 * pkin(4);
t33 = pkin(5) * t46 + t103;
t102 = pkin(9) + t33;
t95 = t46 * t57;
t94 = t47 * t57;
t91 = t53 * t55;
t90 = t53 * t57;
t89 = t56 * t57;
t88 = t57 * t58;
t87 = pkin(2) * t92 + pkin(9) * t93;
t48 = t56 * pkin(4);
t34 = pkin(5) * t47 + t48;
t86 = t101 * pkin(1) + pkin(8) * t80;
t84 = t31 * pkin(2) + t86;
t83 = pkin(9) + t103;
t82 = g(3) * t87;
t22 = t28 * pkin(2);
t79 = t29 * pkin(9) - t22;
t24 = t30 * pkin(2);
t78 = t31 * pkin(9) - t24;
t77 = -t100 * pkin(1) + pkin(8) * t81;
t76 = pkin(3) * t57 + pkin(10) * t54;
t16 = t29 * t54 + t57 * t81;
t20 = t31 * t54 - t57 * t80;
t75 = -g(1) * t16 + g(2) * t20;
t74 = g(1) * t28 - g(2) * t30;
t32 = pkin(3) + t34;
t50 = -qJ(6) + t59;
t71 = t32 * t57 - t50 * t54;
t45 = t48 + pkin(3);
t70 = t45 * t57 - t54 * t59;
t69 = t30 * pkin(9) + t84;
t68 = -t29 * pkin(2) + t77;
t26 = t54 * t93 - t85 * t57;
t67 = g(1) * t20 + g(2) * t16 + g(3) * t26;
t66 = g(1) * t21 + g(2) * t17 + g(3) * t27;
t65 = g(1) * t101 + g(2) * t100;
t64 = -t28 * pkin(9) + t68;
t63 = -g(1) * t30 - g(2) * t28 + g(3) * t92;
t62 = g(1) * t31 + g(2) * t29 + g(3) * t93;
t15 = t63 * t54;
t14 = t21 * t56 + t30 * t53;
t12 = t21 * t47 + t30 * t46;
t8 = t67 * t47;
t7 = t67 * t46;
t6 = g(1) * t110 - g(2) * t12;
t5 = -g(1) * t111 - g(2) * t11;
t4 = -g(1) * (-t30 * t94 + t31 * t46) - g(2) * (-t28 * t94 + t29 * t46) - (t46 * t55 + t47 * t88) * t104;
t3 = -g(1) * (t30 * t95 + t31 * t47) - g(2) * (t28 * t95 + t29 * t47) - (-t46 * t88 + t47 * t55) * t104;
t2 = g(1) * t12 + g(2) * t110 - g(3) * (-t27 * t47 + t46 * t92);
t9 = [0, 0, 0, 0, 0, 0, g(1) * t100 - g(2) * t101, t65, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t29 - g(2) * t31, -t74, -t65 * t52, -g(1) * t77 - g(2) * t86, 0, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t21, t75, t74, -g(1) * t64 - g(2) * t69, 0, 0, 0, 0, 0, 0, g(1) * t108 - g(2) * t14, -g(1) * t109 - g(2) * t13, -t75, -g(1) * (-pkin(3) * t17 - pkin(10) * t16 + t64) - g(2) * (t21 * pkin(3) + t20 * pkin(10) + t69) 0, 0, 0, 0, 0, 0, t6, t5, -t75, -g(1) * (t16 * t59 - t17 * t45 - t83 * t28 + t68) - g(2) * (-t20 * t59 + t21 * t45 + t83 * t30 + t84) 0, 0, 0, 0, 0, 0, t6, t5, -t75, -g(1) * (-t102 * t28 + t16 * t50 - t17 * t32 + t68) - g(2) * (t102 * t30 - t20 * t50 + t21 * t32 + t84); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t62, 0, 0, 0, 0, 0, 0, 0, 0, -t63 * t57, t15, -t62, -g(1) * t78 - g(2) * t79 - t82, 0, 0, 0, 0, 0, 0, -g(1) * (-t30 * t89 + t31 * t53) - g(2) * (-t28 * t89 + t29 * t53) - (t56 * t88 + t91) * t104, -g(1) * (t30 * t90 + t31 * t56) - g(2) * (t28 * t90 + t29 * t56) - (-t53 * t88 + t55 * t56) * t104, -t15, -g(1) * (-t76 * t30 + t78) - g(2) * (-t76 * t28 + t79) - g(3) * (t76 * t92 + t87) 0, 0, 0, 0, 0, 0, t4, t3, -t15, -g(1) * (-t70 * t30 + t83 * t31 - t24) - g(2) * (-t70 * t28 + t83 * t29 - t22) - t82 - (pkin(4) * t91 + t70 * t58) * t104, 0, 0, 0, 0, 0, 0, t4, t3, -t15, -g(1) * (t102 * t31 - t71 * t30 - t24) - g(2) * (t102 * t29 - t71 * t28 - t22) - t82 - (t33 * t55 + t71 * t58) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t66, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t56, -t67 * t53, -t66, -g(1) * (-t20 * pkin(3) + t21 * pkin(10)) - g(2) * (-t16 * pkin(3) + t17 * pkin(10)) - g(3) * (-t26 * pkin(3) + t27 * pkin(10)) 0, 0, 0, 0, 0, 0, t8, -t7, -t66, -g(1) * (-t20 * t45 - t21 * t59) - g(2) * (-t16 * t45 - t17 * t59) - g(3) * (-t26 * t45 - t27 * t59) 0, 0, 0, 0, 0, 0, t8, -t7, -t66, -g(1) * (-t20 * t32 - t21 * t50) - g(2) * (-t16 * t32 - t17 * t50) - g(3) * (-t26 * t32 - t27 * t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, g(1) * t14 + g(2) * t108 - g(3) * (-t27 * t56 + t53 * t92) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t107 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t21 * t33 + t30 * t34) - g(2) * (-t17 * t33 + t28 * t34) - g(3) * (-t27 * t33 - t34 * t92); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67;];
taug_reg  = t9;
