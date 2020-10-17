% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 13:11:32
% EndTime: 2019-05-07 13:11:39
% DurationCPUTime: 1.20s
% Computational Cost: add. (719->143), mult. (1988->273), div. (0->0), fcn. (2603->16), ass. (0->86)
t54 = sin(pkin(13));
t60 = sin(qJ(3));
t65 = cos(qJ(3));
t85 = cos(pkin(13));
t48 = -t65 * t54 - t60 * t85;
t55 = sin(pkin(7));
t37 = t48 * t55;
t57 = cos(pkin(7));
t39 = t48 * t57;
t61 = sin(qJ(2));
t62 = sin(qJ(1));
t66 = cos(qJ(2));
t67 = cos(qJ(1));
t86 = cos(pkin(6));
t81 = t67 * t86;
t43 = t62 * t61 - t66 * t81;
t44 = t61 * t81 + t62 * t66;
t74 = -t60 * t54 + t65 * t85;
t56 = sin(pkin(6));
t94 = t56 * t67;
t16 = -t37 * t94 - t43 * t39 - t44 * t74;
t32 = -t43 * t55 + t57 * t94;
t59 = sin(qJ(5));
t64 = cos(qJ(5));
t5 = t16 * t64 + t32 * t59;
t58 = sin(qJ(6));
t63 = cos(qJ(6));
t36 = t74 * t55;
t38 = t74 * t57;
t76 = -t36 * t94 - t43 * t38 + t44 * t48;
t109 = t5 * t58 - t63 * t76;
t108 = t5 * t63 + t58 * t76;
t107 = t16 * t59 - t32 * t64;
t104 = g(1) * t67 + g(2) * t62;
t101 = g(3) * t56;
t99 = t55 * t59;
t98 = t55 * t64;
t97 = t56 * t61;
t96 = t56 * t62;
t95 = t56 * t66;
t93 = t57 * t60;
t92 = t57 * t65;
t91 = t58 * t64;
t90 = t60 * t66;
t89 = t61 * t65;
t88 = t63 * t64;
t87 = pkin(10) + qJ(4);
t84 = t55 * t97;
t83 = t55 * t94;
t82 = t62 * t86;
t80 = t86 * t55;
t79 = g(2) * t44 + g(3) * t97;
t78 = t43 * t57 + t83;
t45 = -t67 * t61 - t66 * t82;
t77 = t45 * t57 + t55 * t96;
t42 = -t55 * t95 + t86 * t57;
t34 = -t45 * t55 + t57 * t96;
t46 = -t61 * t82 + t67 * t66;
t70 = -t37 * t96 - t45 * t39 + t46 * t74;
t6 = t34 * t64 - t59 * t70;
t69 = -t86 * t37 + (-t39 * t66 + t61 * t74) * t56;
t75 = g(1) * t6 + g(2) * t107 + g(3) * (t42 * t64 - t59 * t69);
t73 = t43 * t93 - t44 * t65 + t60 * t83;
t18 = t36 * t96 + t45 * t38 + t46 * t48;
t21 = t86 * t36 + (t38 * t66 + t48 * t61) * t56;
t72 = g(1) * t18 + g(2) * t76 + g(3) * t21;
t71 = g(1) * t46 + t79;
t68 = g(2) * t78 - g(3) * (t57 * t95 + t80);
t53 = t65 * pkin(3) + pkin(2);
t41 = pkin(3) * t93 - t87 * t55;
t31 = (t39 * t61 + t66 * t74) * t56;
t30 = (t38 * t61 - t48 * t66) * t56;
t29 = t46 * t65 + t77 * t60;
t28 = -t46 * t60 + t77 * t65;
t27 = t31 * t64 + t59 * t84;
t26 = t46 * t39 + t45 * t74;
t25 = t46 * t38 - t45 * t48;
t24 = t44 * t39 - t43 * t74;
t23 = t44 * t38 + t43 * t48;
t11 = t26 * t64 + t46 * t99;
t10 = t24 * t64 + t44 * t99;
t9 = t42 * t59 + t64 * t69;
t7 = t34 * t59 + t64 * t70;
t2 = -t18 * t58 + t7 * t63;
t1 = -t18 * t63 - t7 * t58;
t3 = [0, g(1) * t62 - g(2) * t67, t104, 0, 0, 0, 0, 0, g(1) * t44 - g(2) * t46, -g(1) * t43 - g(2) * t45, 0, 0, 0, 0, 0, -g(1) * t73 - g(2) * t29, -g(1) * (t44 * t60 + t78 * t65) - g(2) * t28, -g(1) * t32 - g(2) * t34, -g(1) * (-t62 * pkin(1) + t43 * t41 - t44 * t53) - g(2) * (t67 * pkin(1) + t45 * t41 + t46 * t53) - t104 * t56 * (t55 * t60 * pkin(3) + t87 * t57 + pkin(9)) 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, g(1) * t107 - g(2) * t6, 0, 0, 0, 0, 0, -g(1) * t108 - g(2) * t2, g(1) * t109 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t45 + g(2) * t43 - g(3) * t95, t71, 0, 0, 0, 0, 0, -g(1) * (t45 * t65 - t46 * t93) - g(2) * (-t43 * t65 - t44 * t93) - (-t61 * t93 + t65 * t66) * t101, -g(1) * (-t45 * t60 - t46 * t92) - g(2) * (t43 * t60 - t44 * t92) - (-t57 * t89 - t90) * t101, -t71 * t55, -g(1) * (-t46 * t41 + t45 * t53) - g(2) * (-t44 * t41 - t43 * t53) - (-t41 * t61 + t53 * t66) * t101, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t10 - g(3) * t27, -g(1) * (-t26 * t59 + t46 * t98) - g(2) * (-t24 * t59 + t44 * t98) - g(3) * (-t31 * t59 + t64 * t84) 0, 0, 0, 0, 0, -g(1) * (t11 * t63 + t25 * t58) - g(2) * (t10 * t63 + t23 * t58) - g(3) * (t27 * t63 + t30 * t58) -g(1) * (-t11 * t58 + t25 * t63) - g(2) * (-t10 * t58 + t23 * t63) - g(3) * (-t27 * t58 + t30 * t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t28 + t79 * t60 + t68 * t65, g(1) * t29 - g(2) * t73 - g(3) * (-t60 * t80 + (-t57 * t90 - t89) * t56) 0 (t71 * t60 + (-g(1) * t77 + t68) * t65) * pkin(3), 0, 0, 0, 0, 0, -t72 * t64, t72 * t59, 0, 0, 0, 0, 0, -g(1) * (t18 * t88 + t58 * t70) - g(2) * (-t16 * t58 + t76 * t88) - g(3) * (t21 * t88 + t58 * t69) -g(1) * (-t18 * t91 + t63 * t70) - g(2) * (-t16 * t63 - t76 * t91) - g(3) * (-t21 * t91 + t63 * t69); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t34 + g(2) * t32 - g(3) * t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, g(1) * t7 - g(2) * t5 + g(3) * t9, 0, 0, 0, 0, 0, -t75 * t63, t75 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t109 - g(3) * (-t21 * t63 - t9 * t58) g(1) * t2 - g(2) * t108 - g(3) * (t21 * t58 - t9 * t63);];
taug_reg  = t3;
