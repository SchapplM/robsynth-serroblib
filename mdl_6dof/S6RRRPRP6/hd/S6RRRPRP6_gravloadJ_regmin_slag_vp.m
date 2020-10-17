% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:08:00
% EndTime: 2019-05-07 08:08:02
% DurationCPUTime: 0.52s
% Computational Cost: add. (402->113), mult. (757->185), div. (0->0), fcn. (905->12), ass. (0->72)
t50 = sin(qJ(2));
t51 = sin(qJ(1));
t54 = cos(qJ(2));
t55 = cos(qJ(1));
t72 = cos(pkin(6));
t65 = t55 * t72;
t22 = t51 * t50 - t54 * t65;
t48 = sin(qJ(5));
t23 = t50 * t65 + t51 * t54;
t44 = qJ(3) + pkin(11);
t41 = sin(t44);
t42 = cos(t44);
t45 = sin(pkin(6));
t77 = t45 * t55;
t5 = t23 * t42 - t41 * t77;
t52 = cos(qJ(5));
t94 = -t22 * t52 + t5 * t48;
t89 = t22 * t48;
t93 = t5 * t52 + t89;
t66 = t51 * t72;
t24 = t55 * t50 + t54 * t66;
t25 = -t50 * t66 + t55 * t54;
t80 = t45 * t51;
t9 = t25 * t42 + t41 * t80;
t1 = t24 * t52 - t9 * t48;
t81 = t45 * t50;
t17 = t72 * t41 + t42 * t81;
t75 = t52 * t54;
t92 = g(2) * t94 - g(3) * (-t17 * t48 - t45 * t75) - g(1) * t1;
t90 = g(3) * t45;
t87 = t23 * t48;
t86 = t24 * t48;
t85 = t25 * t48;
t49 = sin(qJ(3));
t84 = t25 * t49;
t83 = t42 * t48;
t82 = t42 * t52;
t53 = cos(qJ(3));
t79 = t45 * t53;
t78 = t45 * t54;
t76 = t48 * t54;
t40 = t53 * pkin(3) + pkin(2);
t47 = -qJ(4) - pkin(9);
t74 = -t22 * t40 - t23 * t47;
t73 = -t24 * t40 - t25 * t47;
t71 = t49 * t81;
t70 = t49 * t80;
t69 = t51 * t79;
t33 = t49 * t77;
t68 = t53 * t77;
t67 = t23 * t53 - t33;
t64 = t72 * t53;
t63 = t55 * pkin(1) + pkin(3) * t70 + pkin(8) * t80 - t24 * t47 + t25 * t40;
t62 = g(1) * t22 - g(2) * t24;
t39 = t52 * pkin(5) + pkin(4);
t46 = -qJ(6) - pkin(10);
t61 = t39 * t42 - t41 * t46;
t4 = t23 * t41 + t42 * t77;
t60 = t23 * t49 + t68;
t59 = -t51 * pkin(1) + pkin(3) * t33 + pkin(8) * t77 + t22 * t47 - t23 * t40;
t16 = t41 * t81 - t72 * t42;
t8 = t25 * t41 - t42 * t80;
t58 = g(1) * t8 + g(2) * t4 + g(3) * t16;
t3 = -g(1) * t24 - g(2) * t22 + g(3) * t78;
t57 = g(1) * t25 + g(2) * t23 + g(3) * t81;
t38 = pkin(3) * t64;
t30 = pkin(3) * t69;
t26 = t40 * t78;
t11 = t25 * t53 + t70;
t10 = t69 - t84;
t2 = t9 * t52 + t86;
t6 = [0, g(1) * t51 - g(2) * t55, g(1) * t55 + g(2) * t51, 0, 0, 0, 0, 0, g(1) * t23 - g(2) * t25, -t62, 0, 0, 0, 0, 0, g(1) * t67 - g(2) * t11, -g(1) * t60 - g(2) * t10, t62, -g(1) * t59 - g(2) * t63, 0, 0, 0, 0, 0, g(1) * t93 - g(2) * t2, -g(1) * t94 - g(2) * t1, g(1) * t4 - g(2) * t8, -g(1) * (-pkin(5) * t89 - t39 * t5 + t4 * t46 + t59) - g(2) * (pkin(5) * t86 + t9 * t39 - t8 * t46 + t63); 0, 0, 0, 0, 0, 0, 0, 0, -t3, t57, 0, 0, 0, 0, 0, -t3 * t53, t3 * t49, -t57, -g(1) * t73 - g(2) * t74 - g(3) * (-t47 * t81 + t26) 0, 0, 0, 0, 0, -g(1) * (-t24 * t82 + t85) - g(2) * (-t22 * t82 + t87) - (t42 * t75 + t48 * t50) * t90, -g(1) * (t24 * t83 + t25 * t52) - g(2) * (t22 * t83 + t23 * t52) - (-t42 * t76 + t50 * t52) * t90, -t3 * t41, -g(1) * (pkin(5) * t85 - t24 * t61 + t73) - g(2) * (pkin(5) * t87 - t22 * t61 + t74) - g(3) * t26 - (t61 * t54 + (pkin(5) * t48 - t47) * t50) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t60 - g(3) * (t64 - t71) g(1) * t11 + g(2) * t67 - g(3) * (-t49 * t72 - t50 * t79) 0, -g(1) * t30 - g(3) * t38 + (g(2) * t68 + t49 * t57) * pkin(3), 0, 0, 0, 0, 0, t58 * t52, -t58 * t48, -g(1) * t9 - g(2) * t5 - g(3) * t17, -g(1) * (-pkin(3) * t84 - t8 * t39 - t9 * t46 + t30) - g(2) * (-pkin(3) * t60 - t4 * t39 - t5 * t46) - g(3) * (-pkin(3) * t71 - t16 * t39 - t17 * t46 + t38); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, g(1) * t2 + g(2) * t93 - g(3) * (-t17 * t52 + t45 * t76) 0, t92 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58;];
taug_reg  = t6;
