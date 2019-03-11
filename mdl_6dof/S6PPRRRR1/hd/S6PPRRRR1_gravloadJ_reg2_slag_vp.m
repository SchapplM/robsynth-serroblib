% Calculate inertial parameters regressor of gravitation load for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t82 = sin(pkin(7));
t83 = sin(pkin(6));
t84 = cos(pkin(13));
t86 = cos(pkin(7));
t87 = cos(pkin(6));
t100 = t83 * t84 * t86 + t82 * t87;
t50 = sin(qJ(3));
t80 = sin(pkin(13));
t67 = t83 * t80;
t99 = cos(qJ(3));
t30 = t100 * t50 + t67 * t99;
t68 = t83 * t82;
t37 = -t68 * t84 + t86 * t87;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t105 = -t30 * t49 + t37 * t52;
t81 = sin(pkin(12));
t66 = t81 * t84;
t85 = cos(pkin(12));
t69 = t85 * t80;
t57 = t66 * t87 + t69;
t65 = t81 * t83;
t101 = t57 * t86 - t65 * t82;
t64 = t81 * t80;
t71 = t85 * t84;
t39 = -t64 * t87 + t71;
t24 = -t101 * t50 + t39 * t99;
t32 = t57 * t82 + t65 * t86;
t104 = -t24 * t49 + t32 * t52;
t56 = -t71 * t87 + t64;
t102 = t56 * t86 + t68 * t85;
t38 = t69 * t87 + t66;
t22 = -t102 * t50 + t38 * t99;
t70 = t85 * t83;
t31 = t56 * t82 - t70 * t86;
t103 = -t22 * t49 + t31 * t52;
t47 = qJ(4) + qJ(5);
t46 = cos(t47);
t48 = sin(qJ(6));
t92 = t46 * t48;
t51 = cos(qJ(6));
t91 = t46 * t51;
t21 = t102 * t99 + t38 * t50;
t44 = pkin(4) * t52 + pkin(3);
t53 = -pkin(10) - pkin(9);
t90 = -t21 * t44 - t22 * t53;
t23 = t101 * t99 + t39 * t50;
t89 = -t23 * t44 - t24 * t53;
t29 = -t100 * t99 + t50 * t67;
t88 = -t29 * t44 - t30 * t53;
t45 = sin(t47);
t10 = -t22 * t45 + t31 * t46;
t11 = t22 * t46 + t31 * t45;
t79 = pkin(5) * t10 + pkin(11) * t11;
t12 = -t24 * t45 + t32 * t46;
t13 = t24 * t46 + t32 * t45;
t78 = pkin(5) * t12 + pkin(11) * t13;
t17 = -t30 * t45 + t37 * t46;
t18 = t30 * t46 + t37 * t45;
t77 = pkin(5) * t17 + pkin(11) * t18;
t76 = t103 * pkin(4);
t75 = t104 * pkin(4);
t74 = t105 * pkin(4);
t73 = -pkin(5) * t46 - pkin(11) * t45;
t63 = g(1) * t12 + g(2) * t10 + g(3) * t17;
t5 = g(1) * t13 + g(2) * t11 + g(3) * t18;
t62 = g(1) * t23 + g(2) * t21 + g(3) * t29;
t61 = g(1) * t24 + g(2) * t22 + g(3) * t30;
t36 = -g(1) * t65 + g(2) * t70 - g(3) * t87;
t6 = t62 * t45;
t2 = t63 * t51;
t1 = t63 * t48;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t61, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t52, -t62 * t49, -t61, -g(1) * (-pkin(3) * t23 + pkin(9) * t24) - g(2) * (-pkin(3) * t21 + pkin(9) * t22) - g(3) * (-pkin(3) * t29 + pkin(9) * t30) 0, 0, 0, 0, 0, 0, t62 * t46, -t6, -t61, -g(1) * t89 - g(2) * t90 - g(3) * t88, 0, 0, 0, 0, 0, 0, -g(1) * (-t23 * t91 + t24 * t48) - g(2) * (-t21 * t91 + t22 * t48) - g(3) * (-t29 * t91 + t30 * t48) -g(1) * (t23 * t92 + t24 * t51) - g(2) * (t21 * t92 + t22 * t51) - g(3) * (t29 * t92 + t30 * t51) t6, -g(1) * (t23 * t73 + t89) - g(2) * (t21 * t73 + t90) - g(3) * (t29 * t73 + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t103 - g(3) * t105, -g(1) * (-t24 * t52 - t32 * t49) - g(2) * (-t22 * t52 - t31 * t49) - g(3) * (-t30 * t52 - t37 * t49) 0, 0, 0, 0, 0, 0, 0, 0, -t63, t5, 0, -g(1) * t75 - g(2) * t76 - g(3) * t74, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (t75 + t78) - g(2) * (t76 + t79) - g(3) * (t74 + t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * t78 - g(2) * t79 - g(3) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t13 * t48 + t23 * t51) - g(2) * (-t11 * t48 + t21 * t51) - g(3) * (-t18 * t48 + t29 * t51) -g(1) * (-t13 * t51 - t23 * t48) - g(2) * (-t11 * t51 - t21 * t48) - g(3) * (-t18 * t51 - t29 * t48) 0, 0;];
taug_reg  = t3;
