% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t83 = cos(pkin(12));
t85 = cos(pkin(6));
t72 = t85 * t83;
t80 = sin(pkin(12));
t90 = sin(qJ(2));
t92 = cos(qJ(2));
t54 = -t92 * t72 + t80 * t90;
t81 = sin(pkin(7));
t82 = sin(pkin(6));
t68 = t82 * t81;
t84 = cos(pkin(7));
t98 = t54 * t84 + t83 * t68;
t70 = t85 * t80;
t55 = t92 * t70 + t83 * t90;
t67 = t82 * t80;
t97 = t55 * t84 - t81 * t67;
t69 = t84 * t82;
t96 = t92 * t69 + t85 * t81;
t33 = t90 * t72 + t80 * t92;
t47 = sin(qJ(3));
t91 = cos(qJ(3));
t10 = t33 * t47 + t98 * t91;
t34 = -t90 * t70 + t83 * t92;
t12 = t34 * t47 + t97 * t91;
t73 = t82 * t90;
t24 = t47 * t73 - t96 * t91;
t62 = g(1) * t12 + g(2) * t10 + g(3) * t24;
t43 = pkin(13) + qJ(6);
t41 = sin(t43);
t48 = cos(qJ(4));
t89 = t41 * t48;
t42 = cos(t43);
t88 = t42 * t48;
t44 = sin(pkin(13));
t87 = t44 * t48;
t45 = cos(pkin(13));
t86 = t45 * t48;
t79 = pkin(9) * t81;
t46 = sin(qJ(4));
t78 = t46 * t81;
t77 = t47 * t84;
t76 = t48 * t81;
t75 = t84 * t91;
t74 = t92 * t82;
t25 = t96 * t47 + t91 * t73;
t53 = -t92 * t68 + t85 * t84;
t14 = t25 * t46 - t53 * t48;
t11 = t33 * t91 - t98 * t47;
t49 = t54 * t81 - t83 * t69;
t2 = t11 * t46 - t49 * t48;
t13 = t34 * t91 - t97 * t47;
t50 = t55 * t81 + t84 * t67;
t4 = t13 * t46 - t50 * t48;
t65 = g(1) * t4 + g(2) * t2 + g(3) * t14;
t15 = t25 * t48 + t53 * t46;
t3 = t11 * t48 + t49 * t46;
t5 = t13 * t48 + t50 * t46;
t64 = g(1) * t5 + g(2) * t3 + g(3) * t15;
t59 = t90 * t69;
t31 = -t47 * t59 + t91 * t74;
t58 = t90 * t68;
t20 = t31 * t46 - t48 * t58;
t17 = -t33 * t77 - t54 * t91;
t6 = t17 * t46 - t33 * t76;
t19 = -t34 * t77 - t55 * t91;
t8 = t19 * t46 - t34 * t76;
t63 = g(1) * t8 + g(2) * t6 + g(3) * t20;
t61 = g(1) * t13 + g(2) * t11 + g(3) * t25;
t30 = t47 * t74 + t91 * t59;
t21 = t31 * t48 + t46 * t58;
t18 = t34 * t75 - t55 * t47;
t16 = t33 * t75 - t54 * t47;
t9 = t19 * t48 + t34 * t78;
t7 = t17 * t48 + t33 * t78;
t1 = t62 * t46;
t22 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t55 + g(2) * t54 - g(3) * t74, g(1) * t34 + g(2) * t33 + g(3) * t73, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t17 - g(3) * t31, g(1) * t18 + g(2) * t16 + g(3) * t30, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7 - g(3) * t21, t63, -g(1) * (t18 * t44 + t9 * t45) - g(2) * (t16 * t44 + t7 * t45) - g(3) * (t21 * t45 + t30 * t44) -g(1) * (t18 * t45 - t9 * t44) - g(2) * (t16 * t45 - t7 * t44) - g(3) * (-t21 * t44 + t30 * t45) -t63, -g(1) * (-t55 * pkin(2) + t19 * pkin(3) + t9 * pkin(4) + t18 * pkin(10) + t8 * qJ(5) + t34 * t79) - g(2) * (-t54 * pkin(2) + t17 * pkin(3) + t7 * pkin(4) + t16 * pkin(10) + t6 * qJ(5) + t33 * t79) - g(3) * (pkin(2) * t74 + t31 * pkin(3) + t21 * pkin(4) + pkin(9) * t58 + t30 * pkin(10) + t20 * qJ(5)) 0, 0, 0, 0, 0, -g(1) * (t18 * t41 + t9 * t42) - g(2) * (t16 * t41 + t7 * t42) - g(3) * (t21 * t42 + t30 * t41) -g(1) * (t18 * t42 - t9 * t41) - g(2) * (t16 * t42 - t7 * t41) - g(3) * (-t21 * t41 + t30 * t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t61, 0, 0, 0, 0, 0, t62 * t48, -t1, -g(1) * (-t12 * t86 + t13 * t44) - g(2) * (-t10 * t86 + t11 * t44) - g(3) * (-t24 * t86 + t25 * t44) -g(1) * (t12 * t87 + t13 * t45) - g(2) * (t10 * t87 + t11 * t45) - g(3) * (t24 * t87 + t25 * t45) t1, -t61 * pkin(10) + t62 * (pkin(4) * t48 + qJ(5) * t46 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (-t12 * t88 + t13 * t41) - g(2) * (-t10 * t88 + t11 * t41) - g(3) * (-t24 * t88 + t25 * t41) -g(1) * (t12 * t89 + t13 * t42) - g(2) * (t10 * t89 + t11 * t42) - g(3) * (t24 * t89 + t25 * t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t64, t65 * t45, -t65 * t44, -t64, -g(1) * (-t4 * pkin(4) + t5 * qJ(5)) - g(2) * (-t2 * pkin(4) + t3 * qJ(5)) - g(3) * (-t14 * pkin(4) + t15 * qJ(5)) 0, 0, 0, 0, 0, t65 * t42, -t65 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t42 - t5 * t41) - g(2) * (t10 * t42 - t3 * t41) - g(3) * (-t15 * t41 + t24 * t42) -g(1) * (-t12 * t41 - t5 * t42) - g(2) * (-t10 * t41 - t3 * t42) - g(3) * (-t15 * t42 - t24 * t41);];
taug_reg  = t22;
