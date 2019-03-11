% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t46 = cos(qJ(1));
t78 = sin(pkin(12));
t82 = cos(pkin(6));
t67 = t82 * t78;
t80 = cos(pkin(12));
t87 = sin(qJ(1));
t30 = t46 * t67 + t87 * t80;
t43 = sin(qJ(3));
t88 = cos(qJ(3));
t68 = t82 * t80;
t56 = -t46 * t68 + t87 * t78;
t40 = sin(pkin(6));
t79 = sin(pkin(7));
t76 = t40 * t79;
t81 = cos(pkin(7));
t97 = t46 * t76 + t56 * t81;
t14 = t30 * t43 + t97 * t88;
t41 = sin(qJ(6));
t44 = cos(qJ(6));
t17 = -t30 * t88 + t97 * t43;
t42 = sin(qJ(4));
t45 = cos(qJ(4));
t77 = t40 * t81;
t92 = -t46 * t77 + t56 * t79;
t8 = t17 * t42 + t45 * t92;
t101 = -t14 * t44 + t8 * t41;
t100 = t14 * t41 + t8 * t44;
t9 = t17 * t45 - t42 * t92;
t52 = t46 * t78 + t87 * t68;
t93 = t52 * t81 - t87 * t76;
t31 = t46 * t80 - t87 * t67;
t18 = t31 * t43 + t93 * t88;
t65 = t79 * t82;
t66 = t81 * t80;
t23 = -t88 * t65 + (t43 * t78 - t66 * t88) * t40;
t58 = g(1) * t18 + g(2) * t14 + g(3) * t23;
t86 = t41 * t42;
t85 = t42 * t44;
t83 = qJ(2) * t40;
t84 = t46 * pkin(1) + t87 * t83;
t19 = t31 * t88 - t93 * t43;
t47 = -t52 * t79 - t87 * t77;
t10 = t19 * t42 + t47 * t45;
t73 = g(1) * t8 + g(2) * t10;
t11 = t19 * t45 - t47 * t42;
t72 = g(1) * t9 + g(2) * t11;
t71 = -t87 * pkin(1) + t46 * t83;
t70 = -g(1) * t14 + g(2) * t18;
t62 = g(1) * t87 - g(2) * t46;
t61 = -g(1) * t46 - g(2) * t87;
t24 = t43 * t65 + (t43 * t66 + t88 * t78) * t40;
t51 = -t80 * t76 + t82 * t81;
t12 = t24 * t42 - t51 * t45;
t60 = g(1) * t10 - g(2) * t8 + g(3) * t12;
t13 = t24 * t45 + t51 * t42;
t59 = g(1) * t11 - g(2) * t9 + g(3) * t13;
t57 = g(1) * t19 - g(2) * t17 + g(3) * t24;
t29 = -g(3) * t82 - t62 * t40;
t5 = t10 * t41 + t18 * t44;
t4 = t10 * t44 - t18 * t41;
t3 = t58 * t45;
t2 = t58 * t42;
t1 = [0, t62, -t61, g(1) * t30 - g(2) * t31, -g(1) * t56 + g(2) * t52, t61 * t40, -g(1) * t71 - g(2) * t84, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t19, t70, 0, 0, 0, 0, 0, -t72, t73, -t70, t72, -t73, -g(1) * (-t30 * pkin(2) + t17 * pkin(3) + t9 * pkin(4) - pkin(10) * t14 + t8 * qJ(5) + t71) - g(2) * (t31 * pkin(2) + t19 * pkin(3) + t11 * pkin(4) + t18 * pkin(10) + t10 * qJ(5) + t84) + (g(1) * t92 + g(2) * t47) * pkin(9), 0, 0, 0, 0, 0, -g(1) * t101 - g(2) * t5, -g(1) * t100 - g(2) * t4; 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t57, 0, 0, 0, 0, 0, t3, -t2, -t57, -t3, t2, -t57 * pkin(10) + t58 * (pkin(4) * t45 + qJ(5) * t42 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (-t18 * t86 + t19 * t44) - g(2) * (-t14 * t86 - t17 * t44) - g(3) * (-t23 * t86 + t24 * t44) -g(1) * (-t18 * t85 - t19 * t41) - g(2) * (-t14 * t85 + t17 * t41) - g(3) * (-t23 * t85 - t24 * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t59, 0, -t60, -t59, -g(1) * (-t10 * pkin(4) + t11 * qJ(5)) - g(2) * (pkin(4) * t8 - qJ(5) * t9) - g(3) * (-t12 * pkin(4) + t13 * qJ(5)) 0, 0, 0, 0, 0, -t59 * t41, -t59 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 + g(2) * t100 - g(3) * (t12 * t44 - t23 * t41) g(1) * t5 - g(2) * t101 - g(3) * (-t12 * t41 - t23 * t44);];
taug_reg  = t1;
