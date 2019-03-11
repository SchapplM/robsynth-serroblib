% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t45 = cos(qJ(3));
t41 = sin(qJ(3));
t82 = g(3) * t41;
t46 = cos(qJ(1));
t85 = g(2) * t46;
t42 = sin(qJ(1));
t86 = g(1) * t42;
t95 = t85 - t86;
t100 = t95 * t45 + t82;
t44 = cos(qJ(4));
t71 = t46 * t44;
t40 = sin(qJ(4));
t77 = t42 * t40;
t16 = t41 * t77 - t71;
t72 = t46 * t40;
t76 = t42 * t44;
t17 = t41 * t76 + t72;
t39 = sin(qJ(6));
t43 = cos(qJ(6));
t3 = t16 * t43 - t17 * t39;
t54 = t39 * t44 - t40 * t43;
t18 = t41 * t72 + t76;
t19 = t41 * t71 - t77;
t62 = -t18 * t43 + t19 * t39;
t98 = g(3) * t45;
t101 = -g(1) * t3 - g(2) * t62 + t54 * t98;
t97 = t19 * pkin(4) + t18 * qJ(5);
t4 = t16 * t39 + t17 * t43;
t53 = t39 * t40 + t43 * t44;
t55 = t18 * t39 + t19 * t43;
t93 = -g(1) * t4 + g(2) * t55 - t53 * t98;
t90 = -pkin(1) - pkin(7);
t89 = -pkin(4) - pkin(5);
t88 = -pkin(8) + pkin(9);
t87 = pkin(8) * t41;
t80 = t40 * t45;
t79 = t41 * t42;
t78 = t41 * t46;
t75 = t42 * t45;
t74 = t44 * t45;
t73 = t45 * t46;
t70 = pkin(3) * t75 + pkin(8) * t79;
t69 = t46 * pkin(1) + t42 * qJ(2);
t68 = qJ(5) * t45;
t66 = pkin(8) * t75;
t65 = t42 * t74;
t64 = t46 * pkin(7) + t69;
t63 = t90 * t42;
t61 = -qJ(5) * t40 - pkin(3);
t60 = -t16 * pkin(4) + t17 * qJ(5);
t59 = t18 * pkin(4) - t19 * qJ(5);
t58 = pkin(4) * t65 + t68 * t77 + t70;
t57 = pkin(3) * t79 + t64;
t56 = g(1) * t18 + g(2) * t16;
t24 = g(1) * t46 + g(2) * t42;
t34 = t46 * qJ(2);
t52 = pkin(3) * t78 - pkin(8) * t73 + t34;
t51 = -pkin(4) * t44 + t61;
t50 = t17 * pkin(4) + t16 * qJ(5) + t57;
t2 = g(1) * t16 - g(2) * t18 + g(3) * t80;
t49 = g(1) * t17 - g(2) * t19 + g(3) * t74;
t48 = t63 + t52;
t35 = t45 * pkin(8);
t25 = t44 * t68;
t20 = t24 * t45;
t8 = g(1) * t79 - g(2) * t78 + t98;
t7 = t100 * t44;
t6 = t100 * t40;
t5 = -g(1) * t19 - g(2) * t17;
t1 = [0, 0, 0, 0, 0, 0, -t95, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t24, -g(1) * (-t42 * pkin(1) + t34) - g(2) * t69, 0, 0, 0, 0, 0, 0, -t24 * t41, -t20, -t95, -g(1) * (t34 + t63) - g(2) * t64, 0, 0, 0, 0, 0, 0, t5, t56, t20, -g(1) * t48 - g(2) * (t57 - t66) 0, 0, 0, 0, 0, 0, t5, t20, -t56, -g(1) * (t48 + t97) - g(2) * (t50 - t66) 0, 0, 0, 0, 0, 0, -g(1) * t55 - g(2) * t4, g(1) * t62 - g(2) * t3, -t20, -g(1) * (t19 * pkin(5) + pkin(9) * t73 + t52 + t97) - g(2) * (t17 * pkin(5) + t50) + (-g(2) * t88 * t45 - g(1) * t90) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t8, -g(1) * t70 - g(3) * (-t41 * pkin(3) + t35) - (-pkin(3) * t45 - t87) * t85, 0, 0, 0, 0, 0, 0, t7, -t8, t6, -g(1) * t58 - g(3) * t35 - t51 * t82 - (t51 * t45 - t87) * t85, 0, 0, 0, 0, 0, 0, t100 * t53, -t100 * t54, t8, -g(1) * (pkin(5) * t65 + t58) - g(3) * (-t45 * pkin(9) + t35) + (pkin(9) * t86 - g(3) * (-pkin(5) * t44 + t51)) * t41 - (t88 * t41 + (t89 * t44 + t61) * t45) * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t49, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t49, -g(1) * t60 - g(2) * t59 - g(3) * (-pkin(4) * t80 + t25) 0, 0, 0, 0, 0, 0, -t101, t93, 0, -g(1) * (-t16 * pkin(5) + t60) - g(2) * (t18 * pkin(5) + t59) - g(3) * (t89 * t80 + t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t93, 0, 0;];
taug_reg  = t1;
