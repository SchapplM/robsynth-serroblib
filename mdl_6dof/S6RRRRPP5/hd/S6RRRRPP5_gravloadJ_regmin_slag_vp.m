% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t48 = cos(qJ(3));
t38 = t48 * pkin(3) + pkin(2);
t49 = cos(qJ(2));
t29 = t49 * t38;
t90 = pkin(1) + t29;
t47 = sin(qJ(1));
t45 = sin(qJ(3));
t50 = cos(qJ(1));
t72 = t50 * t45;
t23 = t47 * t48 - t49 * t72;
t86 = g(1) * t47;
t59 = -g(2) * t50 + t86;
t84 = g(2) * t47;
t60 = g(1) * t50 + t84;
t46 = sin(qJ(2));
t19 = -g(3) * t49 + t60 * t46;
t89 = -pkin(4) - pkin(5);
t88 = pkin(3) * t45;
t44 = qJ(3) + qJ(4);
t40 = cos(t44);
t87 = pkin(4) * t40;
t82 = g(3) * t46;
t39 = sin(t44);
t80 = t39 * t46;
t79 = t40 * t46;
t51 = -pkin(9) - pkin(8);
t78 = t46 * t51;
t77 = t47 * t45;
t75 = t47 * t49;
t74 = t50 * t39;
t73 = t50 * t40;
t71 = t50 * t48;
t70 = qJ(5) * t39;
t69 = qJ(6) + t51;
t67 = t89 * t39;
t65 = t69 * t50;
t15 = t39 * t75 + t73;
t16 = t40 * t75 - t74;
t64 = -t15 * pkin(4) + t16 * qJ(5);
t17 = -t47 * t40 + t49 * t74;
t18 = t47 * t39 + t49 * t73;
t63 = -t17 * pkin(4) + t18 * qJ(5);
t62 = -t38 - t70;
t61 = g(3) * (t29 + (t70 + t87) * t49);
t5 = g(1) * t15 - g(2) * t17;
t21 = t45 * t75 + t71;
t57 = t60 * t49;
t56 = pkin(3) * t72 - t16 * pkin(4) + t50 * pkin(7) - t15 * qJ(5) + t47 * t78;
t54 = pkin(3) * t77 + t18 * pkin(4) + t47 * pkin(7) + t17 * qJ(5) + t90 * t50;
t2 = g(1) * t17 + g(2) * t15 + g(3) * t80;
t4 = g(1) * t18 + g(2) * t16 + g(3) * t79;
t53 = t23 * pkin(3) + t63;
t52 = -t21 * pkin(3) + t64;
t27 = qJ(5) * t79;
t25 = t59 * t46;
t24 = t49 * t71 + t77;
t22 = -t48 * t75 + t72;
t20 = t57 + t82;
t12 = t17 * pkin(5);
t9 = t15 * pkin(5);
t8 = t19 * t40;
t7 = t19 * t39;
t6 = g(1) * t16 - g(2) * t18;
t1 = [0, t59, t60, 0, 0, 0, 0, 0, t59 * t49, -t25, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t24, -g(1) * t21 - g(2) * t23, 0, 0, 0, 0, 0, t6, -t5, t6, t25, t5, -g(1) * (-t47 * t90 + t56) - g(2) * (-t50 * t78 + t54) t6, t5, -t25, -g(1) * (-t16 * pkin(5) + t56) - g(2) * (t18 * pkin(5) - t46 * t65 + t54) - (t46 * qJ(6) - t90) * t86; 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, t19 * t48, -t19 * t45, 0, 0, 0, 0, 0, t8, -t7, t8, -t20, t7, -t61 + t51 * t57 + (g(3) * t51 + t60 * (-t62 + t87)) * t46, t8, t7, t20, -t61 + (-g(3) * pkin(5) * t40 + g(1) * t65 + t69 * t84) * t49 + (g(3) * t69 + t60 * (-t89 * t40 - t62)) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t23 + g(2) * t21 + t45 * t82, g(1) * t24 - g(2) * t22 + t48 * t82, 0, 0, 0, 0, 0, t2, t4, t2, 0, -t4, -g(1) * t53 - g(2) * t52 - g(3) * (t27 + (-pkin(4) * t39 - t88) * t46) t2, -t4, 0, -g(1) * (-t12 + t53) - g(2) * (-t9 + t52) - g(3) * t27 - (t67 - t88) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t4, t2, 0, -t4, -g(1) * t63 - g(2) * t64 - g(3) * (-pkin(4) * t80 + t27) t2, -t4, 0, -g(1) * (-t12 + t63) - g(2) * (t64 - t9) - g(3) * (t46 * t67 + t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19;];
taug_reg  = t1;
