% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPP6
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t44 = cos(qJ(3));
t34 = t44 * pkin(3) + pkin(2);
t45 = cos(qJ(2));
t26 = t45 * t34;
t88 = pkin(1) + t26;
t43 = sin(qJ(1));
t41 = sin(qJ(3));
t46 = cos(qJ(1));
t71 = t46 * t41;
t20 = t43 * t44 - t45 * t71;
t83 = g(2) * t43;
t57 = g(1) * t46 + t83;
t42 = sin(qJ(2));
t49 = -g(3) * t45 + t57 * t42;
t87 = pkin(3) * t41;
t40 = qJ(3) + qJ(4);
t36 = cos(t40);
t86 = pkin(4) * t36;
t85 = g(1) * t43;
t82 = g(3) * t42;
t47 = -pkin(9) - pkin(8);
t80 = pkin(5) - t47;
t35 = sin(t40);
t79 = t35 * t42;
t78 = t36 * t42;
t77 = t42 * t47;
t76 = t43 * t41;
t74 = t43 * t45;
t73 = t46 * t35;
t72 = t46 * t36;
t70 = t46 * t44;
t69 = -pkin(4) - qJ(6);
t68 = qJ(5) * t35;
t13 = t35 * t74 + t72;
t67 = t13 * qJ(6);
t15 = -t43 * t36 + t45 * t73;
t66 = t15 * qJ(6);
t64 = t80 * t46;
t62 = t69 * t35;
t14 = t36 * t74 - t73;
t61 = -t13 * pkin(4) + t14 * qJ(5);
t16 = t43 * t35 + t45 * t72;
t60 = -t15 * pkin(4) + t16 * qJ(5);
t59 = -t34 - t68;
t58 = g(3) * (t26 + (t68 + t86) * t45);
t5 = g(1) * t13 - g(2) * t15;
t6 = g(1) * t14 - g(2) * t16;
t56 = -g(2) * t46 + t85;
t18 = t41 * t74 + t70;
t54 = t57 * t45;
t53 = pkin(3) * t71 - t14 * pkin(4) + t46 * pkin(7) - t13 * qJ(5) + t43 * t77;
t52 = pkin(3) * t76 + t16 * pkin(4) + t43 * pkin(7) + t15 * qJ(5) + t88 * t46;
t2 = g(1) * t15 + g(2) * t13 + g(3) * t79;
t4 = g(1) * t16 + g(2) * t14 + g(3) * t78;
t50 = t20 * pkin(3) + t60;
t48 = -t18 * pkin(3) + t61;
t24 = qJ(5) * t78;
t22 = t56 * t42;
t21 = t45 * t70 + t76;
t19 = -t44 * t74 + t71;
t17 = t54 + t82;
t8 = t49 * t36;
t7 = t49 * t35;
t1 = [0, t56, t57, 0, 0, 0, 0, 0, t56 * t45, -t22, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t21, -g(1) * t18 - g(2) * t20, 0, 0, 0, 0, 0, t6, -t5, t22, -t6, t5, -g(1) * (-t43 * t88 + t53) - g(2) * (-t46 * t77 + t52) t22, t5, t6, -g(1) * (-t14 * qJ(6) + t53) - g(2) * (t16 * qJ(6) + t42 * t64 + t52) - (-t42 * pkin(5) - t88) * t85; 0, 0, 0, 0, 0, 0, 0, 0, t49, t17, 0, 0, 0, 0, 0, t49 * t44, -t49 * t41, 0, 0, 0, 0, 0, t8, -t7, -t17, -t8, t7, -t58 + t47 * t54 + (g(3) * t47 + t57 * (-t59 + t86)) * t42, -t17, t7, t8, -t58 + (-g(3) * qJ(6) * t36 - g(1) * t64 - t80 * t83) * t45 + (-g(3) * t80 + t57 * (-t69 * t36 - t59)) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t20 + g(2) * t18 + t41 * t82, g(1) * t21 - g(2) * t19 + t44 * t82, 0, 0, 0, 0, 0, t2, t4, 0, -t2, -t4, -g(1) * t50 - g(2) * t48 - g(3) * (t24 + (-pkin(4) * t35 - t87) * t42) 0, -t4, t2, -g(1) * (t50 - t66) - g(2) * (t48 - t67) - g(3) * t24 - (t62 - t87) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t4, 0, -t2, -t4, -g(1) * t60 - g(2) * t61 - g(3) * (-pkin(4) * t79 + t24) 0, -t4, t2, -g(1) * (t60 - t66) - g(2) * (t61 - t67) - g(3) * (t42 * t62 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4;];
taug_reg  = t1;
