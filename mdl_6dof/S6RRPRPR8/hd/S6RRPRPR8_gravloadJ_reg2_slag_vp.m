% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t41 = pkin(10) + qJ(4);
t35 = sin(t41);
t36 = cos(t41);
t50 = cos(qJ(1));
t80 = t50 * t36;
t47 = sin(qJ(1));
t49 = cos(qJ(2));
t82 = t47 * t49;
t15 = t35 * t82 + t80;
t81 = t50 * t35;
t16 = t36 * t82 - t81;
t45 = sin(qJ(6));
t48 = cos(qJ(6));
t105 = t15 * t48 - t16 * t45;
t46 = sin(qJ(2));
t106 = g(3) * t46;
t17 = -t47 * t36 + t49 * t81;
t18 = t47 * t35 + t49 * t80;
t3 = t17 * t48 - t18 * t45;
t58 = t35 * t48 - t36 * t45;
t108 = g(1) * t3 + g(2) * t105 + t58 * t106;
t92 = g(2) * t47;
t26 = g(1) * t50 + t92;
t107 = -g(3) * t49 + t26 * t46;
t103 = -t16 * pkin(4) - t15 * qJ(5);
t4 = t17 * t45 + t18 * t48;
t57 = t35 * t45 + t36 * t48;
t60 = t15 * t45 + t16 * t48;
t101 = g(1) * t4 + g(2) * t60 + t57 * t106;
t96 = -pkin(4) - pkin(5);
t95 = pkin(4) * t36;
t94 = g(1) * t47;
t44 = -pkin(8) - qJ(3);
t89 = pkin(9) + t44;
t87 = t35 * t46;
t86 = t36 * t46;
t85 = t46 * t47;
t84 = t46 * t50;
t42 = sin(pkin(10));
t83 = t47 * t42;
t43 = cos(pkin(10));
t34 = t43 * pkin(3) + pkin(2);
t25 = t49 * t34;
t79 = t50 * t42;
t78 = t50 * t43;
t77 = t50 * pkin(1) + t47 * pkin(7);
t76 = qJ(5) * t35;
t74 = t44 * t84;
t38 = t50 * pkin(7);
t73 = pkin(3) * t79 + t44 * t85 + t38;
t72 = t89 * t50;
t71 = -pkin(1) - t25;
t69 = -t15 * pkin(4) + t16 * qJ(5);
t68 = -t17 * pkin(4) + t18 * qJ(5);
t67 = -t34 - t76;
t66 = pkin(3) * t83 + t50 * t25 + t77;
t65 = g(3) * (t25 + (t76 + t95) * t49);
t64 = g(1) * t15 - g(2) * t17;
t63 = -g(2) * t50 + t94;
t62 = t49 * pkin(2) + t46 * qJ(3);
t55 = t26 * t49;
t53 = t18 * pkin(4) + t17 * qJ(5) + t66;
t2 = g(1) * t17 + g(2) * t15 + g(3) * t87;
t52 = g(1) * t18 + g(2) * t16 + g(3) * t86;
t51 = t71 * t47 + t73;
t23 = qJ(5) * t86;
t21 = g(1) * t85 - g(2) * t84;
t20 = t55 + t106;
t7 = t107 * t36;
t6 = t107 * t35;
t5 = g(1) * t16 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, t63, t26, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t49, -t21, -t26, -g(1) * (-t47 * pkin(1) + t38) - g(2) * t77, 0, 0, 0, 0, 0, 0, -g(1) * (-t43 * t82 + t79) - g(2) * (t49 * t78 + t83) -g(1) * (t42 * t82 + t78) - g(2) * (t47 * t43 - t49 * t79) t21, -g(1) * t38 - g(2) * (t62 * t50 + t77) - (-pkin(1) - t62) * t94, 0, 0, 0, 0, 0, 0, t5, -t64, t21, -g(1) * t51 - g(2) * (t66 - t74) 0, 0, 0, 0, 0, 0, t5, t21, t64, -g(1) * (t51 + t103) - g(2) * (t53 - t74) 0, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t4, g(1) * t105 - g(2) * t3, -t21, -g(1) * (-t16 * pkin(5) + t103 + t73) - g(2) * (t18 * pkin(5) - t46 * t72 + t53) - (t46 * pkin(9) + t71) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t20, 0, 0, 0, 0, 0, 0, 0, 0, t107 * t43, -t107 * t42, -t20, -g(3) * t62 + t26 * (pkin(2) * t46 - qJ(3) * t49) 0, 0, 0, 0, 0, 0, t7, -t6, -t20, -g(3) * (-t46 * t44 + t25) + t26 * (t34 * t46 + t44 * t49) 0, 0, 0, 0, 0, 0, t7, -t20, t6, -t65 + t44 * t55 + (g(3) * t44 + t26 * (-t67 + t95)) * t46, 0, 0, 0, 0, 0, 0, t107 * t57, t107 * t58, t20, -t65 + (-g(3) * pkin(5) * t36 + g(1) * t72 + t89 * t92) * t49 + (g(3) * t89 + t26 * (-t96 * t36 - t67)) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t52, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t52, -g(1) * t68 - g(2) * t69 - g(3) * (-pkin(4) * t87 + t23) 0, 0, 0, 0, 0, 0, t108, -t101, 0, -g(1) * (-t17 * pkin(5) + t68) - g(2) * (-t15 * pkin(5) + t69) - g(3) * (t96 * t87 + t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, t101, 0, 0;];
taug_reg  = t1;
