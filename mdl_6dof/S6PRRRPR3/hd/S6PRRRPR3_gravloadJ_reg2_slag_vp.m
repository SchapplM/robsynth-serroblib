% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t48 = sin(qJ(3));
t51 = cos(qJ(3));
t76 = cos(pkin(6));
t46 = sin(pkin(6));
t49 = sin(qJ(2));
t84 = t46 * t49;
t97 = -t48 * t84 + t76 * t51;
t52 = cos(qJ(2));
t45 = sin(pkin(11));
t72 = t45 * t76;
t75 = cos(pkin(11));
t30 = -t49 * t72 + t75 * t52;
t83 = t46 * t51;
t96 = -t30 * t48 + t45 * t83;
t44 = qJ(3) + qJ(4);
t43 = cos(t44);
t42 = sin(t44);
t77 = qJ(5) * t42;
t82 = t46 * t52;
t95 = (pkin(4) * t43 + t77) * t82;
t94 = g(3) * t46;
t61 = t76 * t75;
t28 = t45 * t52 + t49 * t61;
t71 = t46 * t75;
t12 = t28 * t42 + t43 * t71;
t93 = t12 * pkin(10);
t85 = t45 * t46;
t14 = t30 * t42 - t43 * t85;
t92 = t14 * pkin(10);
t23 = t42 * t84 - t76 * t43;
t91 = t23 * pkin(10);
t27 = t45 * t49 - t52 * t61;
t90 = t27 * t43;
t29 = t75 * t49 + t52 * t72;
t89 = t29 * t43;
t47 = sin(qJ(6));
t87 = t42 * t47;
t50 = cos(qJ(6));
t86 = t42 * t50;
t81 = t47 * t52;
t80 = t50 * t52;
t41 = t51 * pkin(3) + pkin(2);
t53 = -pkin(9) - pkin(8);
t79 = -t27 * t41 - t28 * t53;
t78 = -t29 * t41 - t30 * t53;
t13 = t28 * t43 - t42 * t71;
t69 = -t12 * pkin(4) + t13 * qJ(5);
t15 = t30 * t43 + t42 * t85;
t68 = -t14 * pkin(4) + t15 * qJ(5);
t24 = t76 * t42 + t43 * t84;
t67 = -t23 * pkin(4) + t24 * qJ(5);
t66 = -pkin(4) * t90 - t27 * t77 + t79;
t65 = -pkin(4) * t89 - t29 * t77 + t78;
t64 = t96 * pkin(3);
t32 = t41 * t82;
t63 = -t53 * t84 + t32;
t62 = t97 * pkin(3);
t4 = g(1) * t14 + g(2) * t12 + g(3) * t23;
t6 = g(1) * t15 + g(2) * t13 + g(3) * t24;
t60 = -t28 * t48 - t51 * t71;
t59 = t64 + t68;
t58 = -g(1) * t29 - g(2) * t27 + g(3) * t82;
t57 = g(1) * t30 + g(2) * t28 + g(3) * t84;
t56 = t62 + t67;
t55 = t60 * pkin(3);
t54 = t55 + t69;
t8 = t58 * t43;
t7 = t58 * t42;
t2 = t6 * t50;
t1 = t6 * t47;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t57, 0, 0, 0, 0, 0, 0, 0, 0, -t58 * t51, t58 * t48, -t57, -g(1) * (-t29 * pkin(2) + t30 * pkin(8)) - g(2) * (-t27 * pkin(2) + t28 * pkin(8)) - (pkin(2) * t52 + pkin(8) * t49) * t94, 0, 0, 0, 0, 0, 0, -t8, t7, -t57, -g(1) * t78 - g(2) * t79 - g(3) * t63, 0, 0, 0, 0, 0, 0, -t57, t8, -t7, -g(1) * t65 - g(2) * t66 - g(3) * (t63 + t95) 0, 0, 0, 0, 0, 0, -g(1) * (-t29 * t87 + t30 * t50) - g(2) * (-t27 * t87 + t28 * t50) - (t42 * t81 + t49 * t50) * t94, -g(1) * (-t29 * t86 - t30 * t47) - g(2) * (-t27 * t86 - t28 * t47) - (t42 * t80 - t47 * t49) * t94, -t8, -g(1) * (t30 * pkin(5) - pkin(10) * t89 + t65) - g(2) * (t28 * pkin(5) - pkin(10) * t90 + t66) - g(3) * (t32 + t95) - (pkin(10) * t43 * t52 + (pkin(5) - t53) * t49) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t60 - g(3) * t97, -g(1) * (-t30 * t51 - t48 * t85) - g(2) * (-t28 * t51 + t48 * t71) - g(3) * (-t76 * t48 - t49 * t83) 0, 0, 0, 0, 0, 0, 0, 0, t4, t6, 0, -g(1) * t64 - g(2) * t55 - g(3) * t62, 0, 0, 0, 0, 0, 0, 0, -t4, -t6, -g(1) * t59 - g(2) * t54 - g(3) * t56, 0, 0, 0, 0, 0, 0, -t1, -t2, t4, -g(1) * (t59 - t92) - g(2) * (t54 - t93) - g(3) * (t56 - t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t6, -g(1) * t68 - g(2) * t69 - g(3) * t67, 0, 0, 0, 0, 0, 0, -t1, -t2, t4, -g(1) * (t68 - t92) - g(2) * (t69 - t93) - g(3) * (t67 - t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t50 - t29 * t47) - g(2) * (t12 * t50 - t27 * t47) - g(3) * (t23 * t50 + t46 * t81) -g(1) * (-t14 * t47 - t29 * t50) - g(2) * (-t12 * t47 - t27 * t50) - g(3) * (-t23 * t47 + t46 * t80) 0, 0;];
taug_reg  = t3;
