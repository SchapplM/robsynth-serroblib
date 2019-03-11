% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t69 = cos(pkin(10));
t70 = cos(pkin(6));
t55 = t70 * t69;
t68 = sin(pkin(10));
t26 = t46 * t68 - t49 * t55;
t45 = sin(qJ(3));
t71 = qJ(4) * t45;
t48 = cos(qJ(3));
t81 = t26 * t48;
t86 = -pkin(3) * t81 - t26 * t71;
t54 = t70 * t68;
t28 = t46 * t69 + t49 * t54;
t80 = t28 * t48;
t85 = -pkin(3) * t80 - t28 * t71;
t84 = pkin(4) + pkin(8);
t42 = sin(pkin(6));
t83 = g(3) * t42;
t47 = cos(qJ(5));
t41 = t47 * pkin(5) + pkin(4);
t82 = pkin(8) + t41;
t79 = t42 * t46;
t78 = t42 * t49;
t44 = sin(qJ(5));
t77 = t44 * t45;
t76 = t44 * t49;
t75 = t45 * t47;
t74 = t47 * t49;
t73 = t48 * t49;
t72 = pkin(2) * t78 + pkin(8) * t79;
t23 = t26 * pkin(2);
t67 = -t23 + t86;
t24 = t28 * pkin(2);
t66 = -t24 + t85;
t27 = t46 * t55 + t49 * t68;
t65 = t27 * pkin(8) - t23;
t29 = -t46 * t54 + t49 * t69;
t64 = t29 * pkin(8) - t24;
t63 = pkin(5) * t44 + qJ(4);
t62 = t42 * t69;
t61 = t42 * t68;
t15 = t27 * t45 + t48 * t62;
t13 = t15 * pkin(3);
t16 = t27 * t48 - t45 * t62;
t60 = t16 * qJ(4) - t13;
t17 = t29 * t45 - t48 * t61;
t14 = t17 * pkin(3);
t18 = t29 * t48 + t45 * t61;
t59 = t18 * qJ(4) - t14;
t30 = t45 * t79 - t48 * t70;
t25 = t30 * pkin(3);
t31 = t45 * t70 + t48 * t79;
t58 = t31 * qJ(4) - t25;
t57 = t42 * pkin(3) * t73 + t71 * t78 + t72;
t56 = g(3) * t57;
t43 = -qJ(6) - pkin(9);
t53 = pkin(5) * t77 - t43 * t48;
t8 = g(1) * t17 + g(2) * t15 + g(3) * t30;
t52 = g(1) * t18 + g(2) * t16 + g(3) * t31;
t51 = -g(1) * t28 - g(2) * t26 + g(3) * t78;
t50 = g(1) * t29 + g(2) * t27 + g(3) * t79;
t1 = -g(1) * (t17 * t47 - t28 * t44) - g(2) * (t15 * t47 - t26 * t44) - g(3) * (t30 * t47 + t42 * t76);
t11 = t51 * t48;
t10 = t51 * t45;
t6 = t52 * t47;
t5 = t52 * t44;
t4 = -g(1) * (-t28 * t77 + t29 * t47) - g(2) * (-t26 * t77 + t27 * t47) - (t45 * t76 + t46 * t47) * t83;
t3 = -g(1) * (-t28 * t75 - t29 * t44) - g(2) * (-t26 * t75 - t27 * t44) - (-t44 * t46 + t45 * t74) * t83;
t2 = -g(1) * (-t17 * t44 - t28 * t47) - g(2) * (-t15 * t44 - t26 * t47) - g(3) * (-t30 * t44 + t42 * t74);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, -t50, -g(1) * t64 - g(2) * t65 - g(3) * t72, 0, 0, 0, 0, 0, 0, -t50, t11, -t10, -g(1) * (t64 + t85) - g(2) * (t65 + t86) - t56, 0, 0, 0, 0, 0, 0, t4, t3, -t11, -g(1) * (-pkin(9) * t80 + t29 * t84 + t66) - g(2) * (-pkin(9) * t81 + t27 * t84 + t67) - g(3) * ((pkin(4) * t46 + pkin(9) * t73) * t42 + t57) 0, 0, 0, 0, 0, 0, t4, t3, -t11, -g(1) * (-t28 * t53 + t29 * t82 + t66) - g(2) * (-t26 * t53 + t27 * t82 + t67) - t56 - (t41 * t46 + t49 * t53) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t52, -g(1) * t59 - g(2) * t60 - g(3) * t58, 0, 0, 0, 0, 0, 0, -t5, -t6, t8, -g(1) * (-t17 * pkin(9) + t59) - g(2) * (-t15 * pkin(9) + t60) - g(3) * (-t30 * pkin(9) + t58) 0, 0, 0, 0, 0, 0, -t5, -t6, t8, -g(1) * (t17 * t43 + t18 * t63 - t14) - g(2) * (t15 * t43 + t16 * t63 - t13) - g(3) * (t30 * t43 + t31 * t63 - t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52;];
taug_reg  = t7;
