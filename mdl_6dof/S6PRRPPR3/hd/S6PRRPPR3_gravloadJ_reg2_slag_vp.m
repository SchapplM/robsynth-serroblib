% Calculate inertial parameters regressor of gravitation load for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t43 = cos(pkin(10));
t67 = cos(pkin(6));
t61 = t43 * t67;
t66 = sin(pkin(10));
t27 = -t66 * t46 + t49 * t61;
t45 = sin(qJ(3));
t69 = qJ(4) * t45;
t48 = cos(qJ(3));
t81 = t27 * t48;
t84 = pkin(3) * t81 + t27 * t69;
t54 = t67 * t66;
t29 = -t43 * t46 - t49 * t54;
t80 = t29 * t48;
t83 = pkin(3) * t80 + t29 * t69;
t42 = sin(pkin(6));
t82 = g(3) * t42;
t79 = t42 * t46;
t78 = t42 * t48;
t77 = t42 * t49;
t44 = sin(qJ(6));
t76 = t44 * t45;
t75 = t44 * t49;
t47 = cos(qJ(6));
t74 = t45 * t47;
t73 = t47 * t49;
t72 = pkin(8) - qJ(5);
t71 = qJ(4) + pkin(5);
t70 = pkin(2) * t77 + pkin(8) * t79;
t68 = qJ(5) * t46;
t65 = t48 * t77;
t23 = t27 * pkin(2);
t28 = t46 * t61 + t66 * t49;
t64 = t28 * pkin(8) + t23;
t24 = t29 * pkin(2);
t30 = t43 * t49 - t46 * t54;
t63 = t30 * pkin(8) + t24;
t62 = t42 * t66;
t13 = t28 * t45 + t43 * t78;
t10 = t13 * pkin(3);
t14 = -t43 * t42 * t45 + t28 * t48;
t60 = t14 * qJ(4) - t10;
t15 = t30 * t45 - t48 * t62;
t12 = t15 * pkin(3);
t16 = t30 * t48 + t45 * t62;
t59 = t16 * qJ(4) - t12;
t31 = t45 * t79 - t67 * t48;
t26 = t31 * pkin(3);
t32 = t67 * t45 + t46 * t78;
t58 = t32 * qJ(4) - t26;
t57 = pkin(3) * t65 + t69 * t77 + t70;
t56 = pkin(4) * t65 + t57;
t55 = pkin(5) * t45 + pkin(9) * t48;
t2 = g(1) * t15 + g(2) * t13 + g(3) * t31;
t53 = g(1) * t16 + g(2) * t14 + g(3) * t32;
t52 = pkin(4) * t81 + t72 * t28 + t23 + t84;
t51 = pkin(4) * t80 + t72 * t30 + t24 + t83;
t50 = g(1) * t29 + g(2) * t27 + g(3) * t77;
t8 = g(1) * t30 + g(2) * t28 + g(3) * t79;
t25 = t31 * pkin(4);
t11 = t15 * pkin(4);
t9 = t13 * pkin(4);
t5 = t50 * t48;
t4 = t50 * t45;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t8, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, -t8, -g(1) * t63 - g(2) * t64 - g(3) * t70, 0, 0, 0, 0, 0, 0, -t5, -t8, -t4, -g(1) * (t63 + t83) - g(2) * (t64 + t84) - g(3) * t57, 0, 0, 0, 0, 0, 0, -t4, t5, t8, -g(1) * t51 - g(2) * t52 - g(3) * (-t42 * t68 + t56) 0, 0, 0, 0, 0, 0, -g(1) * (t29 * t74 - t30 * t44) - g(2) * (t27 * t74 - t28 * t44) - (-t44 * t46 + t45 * t73) * t82, -g(1) * (-t29 * t76 - t30 * t47) - g(2) * (-t27 * t76 - t28 * t47) - (-t45 * t75 - t46 * t47) * t82, -t5, -g(1) * (t55 * t29 + t51) - g(2) * (t55 * t27 + t52) - g(3) * t56 - (t55 * t49 - t68) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t53, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t53, -g(1) * t59 - g(2) * t60 - g(3) * t58, 0, 0, 0, 0, 0, 0, -t53, -t2, 0, -g(1) * (-t11 + t59) - g(2) * (t60 - t9) - g(3) * (-t25 + t58) 0, 0, 0, 0, 0, 0, -t53 * t47, t53 * t44, t2, -g(1) * (-t15 * pkin(9) + t71 * t16 - t11 - t12) - g(2) * (-t13 * pkin(9) + t71 * t14 - t10 - t9) - g(3) * (-t31 * pkin(9) + t71 * t32 - t25 - t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t44 + t29 * t47) - g(2) * (-t13 * t44 + t27 * t47) - g(3) * (-t31 * t44 + t42 * t73) -g(1) * (-t15 * t47 - t29 * t44) - g(2) * (-t13 * t47 - t27 * t44) - g(3) * (-t31 * t47 - t42 * t75) 0, 0;];
taug_reg  = t1;
