% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t45 = sin(qJ(2));
t46 = sin(qJ(1));
t49 = cos(qJ(2));
t72 = cos(pkin(5));
t85 = cos(qJ(1));
t59 = t72 * t85;
t27 = t45 * t59 + t46 * t49;
t44 = sin(qJ(3));
t48 = cos(qJ(3));
t41 = sin(pkin(5));
t68 = t41 * t85;
t15 = t27 * t48 - t44 * t68;
t26 = t45 * t46 - t49 * t59;
t43 = sin(qJ(4));
t47 = cos(qJ(4));
t89 = t15 * t43 - t26 * t47;
t88 = t15 * t47 + t26 * t43;
t64 = t46 * t72;
t29 = -t45 * t64 + t49 * t85;
t81 = t41 * t46;
t19 = t29 * t48 + t44 * t81;
t28 = t45 * t85 + t49 * t64;
t11 = -t19 * t43 + t28 * t47;
t80 = t41 * t48;
t25 = t44 * t72 + t45 * t80;
t79 = t41 * t49;
t1 = g(2) * t89 - g(3) * (-t25 * t43 - t47 * t79) - g(1) * t11;
t86 = g(3) * t41;
t82 = t41 * t45;
t78 = t43 * t45;
t77 = t43 * t48;
t76 = t47 * t48;
t75 = t48 * t49;
t74 = pkin(2) * t79 + pkin(8) * t82;
t73 = t85 * pkin(1) + pkin(7) * t81;
t71 = pkin(2) * t29 + t73;
t70 = pkin(4) * t43 + pkin(8);
t69 = g(3) * t74;
t67 = -t46 * pkin(1) + pkin(7) * t68;
t20 = t26 * pkin(2);
t66 = t27 * pkin(8) - t20;
t22 = t28 * pkin(2);
t65 = t29 * pkin(8) - t22;
t63 = -pkin(2) * t27 + t67;
t62 = pkin(3) * t48 + pkin(9) * t44;
t14 = t27 * t44 + t48 * t68;
t18 = t29 * t44 - t46 * t80;
t61 = -g(1) * t14 + g(2) * t18;
t60 = g(1) * t26 - g(2) * t28;
t39 = pkin(4) * t47 + pkin(3);
t42 = -qJ(5) - pkin(9);
t58 = t39 * t48 - t42 * t44;
t57 = pkin(8) * t28 + t71;
t56 = g(1) * t85 + g(2) * t46;
t55 = -pkin(8) * t26 + t63;
t24 = t44 * t82 - t48 * t72;
t54 = g(1) * t18 + g(2) * t14 + g(3) * t24;
t53 = g(1) * t19 + g(2) * t15 + g(3) * t25;
t52 = -g(1) * t28 - g(2) * t26 + g(3) * t79;
t51 = g(1) * t29 + g(2) * t27 + g(3) * t82;
t13 = t52 * t44;
t12 = t19 * t47 + t28 * t43;
t8 = t54 * t47;
t7 = t54 * t43;
t6 = g(1) * t88 - g(2) * t12;
t5 = -g(1) * t89 - g(2) * t11;
t4 = -g(1) * (-t28 * t76 + t29 * t43) - g(2) * (-t26 * t76 + t27 * t43) - (t47 * t75 + t78) * t86;
t3 = -g(1) * (t28 * t77 + t29 * t47) - g(2) * (t26 * t77 + t27 * t47) - (-t43 * t75 + t45 * t47) * t86;
t2 = g(1) * t12 + g(2) * t88 - g(3) * (-t25 * t47 + t43 * t79);
t9 = [0, 0, 0, 0, 0, 0, g(1) * t46 - g(2) * t85, t56, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t27 - g(2) * t29, -t60, -t56 * t41, -g(1) * t67 - g(2) * t73, 0, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t19, t61, t60, -g(1) * t55 - g(2) * t57, 0, 0, 0, 0, 0, 0, t6, t5, -t61, -g(1) * (-pkin(3) * t15 - pkin(9) * t14 + t55) - g(2) * (pkin(3) * t19 + pkin(9) * t18 + t57), 0, 0, 0, 0, 0, 0, t6, t5, -t61, -g(1) * (t14 * t42 - t15 * t39 - t26 * t70 + t63) - g(2) * (-t18 * t42 + t19 * t39 + t28 * t70 + t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t51, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * t48, t13, -t51, -g(1) * t65 - g(2) * t66 - t69, 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(1) * (-t28 * t62 + t65) - g(2) * (-t26 * t62 + t66) - g(3) * (t62 * t79 + t74), 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(1) * (-t28 * t58 + t29 * t70 - t22) - g(2) * (-t26 * t58 + t27 * t70 - t20) - t69 - (pkin(4) * t78 + t49 * t58) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t53, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t53, -g(1) * (-pkin(3) * t18 + pkin(9) * t19) - g(2) * (-pkin(3) * t14 + pkin(9) * t15) - g(3) * (-pkin(3) * t24 + pkin(9) * t25), 0, 0, 0, 0, 0, 0, t8, -t7, -t53, -g(1) * (-t18 * t39 - t19 * t42) - g(2) * (-t14 * t39 - t15 * t42) - g(3) * (-t24 * t39 - t25 * t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54;];
taug_reg = t9;
