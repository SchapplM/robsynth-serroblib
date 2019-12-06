% Calculate inertial parameters regressor of gravitation load for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t48 = sin(qJ(3));
t51 = cos(qJ(3));
t88 = -pkin(3) * t51 - qJ(4) * t48;
t44 = sin(pkin(10));
t86 = t44 * t51;
t45 = sin(pkin(5));
t49 = sin(qJ(2));
t85 = t45 * t49;
t52 = cos(qJ(2));
t84 = t45 * t52;
t46 = cos(pkin(10));
t47 = sin(qJ(5));
t83 = t46 * t47;
t50 = cos(qJ(5));
t82 = t46 * t50;
t81 = t46 * t51;
t80 = t47 * t48;
t79 = t48 * t50;
t78 = t51 * t52;
t77 = pkin(2) * t84 + pkin(7) * t85;
t75 = cos(pkin(5));
t74 = cos(pkin(9));
t73 = sin(pkin(9));
t72 = t48 * t84;
t71 = t44 * t84;
t59 = t75 * t74;
t28 = t73 * t49 - t52 * t59;
t29 = t49 * t59 + t73 * t52;
t70 = -t28 * pkin(2) + t29 * pkin(7);
t58 = t75 * t73;
t30 = t74 * t49 + t52 * t58;
t31 = -t49 * t58 + t74 * t52;
t69 = -t30 * pkin(2) + t31 * pkin(7);
t68 = t45 * t74;
t67 = t45 * t73;
t13 = t29 * t48 + t51 * t68;
t14 = t29 * t51 - t48 * t68;
t66 = -t13 * pkin(3) + t14 * qJ(4);
t15 = t31 * t48 - t51 * t67;
t16 = t31 * t51 + t48 * t67;
t65 = -t15 * pkin(3) + t16 * qJ(4);
t32 = t48 * t85 - t75 * t51;
t33 = t75 * t48 + t51 * t85;
t64 = -t32 * pkin(3) + t33 * qJ(4);
t63 = t45 * pkin(3) * t78 + qJ(4) * t72 + t77;
t62 = -pkin(4) * t46 - pkin(8) * t44;
t61 = t88 * t28 + t70;
t60 = t88 * t30 + t69;
t17 = -t46 * t85 + t51 * t71;
t6 = -t28 * t86 - t29 * t46;
t8 = -t30 * t86 - t31 * t46;
t57 = g(1) * t8 + g(2) * t6 + g(3) * t17;
t56 = g(1) * t15 + g(2) * t13 + g(3) * t32;
t55 = g(1) * t16 + g(2) * t14 + g(3) * t33;
t54 = -g(1) * t30 - g(2) * t28 + g(3) * t84;
t53 = g(1) * t31 + g(2) * t29 + g(3) * t85;
t18 = (t44 * t49 + t46 * t78) * t45;
t12 = t33 * t46 - t71;
t9 = -t30 * t81 + t31 * t44;
t7 = -t28 * t81 + t29 * t44;
t5 = t54 * t48;
t4 = t16 * t46 + t30 * t44;
t3 = t14 * t46 + t28 * t44;
t1 = t56 * t44;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t53, 0, 0, 0, 0, 0, 0, 0, 0, -t54 * t51, t5, -t53, -g(1) * t69 - g(2) * t70 - g(3) * t77, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7 - g(3) * t18, t57, -t5, -g(1) * t60 - g(2) * t61 - g(3) * t63, 0, 0, 0, 0, 0, 0, -g(1) * (-t30 * t80 + t9 * t50) - g(2) * (-t28 * t80 + t7 * t50) - g(3) * (t18 * t50 + t47 * t72), -g(1) * (-t30 * t79 - t9 * t47) - g(2) * (-t28 * t79 - t7 * t47) - g(3) * (-t18 * t47 + t50 * t72), -t57, -g(1) * (t9 * pkin(4) + t8 * pkin(8) + t60) - g(2) * (t7 * pkin(4) + t6 * pkin(8) + t61) - g(3) * (t18 * pkin(4) + t17 * pkin(8) + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t55, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t46, -t1, -t55, -g(1) * t65 - g(2) * t66 - g(3) * t64, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t82 + t16 * t47) - g(2) * (-t13 * t82 + t14 * t47) - g(3) * (-t32 * t82 + t33 * t47), -g(1) * (t15 * t83 + t16 * t50) - g(2) * (t13 * t83 + t14 * t50) - g(3) * (t32 * t83 + t33 * t50), t1, -g(1) * (t62 * t15 + t65) - g(2) * (t62 * t13 + t66) - g(3) * (t62 * t32 + t64); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t15 * t50 - t4 * t47) - g(2) * (t13 * t50 - t3 * t47) - g(3) * (-t12 * t47 + t32 * t50), -g(1) * (-t15 * t47 - t4 * t50) - g(2) * (-t13 * t47 - t3 * t50) - g(3) * (-t12 * t50 - t32 * t47), 0, 0;];
taug_reg = t2;
