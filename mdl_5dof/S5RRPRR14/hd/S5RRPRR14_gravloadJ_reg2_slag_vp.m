% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR14_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t47 = sin(qJ(2));
t48 = sin(qJ(1));
t50 = cos(qJ(2));
t66 = cos(pkin(5));
t80 = cos(qJ(1));
t56 = t66 * t80;
t23 = t48 * t47 - t50 * t56;
t46 = sin(qJ(5));
t49 = cos(qJ(5));
t24 = t47 * t56 + t48 * t50;
t41 = pkin(10) + qJ(4);
t38 = sin(t41);
t39 = cos(t41);
t43 = sin(pkin(5));
t64 = t43 * t80;
t8 = t24 * t39 - t38 * t64;
t83 = -t23 * t49 + t8 * t46;
t82 = t23 * t46 + t8 * t49;
t81 = g(3) * t43;
t77 = t39 * t46;
t76 = t39 * t49;
t75 = t43 * t47;
t74 = t43 * t48;
t73 = t43 * t50;
t45 = -pkin(8) - qJ(3);
t72 = t45 * t47;
t71 = t46 * t50;
t70 = t49 * t50;
t44 = cos(pkin(10));
t37 = t44 * pkin(3) + pkin(2);
t69 = -t23 * t37 - t24 * t45;
t61 = t48 * t66;
t25 = t80 * t47 + t50 * t61;
t26 = -t47 * t61 + t80 * t50;
t68 = -t25 * t37 - t26 * t45;
t67 = t80 * pkin(1) + pkin(7) * t74;
t42 = sin(pkin(10));
t65 = t42 * t74;
t63 = -t48 * pkin(1) + pkin(7) * t64;
t62 = -t24 * t38 - t39 * t64;
t60 = t42 * t64;
t11 = t26 * t38 - t39 * t74;
t59 = g(1) * t62 + g(2) * t11;
t58 = pkin(3) * t65 - t25 * t45 + t26 * t37 + t67;
t57 = pkin(4) * t39 + pkin(9) * t38;
t6 = g(1) * t23 - g(2) * t25;
t55 = g(1) * t80 + g(2) * t48;
t54 = pkin(3) * t60 + t23 * t45 - t24 * t37 + t63;
t17 = -t38 * t75 + t66 * t39;
t53 = g(1) * t11 - g(2) * t62 - g(3) * t17;
t12 = t26 * t39 + t38 * t74;
t18 = t66 * t38 + t39 * t75;
t52 = g(1) * t12 + g(2) * t8 + g(3) * t18;
t4 = -g(1) * t25 - g(2) * t23 + g(3) * t73;
t51 = g(1) * t26 + g(2) * t24 + g(3) * t75;
t27 = t37 * t73;
t3 = t12 * t49 + t25 * t46;
t2 = -t12 * t46 + t25 * t49;
t1 = t4 * t38;
t5 = [0, 0, 0, 0, 0, 0, g(1) * t48 - g(2) * t80, t55, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t24 - g(2) * t26, -t6, -t55 * t43, -g(1) * t63 - g(2) * t67, 0, 0, 0, 0, 0, 0, -g(1) * (-t24 * t44 + t60) - g(2) * (t26 * t44 + t65), -g(1) * (t24 * t42 + t44 * t64) - g(2) * (-t26 * t42 + t44 * t74), t6, -g(1) * (-t24 * pkin(2) - t23 * qJ(3) + t63) - g(2) * (t26 * pkin(2) + t25 * qJ(3) + t67), 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t12, t59, t6, -g(1) * t54 - g(2) * t58, 0, 0, 0, 0, 0, 0, g(1) * t82 - g(2) * t3, -g(1) * t83 - g(2) * t2, -t59, -g(1) * (-pkin(4) * t8 + pkin(9) * t62 + t54) - g(2) * (t12 * pkin(4) + t11 * pkin(9) + t58); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t51, 0, 0, 0, 0, 0, 0, 0, 0, -t4 * t44, t4 * t42, -t51, -g(1) * (-t25 * pkin(2) + t26 * qJ(3)) - g(2) * (-t23 * pkin(2) + t24 * qJ(3)) - (pkin(2) * t50 + qJ(3) * t47) * t81, 0, 0, 0, 0, 0, 0, -t4 * t39, t1, -t51, -g(1) * t68 - g(2) * t69 - g(3) * (-t43 * t72 + t27), 0, 0, 0, 0, 0, 0, -g(1) * (-t25 * t76 + t26 * t46) - g(2) * (-t23 * t76 + t24 * t46) - (t39 * t70 + t46 * t47) * t81, -g(1) * (t25 * t77 + t26 * t49) - g(2) * (t23 * t77 + t24 * t49) - (-t39 * t71 + t47 * t49) * t81, -t1, -g(1) * (-t57 * t25 + t68) - g(2) * (-t57 * t23 + t69) - g(3) * t27 - (t57 * t50 - t72) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t52, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t49, -t53 * t46, -t52, -g(1) * (-t11 * pkin(4) + t12 * pkin(9)) - g(2) * (pkin(4) * t62 + t8 * pkin(9)) - g(3) * (t17 * pkin(4) + t18 * pkin(9)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t83 - g(3) * (-t18 * t46 - t43 * t70), g(1) * t3 + g(2) * t82 - g(3) * (-t18 * t49 + t43 * t71), 0, 0;];
taug_reg = t5;
