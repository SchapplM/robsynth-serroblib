% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP11
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
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t48 = sin(qJ(2));
t49 = sin(qJ(1));
t52 = cos(qJ(2));
t67 = cos(pkin(5));
t77 = cos(qJ(1));
t60 = t67 * t77;
t35 = t48 * t60 + t49 * t52;
t47 = sin(qJ(3));
t51 = cos(qJ(3));
t45 = sin(pkin(5));
t65 = t45 * t77;
t20 = t35 * t51 - t47 * t65;
t34 = t49 * t48 - t52 * t60;
t46 = sin(qJ(4));
t50 = cos(qJ(4));
t6 = t20 * t46 - t34 * t50;
t7 = t20 * t50 + t34 * t46;
t63 = t49 * t67;
t36 = t77 * t48 + t52 * t63;
t83 = g(1) * t36 + g(2) * t34;
t37 = -t48 * t63 + t77 * t52;
t72 = t45 * t51;
t23 = t37 * t47 - t49 * t72;
t64 = -t35 * t47 - t51 * t65;
t74 = t45 * t48;
t55 = -g(3) * (-t47 * t74 + t67 * t51) - g(2) * t64 + g(1) * t23;
t73 = t45 * t49;
t71 = t45 * t52;
t70 = t46 * t51;
t69 = t50 * t51;
t68 = t50 * t52;
t66 = t46 * t71;
t24 = t37 * t51 + t47 * t73;
t10 = t24 * t46 - t36 * t50;
t62 = -g(1) * t6 + g(2) * t10;
t61 = g(1) * t64 + g(2) * t23;
t59 = pkin(3) * t51 + pkin(9) * t47 + pkin(2);
t33 = t67 * t47 + t48 * t72;
t17 = t33 * t46 + t45 * t68;
t1 = g(1) * t10 + g(2) * t6 + g(3) * t17;
t11 = t24 * t50 + t36 * t46;
t18 = t33 * t50 - t66;
t57 = g(1) * t11 + g(2) * t7 + g(3) * t18;
t13 = -t34 * t70 - t35 * t50;
t15 = -t36 * t70 - t37 * t50;
t25 = -t50 * t74 + t51 * t66;
t56 = g(1) * t15 + g(2) * t13 + g(3) * t25;
t54 = g(1) * t24 + g(2) * t20 + g(3) * t33;
t53 = g(3) * t71 - t83;
t26 = (t46 * t48 + t51 * t68) * t45;
t16 = -t36 * t69 + t37 * t46;
t14 = -t34 * t69 + t35 * t46;
t12 = t53 * t47;
t5 = t55 * t50;
t4 = t55 * t46;
t3 = g(1) * t7 - g(2) * t11;
t2 = -g(1) * t16 - g(2) * t14 - g(3) * t26;
t8 = [0, g(1) * t49 - g(2) * t77, g(1) * t77 + g(2) * t49, 0, 0, 0, 0, 0, g(1) * t35 - g(2) * t37, -g(1) * t34 + g(2) * t36, 0, 0, 0, 0, 0, g(1) * t20 - g(2) * t24, t61, 0, 0, 0, 0, 0, t3, t62, t3, -t61, -t62, -g(1) * (-t49 * pkin(1) - t35 * pkin(2) - pkin(3) * t20 - pkin(4) * t7 + pkin(7) * t65 - t34 * pkin(8) + pkin(9) * t64 - qJ(5) * t6) - g(2) * (t77 * pkin(1) + t37 * pkin(2) + t24 * pkin(3) + t11 * pkin(4) + pkin(7) * t73 + t36 * pkin(8) + t23 * pkin(9) + t10 * qJ(5)); 0, 0, 0, 0, 0, 0, 0, 0, -t53, g(1) * t37 + g(2) * t35 + g(3) * t74, 0, 0, 0, 0, 0, -t53 * t51, t12, 0, 0, 0, 0, 0, t2, t56, t2, -t12, -t56, -g(1) * (t16 * pkin(4) + t37 * pkin(8) + t15 * qJ(5)) - g(2) * (t14 * pkin(4) + t35 * pkin(8) + t13 * qJ(5)) + t83 * t59 + (-t26 * pkin(4) - t25 * qJ(5) - (pkin(8) * t48 + t59 * t52) * t45) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t54, 0, 0, 0, 0, 0, t5, -t4, t5, -t54, t4, -t54 * pkin(9) + t55 * (pkin(4) * t50 + qJ(5) * t46 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t57, t1, 0, -t57, -g(1) * (-t10 * pkin(4) + t11 * qJ(5)) - g(2) * (-t6 * pkin(4) + t7 * qJ(5)) - g(3) * (-t17 * pkin(4) + t18 * qJ(5)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t8;
