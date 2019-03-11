% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t42 = sin(qJ(5));
t45 = cos(qJ(5));
t75 = pkin(5) * t45 + qJ(6) * t42;
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t25 = g(1) * t47 + g(2) * t44;
t40 = qJ(2) + pkin(10);
t36 = qJ(4) + t40;
t29 = sin(t36);
t52 = t25 * t29;
t30 = cos(t36);
t74 = t75 * t30;
t61 = t30 * pkin(4) + t29 * pkin(9);
t70 = g(3) * t29;
t69 = g(3) * t42;
t68 = t29 * t47;
t67 = t30 * t47;
t66 = t44 * t42;
t65 = t44 * t45;
t41 = -qJ(3) - pkin(7);
t39 = -pkin(8) + t41;
t64 = t47 * t39;
t63 = t47 * t42;
t62 = t47 * t45;
t35 = cos(t40);
t46 = cos(qJ(2));
t37 = t46 * pkin(2);
t60 = pkin(3) * t35 + t37;
t15 = pkin(1) + t60;
t12 = t47 * t15;
t58 = pkin(4) * t67 + pkin(9) * t68 + t12;
t57 = t60 + t61;
t18 = t44 * t30 * pkin(9);
t56 = -t44 * t29 * pkin(4) + t18;
t21 = pkin(9) * t67;
t55 = -pkin(4) * t68 + t21;
t10 = t30 * t63 - t65;
t8 = t30 * t66 + t62;
t54 = g(1) * t8 - g(2) * t10;
t24 = g(1) * t44 - g(2) * t47;
t1 = g(1) * t10 + g(2) * t8 + t29 * t69;
t11 = t30 * t62 + t66;
t9 = t30 * t65 - t63;
t51 = g(1) * t11 + g(2) * t9 + t45 * t70;
t5 = -g(3) * t30 + t52;
t43 = sin(qJ(2));
t50 = -g(3) * t46 + t25 * t43;
t49 = (-g(1) * (-t15 - t61) + g(2) * t39) * t44;
t48 = (pkin(4) + t75) * t52;
t34 = sin(t40);
t33 = t37 + pkin(1);
t16 = -t43 * pkin(2) - pkin(3) * t34;
t14 = t47 * t16;
t13 = t44 * t16;
t7 = t24 * t29;
t6 = t25 * t30 + t70;
t4 = t5 * t45;
t3 = -t30 * t69 + t42 * t52;
t2 = g(1) * t9 - g(2) * t11;
t17 = [0, 0, 0, 0, 0, 0, t24, t25, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t46, -t24 * t43, -t25, -g(1) * (-t44 * pkin(1) + t47 * pkin(7)) - g(2) * (t47 * pkin(1) + t44 * pkin(7)) 0, 0, 0, 0, 0, 0, t24 * t35, -t24 * t34, -t25, -g(1) * (-t44 * t33 - t47 * t41) - g(2) * (t47 * t33 - t44 * t41) 0, 0, 0, 0, 0, 0, t24 * t30, -t7, -t25, -g(1) * (-t44 * t15 - t64) - g(2) * (-t44 * t39 + t12) 0, 0, 0, 0, 0, 0, t2, -t54, t7, g(1) * t64 - g(2) * t58 + t49, 0, 0, 0, 0, 0, 0, t2, t7, t54, -g(1) * (-t9 * pkin(5) - t8 * qJ(6) - t64) - g(2) * (t11 * pkin(5) + t10 * qJ(6) + t58) + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, g(3) * t43 + t25 * t46, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t35 + t25 * t34, g(3) * t34 + t25 * t35, 0, t50 * pkin(2), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t60 - t25 * t16, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t14 + t55) - g(2) * (t13 + t56) - g(3) * t57, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (t14 + t21) - g(2) * (t13 + t18) - g(3) * (t57 + t74) + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t55 - g(2) * t56 - g(3) * t61, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t21 - g(2) * t18 - g(3) * (t61 + t74) + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t51, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t51, -g(1) * (-t10 * pkin(5) + t11 * qJ(6)) - g(2) * (-t8 * pkin(5) + t9 * qJ(6)) - (-pkin(5) * t42 + qJ(6) * t45) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t17;
