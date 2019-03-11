% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t33 = sin(qJ(2));
t34 = sin(qJ(1));
t37 = cos(qJ(2));
t48 = cos(pkin(6));
t62 = cos(qJ(1));
t43 = t48 * t62;
t16 = t34 * t33 - t37 * t43;
t28 = qJ(4) + pkin(12) + qJ(6);
t25 = sin(t28);
t26 = cos(t28);
t17 = t33 * t43 + t34 * t37;
t32 = sin(qJ(3));
t36 = cos(qJ(3));
t29 = sin(pkin(6));
t46 = t29 * t62;
t9 = t17 * t36 - t32 * t46;
t71 = -t16 * t26 + t9 * t25;
t70 = t16 * t25 + t9 * t26;
t31 = sin(qJ(4));
t35 = cos(qJ(4));
t69 = -t16 * t35 + t9 * t31;
t68 = t16 * t31 + t9 * t35;
t53 = t29 * t36;
t15 = t48 * t32 + t33 * t53;
t45 = t34 * t48;
t19 = -t33 * t45 + t62 * t37;
t54 = t29 * t34;
t13 = t19 * t36 + t32 * t54;
t18 = t62 * t33 + t37 * t45;
t5 = -t13 * t31 + t18 * t35;
t52 = t29 * t37;
t67 = g(2) * t69 - g(3) * (-t15 * t31 - t35 * t52) - g(1) * t5;
t65 = g(2) * t16;
t64 = g(2) * t17;
t63 = g(3) * t29;
t57 = t25 * t36;
t56 = t26 * t36;
t55 = t29 * t33;
t51 = t31 * t36;
t50 = t35 * t36;
t49 = t36 * t37;
t47 = pkin(4) * t31 + pkin(9);
t12 = t19 * t32 - t34 * t53;
t8 = t17 * t32 + t36 * t46;
t44 = -g(1) * t8 + g(2) * t12;
t27 = t35 * pkin(4) + pkin(3);
t30 = -qJ(5) - pkin(10);
t42 = t27 * t36 - t30 * t32 + pkin(2);
t14 = t32 * t55 - t48 * t36;
t41 = g(1) * t12 + g(2) * t8 + g(3) * t14;
t40 = g(1) * t13 + g(2) * t9 + g(3) * t15;
t39 = -g(1) * t18 + g(3) * t52 - t65;
t7 = t39 * t32;
t6 = t13 * t35 + t18 * t31;
t4 = t13 * t26 + t18 * t25;
t3 = -t13 * t25 + t18 * t26;
t2 = g(1) * t4 + g(2) * t70 - g(3) * (-t15 * t26 + t25 * t52);
t1 = -g(1) * t3 + g(2) * t71 - g(3) * (-t15 * t25 - t26 * t52);
t10 = [0, g(1) * t34 - g(2) * t62, g(1) * t62 + g(2) * t34, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t19, -g(1) * t16 + g(2) * t18, 0, 0, 0, 0, 0, g(1) * t9 - g(2) * t13, t44, 0, 0, 0, 0, 0, g(1) * t68 - g(2) * t6, -g(1) * t69 - g(2) * t5, -t44, -g(1) * (-t34 * pkin(1) - t17 * pkin(2) + pkin(8) * t46 - t47 * t16 - t27 * t9 + t30 * t8) - g(2) * (t62 * pkin(1) + t19 * pkin(2) + pkin(8) * t54 - t12 * t30 + t13 * t27 + t47 * t18) 0, 0, 0, 0, 0, g(1) * t70 - g(2) * t4, -g(1) * t71 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, -t39, g(1) * t19 + g(3) * t55 + t64, 0, 0, 0, 0, 0, -t39 * t36, t7, 0, 0, 0, 0, 0, -g(1) * (-t18 * t50 + t19 * t31) - g(2) * (-t16 * t50 + t17 * t31) - (t31 * t33 + t35 * t49) * t63, -g(1) * (t18 * t51 + t19 * t35) - g(2) * (t16 * t51 + t17 * t35) - (-t31 * t49 + t33 * t35) * t63, -t7, -g(1) * (-t42 * t18 + t47 * t19) - t47 * t64 + t42 * t65 - (t47 * t33 + t42 * t37) * t63, 0, 0, 0, 0, 0, -g(1) * (-t18 * t56 + t19 * t25) - g(2) * (-t16 * t56 + t17 * t25) - (t25 * t33 + t26 * t49) * t63, -g(1) * (t18 * t57 + t19 * t26) - g(2) * (t16 * t57 + t17 * t26) - (-t25 * t49 + t26 * t33) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, 0, 0, 0, 0, 0, t41 * t35, -t41 * t31, -t40, -g(1) * (-t12 * t27 - t13 * t30) - g(2) * (-t8 * t27 - t9 * t30) - g(3) * (-t14 * t27 - t15 * t30) 0, 0, 0, 0, 0, t41 * t26, -t41 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, g(1) * t6 + g(2) * t68 - g(3) * (-t15 * t35 + t31 * t52) 0, t67 * pkin(4), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t10;
