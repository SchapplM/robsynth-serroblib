% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t44 = cos(qJ(1));
t65 = sin(pkin(13));
t69 = cos(pkin(6));
t58 = t69 * t65;
t67 = cos(pkin(13));
t75 = sin(qJ(1));
t28 = t44 * t58 + t75 * t67;
t41 = sin(qJ(3));
t76 = cos(qJ(3));
t59 = t69 * t67;
t50 = -t44 * t59 + t75 * t65;
t38 = sin(pkin(6));
t66 = sin(pkin(7));
t63 = t38 * t66;
t68 = cos(pkin(7));
t82 = t44 * t63 + t50 * t68;
t14 = t28 * t41 + t82 * t76;
t37 = qJ(5) + qJ(6);
t35 = sin(t37);
t36 = cos(t37);
t17 = -t28 * t76 + t82 * t41;
t64 = t38 * t68;
t23 = -t44 * t64 + t50 * t66;
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t9 = t17 * t43 - t23 * t40;
t90 = t14 * t36 + t9 * t35;
t89 = -t14 * t35 + t9 * t36;
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t88 = t14 * t42 + t9 * t39;
t87 = -t14 * t39 + t9 * t42;
t81 = t17 * t40 + t23 * t43;
t46 = t44 * t65 + t75 * t59;
t77 = t46 * t68 - t75 * t63;
t74 = t35 * t43;
t73 = t36 * t43;
t72 = t39 * t43;
t71 = t42 * t43;
t70 = qJ(2) * t38;
t57 = t68 * t67;
t56 = t66 * t69;
t54 = g(1) * t75 - g(2) * t44;
t53 = -g(1) * t44 - g(2) * t75;
t29 = t44 * t67 - t75 * t58;
t19 = t29 * t76 - t77 * t41;
t24 = t46 * t66 + t75 * t64;
t10 = -t19 * t40 + t24 * t43;
t21 = t41 * t56 + (t41 * t57 + t76 * t65) * t38;
t27 = -t67 * t63 + t69 * t68;
t52 = g(1) * t10 + g(2) * t81 + g(3) * (-t21 * t40 + t27 * t43);
t18 = t29 * t41 + t77 * t76;
t20 = -t76 * t56 + (t41 * t65 - t57 * t76) * t38;
t51 = g(1) * t18 + g(2) * t14 + g(3) * t20;
t13 = t21 * t43 + t27 * t40;
t11 = t19 * t43 + t24 * t40;
t6 = t11 * t42 + t18 * t39;
t5 = -t11 * t39 + t18 * t42;
t4 = t11 * t36 + t18 * t35;
t3 = -t11 * t35 + t18 * t36;
t2 = g(1) * t4 - g(2) * t89 - g(3) * (-t13 * t36 - t20 * t35);
t1 = -g(1) * t3 - g(2) * t90 - g(3) * (-t13 * t35 + t20 * t36);
t7 = [0, t54, -t53, g(1) * t28 - g(2) * t29, -g(1) * t50 + g(2) * t46, t53 * t38, -g(1) * (-t75 * pkin(1) + t44 * t70) - g(2) * (t44 * pkin(1) + t75 * t70) 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t19, -g(1) * t14 + g(2) * t18, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, g(1) * t81 - g(2) * t10, 0, 0, 0, 0, 0, -g(1) * t87 - g(2) * t6, g(1) * t88 - g(2) * t5, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t4, g(1) * t90 - g(2) * t3; 0, 0, 0, 0, 0, 0, -g(3) * t69 - t54 * t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, g(1) * t19 - g(2) * t17 + g(3) * t21, 0, 0, 0, 0, 0, t51 * t43, -t51 * t40, 0, 0, 0, 0, 0, -g(1) * (-t18 * t71 + t19 * t39) - g(2) * (-t14 * t71 - t17 * t39) - g(3) * (-t20 * t71 + t21 * t39) -g(1) * (t18 * t72 + t19 * t42) - g(2) * (t14 * t72 - t17 * t42) - g(3) * (t20 * t72 + t21 * t42) 0, 0, 0, 0, 0, -g(1) * (-t18 * t73 + t19 * t35) - g(2) * (-t14 * t73 - t17 * t35) - g(3) * (-t20 * t73 + t21 * t35) -g(1) * (t18 * t74 + t19 * t36) - g(2) * (t14 * t74 - t17 * t36) - g(3) * (t20 * t74 + t21 * t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, g(1) * t11 - g(2) * t9 + g(3) * t13, 0, 0, 0, 0, 0, -t52 * t42, t52 * t39, 0, 0, 0, 0, 0, -t52 * t36, t52 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t88 - g(3) * (-t13 * t39 + t20 * t42) g(1) * t6 - g(2) * t87 - g(3) * (-t13 * t42 - t20 * t39) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t7;
