% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t46 = cos(qJ(2));
t62 = cos(pkin(6));
t69 = cos(qJ(1));
t55 = t62 * t69;
t28 = t42 * t41 - t46 * t55;
t38 = sin(qJ(6));
t43 = cos(qJ(6));
t29 = t41 * t55 + t42 * t46;
t40 = sin(qJ(3));
t45 = cos(qJ(3));
t37 = sin(pkin(6));
t60 = t37 * t69;
t19 = t29 * t40 + t45 * t60;
t20 = t29 * t45 - t40 * t60;
t39 = sin(qJ(5));
t44 = cos(qJ(5));
t6 = t19 * t39 + t20 * t44;
t80 = t28 * t43 + t6 * t38;
t79 = -t28 * t38 + t6 * t43;
t66 = t37 * t41;
t26 = t40 * t66 - t62 * t45;
t64 = t37 * t45;
t27 = t62 * t40 + t41 * t64;
t76 = t19 * t44 - t20 * t39;
t59 = t42 * t62;
t31 = -t41 * t59 + t69 * t46;
t23 = t31 * t40 - t42 * t64;
t65 = t37 * t42;
t24 = t31 * t45 + t40 * t65;
t9 = t23 * t44 - t24 * t39;
t50 = g(1) * t9 + g(2) * t76 + g(3) * (t26 * t44 - t27 * t39);
t78 = t50 * t38;
t77 = t50 * t43;
t63 = t37 * t46;
t71 = g(2) * t28;
t30 = t69 * t41 + t46 * t59;
t72 = g(1) * t30;
t48 = g(3) * t63 - t71 - t72;
t10 = t23 * t39 + t24 * t44;
t18 = t26 * t39 + t27 * t44;
t75 = g(1) * t10 + g(2) * t6 + g(3) * t18;
t58 = -g(1) * t19 + g(2) * t23;
t57 = g(1) * t28 - g(2) * t30;
t56 = -g(1) * t31 - g(2) * t29;
t53 = t39 * t40 + t44 * t45;
t52 = pkin(3) * t45 + qJ(4) * t40 + pkin(2);
t3 = g(1) * t23 + g(2) * t19 + g(3) * t26;
t49 = g(1) * t24 + g(2) * t20 + g(3) * t27;
t47 = g(3) * t66 - t56;
t25 = t53 * t63;
t15 = t53 * t30;
t14 = t53 * t28;
t13 = t48 * t45;
t12 = t48 * t40;
t11 = g(1) * t20 - g(2) * t24;
t2 = t10 * t43 - t30 * t38;
t1 = -t10 * t38 - t30 * t43;
t4 = [0, g(1) * t42 - g(2) * t69, g(1) * t69 + g(2) * t42, 0, 0, 0, 0, 0, g(1) * t29 - g(2) * t31, -t57, 0, 0, 0, 0, 0, t11, t58, t11, t57, -t58, -g(1) * (-t42 * pkin(1) - t29 * pkin(2) - pkin(3) * t20 + pkin(8) * t60 - t28 * pkin(9) - qJ(4) * t19) - g(2) * (t69 * pkin(1) + t31 * pkin(2) + t24 * pkin(3) + pkin(8) * t65 + t30 * pkin(9) + t23 * qJ(4)) 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t10, g(1) * t76 - g(2) * t9, 0, 0, 0, 0, 0, g(1) * t79 - g(2) * t2, -g(1) * t80 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, 0, 0, 0, 0, -t13, t12, -t13, -t47, -t12, t56 * pkin(9) + t52 * t72 + t52 * t71 - g(3) * (pkin(9) * t41 + t52 * t46) * t37, 0, 0, 0, 0, 0, g(1) * t15 + g(2) * t14 - g(3) * t25, t48 * (t39 * t45 - t40 * t44) 0, 0, 0, 0, 0, -g(1) * (-t15 * t43 - t31 * t38) - g(2) * (-t14 * t43 - t29 * t38) - g(3) * (t25 * t43 - t38 * t66) -g(1) * (t15 * t38 - t31 * t43) - g(2) * (t14 * t38 - t29 * t43) - g(3) * (-t25 * t38 - t43 * t66); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t49, t3, 0, -t49, -g(1) * (-t23 * pkin(3) + t24 * qJ(4)) - g(2) * (-t19 * pkin(3) + t20 * qJ(4)) - g(3) * (-t26 * pkin(3) + t27 * qJ(4)) 0, 0, 0, 0, 0, t50, -t75, 0, 0, 0, 0, 0, t77, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t75, 0, 0, 0, 0, 0, -t77, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t80 - g(3) * (-t18 * t38 + t43 * t63) g(1) * t2 + g(2) * t79 - g(3) * (-t18 * t43 - t38 * t63);];
taug_reg  = t4;
