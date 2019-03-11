% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR10
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t46 = cos(qJ(1));
t67 = sin(pkin(13));
t71 = cos(pkin(6));
t60 = t71 * t67;
t69 = cos(pkin(13));
t75 = sin(qJ(1));
t30 = t46 * t60 + t75 * t69;
t43 = sin(qJ(3));
t76 = cos(qJ(3));
t61 = t71 * t69;
t52 = -t46 * t61 + t75 * t67;
t40 = sin(pkin(6));
t68 = sin(pkin(7));
t65 = t40 * t68;
t70 = cos(pkin(7));
t86 = t46 * t65 + t52 * t70;
t16 = t30 * t43 + t76 * t86;
t41 = sin(qJ(6));
t44 = cos(qJ(6));
t19 = -t30 * t76 + t43 * t86;
t66 = t40 * t70;
t25 = -t46 * t66 + t52 * t68;
t39 = qJ(4) + qJ(5);
t37 = sin(t39);
t38 = cos(t39);
t9 = t19 * t38 - t25 * t37;
t90 = t16 * t44 + t9 * t41;
t89 = -t16 * t41 + t9 * t44;
t85 = t19 * t37 + t25 * t38;
t42 = sin(qJ(4));
t45 = cos(qJ(4));
t84 = t19 * t42 + t25 * t45;
t83 = t19 * t45 - t25 * t42;
t48 = t46 * t67 + t75 * t61;
t77 = t48 * t70 - t75 * t65;
t74 = t38 * t41;
t73 = t38 * t44;
t72 = qJ(2) * t40;
t59 = t70 * t69;
t58 = t68 * t71;
t56 = g(1) * t75 - g(2) * t46;
t55 = -g(1) * t46 - g(2) * t75;
t31 = t46 * t69 - t75 * t60;
t21 = t31 * t76 - t77 * t43;
t26 = t48 * t68 + t75 * t66;
t10 = -t21 * t37 + t26 * t38;
t23 = t43 * t58 + (t43 * t59 + t76 * t67) * t40;
t29 = -t69 * t65 + t71 * t70;
t54 = g(1) * t10 + g(2) * t85 + g(3) * (-t23 * t37 + t29 * t38);
t20 = t31 * t43 + t77 * t76;
t22 = -t76 * t58 + (t43 * t67 - t59 * t76) * t40;
t53 = g(1) * t20 + g(2) * t16 + g(3) * t22;
t15 = t23 * t38 + t29 * t37;
t13 = t21 * t45 + t26 * t42;
t12 = -t21 * t42 + t26 * t45;
t11 = t21 * t38 + t26 * t37;
t6 = t11 * t44 + t20 * t41;
t5 = -t11 * t41 + t20 * t44;
t4 = g(1) * t11 - g(2) * t9 + g(3) * t15;
t2 = t54 * t44;
t1 = t54 * t41;
t3 = [0, t56, -t55, g(1) * t30 - g(2) * t31, -g(1) * t52 + g(2) * t48, t55 * t40, -g(1) * (-t75 * pkin(1) + t46 * t72) - g(2) * (t46 * pkin(1) + t75 * t72) 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t21, -g(1) * t16 + g(2) * t20, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t13, g(1) * t84 - g(2) * t12, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, g(1) * t85 - g(2) * t10, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t6, g(1) * t90 - g(2) * t5; 0, 0, 0, 0, 0, 0, -g(3) * t71 - t56 * t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, g(1) * t21 - g(2) * t19 + g(3) * t23, 0, 0, 0, 0, 0, t53 * t45, -t53 * t42, 0, 0, 0, 0, 0, t53 * t38, -t53 * t37, 0, 0, 0, 0, 0, -g(1) * (-t20 * t73 + t21 * t41) - g(2) * (-t16 * t73 - t19 * t41) - g(3) * (-t22 * t73 + t23 * t41) -g(1) * (t20 * t74 + t21 * t44) - g(2) * (t16 * t74 - t19 * t44) - g(3) * (t22 * t74 + t23 * t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t84 - g(3) * (-t23 * t42 + t29 * t45) g(1) * t13 - g(2) * t83 - g(3) * (-t23 * t45 - t29 * t42) 0, 0, 0, 0, 0, -t54, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t90 - g(3) * (-t15 * t41 + t22 * t44) g(1) * t6 - g(2) * t89 - g(3) * (-t15 * t44 - t22 * t41);];
taug_reg  = t3;
