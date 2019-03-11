% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t38 = sin(pkin(11));
t43 = sin(qJ(2));
t47 = cos(qJ(2));
t62 = cos(pkin(11));
t30 = -t47 * t38 - t43 * t62;
t52 = -t43 * t38 + t47 * t62;
t40 = cos(pkin(6));
t48 = cos(qJ(1));
t65 = t48 * t43;
t44 = sin(qJ(1));
t68 = t44 * t47;
t26 = -t40 * t68 - t65;
t39 = sin(pkin(6));
t75 = t39 * t47;
t82 = -g(1) * t26 - g(3) * t75;
t63 = t30 * t40;
t10 = -t44 * t52 + t48 * t63;
t41 = sin(qJ(6));
t45 = cos(qJ(6));
t49 = t52 * t40;
t11 = t44 * t30 + t48 * t49;
t42 = sin(qJ(5));
t46 = cos(qJ(5));
t74 = t39 * t48;
t53 = t11 * t42 + t46 * t74;
t81 = -t10 * t45 + t41 * t53;
t15 = t44 * t63 + t48 * t52;
t80 = t10 * t41 + t45 * t53;
t76 = t39 * t44;
t73 = t41 * t42;
t72 = t42 * t45;
t69 = t44 * t43;
t64 = t48 * t47;
t60 = t40 * t64;
t23 = t40 * t43 * pkin(2) + (-pkin(8) - qJ(3)) * t39;
t37 = t47 * pkin(2) + pkin(1);
t59 = -t44 * t23 + t48 * t37;
t56 = g(1) * t48 + g(2) * t44;
t55 = g(1) * t44 - g(2) * t48;
t54 = -t48 * t23 - t44 * t37;
t6 = -t11 * t46 + t42 * t74;
t21 = t52 * t39;
t14 = t48 * t30 - t44 * t49;
t4 = -t14 * t46 - t42 * t76;
t51 = g(1) * t4 + g(2) * t6 + g(3) * (-t21 * t46 - t40 * t42);
t22 = t30 * t39;
t50 = -g(1) * t15 + g(2) * t10 + g(3) * t22;
t31 = pkin(2) * t60;
t28 = t56 * t39;
t27 = -t40 * t69 + t64;
t25 = -t40 * t65 - t68;
t24 = -t60 + t69;
t20 = -g(3) * t40 - t55 * t39;
t17 = -t21 * t42 + t40 * t46;
t5 = -t14 * t42 + t46 * t76;
t3 = t15 * t41 + t5 * t45;
t2 = t15 * t45 - t5 * t41;
t1 = g(1) * t14 + g(2) * t11 + g(3) * t21;
t7 = [0, t55, t56, 0, 0, 0, 0, 0, -g(1) * t25 - g(2) * t27, -g(1) * t24 - g(2) * t26, -t28, -g(1) * t54 - g(2) * t59, -t28, g(1) * t10 + g(2) * t15, -g(1) * t11 + g(2) * t14, -g(1) * (pkin(3) * t10 + t11 * qJ(4) + t54) - g(2) * (t15 * pkin(3) - t14 * qJ(4) + t59) 0, 0, 0, 0, 0, -g(1) * t53 - g(2) * t5, g(1) * t6 - g(2) * t4, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t3, g(1) * t81 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t24 + t82, g(3) * t39 * t43 + g(1) * t27 - g(2) * t25, 0, -g(2) * t31 + (g(2) * t69 + t82) * pkin(2), 0, t1, t50, -g(1) * (t26 * pkin(2) + t14 * pkin(3) + qJ(4) * t15) - g(2) * (-pkin(2) * t69 + t11 * pkin(3) - t10 * qJ(4) + t31) - g(3) * (pkin(2) * t75 + t21 * pkin(3) - t22 * qJ(4)) 0, 0, 0, 0, 0, t50 * t42, t50 * t46, 0, 0, 0, 0, 0, -g(1) * (t14 * t41 + t15 * t72) - g(2) * (-t10 * t72 + t11 * t41) - g(3) * (t21 * t41 - t22 * t72) -g(1) * (t14 * t45 - t15 * t73) - g(2) * (t10 * t73 + t11 * t45) - g(3) * (t21 * t45 + t22 * t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, g(1) * t5 - g(2) * t53 + g(3) * t17, 0, 0, 0, 0, 0, -t51 * t45, t51 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t81 - g(3) * (-t17 * t41 - t22 * t45) g(1) * t3 - g(2) * t80 - g(3) * (-t17 * t45 + t22 * t41);];
taug_reg  = t7;
