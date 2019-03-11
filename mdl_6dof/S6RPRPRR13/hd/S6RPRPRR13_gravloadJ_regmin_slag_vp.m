% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t42 = cos(pkin(6));
t40 = cos(pkin(12));
t49 = cos(qJ(1));
t66 = t49 * t40;
t37 = sin(pkin(12));
t46 = sin(qJ(1));
t70 = t46 * t37;
t27 = -t42 * t66 + t70;
t67 = t49 * t37;
t68 = t46 * t40;
t28 = t42 * t67 + t68;
t41 = cos(pkin(7));
t45 = sin(qJ(3));
t73 = t41 * t45;
t39 = sin(pkin(6));
t74 = t39 * t49;
t38 = sin(pkin(7));
t76 = t38 * t45;
t77 = cos(qJ(3));
t12 = -t27 * t73 + t28 * t77 - t74 * t76;
t43 = sin(qJ(6));
t47 = cos(qJ(6));
t63 = t38 * t77;
t60 = t39 * t63;
t62 = t41 * t77;
t11 = t27 * t62 + t28 * t45 + t49 * t60;
t19 = -t27 * t38 + t41 * t74;
t44 = sin(qJ(5));
t48 = cos(qJ(5));
t6 = t11 * t44 - t19 * t48;
t82 = -t12 * t47 + t6 * t43;
t81 = t12 * t43 + t6 * t47;
t78 = t11 * t48 + t19 * t44;
t75 = t39 * t40;
t72 = t43 * t44;
t71 = t44 * t47;
t69 = t46 * t39;
t64 = qJ(2) * t39;
t65 = t49 * pkin(1) + t46 * t64;
t61 = -t46 * pkin(1) + t49 * t64;
t29 = -t42 * t70 + t66;
t54 = t42 * t68 + t67;
t50 = t54 * t41;
t15 = t29 * t45 - t46 * t60 + t77 * t50;
t59 = -g(1) * t11 + g(2) * t15;
t16 = t29 * t77 + (t38 * t69 - t50) * t45;
t58 = -g(1) * t12 + g(2) * t16;
t21 = t54 * t38 + t41 * t69;
t57 = -g(1) * t19 - g(2) * t21;
t56 = g(1) * t49 + g(2) * t46;
t55 = g(1) * t46 - g(2) * t49;
t17 = t39 * t37 * t45 - t42 * t63 - t62 * t75;
t26 = -t38 * t75 + t42 * t41;
t7 = t15 * t48 - t21 * t44;
t53 = g(1) * t7 + g(2) * t78 + g(3) * (t17 * t48 - t26 * t44);
t52 = g(1) * t15 + g(2) * t11 + g(3) * t17;
t18 = t42 * t76 + (t77 * t37 + t40 * t73) * t39;
t51 = g(1) * t16 + g(2) * t12 + g(3) * t18;
t25 = -g(3) * t42 - t55 * t39;
t10 = t17 * t44 + t26 * t48;
t8 = t15 * t44 + t21 * t48;
t2 = t16 * t43 + t8 * t47;
t1 = t16 * t47 - t8 * t43;
t3 = [0, t55, t56, g(1) * t28 - g(2) * t29, -g(1) * t27 + g(2) * t54, -t56 * t39, -g(1) * t61 - g(2) * t65, 0, 0, 0, 0, 0, -t58, t59, t57, t58, -t59, -g(1) * (-t28 * pkin(2) - pkin(3) * t12 - qJ(4) * t11 + t61) - g(2) * (t29 * pkin(2) + t16 * pkin(3) + t15 * qJ(4) + t65) + t57 * pkin(9), 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t8, g(1) * t78 - g(2) * t7, 0, 0, 0, 0, 0, g(1) * t81 - g(2) * t2, -g(1) * t82 - g(2) * t1; 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t51, 0, -t52, -t51, -g(1) * (-t15 * pkin(3) + t16 * qJ(4)) - g(2) * (-t11 * pkin(3) + t12 * qJ(4)) - g(3) * (-t17 * pkin(3) + t18 * qJ(4)) 0, 0, 0, 0, 0, -t51 * t44, -t51 * t48, 0, 0, 0, 0, 0, -g(1) * (-t15 * t43 + t16 * t71) - g(2) * (-t11 * t43 + t12 * t71) - g(3) * (-t17 * t43 + t18 * t71) -g(1) * (-t15 * t47 - t16 * t72) - g(2) * (-t11 * t47 - t12 * t72) - g(3) * (-t17 * t47 - t18 * t72); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, g(1) * t8 + g(2) * t6 + g(3) * t10, 0, 0, 0, 0, 0, -t53 * t47, t53 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t82 - g(3) * (-t10 * t43 + t18 * t47) g(1) * t2 + g(2) * t81 - g(3) * (-t10 * t47 - t18 * t43);];
taug_reg  = t3;
