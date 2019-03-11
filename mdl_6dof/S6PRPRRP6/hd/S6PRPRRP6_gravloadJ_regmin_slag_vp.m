% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t43 = sin(pkin(10));
t45 = cos(pkin(10));
t52 = cos(qJ(2));
t46 = cos(pkin(6));
t49 = sin(qJ(2));
t67 = t46 * t49;
t31 = t43 * t52 + t45 * t67;
t33 = -t43 * t67 + t45 * t52;
t77 = -g(1) * t33 - g(2) * t31;
t66 = t46 * t52;
t30 = t43 * t49 - t45 * t66;
t32 = t43 * t66 + t45 * t49;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t44 = sin(pkin(6));
t68 = t44 * t52;
t71 = t44 * t48;
t55 = g(3) * (-t46 * t48 - t51 * t68) + g(2) * (t30 * t51 + t45 * t71) + g(1) * (t32 * t51 - t43 * t71);
t70 = t44 * t49;
t69 = t44 * t51;
t47 = sin(qJ(5));
t65 = t47 * t48;
t50 = cos(qJ(5));
t64 = t48 * t50;
t63 = t49 * t50;
t62 = pkin(2) * t68 + qJ(3) * t70;
t61 = t47 * t70;
t60 = pkin(4) * t48 - pkin(9) * t51;
t35 = t46 * t51 - t48 * t68;
t19 = t35 * t47 - t44 * t63;
t16 = t32 * t48 + t43 * t69;
t5 = t16 * t47 - t33 * t50;
t18 = -t30 * t48 + t45 * t69;
t7 = -t18 * t47 - t31 * t50;
t1 = g(1) * t5 + g(2) * t7 + g(3) * t19;
t20 = t35 * t50 + t61;
t6 = t16 * t50 + t33 * t47;
t8 = -t18 * t50 + t31 * t47;
t57 = g(1) * t6 + g(2) * t8 + g(3) * t20;
t11 = t30 * t50 + t31 * t65;
t13 = t32 * t50 + t33 * t65;
t21 = t48 * t61 - t50 * t68;
t56 = g(1) * t13 + g(2) * t11 + g(3) * t21;
t54 = g(1) * t16 - g(2) * t18 + g(3) * t35;
t10 = -g(1) * t32 - g(2) * t30 + g(3) * t68;
t53 = g(3) * t70 - t77;
t29 = t32 * pkin(2);
t28 = t30 * pkin(2);
t22 = (t47 * t52 + t48 * t63) * t44;
t14 = -t32 * t47 + t33 * t64;
t12 = -t30 * t47 + t31 * t64;
t9 = t53 * t51;
t4 = t55 * t50;
t3 = t55 * t47;
t2 = -g(1) * t14 - g(2) * t12 - g(3) * t22;
t15 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t10, t53, t10, -t53, -g(1) * (t33 * qJ(3) - t29) - g(2) * (t31 * qJ(3) - t28) - g(3) * t62, 0, 0, 0, 0, 0, -t53 * t48, -t9, 0, 0, 0, 0, 0, t2, t56, t2, t9, -t56, -g(1) * (t14 * pkin(5) - t32 * pkin(8) + t13 * qJ(6) - t29) - g(2) * (t12 * pkin(5) - t30 * pkin(8) + t11 * qJ(6) - t28) + t77 * (qJ(3) + t60) + (-t22 * pkin(5) - t21 * qJ(6) - t62 - (pkin(8) * t52 + t60 * t49) * t44) * g(3); 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t54, 0, 0, 0, 0, 0, -t4, t3, -t4, -t54, -t3, -t54 * pkin(9) - t55 * (pkin(5) * t50 + qJ(6) * t47 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t57, t1, 0, -t57, -g(1) * (-t5 * pkin(5) + t6 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - g(3) * (-t19 * pkin(5) + t20 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t15;
