% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:20:06
% EndTime: 2019-05-05 06:20:07
% DurationCPUTime: 0.43s
% Computational Cost: add. (460->115), mult. (1314->203), div. (0->0), fcn. (1698->14), ass. (0->67)
t39 = sin(pkin(7));
t79 = pkin(9) * t39;
t78 = cos(qJ(3));
t44 = sin(qJ(5));
t77 = t39 * t44;
t48 = cos(qJ(5));
t76 = t39 * t48;
t40 = sin(pkin(6));
t42 = cos(pkin(7));
t75 = t40 * t42;
t46 = sin(qJ(2));
t74 = t40 * t46;
t49 = cos(qJ(2));
t73 = t40 * t49;
t45 = sin(qJ(3));
t72 = t42 * t45;
t43 = sin(qJ(6));
t71 = t43 * t44;
t47 = cos(qJ(6));
t70 = t44 * t47;
t69 = t45 * t46;
t68 = t45 * t49;
t67 = cos(pkin(6));
t66 = sin(pkin(12));
t41 = cos(pkin(12));
t65 = t39 * t40 * t41;
t64 = t39 * t74;
t63 = t42 * t78;
t62 = t78 * t46;
t61 = t78 * t49;
t60 = t40 * t66;
t59 = t41 * t67;
t58 = t67 * t39;
t57 = t39 * t60;
t56 = t67 * t66;
t32 = -t41 * t46 - t49 * t56;
t33 = t41 * t49 - t46 * t56;
t10 = -t32 * t63 + t33 * t45 - t78 * t57;
t19 = t40 * t69 - t78 * t58 - t61 * t75;
t30 = -t66 * t46 + t49 * t59;
t21 = -t30 * t39 - t41 * t75;
t22 = -t32 * t39 + t42 * t60;
t29 = -t39 * t73 + t67 * t42;
t31 = t46 * t59 + t66 * t49;
t8 = -t30 * t63 + t31 * t45 + t78 * t65;
t55 = g(1) * (t10 * t48 - t22 * t44) + g(2) * (-t21 * t44 + t8 * t48) + g(3) * (t19 * t48 - t29 * t44);
t54 = g(1) * t10 + g(2) * t8 + g(3) * t19;
t11 = t33 * t78 + (t32 * t42 + t57) * t45;
t20 = t45 * t58 + (t42 * t68 + t62) * t40;
t9 = t30 * t72 + t31 * t78 - t45 * t65;
t53 = g(1) * t11 + g(2) * t9 + g(3) * t20;
t14 = t30 * t45 + t31 * t63;
t16 = t32 * t45 + t33 * t63;
t27 = (t42 * t62 + t68) * t40;
t52 = g(1) * t16 + g(2) * t14 + g(3) * t27;
t15 = t30 * t78 - t31 * t72;
t17 = t32 * t78 - t33 * t72;
t28 = (-t42 * t69 + t61) * t40;
t51 = g(1) * t17 + g(2) * t15 + g(3) * t28;
t50 = g(1) * t33 + g(2) * t31 + g(3) * t74;
t18 = t27 * t44 + t48 * t64;
t13 = t19 * t44 + t29 * t48;
t7 = t16 * t44 + t33 * t76;
t6 = t14 * t44 + t31 * t76;
t5 = t10 * t44 + t22 * t48;
t3 = t21 * t48 + t8 * t44;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(1) * t32 - g(2) * t30 - g(3) * t73, t50, 0, 0, 0, 0, 0, -t51, t52, -t50 * t39, t51, -t52, -g(1) * (t32 * pkin(2) + t17 * pkin(3) + t16 * qJ(4) + t33 * t79) - g(2) * (t30 * pkin(2) + t15 * pkin(3) + t14 * qJ(4) + t31 * t79) - g(3) * (t28 * pkin(3) + t27 * qJ(4) + (pkin(2) * t49 + t46 * t79) * t40) 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t6 - g(3) * t18, -g(1) * (t16 * t48 - t33 * t77) - g(2) * (t14 * t48 - t31 * t77) - g(3) * (t27 * t48 - t44 * t64) 0, 0, 0, 0, 0, -g(1) * (t17 * t43 + t7 * t47) - g(2) * (t15 * t43 + t6 * t47) - g(3) * (t18 * t47 + t28 * t43) -g(1) * (t17 * t47 - t7 * t43) - g(2) * (t15 * t47 - t6 * t43) - g(3) * (-t18 * t43 + t28 * t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t53, 0, -t54, -t53, -g(1) * (-t10 * pkin(3) + t11 * qJ(4)) - g(2) * (-t8 * pkin(3) + t9 * qJ(4)) - g(3) * (-t19 * pkin(3) + t20 * qJ(4)) 0, 0, 0, 0, 0, -t53 * t44, -t53 * t48, 0, 0, 0, 0, 0, -g(1) * (-t10 * t43 + t11 * t70) - g(2) * (-t8 * t43 + t9 * t70) - g(3) * (-t19 * t43 + t20 * t70) -g(1) * (-t10 * t47 - t11 * t71) - g(2) * (-t8 * t47 - t9 * t71) - g(3) * (-t19 * t47 - t20 * t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, g(1) * t5 + g(2) * t3 + g(3) * t13, 0, 0, 0, 0, 0, -t55 * t47, t55 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t47 - t5 * t43) - g(2) * (-t3 * t43 + t9 * t47) - g(3) * (-t13 * t43 + t20 * t47) -g(1) * (-t11 * t43 - t5 * t47) - g(2) * (-t3 * t47 - t9 * t43) - g(3) * (-t13 * t47 - t20 * t43);];
taug_reg  = t1;
