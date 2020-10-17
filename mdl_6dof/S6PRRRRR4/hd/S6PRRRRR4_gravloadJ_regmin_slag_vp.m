% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:28:45
% EndTime: 2019-05-05 11:28:47
% DurationCPUTime: 0.49s
% Computational Cost: add. (615->113), mult. (1492->213), div. (0->0), fcn. (1950->16), ass. (0->66)
t40 = sin(pkin(7));
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t42 = cos(pkin(13));
t67 = cos(pkin(6));
t61 = t42 * t67;
t65 = sin(pkin(13));
t52 = t46 * t65 - t49 * t61;
t66 = cos(pkin(7));
t41 = sin(pkin(6));
t70 = t41 * t42;
t79 = t40 * t70 + t52 * t66;
t56 = t67 * t65;
t53 = t42 * t46 + t49 * t56;
t62 = t41 * t65;
t78 = -t40 * t62 + t53 * t66;
t77 = cos(qJ(3));
t39 = qJ(4) + qJ(5);
t37 = sin(t39);
t76 = t37 * t40;
t38 = cos(t39);
t75 = t38 * t40;
t43 = sin(qJ(6));
t74 = t38 * t43;
t47 = cos(qJ(6));
t73 = t38 * t47;
t44 = sin(qJ(4));
t72 = t40 * t44;
t48 = cos(qJ(4));
t71 = t40 * t48;
t69 = t41 * t46;
t68 = t41 * t49;
t64 = t40 * t69;
t45 = sin(qJ(3));
t60 = t45 * t66;
t59 = t67 * t40;
t57 = t66 * t77;
t31 = t46 * t61 + t49 * t65;
t14 = t31 * t77 - t79 * t45;
t32 = t42 * t49 - t46 * t56;
t16 = t32 * t77 - t78 * t45;
t23 = t45 * t59 + (t46 * t77 + t49 * t60) * t41;
t24 = t40 * t52 - t66 * t70;
t25 = t40 * t53 + t62 * t66;
t30 = -t40 * t68 + t67 * t66;
t55 = g(1) * (-t16 * t37 + t25 * t38) + g(2) * (-t14 * t37 + t24 * t38) + g(3) * (-t23 * t37 + t30 * t38);
t13 = t31 * t45 + t79 * t77;
t15 = t32 * t45 + t78 * t77;
t22 = t45 * t69 - t57 * t68 - t59 * t77;
t54 = g(1) * t15 + g(2) * t13 + g(3) * t22;
t29 = (-t46 * t60 + t49 * t77) * t41;
t28 = (t45 * t49 + t46 * t57) * t41;
t21 = t29 * t38 + t37 * t64;
t20 = -t32 * t60 - t53 * t77;
t19 = t32 * t57 - t45 * t53;
t18 = -t31 * t60 - t52 * t77;
t17 = t31 * t57 - t45 * t52;
t12 = t23 * t38 + t30 * t37;
t10 = t20 * t38 + t32 * t76;
t9 = t18 * t38 + t31 * t76;
t8 = t16 * t38 + t25 * t37;
t6 = t14 * t38 + t24 * t37;
t4 = g(1) * t8 + g(2) * t6 + g(3) * t12;
t2 = t55 * t47;
t1 = t55 * t43;
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t53 + g(2) * t52 - g(3) * t68, g(1) * t32 + g(2) * t31 + g(3) * t69, 0, 0, 0, 0, 0, -g(1) * t20 - g(2) * t18 - g(3) * t29, g(1) * t19 + g(2) * t17 + g(3) * t28, 0, 0, 0, 0, 0, -g(1) * (t20 * t48 + t32 * t72) - g(2) * (t18 * t48 + t31 * t72) - g(3) * (t29 * t48 + t44 * t64) -g(1) * (-t20 * t44 + t32 * t71) - g(2) * (-t18 * t44 + t31 * t71) - g(3) * (-t29 * t44 + t48 * t64) 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t9 - g(3) * t21, -g(1) * (-t20 * t37 + t32 * t75) - g(2) * (-t18 * t37 + t31 * t75) - g(3) * (-t29 * t37 + t38 * t64) 0, 0, 0, 0, 0, -g(1) * (t10 * t47 + t19 * t43) - g(2) * (t17 * t43 + t9 * t47) - g(3) * (t21 * t47 + t28 * t43) -g(1) * (-t10 * t43 + t19 * t47) - g(2) * (t17 * t47 - t9 * t43) - g(3) * (-t21 * t43 + t28 * t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, g(1) * t16 + g(2) * t14 + g(3) * t23, 0, 0, 0, 0, 0, t54 * t48, -t54 * t44, 0, 0, 0, 0, 0, t54 * t38, -t54 * t37, 0, 0, 0, 0, 0, -g(1) * (-t15 * t73 + t16 * t43) - g(2) * (-t13 * t73 + t14 * t43) - g(3) * (-t22 * t73 + t23 * t43) -g(1) * (t15 * t74 + t16 * t47) - g(2) * (t13 * t74 + t14 * t47) - g(3) * (t22 * t74 + t23 * t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t44 + t25 * t48) - g(2) * (-t14 * t44 + t24 * t48) - g(3) * (-t23 * t44 + t30 * t48) -g(1) * (-t16 * t48 - t25 * t44) - g(2) * (-t14 * t48 - t24 * t44) - g(3) * (-t23 * t48 - t30 * t44) 0, 0, 0, 0, 0, -t55, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t15 * t47 - t8 * t43) - g(2) * (t13 * t47 - t6 * t43) - g(3) * (-t12 * t43 + t22 * t47) -g(1) * (-t15 * t43 - t8 * t47) - g(2) * (-t13 * t43 - t6 * t47) - g(3) * (-t12 * t47 - t22 * t43);];
taug_reg  = t3;
