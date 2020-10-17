% Calculate minimal parameter regressor of gravitation load for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:17:33
% EndTime: 2019-05-04 20:17:34
% DurationCPUTime: 0.28s
% Computational Cost: add. (439->71), mult. (1241->119), div. (0->0), fcn. (1615->14), ass. (0->54)
t57 = sin(pkin(12));
t58 = sin(pkin(11));
t48 = t58 * t57;
t61 = cos(pkin(12));
t62 = cos(pkin(11));
t55 = t62 * t61;
t64 = cos(pkin(6));
t38 = -t64 * t55 + t48;
t59 = sin(pkin(7));
t60 = sin(pkin(6));
t52 = t60 * t59;
t63 = cos(pkin(7));
t73 = t38 * t63 + t62 * t52;
t50 = t58 * t61;
t53 = t62 * t57;
t39 = t64 * t50 + t53;
t49 = t58 * t60;
t72 = t39 * t63 - t59 * t49;
t71 = t61 * t63 * t60 + t64 * t59;
t23 = -t64 * t48 + t55;
t30 = sin(qJ(3));
t67 = cos(qJ(3));
t10 = t23 * t30 + t72 * t67;
t51 = t60 * t57;
t16 = t30 * t51 - t71 * t67;
t22 = t64 * t53 + t50;
t8 = t22 * t30 + t73 * t67;
t44 = g(1) * t10 + g(2) * t8 + g(3) * t16;
t28 = sin(qJ(6));
t29 = sin(qJ(4));
t66 = t28 * t29;
t31 = cos(qJ(6));
t65 = t29 * t31;
t54 = t62 * t60;
t17 = t71 * t30 + t67 * t51;
t32 = cos(qJ(4));
t37 = -t61 * t52 + t64 * t63;
t12 = t17 * t29 - t37 * t32;
t33 = t38 * t59 - t63 * t54;
t9 = t22 * t67 - t73 * t30;
t4 = t9 * t29 - t33 * t32;
t11 = t23 * t67 - t72 * t30;
t34 = t39 * t59 + t63 * t49;
t6 = t11 * t29 - t34 * t32;
t46 = g(1) * t6 + g(2) * t4 + g(3) * t12;
t13 = t17 * t32 + t37 * t29;
t5 = t33 * t29 + t9 * t32;
t7 = t11 * t32 + t34 * t29;
t45 = g(1) * t7 + g(2) * t5 + g(3) * t13;
t43 = g(1) * t11 + g(2) * t9 + g(3) * t17;
t20 = -g(1) * t49 + g(2) * t54 - g(3) * t64;
t3 = t44 * t32;
t2 = t44 * t29;
t1 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t44, t43, 0, 0, 0, 0, 0, t3, -t2, -t43, -t3, t2, -t43 * pkin(9) + t44 * (pkin(4) * t32 + qJ(5) * t29 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (-t10 * t66 + t11 * t31) - g(2) * (t9 * t31 - t8 * t66) - g(3) * (-t16 * t66 + t17 * t31) -g(1) * (-t10 * t65 - t11 * t28) - g(2) * (-t9 * t28 - t8 * t65) - g(3) * (-t16 * t65 - t17 * t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t45, 0, -t46, -t45, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (-t4 * pkin(4) + t5 * qJ(5)) - g(3) * (-t12 * pkin(4) + t13 * qJ(5)) 0, 0, 0, 0, 0, -t45 * t28, -t45 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t10 * t28 + t6 * t31) - g(2) * (-t8 * t28 + t4 * t31) - g(3) * (t12 * t31 - t16 * t28) -g(1) * (-t10 * t31 - t6 * t28) - g(2) * (-t4 * t28 - t8 * t31) - g(3) * (-t12 * t28 - t16 * t31);];
taug_reg  = t1;
