% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:03:57
% EndTime: 2019-05-06 10:03:59
% DurationCPUTime: 0.51s
% Computational Cost: add. (408->92), mult. (898->167), div. (0->0), fcn. (1148->14), ass. (0->60)
t39 = sin(pkin(11));
t44 = sin(qJ(2));
t47 = cos(qJ(2));
t61 = cos(pkin(11));
t28 = -t39 * t47 - t44 * t61;
t42 = cos(pkin(6));
t21 = t28 * t42;
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t52 = -t39 * t44 + t47 * t61;
t10 = -t21 * t48 + t45 * t52;
t37 = pkin(12) + qJ(5);
t35 = sin(t37);
t36 = cos(t37);
t40 = sin(pkin(6));
t67 = t40 * t48;
t4 = t10 * t36 - t35 * t67;
t43 = sin(qJ(6));
t46 = cos(qJ(6));
t49 = t52 * t42;
t9 = t45 * t28 + t48 * t49;
t77 = t4 * t43 + t46 * t9;
t76 = t4 * t46 - t43 * t9;
t63 = t48 * t44;
t64 = t45 * t47;
t25 = -t42 * t64 - t63;
t68 = t40 * t47;
t75 = -g(1) * t25 - g(3) * t68;
t71 = t36 * t43;
t70 = t36 * t46;
t69 = t40 * t45;
t65 = t45 * t44;
t62 = t48 * t47;
t59 = t42 * t62;
t22 = t42 * t44 * pkin(2) + (-pkin(8) - qJ(3)) * t40;
t34 = pkin(2) * t47 + pkin(1);
t58 = -t22 * t45 + t34 * t48;
t56 = g(1) * t48 + g(2) * t45;
t55 = g(1) * t45 - g(2) * t48;
t11 = -t21 * t45 - t48 * t52;
t54 = -t22 * t48 - t34 * t45;
t53 = t10 * t35 + t36 * t67;
t20 = t28 * t40;
t6 = t11 * t35 + t36 * t69;
t51 = g(1) * t6 - g(2) * t53 + g(3) * (t20 * t35 + t36 * t42);
t12 = t28 * t48 - t45 * t49;
t19 = t52 * t40;
t50 = g(1) * t12 + g(2) * t9 + g(3) * t19;
t41 = cos(pkin(12));
t38 = sin(pkin(12));
t30 = pkin(2) * t59;
t26 = -t42 * t65 + t62;
t24 = -t42 * t63 - t64;
t23 = -t59 + t65;
t18 = -g(3) * t42 - t40 * t55;
t15 = -t20 * t36 + t35 * t42;
t7 = -t11 * t36 + t35 * t69;
t2 = -t12 * t43 + t46 * t7;
t1 = -t12 * t46 - t43 * t7;
t3 = [0, t55, t56, 0, 0, 0, 0, 0, -g(1) * t24 - g(2) * t26, -g(1) * t23 - g(2) * t25, -t56 * t40, -g(1) * t54 - g(2) * t58, -g(1) * (-t10 * t41 + t38 * t67) - g(2) * (-t11 * t41 + t38 * t69) -g(1) * (t10 * t38 + t41 * t67) - g(2) * (t11 * t38 + t41 * t69) -g(1) * t9 + g(2) * t12, -g(1) * (-pkin(3) * t10 + qJ(4) * t9 + t54) - g(2) * (-pkin(3) * t11 - qJ(4) * t12 + t58) 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t7, -g(1) * t53 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t76 - g(2) * t2, -g(1) * t77 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t23 + t75, g(3) * t40 * t44 + g(1) * t26 - g(2) * t24, 0, -g(2) * t30 + (g(2) * t65 + t75) * pkin(2), -t50 * t41, t50 * t38, g(1) * t11 - g(2) * t10 + g(3) * t20, -g(1) * (pkin(2) * t25 + t12 * pkin(3) - t11 * qJ(4)) - g(2) * (-pkin(2) * t65 + pkin(3) * t9 + qJ(4) * t10 + t30) - g(3) * (pkin(2) * t68 + pkin(3) * t19 - qJ(4) * t20) 0, 0, 0, 0, 0, -t50 * t36, t50 * t35, 0, 0, 0, 0, 0, -g(1) * (-t11 * t43 + t12 * t70) - g(2) * (t10 * t43 + t70 * t9) - g(3) * (t19 * t70 - t20 * t43) -g(1) * (-t11 * t46 - t12 * t71) - g(2) * (t10 * t46 - t71 * t9) - g(3) * (-t19 * t71 - t20 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, g(1) * t7 + g(2) * t4 + g(3) * t15, 0, 0, 0, 0, 0, -t51 * t46, t51 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t77 - g(3) * (-t15 * t43 - t19 * t46) g(1) * t2 + g(2) * t76 - g(3) * (-t15 * t46 + t19 * t43);];
taug_reg  = t3;
