% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:24:37
% EndTime: 2019-05-07 12:24:39
% DurationCPUTime: 0.53s
% Computational Cost: add. (408->98), mult. (710->182), div. (0->0), fcn. (888->14), ass. (0->61)
t34 = sin(qJ(2));
t35 = sin(qJ(1));
t38 = cos(qJ(2));
t39 = cos(qJ(1));
t51 = cos(pkin(6));
t47 = t39 * t51;
t16 = t35 * t34 - t38 * t47;
t29 = qJ(5) + qJ(6);
t26 = sin(t29);
t27 = cos(t29);
t17 = t34 * t47 + t35 * t38;
t28 = qJ(3) + pkin(12);
t24 = sin(t28);
t25 = cos(t28);
t30 = sin(pkin(6));
t54 = t30 * t39;
t9 = -t17 * t25 + t24 * t54;
t77 = t16 * t27 + t9 * t26;
t76 = -t16 * t26 + t9 * t27;
t32 = sin(qJ(5));
t36 = cos(qJ(5));
t75 = t16 * t36 + t9 * t32;
t74 = -t16 * t32 + t9 * t36;
t73 = g(1) * t39 + g(2) * t35;
t48 = t35 * t51;
t19 = -t34 * t48 + t39 * t38;
t33 = sin(qJ(3));
t37 = cos(qJ(3));
t56 = t30 * t37;
t12 = -t19 * t33 + t35 * t56;
t45 = t17 * t33 + t37 * t54;
t58 = t30 * t34;
t72 = -g(3) * (-t33 * t58 + t51 * t37) + g(2) * t45 - g(1) * t12;
t68 = g(3) * t30;
t63 = t25 * t26;
t62 = t25 * t27;
t61 = t25 * t32;
t60 = t25 * t36;
t59 = t25 * t38;
t57 = t30 * t35;
t55 = t30 * t38;
t53 = t32 * t38;
t52 = t36 * t38;
t49 = t17 * t37 - t33 * t54;
t18 = t39 * t34 + t38 * t48;
t46 = g(1) * t16 - g(2) * t18;
t44 = g(1) * (-t19 * t24 + t25 * t57) + g(2) * (-t17 * t24 - t25 * t54) + g(3) * (-t24 * t58 + t51 * t25);
t42 = -g(1) * t18 - g(2) * t16 + g(3) * t55;
t41 = g(1) * t19 + g(2) * t17 + g(3) * t58;
t31 = -qJ(4) - pkin(9);
t23 = t37 * pkin(3) + pkin(2);
t15 = t51 * t24 + t25 * t58;
t13 = t19 * t37 + t33 * t57;
t11 = t19 * t25 + t24 * t57;
t6 = t11 * t36 + t18 * t32;
t5 = -t11 * t32 + t18 * t36;
t4 = t11 * t27 + t18 * t26;
t3 = -t11 * t26 + t18 * t27;
t2 = g(1) * t4 - g(2) * t76 - g(3) * (-t15 * t27 + t26 * t55);
t1 = -g(1) * t3 - g(2) * t77 - g(3) * (-t15 * t26 - t27 * t55);
t7 = [0, g(1) * t35 - g(2) * t39, t73, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t19, -t46, 0, 0, 0, 0, 0, g(1) * t49 - g(2) * t13, -g(1) * t45 - g(2) * t12, t46, -g(1) * (-t35 * pkin(1) + t16 * t31 - t17 * t23) - g(2) * (t39 * pkin(1) - t18 * t31 + t19 * t23) - t73 * t30 * (pkin(3) * t33 + pkin(8)) 0, 0, 0, 0, 0, -g(1) * t74 - g(2) * t6, g(1) * t75 - g(2) * t5, 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t4, g(1) * t77 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, -t42, t41, 0, 0, 0, 0, 0, -t42 * t37, t42 * t33, -t41, -g(1) * (-t18 * t23 - t19 * t31) - g(2) * (-t16 * t23 - t17 * t31) - (t23 * t38 - t31 * t34) * t68, 0, 0, 0, 0, 0, -g(1) * (-t18 * t60 + t19 * t32) - g(2) * (-t16 * t60 + t17 * t32) - (t25 * t52 + t32 * t34) * t68, -g(1) * (t18 * t61 + t19 * t36) - g(2) * (t16 * t61 + t17 * t36) - (-t25 * t53 + t34 * t36) * t68, 0, 0, 0, 0, 0, -g(1) * (-t18 * t62 + t19 * t26) - g(2) * (-t16 * t62 + t17 * t26) - (t26 * t34 + t27 * t59) * t68, -g(1) * (t18 * t63 + t19 * t27) - g(2) * (t16 * t63 + t17 * t27) - (-t26 * t59 + t27 * t34) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, g(1) * t13 + g(2) * t49 - g(3) * (-t51 * t33 - t34 * t56) 0, t72 * pkin(3), 0, 0, 0, 0, 0, -t44 * t36, t44 * t32, 0, 0, 0, 0, 0, -t44 * t27, t44 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t75 - g(3) * (-t15 * t32 - t30 * t52) g(1) * t6 - g(2) * t74 - g(3) * (-t15 * t36 + t30 * t53) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t7;
