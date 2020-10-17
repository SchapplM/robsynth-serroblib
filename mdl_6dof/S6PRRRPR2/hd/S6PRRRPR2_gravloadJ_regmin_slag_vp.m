% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:16:23
% EndTime: 2019-05-05 07:16:24
% DurationCPUTime: 0.33s
% Computational Cost: add. (461->86), mult. (732->142), div. (0->0), fcn. (900->14), ass. (0->52)
t35 = sin(pkin(11));
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t51 = cos(pkin(11));
t52 = cos(pkin(6));
t47 = t52 * t51;
t18 = t35 * t39 - t41 * t47;
t50 = t35 * t52;
t20 = t51 * t39 + t41 * t50;
t65 = -g(1) * t20 - g(2) * t18;
t36 = sin(pkin(6));
t62 = g(3) * t36;
t32 = pkin(12) + qJ(6);
t28 = sin(t32);
t33 = qJ(3) + qJ(4);
t31 = cos(t33);
t61 = t28 * t31;
t29 = cos(t32);
t60 = t29 * t31;
t34 = sin(pkin(12));
t59 = t31 * t34;
t37 = cos(pkin(12));
t58 = t31 * t37;
t57 = t31 * t41;
t56 = t35 * t36;
t55 = t36 * t39;
t40 = cos(qJ(3));
t54 = t36 * t40;
t53 = t36 * t41;
t49 = t36 * t51;
t19 = t35 * t41 + t39 * t47;
t21 = -t39 * t50 + t51 * t41;
t48 = g(1) * t21 + g(2) * t19;
t30 = sin(t33);
t11 = t19 * t30 + t31 * t49;
t13 = t21 * t30 - t31 * t56;
t16 = t30 * t55 - t52 * t31;
t5 = g(1) * t13 + g(2) * t11 + g(3) * t16;
t12 = t19 * t31 - t30 * t49;
t14 = t21 * t31 + t30 * t56;
t17 = t52 * t30 + t31 * t55;
t7 = g(1) * t14 + g(2) * t12 + g(3) * t17;
t45 = g(3) * t53 + t65;
t44 = -g(1) * (-t13 * pkin(4) + t14 * qJ(5)) - g(2) * (-t11 * pkin(4) + t12 * qJ(5)) - g(3) * (-t16 * pkin(4) + t17 * qJ(5));
t38 = sin(qJ(3));
t43 = -g(1) * (-t21 * t38 + t35 * t54) - g(2) * (-t19 * t38 - t40 * t49) - g(3) * (-t38 * t55 + t52 * t40);
t8 = t45 * t30;
t4 = t5 * t37;
t3 = t5 * t34;
t2 = t5 * t29;
t1 = t5 * t28;
t6 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t45, g(3) * t55 + t48, 0, 0, 0, 0, 0, -t45 * t40, t45 * t38, 0, 0, 0, 0, 0, -t45 * t31, t8, -g(1) * (-t20 * t58 + t21 * t34) - g(2) * (-t18 * t58 + t19 * t34) - (t34 * t39 + t37 * t57) * t62, -g(1) * (t20 * t59 + t21 * t37) - g(2) * (t18 * t59 + t19 * t37) - (-t34 * t57 + t37 * t39) * t62, -t8 (t39 * t62 + t48) * (-pkin(9) - pkin(8)) + (-t41 * t62 - t65) * (t40 * pkin(3) + pkin(4) * t31 + qJ(5) * t30 + pkin(2)) 0, 0, 0, 0, 0, -g(1) * (-t20 * t60 + t21 * t28) - g(2) * (-t18 * t60 + t19 * t28) - (t28 * t39 + t29 * t57) * t62, -g(1) * (t20 * t61 + t21 * t29) - g(2) * (t18 * t61 + t19 * t29) - (-t28 * t57 + t29 * t39) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -g(1) * (-t21 * t40 - t38 * t56) - g(2) * (-t19 * t40 + t38 * t49) - g(3) * (-t52 * t38 - t39 * t54) 0, 0, 0, 0, 0, t5, t7, t4, -t3, -t7, t43 * pkin(3) + t44, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t7, t4, -t3, -t7, t44, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t28 + t20 * t29) - g(2) * (-t12 * t28 + t18 * t29) - g(3) * (-t17 * t28 - t29 * t53) -g(1) * (-t14 * t29 - t20 * t28) - g(2) * (-t12 * t29 - t18 * t28) - g(3) * (-t17 * t29 + t28 * t53);];
taug_reg  = t6;
