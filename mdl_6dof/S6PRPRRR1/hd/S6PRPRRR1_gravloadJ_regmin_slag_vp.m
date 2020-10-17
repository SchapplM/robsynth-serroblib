% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:07:45
% EndTime: 2019-05-05 00:07:46
% DurationCPUTime: 0.27s
% Computational Cost: add. (326->67), mult. (684->125), div. (0->0), fcn. (887->14), ass. (0->44)
t28 = sin(pkin(12));
t35 = sin(qJ(2));
t38 = cos(qJ(2));
t47 = cos(pkin(12));
t43 = -t35 * t28 + t38 * t47;
t30 = sin(pkin(6));
t57 = g(3) * t30;
t27 = qJ(4) + qJ(5);
t26 = cos(t27);
t33 = sin(qJ(6));
t56 = t26 * t33;
t36 = cos(qJ(6));
t55 = t26 * t36;
t29 = sin(pkin(11));
t54 = t29 * t30;
t31 = cos(pkin(11));
t53 = t30 * t31;
t34 = sin(qJ(4));
t52 = t30 * t34;
t37 = cos(qJ(4));
t51 = t30 * t37;
t32 = cos(pkin(6));
t50 = t32 * t35;
t49 = t32 * t38;
t22 = -t38 * t28 - t35 * t47;
t20 = t22 * t32;
t45 = -t31 * t20 + t29 * t43;
t44 = t29 * t20 + t31 * t43;
t19 = t22 * t30;
t25 = sin(t27);
t42 = g(1) * (-t25 * t44 + t26 * t54) + g(2) * (-t25 * t45 - t26 * t53) + g(3) * (t19 * t25 + t32 * t26);
t40 = t43 * t32;
t10 = t29 * t22 + t31 * t40;
t13 = t31 * t22 - t29 * t40;
t18 = t43 * t30;
t41 = g(1) * t13 + g(2) * t10 + g(3) * t18;
t39 = -g(1) * (-t29 * t49 - t31 * t35) - g(2) * (-t29 * t35 + t31 * t49) - t38 * t57;
t16 = -t19 * t26 + t32 * t25;
t8 = t25 * t54 + t26 * t44;
t6 = -t25 * t53 + t26 * t45;
t4 = g(1) * t8 + g(2) * t6 + g(3) * t16;
t2 = t42 * t36;
t1 = t42 * t33;
t3 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t39, -g(1) * (t29 * t50 - t31 * t38) - g(2) * (-t29 * t38 - t31 * t50) + t35 * t57, t39 * pkin(2), 0, 0, 0, 0, 0, -t41 * t37, t41 * t34, 0, 0, 0, 0, 0, -t41 * t26, t41 * t25, 0, 0, 0, 0, 0, -g(1) * (t13 * t55 + t33 * t44) - g(2) * (t10 * t55 + t33 * t45) - g(3) * (t18 * t55 - t19 * t33) -g(1) * (-t13 * t56 + t36 * t44) - g(2) * (-t10 * t56 + t36 * t45) - g(3) * (-t18 * t56 - t19 * t36); 0, 0, 0, 0, -g(3) * t32 + (-g(1) * t29 + g(2) * t31) * t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t29 * t51 - t34 * t44) - g(2) * (-t31 * t51 - t34 * t45) - g(3) * (t19 * t34 + t32 * t37) -g(1) * (-t29 * t52 - t37 * t44) - g(2) * (t31 * t52 - t37 * t45) - g(3) * (t19 * t37 - t32 * t34) 0, 0, 0, 0, 0, -t42, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t13 * t36 - t8 * t33) - g(2) * (-t10 * t36 - t6 * t33) - g(3) * (-t16 * t33 - t18 * t36) -g(1) * (t13 * t33 - t8 * t36) - g(2) * (t10 * t33 - t6 * t36) - g(3) * (-t16 * t36 + t18 * t33);];
taug_reg  = t3;
