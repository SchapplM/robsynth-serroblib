% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t34 = sin(qJ(2));
t35 = sin(qJ(1));
t37 = cos(qJ(2));
t38 = cos(qJ(1));
t48 = cos(pkin(6));
t45 = t38 * t48;
t16 = t35 * t34 - t37 * t45;
t33 = sin(qJ(6));
t36 = cos(qJ(6));
t17 = t34 * t45 + t35 * t37;
t29 = pkin(12) + qJ(4);
t28 = qJ(5) + t29;
t24 = sin(t28);
t25 = cos(t28);
t31 = sin(pkin(6));
t51 = t31 * t38;
t8 = t17 * t25 - t24 * t51;
t60 = -t16 * t36 + t8 * t33;
t59 = t16 * t33 + t8 * t36;
t58 = g(3) * t31;
t55 = t25 * t33;
t54 = t25 * t36;
t53 = t31 * t34;
t52 = t31 * t35;
t50 = t33 * t37;
t49 = t36 * t37;
t26 = sin(t29);
t27 = cos(t29);
t47 = t17 * t27 - t26 * t51;
t46 = t35 * t48;
t18 = t38 * t34 + t37 * t46;
t44 = g(1) * t16 - g(2) * t18;
t43 = t17 * t24 + t25 * t51;
t42 = t17 * t26 + t27 * t51;
t19 = -t34 * t46 + t38 * t37;
t10 = -t19 * t24 + t25 * t52;
t41 = g(1) * t10 - g(2) * t43 + g(3) * (-t24 * t53 + t48 * t25);
t40 = -g(1) * t18 - g(2) * t16 + t37 * t58;
t39 = g(1) * t19 + g(2) * t17 + g(3) * t53;
t32 = cos(pkin(12));
t30 = sin(pkin(12));
t15 = t48 * t24 + t25 * t53;
t13 = t19 * t27 + t26 * t52;
t12 = -t19 * t26 + t27 * t52;
t11 = t19 * t25 + t24 * t52;
t6 = t11 * t36 + t18 * t33;
t5 = -t11 * t33 + t18 * t36;
t4 = g(1) * t11 + g(2) * t8 + g(3) * t15;
t2 = t41 * t36;
t1 = t41 * t33;
t3 = [0, g(1) * t35 - g(2) * t38, g(1) * t38 + g(2) * t35, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t19, -t44, -g(1) * (-t17 * t32 + t30 * t51) - g(2) * (t19 * t32 + t30 * t52) -g(1) * (t17 * t30 + t32 * t51) - g(2) * (-t19 * t30 + t32 * t52) t44, -g(1) * (-t35 * pkin(1) - t17 * pkin(2) + pkin(8) * t51 - t16 * qJ(3)) - g(2) * (t38 * pkin(1) + t19 * pkin(2) + pkin(8) * t52 + t18 * qJ(3)) 0, 0, 0, 0, 0, g(1) * t47 - g(2) * t13, -g(1) * t42 - g(2) * t12, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t11, -g(1) * t43 - g(2) * t10, 0, 0, 0, 0, 0, g(1) * t59 - g(2) * t6, -g(1) * t60 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, -t40, t39, -t40 * t32, t40 * t30, -t39, -g(1) * (-t18 * pkin(2) + t19 * qJ(3)) - g(2) * (-t16 * pkin(2) + t17 * qJ(3)) - (pkin(2) * t37 + qJ(3) * t34) * t58, 0, 0, 0, 0, 0, -t40 * t27, t40 * t26, 0, 0, 0, 0, 0, -t40 * t25, t40 * t24, 0, 0, 0, 0, 0, -g(1) * (-t18 * t54 + t19 * t33) - g(2) * (-t16 * t54 + t17 * t33) - (t25 * t49 + t33 * t34) * t58, -g(1) * (t18 * t55 + t19 * t36) - g(2) * (t16 * t55 + t17 * t36) - (-t25 * t50 + t34 * t36) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 + g(2) * t42 - g(3) * (-t26 * t53 + t48 * t27) g(1) * t13 + g(2) * t47 - g(3) * (-t48 * t26 - t27 * t53) 0, 0, 0, 0, 0, -t41, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t60 - g(3) * (-t15 * t33 - t31 * t49) g(1) * t6 + g(2) * t59 - g(3) * (-t15 * t36 + t31 * t50);];
taug_reg  = t3;
