% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR12
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
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t35 = sin(qJ(2));
t37 = cos(qJ(2));
t48 = cos(pkin(6));
t63 = cos(qJ(1));
t43 = t48 * t63;
t62 = sin(qJ(1));
t16 = t62 * t35 - t37 * t43;
t30 = pkin(12) + qJ(5);
t29 = qJ(6) + t30;
t25 = sin(t29);
t26 = cos(t29);
t17 = t35 * t43 + t62 * t37;
t34 = sin(qJ(3));
t36 = cos(qJ(3));
t32 = sin(pkin(6));
t47 = t32 * t63;
t9 = t17 * t36 - t34 * t47;
t71 = -t16 * t26 + t9 * t25;
t70 = t16 * t25 + t9 * t26;
t27 = sin(t30);
t28 = cos(t30);
t69 = -t16 * t28 + t9 * t27;
t68 = t16 * t27 + t9 * t28;
t42 = t48 * t62;
t18 = t63 * t35 + t37 * t42;
t67 = -g(1) * t18 - g(2) * t16;
t64 = g(3) * t32;
t57 = t25 * t36;
t56 = t26 * t36;
t55 = t27 * t36;
t54 = t28 * t36;
t31 = sin(pkin(12));
t53 = t31 * t36;
t52 = t32 * t35;
t51 = t32 * t37;
t33 = cos(pkin(12));
t50 = t33 * t36;
t49 = t36 * t37;
t46 = t32 * t62;
t19 = -t35 * t42 + t63 * t37;
t12 = t19 * t34 - t36 * t46;
t8 = t17 * t34 + t36 * t47;
t45 = -g(1) * t8 + g(2) * t12;
t44 = -g(1) * t19 - g(2) * t17;
t14 = t34 * t52 - t48 * t36;
t40 = g(1) * t12 + g(2) * t8 + g(3) * t14;
t13 = t19 * t36 + t34 * t46;
t15 = t48 * t34 + t36 * t52;
t39 = g(1) * t13 + g(2) * t9 + g(3) * t15;
t38 = g(3) * t51 + t67;
t7 = t38 * t34;
t6 = t13 * t28 + t18 * t27;
t5 = -t13 * t27 + t18 * t28;
t4 = t13 * t26 + t18 * t25;
t3 = -t13 * t25 + t18 * t26;
t2 = g(1) * t4 + g(2) * t70 - g(3) * (-t15 * t26 + t25 * t51);
t1 = -g(1) * t3 + g(2) * t71 - g(3) * (-t15 * t25 - t26 * t51);
t10 = [0, g(1) * t62 - g(2) * t63, g(1) * t63 + g(2) * t62, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t19, -g(1) * t16 + g(2) * t18, 0, 0, 0, 0, 0, g(1) * t9 - g(2) * t13, t45, -g(1) * (-t16 * t31 - t33 * t9) - g(2) * (t13 * t33 + t18 * t31) -g(1) * (-t16 * t33 + t31 * t9) - g(2) * (-t13 * t31 + t18 * t33) -t45, -g(1) * (-t62 * pkin(1) - t17 * pkin(2) - pkin(3) * t9 + pkin(8) * t47 - t16 * pkin(9) - qJ(4) * t8) - g(2) * (t63 * pkin(1) + t19 * pkin(2) + t13 * pkin(3) + pkin(8) * t46 + t18 * pkin(9) + t12 * qJ(4)) 0, 0, 0, 0, 0, g(1) * t68 - g(2) * t6, -g(1) * t69 - g(2) * t5, 0, 0, 0, 0, 0, g(1) * t70 - g(2) * t4, -g(1) * t71 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, -t38, g(3) * t52 - t44, 0, 0, 0, 0, 0, -t38 * t36, t7, -g(1) * (-t18 * t50 + t19 * t31) - g(2) * (-t16 * t50 + t17 * t31) - (t31 * t35 + t33 * t49) * t64, -g(1) * (t18 * t53 + t19 * t33) - g(2) * (t16 * t53 + t17 * t33) - (-t31 * t49 + t33 * t35) * t64, -t7 (-t35 * t64 + t44) * pkin(9) + (-t37 * t64 - t67) * (pkin(3) * t36 + qJ(4) * t34 + pkin(2)) 0, 0, 0, 0, 0, -g(1) * (-t18 * t54 + t19 * t27) - g(2) * (-t16 * t54 + t17 * t27) - (t27 * t35 + t28 * t49) * t64, -g(1) * (t18 * t55 + t19 * t28) - g(2) * (t16 * t55 + t17 * t28) - (-t27 * t49 + t28 * t35) * t64, 0, 0, 0, 0, 0, -g(1) * (-t18 * t56 + t19 * t25) - g(2) * (-t16 * t56 + t17 * t25) - (t25 * t35 + t26 * t49) * t64, -g(1) * (t18 * t57 + t19 * t26) - g(2) * (t16 * t57 + t17 * t26) - (-t25 * t49 + t26 * t35) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39, t40 * t33, -t40 * t31, -t39, -g(1) * (-t12 * pkin(3) + t13 * qJ(4)) - g(2) * (-t8 * pkin(3) + t9 * qJ(4)) - g(3) * (-t14 * pkin(3) + t15 * qJ(4)) 0, 0, 0, 0, 0, t40 * t28, -t40 * t27, 0, 0, 0, 0, 0, t40 * t26, -t40 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t69 - g(3) * (-t15 * t27 - t28 * t51) g(1) * t6 + g(2) * t68 - g(3) * (-t15 * t28 + t27 * t51) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t10;
