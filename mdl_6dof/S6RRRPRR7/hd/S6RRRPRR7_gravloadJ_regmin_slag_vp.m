% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR7
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
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t32 = sin(qJ(2));
t33 = sin(qJ(1));
t36 = cos(qJ(2));
t37 = cos(qJ(1));
t50 = cos(pkin(6));
t46 = t37 * t50;
t16 = t33 * t32 - t36 * t46;
t30 = sin(qJ(6));
t34 = cos(qJ(6));
t17 = t32 * t46 + t33 * t36;
t27 = qJ(3) + pkin(12) + qJ(5);
t24 = sin(t27);
t25 = cos(t27);
t28 = sin(pkin(6));
t53 = t28 * t37;
t8 = t17 * t25 - t24 * t53;
t68 = -t16 * t34 + t8 * t30;
t67 = t16 * t30 + t8 * t34;
t66 = g(1) * t37 + g(2) * t33;
t47 = t33 * t50;
t19 = -t32 * t47 + t37 * t36;
t31 = sin(qJ(3));
t35 = cos(qJ(3));
t54 = t28 * t35;
t12 = -t19 * t31 + t33 * t54;
t43 = t17 * t31 + t35 * t53;
t56 = t28 * t32;
t65 = -g(3) * (-t31 * t56 + t50 * t35) + g(2) * t43 - g(1) * t12;
t61 = g(3) * t28;
t58 = t25 * t30;
t57 = t25 * t34;
t55 = t28 * t33;
t52 = t30 * t36;
t51 = t34 * t36;
t48 = t17 * t35 - t31 * t53;
t18 = t37 * t32 + t36 * t47;
t45 = g(1) * t16 - g(2) * t18;
t44 = t17 * t24 + t25 * t53;
t10 = -t19 * t24 + t25 * t55;
t42 = g(1) * t10 - g(2) * t44 + g(3) * (-t24 * t56 + t50 * t25);
t40 = -g(1) * t18 - g(2) * t16 + t36 * t61;
t39 = g(1) * t19 + g(2) * t17 + g(3) * t56;
t29 = -qJ(4) - pkin(9);
t26 = t35 * pkin(3) + pkin(2);
t15 = t50 * t24 + t25 * t56;
t13 = t19 * t35 + t31 * t55;
t11 = t19 * t25 + t24 * t55;
t6 = t11 * t34 + t18 * t30;
t5 = -t11 * t30 + t18 * t34;
t4 = g(1) * t11 + g(2) * t8 + g(3) * t15;
t2 = t42 * t34;
t1 = t42 * t30;
t3 = [0, g(1) * t33 - g(2) * t37, t66, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t19, -t45, 0, 0, 0, 0, 0, g(1) * t48 - g(2) * t13, -g(1) * t43 - g(2) * t12, t45, -g(1) * (-t33 * pkin(1) + t16 * t29 - t17 * t26) - g(2) * (t37 * pkin(1) - t18 * t29 + t19 * t26) - t66 * t28 * (pkin(3) * t31 + pkin(8)) 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t11, -g(1) * t44 - g(2) * t10, 0, 0, 0, 0, 0, g(1) * t67 - g(2) * t6, -g(1) * t68 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, -t40, t39, 0, 0, 0, 0, 0, -t40 * t35, t40 * t31, -t39, -g(1) * (-t18 * t26 - t19 * t29) - g(2) * (-t16 * t26 - t17 * t29) - (t26 * t36 - t29 * t32) * t61, 0, 0, 0, 0, 0, -t40 * t25, t40 * t24, 0, 0, 0, 0, 0, -g(1) * (-t18 * t57 + t19 * t30) - g(2) * (-t16 * t57 + t17 * t30) - (t25 * t51 + t30 * t32) * t61, -g(1) * (t18 * t58 + t19 * t34) - g(2) * (t16 * t58 + t17 * t34) - (-t25 * t52 + t32 * t34) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, g(1) * t13 + g(2) * t48 - g(3) * (-t50 * t31 - t32 * t54) 0, t65 * pkin(3), 0, 0, 0, 0, 0, -t42, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t68 - g(3) * (-t15 * t30 - t28 * t51) g(1) * t6 + g(2) * t67 - g(3) * (-t15 * t34 + t28 * t52);];
taug_reg  = t3;
