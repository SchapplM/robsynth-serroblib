% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t32 = qJ(3) + qJ(4);
t30 = sin(t32);
t31 = cos(t32);
t33 = sin(qJ(6));
t37 = cos(qJ(6));
t47 = t30 * t37 - t31 * t33;
t36 = sin(qJ(1));
t39 = cos(qJ(2));
t40 = cos(qJ(1));
t58 = t40 * t30;
t18 = -t36 * t31 + t39 * t58;
t57 = t40 * t31;
t19 = t36 * t30 + t39 * t57;
t6 = t18 * t37 - t19 * t33;
t35 = sin(qJ(2));
t64 = g(3) * t35;
t59 = t36 * t39;
t16 = t30 * t59 + t57;
t17 = t31 * t59 - t58;
t77 = t16 * t37 - t17 * t33;
t2 = g(1) * t6 + g(2) * t77 + t47 * t64;
t50 = g(1) * t40 + g(2) * t36;
t79 = -g(3) * t39 + t50 * t35;
t34 = sin(qJ(3));
t38 = cos(qJ(3));
t55 = t40 * t38;
t21 = t34 * t59 + t55;
t56 = t40 * t34;
t23 = t36 * t38 - t39 * t56;
t75 = -g(1) * t23 + g(2) * t21 + t34 * t64;
t46 = t30 * t33 + t31 * t37;
t48 = t16 * t33 + t17 * t37;
t7 = t18 * t33 + t19 * t37;
t74 = g(1) * t7 + g(2) * t48 + t46 * t64;
t61 = t30 * t35;
t60 = t31 * t35;
t53 = pkin(3) * t34 + pkin(7);
t51 = g(1) * t16 - g(2) * t18;
t49 = g(1) * t36 - g(2) * t40;
t29 = t38 * pkin(3) + pkin(2);
t41 = -pkin(9) - pkin(8);
t45 = t39 * t29 - t35 * t41 + pkin(1);
t44 = pkin(4) * t31 + qJ(5) * t30 + t29;
t3 = g(1) * t18 + g(2) * t16 + g(3) * t61;
t5 = g(1) * t19 + g(2) * t17 + g(3) * t60;
t42 = -g(1) * (-t18 * pkin(4) + t19 * qJ(5)) - g(2) * (-t16 * pkin(4) + t17 * qJ(5)) - g(3) * (-pkin(4) * t61 + qJ(5) * t60);
t25 = t49 * t35;
t24 = t36 * t34 + t39 * t55;
t22 = -t38 * t59 + t56;
t20 = t50 * t39 + t64;
t10 = t79 * t31;
t9 = t79 * t30;
t8 = g(1) * t17 - g(2) * t19;
t1 = [0, t49, t50, 0, 0, 0, 0, 0, t49 * t39, -t25, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t24, -g(1) * t21 - g(2) * t23, 0, 0, 0, 0, 0, t8, -t51, t8, t25, t51, -g(1) * (-t17 * pkin(4) - t16 * qJ(5)) - g(2) * (t19 * pkin(4) + t18 * qJ(5)) + (-g(1) * t53 - g(2) * t45) * t40 + (g(1) * t45 - g(2) * t53) * t36, 0, 0, 0, 0, 0, g(1) * t48 - g(2) * t7, g(1) * t77 - g(2) * t6; 0, 0, 0, 0, 0, 0, 0, 0, t79, t20, 0, 0, 0, 0, 0, t79 * t38, -t79 * t34, 0, 0, 0, 0, 0, t10, -t9, t10, -t20, t9 (-g(3) * t44 + t50 * t41) * t39 + (g(3) * t41 + t44 * t50) * t35, 0, 0, 0, 0, 0, t79 * t46, t79 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, g(1) * t24 - g(2) * t22 + t38 * t64, 0, 0, 0, 0, 0, t3, t5, t3, 0, -t5, pkin(3) * t75 + t42, 0, 0, 0, 0, 0, t2, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, t3, 0, -t5, t42, 0, 0, 0, 0, 0, t2, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t74;];
taug_reg  = t1;
