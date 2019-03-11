% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t69 = pkin(3) * t42 + qJ(4) * t39 + pkin(2);
t36 = sin(pkin(6));
t67 = g(3) * t36;
t41 = cos(qJ(5));
t34 = t41 * pkin(5) + pkin(4);
t66 = pkin(8) + t34;
t40 = sin(qJ(2));
t65 = t36 * t40;
t64 = t36 * t42;
t43 = cos(qJ(2));
t63 = t36 * t43;
t38 = sin(qJ(5));
t62 = t38 * t39;
t61 = t38 * t43;
t60 = t39 * t41;
t59 = t41 * t43;
t57 = cos(pkin(6));
t56 = cos(pkin(10));
t35 = sin(pkin(10));
t49 = t57 * t56;
t19 = t35 * t40 - t43 * t49;
t55 = t69 * t19;
t52 = t35 * t57;
t21 = t56 * t40 + t43 * t52;
t54 = t69 * t21;
t53 = pkin(5) * t38 + qJ(4);
t51 = t36 * t56;
t50 = g(3) * (pkin(8) * t65 + t69 * t63);
t37 = -qJ(6) - pkin(9);
t48 = pkin(5) * t62 - t37 * t42;
t22 = -t40 * t52 + t56 * t43;
t10 = t22 * t39 - t35 * t64;
t23 = t39 * t65 - t57 * t42;
t20 = t35 * t43 + t40 * t49;
t8 = t20 * t39 + t42 * t51;
t2 = g(1) * t10 + g(2) * t8 + g(3) * t23;
t11 = t35 * t36 * t39 + t22 * t42;
t24 = t57 * t39 + t40 * t64;
t9 = t20 * t42 - t39 * t51;
t47 = g(1) * t11 + g(2) * t9 + g(3) * t24;
t46 = -g(1) * t21 - g(2) * t19 + g(3) * t63;
t45 = g(1) * t22 + g(2) * t20 + g(3) * t65;
t44 = -g(1) * (t10 * t41 - t21 * t38) - g(2) * (-t19 * t38 + t8 * t41) - g(3) * (t23 * t41 + t36 * t61);
t18 = t23 * pkin(3);
t7 = t10 * pkin(3);
t6 = t8 * pkin(3);
t5 = t46 * t42;
t4 = t46 * t39;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t46, t45, 0, 0, 0, 0, 0, -t5, t4, -t45, t5, -t4, -g(1) * (t22 * pkin(8) - t54) - g(2) * (t20 * pkin(8) - t55) - t50, 0, 0, 0, 0, 0, -g(1) * (-t21 * t62 + t22 * t41) - g(2) * (-t19 * t62 + t20 * t41) - (t39 * t61 + t40 * t41) * t67, -g(1) * (-t21 * t60 - t22 * t38) - g(2) * (-t19 * t60 - t20 * t38) - (-t38 * t40 + t39 * t59) * t67, -t5, -g(1) * (-t48 * t21 + t66 * t22 - t54) - g(2) * (-t48 * t19 + t66 * t20 - t55) - t50 - (t34 * t40 + t48 * t43) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t47, 0, -t2, -t47, -g(1) * (t11 * qJ(4) - t7) - g(2) * (t9 * qJ(4) - t6) - g(3) * (t24 * qJ(4) - t18) 0, 0, 0, 0, 0, -t47 * t38, -t47 * t41, t2, -g(1) * (t10 * t37 + t53 * t11 - t7) - g(2) * (t8 * t37 + t53 * t9 - t6) - g(3) * (t23 * t37 + t53 * t24 - t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -g(1) * (-t10 * t38 - t21 * t41) - g(2) * (-t19 * t41 - t8 * t38) - g(3) * (-t23 * t38 + t36 * t59) 0, t44 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47;];
taug_reg  = t1;
