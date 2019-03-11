% Calculate minimal parameter regressor of gravitation load for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t72 = pkin(3) * t43 + qJ(4) * t41 + pkin(2);
t71 = pkin(4) + pkin(8);
t39 = sin(pkin(6));
t69 = g(3) * t39;
t36 = pkin(11) + qJ(6);
t34 = sin(t36);
t68 = t34 * t41;
t35 = cos(t36);
t67 = t35 * t41;
t37 = sin(pkin(11));
t66 = t37 * t41;
t42 = sin(qJ(2));
t65 = t39 * t42;
t64 = t39 * t43;
t44 = cos(qJ(2));
t63 = t39 * t44;
t40 = cos(pkin(11));
t62 = t40 * t41;
t61 = t41 * t44;
t59 = qJ(5) * t43;
t58 = cos(pkin(6));
t57 = cos(pkin(10));
t38 = sin(pkin(10));
t48 = t58 * t57;
t19 = t38 * t42 - t44 * t48;
t56 = t72 * t19;
t52 = t38 * t58;
t21 = t57 * t42 + t44 * t52;
t55 = t72 * t21;
t20 = t38 * t44 + t42 * t48;
t51 = t39 * t57;
t8 = t20 * t41 + t43 * t51;
t9 = t20 * t43 - t41 * t51;
t54 = -t8 * pkin(3) + t9 * qJ(4);
t22 = -t42 * t52 + t57 * t44;
t10 = t22 * t41 - t38 * t64;
t11 = t38 * t39 * t41 + t22 * t43;
t53 = -t10 * pkin(3) + t11 * qJ(4);
t23 = t41 * t65 - t58 * t43;
t24 = t58 * t41 + t42 * t64;
t50 = -t23 * pkin(3) + t24 * qJ(4);
t49 = pkin(8) * t65 + t72 * t63;
t2 = g(1) * t10 + g(2) * t8 + g(3) * t23;
t47 = g(1) * t11 + g(2) * t9 + g(3) * t24;
t46 = -g(1) * t21 - g(2) * t19 + g(3) * t63;
t45 = g(1) * t22 + g(2) * t20 + g(3) * t65;
t5 = t46 * t43;
t4 = t46 * t41;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t46, t45, 0, 0, 0, 0, 0, -t5, t4, -t45, t5, -t4, -g(1) * (t22 * pkin(8) - t55) - g(2) * (t20 * pkin(8) - t56) - g(3) * t49, -g(1) * (-t21 * t66 + t22 * t40) - g(2) * (-t19 * t66 + t20 * t40) - (t37 * t61 + t40 * t42) * t69, -g(1) * (-t21 * t62 - t22 * t37) - g(2) * (-t19 * t62 - t20 * t37) - (-t37 * t42 + t40 * t61) * t69, -t5, -g(1) * (-t21 * t59 + t71 * t22 - t55) - g(2) * (-t19 * t59 + t71 * t20 - t56) - g(3) * ((pkin(4) * t42 + t44 * t59) * t39 + t49) 0, 0, 0, 0, 0, -g(1) * (-t21 * t68 + t22 * t35) - g(2) * (-t19 * t68 + t20 * t35) - (t34 * t61 + t35 * t42) * t69, -g(1) * (-t21 * t67 - t22 * t34) - g(2) * (-t19 * t67 - t20 * t34) - (-t34 * t42 + t35 * t61) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t47, 0, -t2, -t47, -g(1) * t53 - g(2) * t54 - g(3) * t50, -t47 * t37, -t47 * t40, t2, -g(1) * (-t10 * qJ(5) + t53) - g(2) * (-t8 * qJ(5) + t54) - g(3) * (-t23 * qJ(5) + t50) 0, 0, 0, 0, 0, -t47 * t34, -t47 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t35 - t21 * t34) - g(2) * (-t19 * t34 + t8 * t35) - g(3) * (t23 * t35 + t34 * t63) -g(1) * (-t10 * t34 - t21 * t35) - g(2) * (-t19 * t35 - t8 * t34) - g(3) * (-t23 * t34 + t35 * t63);];
taug_reg  = t1;
