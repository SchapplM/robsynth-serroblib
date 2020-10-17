% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR5
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
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 21:26:11
% EndTime: 2019-05-06 21:26:12
% DurationCPUTime: 0.50s
% Computational Cost: add. (466->89), mult. (1104->170), div. (0->0), fcn. (1449->14), ass. (0->63)
t40 = cos(pkin(6));
t43 = sin(qJ(2));
t48 = cos(qJ(1));
t63 = t48 * t43;
t44 = sin(qJ(1));
t47 = cos(qJ(2));
t65 = t44 * t47;
t27 = -t40 * t65 - t63;
t62 = t48 * t47;
t66 = t44 * t43;
t53 = t40 * t62 - t66;
t39 = sin(pkin(6));
t77 = g(3) * t39;
t84 = -g(1) * t27 - g(2) * t53 - t47 * t77;
t38 = sin(pkin(12));
t61 = cos(pkin(12));
t30 = -t38 * t47 - t43 * t61;
t52 = -t38 * t43 + t47 * t61;
t49 = t52 * t40;
t13 = t44 * t30 + t48 * t49;
t37 = qJ(5) + qJ(6);
t35 = sin(t37);
t36 = cos(t37);
t42 = sin(qJ(4));
t46 = cos(qJ(4));
t23 = t30 * t40;
t56 = -t23 * t48 + t44 * t52;
t69 = t39 * t48;
t8 = -t42 * t69 + t46 * t56;
t83 = t13 * t36 + t35 * t8;
t82 = -t13 * t35 + t36 * t8;
t41 = sin(qJ(5));
t45 = cos(qJ(5));
t81 = t13 * t45 + t41 * t8;
t80 = -t13 * t41 + t45 * t8;
t72 = t35 * t46;
t71 = t36 * t46;
t70 = t39 * t44;
t68 = t41 * t46;
t64 = t45 * t46;
t58 = g(1) * t48 + g(2) * t44;
t57 = g(1) * t44 - g(2) * t48;
t55 = t23 * t44 + t48 * t52;
t54 = t42 * t56 + t46 * t69;
t10 = -t42 * t55 + t46 * t70;
t22 = t30 * t39;
t51 = g(1) * t10 - g(2) * t54 + g(3) * (t22 * t42 + t40 * t46);
t16 = t30 * t48 - t44 * t49;
t21 = t52 * t39;
t50 = g(1) * t16 + g(2) * t13 + g(3) * t21;
t34 = pkin(2) * t47 + pkin(1);
t28 = -t40 * t66 + t62;
t26 = -t40 * t63 - t65;
t24 = t40 * t43 * pkin(2) + (-pkin(8) - qJ(3)) * t39;
t19 = -t22 * t46 + t40 * t42;
t11 = t42 * t70 + t46 * t55;
t6 = t11 * t45 - t16 * t41;
t5 = -t11 * t41 - t16 * t45;
t4 = t11 * t36 - t16 * t35;
t3 = -t11 * t35 - t16 * t36;
t2 = g(1) * t4 + g(2) * t82 - g(3) * (-t19 * t36 + t21 * t35);
t1 = -g(1) * t3 + g(2) * t83 - g(3) * (-t19 * t35 - t21 * t36);
t7 = [0, t57, t58, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t28, g(1) * t53 - g(2) * t27, -t58 * t39, -g(1) * (-t24 * t48 - t34 * t44) - g(2) * (-t24 * t44 + t34 * t48) 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t11, -g(1) * t54 - g(2) * t10, 0, 0, 0, 0, 0, g(1) * t80 - g(2) * t6, -g(1) * t81 - g(2) * t5, 0, 0, 0, 0, 0, g(1) * t82 - g(2) * t4, -g(1) * t83 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t84, g(1) * t28 - g(2) * t26 + t43 * t77, 0, t84 * pkin(2), 0, 0, 0, 0, 0, -t50 * t46, t50 * t42, 0, 0, 0, 0, 0, -g(1) * (t16 * t64 + t41 * t55) - g(2) * (t13 * t64 + t41 * t56) - g(3) * (t21 * t64 - t22 * t41) -g(1) * (-t16 * t68 + t45 * t55) - g(2) * (-t13 * t68 + t45 * t56) - g(3) * (-t21 * t68 - t22 * t45) 0, 0, 0, 0, 0, -g(1) * (t16 * t71 + t35 * t55) - g(2) * (t13 * t71 + t35 * t56) - g(3) * (t21 * t71 - t22 * t35) -g(1) * (-t16 * t72 + t36 * t55) - g(2) * (-t13 * t72 + t36 * t56) - g(3) * (-t21 * t72 - t22 * t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t40 - t39 * t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, g(1) * t11 + g(2) * t8 + g(3) * t19, 0, 0, 0, 0, 0, -t51 * t45, t51 * t41, 0, 0, 0, 0, 0, -t51 * t36, t51 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t81 - g(3) * (-t19 * t41 - t21 * t45) g(1) * t6 + g(2) * t80 - g(3) * (-t19 * t45 + t21 * t41) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t7;
