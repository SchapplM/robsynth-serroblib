% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR4
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 20:47:43
% EndTime: 2019-05-06 20:47:45
% DurationCPUTime: 0.42s
% Computational Cost: add. (428->81), mult. (908->150), div. (0->0), fcn. (1177->14), ass. (0->62)
t43 = cos(pkin(6));
t46 = sin(qJ(2));
t51 = cos(qJ(1));
t68 = t51 * t46;
t47 = sin(qJ(1));
t50 = cos(qJ(2));
t69 = t47 * t50;
t29 = -t43 * t69 - t68;
t67 = t51 * t50;
t70 = t47 * t46;
t56 = t43 * t67 - t70;
t42 = sin(pkin(6));
t78 = g(3) * t42;
t83 = -g(1) * t29 - g(2) * t56 - t50 * t78;
t41 = sin(pkin(12));
t66 = cos(pkin(12));
t32 = -t50 * t41 - t46 * t66;
t55 = -t46 * t41 + t50 * t66;
t52 = t55 * t43;
t15 = t47 * t32 + t51 * t52;
t44 = sin(qJ(6));
t48 = cos(qJ(6));
t40 = qJ(4) + qJ(5);
t38 = sin(t40);
t39 = cos(t40);
t25 = t32 * t43;
t60 = -t51 * t25 + t47 * t55;
t72 = t42 * t51;
t8 = -t38 * t72 + t39 * t60;
t82 = t15 * t48 + t8 * t44;
t81 = -t15 * t44 + t8 * t48;
t75 = t39 * t44;
t74 = t39 * t48;
t73 = t42 * t47;
t45 = sin(qJ(4));
t49 = cos(qJ(4));
t64 = -t45 * t72 + t49 * t60;
t62 = g(1) * t51 + g(2) * t47;
t61 = g(1) * t47 - g(2) * t51;
t59 = t47 * t25 + t51 * t55;
t58 = t38 * t60 + t39 * t72;
t57 = t45 * t60 + t49 * t72;
t10 = -t38 * t59 + t39 * t73;
t24 = t32 * t42;
t54 = g(1) * t10 - g(2) * t58 + g(3) * (t24 * t38 + t43 * t39);
t18 = t51 * t32 - t47 * t52;
t23 = t55 * t42;
t53 = g(1) * t18 + g(2) * t15 + g(3) * t23;
t37 = t50 * pkin(2) + pkin(1);
t30 = -t43 * t70 + t67;
t28 = -t43 * t68 - t69;
t26 = t43 * t46 * pkin(2) + (-pkin(8) - qJ(3)) * t42;
t21 = -t24 * t39 + t43 * t38;
t13 = t45 * t73 + t49 * t59;
t12 = -t45 * t59 + t49 * t73;
t11 = t38 * t73 + t39 * t59;
t6 = t11 * t48 - t18 * t44;
t5 = -t11 * t44 - t18 * t48;
t4 = g(1) * t11 + g(2) * t8 + g(3) * t21;
t2 = t54 * t48;
t1 = t54 * t44;
t3 = [0, t61, t62, 0, 0, 0, 0, 0, -g(1) * t28 - g(2) * t30, g(1) * t56 - g(2) * t29, -t62 * t42, -g(1) * (-t51 * t26 - t47 * t37) - g(2) * (-t47 * t26 + t51 * t37) 0, 0, 0, 0, 0, g(1) * t64 - g(2) * t13, -g(1) * t57 - g(2) * t12, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t11, -g(1) * t58 - g(2) * t10, 0, 0, 0, 0, 0, g(1) * t81 - g(2) * t6, -g(1) * t82 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t83, g(1) * t30 - g(2) * t28 + t46 * t78, 0, t83 * pkin(2), 0, 0, 0, 0, 0, -t53 * t49, t53 * t45, 0, 0, 0, 0, 0, -t53 * t39, t53 * t38, 0, 0, 0, 0, 0, -g(1) * (t18 * t74 + t44 * t59) - g(2) * (t15 * t74 + t44 * t60) - g(3) * (t23 * t74 - t24 * t44) -g(1) * (-t18 * t75 + t48 * t59) - g(2) * (-t15 * t75 + t48 * t60) - g(3) * (-t23 * t75 - t24 * t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t43 - t61 * t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 + g(2) * t57 - g(3) * (t24 * t45 + t43 * t49) g(1) * t13 + g(2) * t64 - g(3) * (t24 * t49 - t43 * t45) 0, 0, 0, 0, 0, -t54, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t82 - g(3) * (-t21 * t44 - t23 * t48) g(1) * t6 + g(2) * t81 - g(3) * (-t21 * t48 + t23 * t44);];
taug_reg  = t3;
