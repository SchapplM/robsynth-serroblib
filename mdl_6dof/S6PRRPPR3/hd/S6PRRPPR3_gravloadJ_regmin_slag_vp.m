% Calculate minimal parameter regressor of gravitation load for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:00:31
% EndTime: 2019-05-05 03:00:32
% DurationCPUTime: 0.31s
% Computational Cost: add. (245->82), mult. (645->126), div. (0->0), fcn. (784->10), ass. (0->49)
t38 = sin(qJ(3));
t69 = qJ(4) * t38 + pkin(2);
t35 = sin(pkin(6));
t68 = g(3) * t35;
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t36 = cos(pkin(10));
t55 = cos(pkin(6));
t48 = t36 * t55;
t54 = sin(pkin(10));
t21 = -t39 * t54 + t42 * t48;
t41 = cos(qJ(3));
t67 = t21 * t41;
t45 = t55 * t54;
t23 = -t36 * t39 - t42 * t45;
t66 = t23 * t41;
t65 = t35 * t39;
t64 = t35 * t41;
t63 = t35 * t42;
t37 = sin(qJ(6));
t62 = t37 * t38;
t61 = t37 * t42;
t40 = cos(qJ(6));
t60 = t38 * t40;
t59 = t40 * t42;
t58 = t41 * t42;
t57 = pkin(8) - qJ(5);
t53 = pkin(3) * t67 + t21 * t69;
t52 = pkin(3) * t66 + t23 * t69;
t22 = t39 * t48 + t42 * t54;
t10 = t22 * t38 + t36 * t64;
t11 = -t35 * t36 * t38 + t22 * t41;
t51 = -pkin(3) * t10 + qJ(4) * t11;
t24 = t36 * t42 - t39 * t45;
t49 = t35 * t54;
t12 = t24 * t38 - t41 * t49;
t13 = t24 * t41 + t38 * t49;
t50 = -pkin(3) * t12 + qJ(4) * t13;
t25 = t38 * t65 - t41 * t55;
t26 = t38 * t55 + t39 * t64;
t47 = -pkin(3) * t25 + qJ(4) * t26;
t46 = t35 * pkin(3) * t58 + pkin(8) * t65 + t63 * t69;
t2 = g(1) * t12 + g(2) * t10 + g(3) * t25;
t44 = g(1) * t13 + g(2) * t11 + g(3) * t26;
t43 = g(1) * t23 + g(2) * t21 + g(3) * t63;
t7 = g(1) * t24 + g(2) * t22 + g(3) * t65;
t5 = t43 * t41;
t4 = t43 * t38;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t43, t7, 0, 0, 0, 0, 0, -t5, t4, -t5, -t7, -t4, -g(1) * (pkin(8) * t24 + t52) - g(2) * (pkin(8) * t22 + t53) - g(3) * t46, -t4, t5, t7, -g(1) * (pkin(4) * t66 + t24 * t57 + t52) - g(2) * (pkin(4) * t67 + t22 * t57 + t53) - g(3) * ((pkin(4) * t58 - qJ(5) * t39) * t35 + t46) 0, 0, 0, 0, 0, -g(1) * (t23 * t60 - t24 * t37) - g(2) * (t21 * t60 - t22 * t37) - (-t37 * t39 + t38 * t59) * t68, -g(1) * (-t23 * t62 - t24 * t40) - g(2) * (-t21 * t62 - t22 * t40) - (-t38 * t61 - t39 * t40) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t44, t2, 0, -t44, -g(1) * t50 - g(2) * t51 - g(3) * t47, -t44, -t2, 0, -g(1) * (-pkin(4) * t12 + t50) - g(2) * (-pkin(4) * t10 + t51) - g(3) * (-pkin(4) * t25 + t47) 0, 0, 0, 0, 0, -t44 * t40, t44 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t37 + t23 * t40) - g(2) * (-t10 * t37 + t21 * t40) - g(3) * (-t25 * t37 + t35 * t59) -g(1) * (-t12 * t40 - t23 * t37) - g(2) * (-t10 * t40 - t21 * t37) - g(3) * (-t25 * t40 - t35 * t61);];
taug_reg  = t1;
