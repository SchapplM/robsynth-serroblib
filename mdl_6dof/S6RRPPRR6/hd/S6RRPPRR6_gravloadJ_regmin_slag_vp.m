% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:54:47
% EndTime: 2019-05-06 10:54:47
% DurationCPUTime: 0.25s
% Computational Cost: add. (214->63), mult. (367->98), div. (0->0), fcn. (413->10), ass. (0->53)
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t20 = g(1) * t41 + g(2) * t38;
t37 = sin(qJ(2));
t68 = t20 * t37;
t36 = sin(qJ(6));
t53 = pkin(10) + qJ(5);
t27 = sin(t53);
t40 = cos(qJ(2));
t51 = cos(t53);
t46 = t37 * t51;
t65 = -t40 * t27 + t46;
t3 = t65 * t38;
t42 = t37 * t27 + t40 * t51;
t56 = t40 * t41;
t5 = t27 * t56 - t41 * t46;
t44 = g(1) * t5 - g(2) * t3 + g(3) * t42;
t67 = t44 * t36;
t39 = cos(qJ(6));
t66 = t44 * t39;
t28 = t37 * qJ(3);
t55 = t40 * pkin(2) + t28;
t63 = pkin(2) * t37;
t62 = g(1) * t38;
t59 = g(3) * t65;
t58 = t40 * pkin(3);
t54 = qJ(3) * t40;
t52 = pkin(2) * t56 + t38 * pkin(7) + (pkin(1) + t28) * t41;
t19 = -g(2) * t41 + t62;
t4 = t42 * t38;
t50 = t41 * t36 + t4 * t39;
t49 = t4 * t36 - t41 * t39;
t34 = sin(pkin(10));
t35 = cos(pkin(10));
t48 = t40 * t34 - t37 * t35;
t47 = t37 * t34 + t40 * t35;
t45 = -pkin(1) - t55;
t6 = t42 * t41;
t43 = g(1) * t6 + g(2) * t4 + t59;
t31 = t41 * pkin(7);
t25 = t41 * t54;
t23 = t38 * t54;
t16 = t19 * t40;
t15 = t19 * t37;
t12 = t47 * t41;
t11 = t48 * t41;
t10 = t47 * t38;
t9 = t48 * t38;
t8 = g(3) * t37 + t20 * t40;
t7 = -g(3) * t40 + t68;
t2 = -t38 * t36 + t6 * t39;
t1 = -t6 * t36 - t38 * t39;
t13 = [0, t19, t20, 0, 0, 0, 0, 0, t16, -t15, t16, -t20, t15, -g(1) * t31 - g(2) * t52 - t45 * t62, g(1) * t10 - g(2) * t12, -g(1) * t9 + g(2) * t11, t20, -g(1) * (-t41 * qJ(4) + t31) - g(2) * (pkin(3) * t56 + t52) + (-g(1) * (t45 - t58) + g(2) * qJ(4)) * t38, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t6, g(1) * t3 + g(2) * t5, 0, 0, 0, 0, 0, g(1) * t50 - g(2) * t2, -g(1) * t49 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, t7, 0, -t8, -g(1) * (-t41 * t63 + t25) - g(2) * (-t38 * t63 + t23) - g(3) * t55, -g(1) * t11 - g(2) * t9 - g(3) * t47, -g(1) * t12 - g(2) * t10 + g(3) * t48, 0, -g(1) * t25 - g(2) * t23 - g(3) * (t55 + t58) + (pkin(2) + pkin(3)) * t68, 0, 0, 0, 0, 0, -t44, -t43, 0, 0, 0, 0, 0, -t66, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t43, 0, 0, 0, 0, 0, t66, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t49 + t36 * t59, g(1) * t2 + g(2) * t50 + t39 * t59;];
taug_reg  = t13;
