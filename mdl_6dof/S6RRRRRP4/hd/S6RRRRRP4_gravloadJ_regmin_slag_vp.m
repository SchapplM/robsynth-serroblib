% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:51:44
% EndTime: 2019-05-08 04:51:46
% DurationCPUTime: 0.31s
% Computational Cost: add. (471->69), mult. (483->98), div. (0->0), fcn. (515->10), ass. (0->57)
t37 = qJ(2) + qJ(3);
t32 = sin(t37);
t34 = cos(t37);
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t53 = g(1) * t43 + g(2) * t40;
t9 = -g(3) * t34 + t53 * t32;
t36 = qJ(4) + qJ(5);
t31 = sin(t36);
t33 = cos(t36);
t79 = pkin(5) * t33 + qJ(6) * t31;
t42 = cos(qJ(2));
t44 = -pkin(10) - pkin(9);
t78 = t42 * pkin(2) - t32 * t44;
t41 = cos(qJ(4));
t59 = t43 * t41;
t38 = sin(qJ(4));
t64 = t40 * t38;
t18 = t34 * t64 + t59;
t60 = t43 * t38;
t63 = t40 * t41;
t20 = -t34 * t60 + t63;
t71 = g(3) * t32;
t77 = -g(1) * t20 + g(2) * t18 + t38 * t71;
t69 = t31 * t32;
t68 = t32 * t33;
t29 = t41 * pkin(4) + pkin(3);
t24 = t34 * t29;
t66 = t40 * t31;
t65 = t40 * t33;
t62 = t43 * t31;
t61 = t43 * t33;
t56 = t34 * t79 + t24;
t55 = pkin(4) * t38 + pkin(7) + pkin(8);
t13 = t34 * t66 + t61;
t15 = t34 * t62 - t65;
t54 = g(1) * t13 - g(2) * t15;
t52 = g(1) * t40 - g(2) * t43;
t51 = t29 + t79;
t50 = t24 + pkin(1) + t78;
t48 = t53 * t34;
t1 = g(1) * t15 + g(2) * t13 + g(3) * t69;
t14 = t34 * t65 - t62;
t16 = t34 * t61 + t66;
t3 = g(1) * t16 + g(2) * t14 + g(3) * t68;
t46 = -g(1) * (-t15 * pkin(5) + t16 * qJ(6)) - g(2) * (-t13 * pkin(5) + t14 * qJ(6)) - g(3) * (-pkin(5) * t69 + qJ(6) * t68);
t39 = sin(qJ(2));
t21 = t34 * t59 + t64;
t19 = -t34 * t63 + t60;
t17 = t52 * t32;
t10 = t48 + t71;
t8 = t9 * t41;
t7 = t9 * t38;
t6 = t9 * t33;
t5 = t9 * t31;
t4 = g(1) * t14 - g(2) * t16;
t2 = [0, t52, t53, 0, 0, 0, 0, 0, t52 * t42, -t52 * t39, 0, 0, 0, 0, 0, t52 * t34, -t17, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t21, -g(1) * t18 - g(2) * t20, 0, 0, 0, 0, 0, t4, -t54, t4, t17, t54, -g(1) * (-t14 * pkin(5) - t13 * qJ(6)) - g(2) * (t16 * pkin(5) + t15 * qJ(6)) + (-g(1) * t55 - g(2) * t50) * t43 + (g(1) * t50 - g(2) * t55) * t40; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t42 + t53 * t39, g(3) * t39 + t53 * t42, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t6, -t5, t6, -t10, t5, -g(3) * (t56 + t78) + t53 * (pkin(2) * t39 + t51 * t32 + t34 * t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t6, -t5, t6, -t10, t5, -g(3) * t56 + t44 * t48 + (g(3) * t44 + t53 * t51) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, g(1) * t21 - g(2) * t19 + t41 * t71, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t77 * pkin(4) + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
