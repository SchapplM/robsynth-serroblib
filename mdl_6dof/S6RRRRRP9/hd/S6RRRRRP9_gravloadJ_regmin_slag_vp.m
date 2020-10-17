% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:11:28
% EndTime: 2019-05-08 06:11:30
% DurationCPUTime: 0.56s
% Computational Cost: add. (451->110), mult. (927->190), div. (0->0), fcn. (1147->12), ass. (0->60)
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t41 = cos(qJ(2));
t51 = cos(pkin(6));
t65 = cos(qJ(1));
t47 = t51 * t65;
t16 = t38 * t37 - t41 * t47;
t33 = qJ(4) + qJ(5);
t29 = sin(t33);
t30 = cos(t33);
t17 = t37 * t47 + t38 * t41;
t36 = sin(qJ(3));
t40 = cos(qJ(3));
t34 = sin(pkin(6));
t50 = t34 * t65;
t9 = t17 * t40 - t36 * t50;
t74 = -t16 * t30 + t9 * t29;
t73 = t16 * t29 + t9 * t30;
t35 = sin(qJ(4));
t39 = cos(qJ(4));
t72 = -t16 * t39 + t9 * t35;
t71 = t16 * t35 + t9 * t39;
t56 = t34 * t40;
t15 = t51 * t36 + t37 * t56;
t49 = t38 * t51;
t19 = -t37 * t49 + t65 * t41;
t57 = t34 * t38;
t13 = t19 * t40 + t36 * t57;
t18 = t65 * t37 + t41 * t49;
t3 = -t13 * t29 + t18 * t30;
t55 = t34 * t41;
t1 = g(2) * t74 - g(3) * (-t15 * t29 - t30 * t55) - g(1) * t3;
t69 = g(2) * t16;
t68 = g(2) * t17;
t67 = g(3) * t34;
t21 = t35 * pkin(4) + pkin(5) * t29;
t66 = pkin(9) + t21;
t60 = t29 * t40;
t59 = t30 * t40;
t58 = t34 * t37;
t54 = t35 * t40;
t53 = t39 * t40;
t52 = t40 * t41;
t22 = t39 * pkin(4) + pkin(5) * t30;
t12 = t19 * t36 - t38 * t56;
t8 = t17 * t36 + t40 * t50;
t48 = -g(1) * t8 + g(2) * t12;
t20 = pkin(3) + t22;
t32 = -qJ(6) - pkin(11) - pkin(10);
t46 = t20 * t40 - t32 * t36 + pkin(2);
t14 = t36 * t58 - t51 * t40;
t45 = g(1) * t12 + g(2) * t8 + g(3) * t14;
t44 = g(1) * t13 + g(2) * t9 + g(3) * t15;
t43 = -g(1) * t18 + g(3) * t55 - t69;
t7 = t43 * t36;
t6 = t13 * t39 + t18 * t35;
t5 = -t13 * t35 + t18 * t39;
t4 = t13 * t30 + t18 * t29;
t2 = g(1) * t4 + g(2) * t73 - g(3) * (-t15 * t30 + t29 * t55);
t10 = [0, g(1) * t38 - g(2) * t65, g(1) * t65 + g(2) * t38, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t19, -g(1) * t16 + g(2) * t18, 0, 0, 0, 0, 0, g(1) * t9 - g(2) * t13, t48, 0, 0, 0, 0, 0, g(1) * t71 - g(2) * t6, -g(1) * t72 - g(2) * t5, 0, 0, 0, 0, 0, g(1) * t73 - g(2) * t4, -g(1) * t74 - g(2) * t3, -t48, -g(1) * (-t38 * pkin(1) - t17 * pkin(2) + pkin(8) * t50 - t66 * t16 - t20 * t9 + t32 * t8) - g(2) * (t65 * pkin(1) + t19 * pkin(2) + pkin(8) * t57 - t12 * t32 + t13 * t20 + t66 * t18); 0, 0, 0, 0, 0, 0, 0, 0, -t43, g(1) * t19 + g(3) * t58 + t68, 0, 0, 0, 0, 0, -t43 * t40, t7, 0, 0, 0, 0, 0, -g(1) * (-t18 * t53 + t19 * t35) - g(2) * (-t16 * t53 + t17 * t35) - (t35 * t37 + t39 * t52) * t67, -g(1) * (t18 * t54 + t19 * t39) - g(2) * (t16 * t54 + t17 * t39) - (-t35 * t52 + t37 * t39) * t67, 0, 0, 0, 0, 0, -g(1) * (-t18 * t59 + t19 * t29) - g(2) * (-t16 * t59 + t17 * t29) - (t29 * t37 + t30 * t52) * t67, -g(1) * (t18 * t60 + t19 * t30) - g(2) * (t16 * t60 + t17 * t30) - (-t29 * t52 + t30 * t37) * t67, -t7, -g(1) * (-t46 * t18 + t66 * t19) - t66 * t68 + t46 * t69 - (t66 * t37 + t46 * t41) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t44, 0, 0, 0, 0, 0, t45 * t39, -t45 * t35, 0, 0, 0, 0, 0, t45 * t30, -t45 * t29, -t44, -g(1) * (-t12 * t20 - t13 * t32) - g(2) * (-t8 * t20 - t9 * t32) - g(3) * (-t14 * t20 - t15 * t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t72 - g(3) * (-t15 * t35 - t39 * t55) g(1) * t6 + g(2) * t71 - g(3) * (-t15 * t39 + t35 * t55) 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t13 * t21 + t18 * t22) - g(2) * (t16 * t22 - t9 * t21) - g(3) * (-t15 * t21 - t22 * t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45;];
taug_reg  = t10;
