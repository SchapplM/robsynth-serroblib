% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:06:30
% EndTime: 2019-05-05 10:06:32
% DurationCPUTime: 0.48s
% Computational Cost: add. (639->119), mult. (1780->209), div. (0->0), fcn. (2295->14), ass. (0->75)
t83 = cos(pkin(12));
t85 = cos(pkin(6));
t71 = t85 * t83;
t80 = sin(pkin(12));
t88 = sin(qJ(2));
t90 = cos(qJ(2));
t54 = -t71 * t90 + t80 * t88;
t81 = sin(pkin(7));
t82 = sin(pkin(6));
t67 = t82 * t81;
t84 = cos(pkin(7));
t98 = t54 * t84 + t67 * t83;
t69 = t85 * t80;
t55 = t69 * t90 + t83 * t88;
t66 = t82 * t80;
t97 = t55 * t84 - t66 * t81;
t68 = t84 * t82;
t96 = t68 * t90 + t81 * t85;
t33 = t71 * t88 + t80 * t90;
t45 = sin(qJ(3));
t89 = cos(qJ(3));
t10 = t33 * t45 + t89 * t98;
t34 = -t69 * t88 + t83 * t90;
t12 = t34 * t45 + t89 * t97;
t72 = t82 * t88;
t24 = t45 * t72 - t89 * t96;
t61 = g(1) * t12 + g(2) * t10 + g(3) * t24;
t11 = t33 * t89 - t45 * t98;
t13 = t34 * t89 - t45 * t97;
t25 = t45 * t96 + t72 * t89;
t95 = g(1) * t13 + g(2) * t11 + g(3) * t25;
t43 = sin(qJ(5));
t47 = cos(qJ(4));
t87 = t43 * t47;
t46 = cos(qJ(5));
t86 = t46 * t47;
t79 = pkin(5) * t43 + pkin(10);
t78 = pkin(9) * t81;
t44 = sin(qJ(4));
t77 = t44 * t81;
t76 = t45 * t84;
t75 = t47 * t81;
t74 = t84 * t89;
t73 = t90 * t82;
t53 = -t67 * t90 + t84 * t85;
t14 = t25 * t44 - t47 * t53;
t48 = t54 * t81 - t68 * t83;
t2 = t11 * t44 - t47 * t48;
t49 = t55 * t81 + t66 * t84;
t4 = t13 * t44 - t47 * t49;
t64 = g(1) * t4 + g(2) * t2 + g(3) * t14;
t15 = t25 * t47 + t44 * t53;
t3 = t11 * t47 + t44 * t48;
t5 = t13 * t47 + t44 * t49;
t63 = g(1) * t5 + g(2) * t3 + g(3) * t15;
t59 = t88 * t68;
t31 = -t45 * t59 + t73 * t89;
t58 = t88 * t67;
t20 = t31 * t44 - t47 * t58;
t17 = -t33 * t76 - t54 * t89;
t6 = t17 * t44 - t33 * t75;
t19 = -t34 * t76 - t55 * t89;
t8 = t19 * t44 - t34 * t75;
t62 = g(1) * t8 + g(2) * t6 + g(3) * t20;
t50 = -g(1) * (t12 * t46 - t43 * t5) - g(2) * (t10 * t46 - t3 * t43) - g(3) * (-t15 * t43 + t24 * t46);
t42 = -qJ(6) - pkin(11);
t41 = pkin(5) * t46 + pkin(4);
t30 = t45 * t73 + t59 * t89;
t21 = t31 * t47 + t44 * t58;
t18 = t34 * t74 - t45 * t55;
t16 = t33 * t74 - t45 * t54;
t9 = t19 * t47 + t34 * t77;
t7 = t17 * t47 + t33 * t77;
t1 = t61 * t44;
t22 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, g(1) * t55 + g(2) * t54 - g(3) * t73, g(1) * t34 + g(2) * t33 + g(3) * t72, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t17 - g(3) * t31, g(1) * t18 + g(2) * t16 + g(3) * t30, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7 - g(3) * t21, t62, 0, 0, 0, 0, 0, -g(1) * (t18 * t43 + t46 * t9) - g(2) * (t16 * t43 + t46 * t7) - g(3) * (t21 * t46 + t30 * t43) -g(1) * (t18 * t46 - t43 * t9) - g(2) * (t16 * t46 - t43 * t7) - g(3) * (-t21 * t43 + t30 * t46) -t62, -g(1) * (-pkin(2) * t55 + t19 * pkin(3) + t18 * t79 + t34 * t78 + t9 * t41 - t8 * t42) - g(2) * (-pkin(2) * t54 + t17 * pkin(3) + t16 * t79 + t33 * t78 + t7 * t41 - t6 * t42) - g(3) * (pkin(2) * t73 + t31 * pkin(3) + pkin(9) * t58 - t20 * t42 + t21 * t41 + t30 * t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t95, 0, 0, 0, 0, 0, t61 * t47, -t1, 0, 0, 0, 0, 0, -g(1) * (-t12 * t86 + t13 * t43) - g(2) * (-t10 * t86 + t11 * t43) - g(3) * (-t24 * t86 + t25 * t43) -g(1) * (t12 * t87 + t13 * t46) - g(2) * (t10 * t87 + t11 * t46) - g(3) * (t24 * t87 + t25 * t46) t1, -t95 * t79 + t61 * (t41 * t47 - t42 * t44 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t63, 0, 0, 0, 0, 0, t64 * t46, -t64 * t43, -t63, -g(1) * (-t4 * t41 - t42 * t5) - g(2) * (-t2 * t41 - t3 * t42) - g(3) * (-t14 * t41 - t15 * t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -g(1) * (-t12 * t43 - t46 * t5) - g(2) * (-t10 * t43 - t3 * t46) - g(3) * (-t15 * t46 - t24 * t43) 0, t50 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64;];
taug_reg  = t22;
