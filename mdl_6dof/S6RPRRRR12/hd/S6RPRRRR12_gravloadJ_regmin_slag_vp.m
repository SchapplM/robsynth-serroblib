% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_gravloadJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 06:43:26
% EndTime: 2019-05-06 06:43:30
% DurationCPUTime: 1.24s
% Computational Cost: add. (1234->121), mult. (3569->231), div. (0->0), fcn. (4743->18), ass. (0->83)
t51 = sin(pkin(8));
t56 = cos(pkin(6));
t53 = cos(pkin(14));
t65 = cos(qJ(1));
t85 = t65 * t53;
t50 = sin(pkin(14));
t61 = sin(qJ(1));
t89 = t61 * t50;
t46 = -t56 * t85 + t89;
t83 = sin(pkin(7));
t52 = sin(pkin(6));
t55 = cos(pkin(7));
t93 = t52 * t55;
t73 = t46 * t83 - t65 * t93;
t100 = t73 * t51;
t86 = t65 * t50;
t88 = t61 * t53;
t47 = t56 * t86 + t88;
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t81 = t52 * t83;
t79 = t65 * t81;
t35 = (t46 * t55 + t79) * t64 + t47 * t60;
t91 = t55 * t60;
t36 = t46 * t91 - t47 * t64 + t60 * t79;
t59 = sin(qJ(4));
t54 = cos(pkin(8));
t98 = cos(qJ(4));
t82 = t54 * t98;
t14 = t100 * t98 - t35 * t82 + t36 * t59;
t15 = (t35 * t54 - t100) * t59 + t36 * t98;
t25 = t35 * t51 + t54 * t73;
t58 = sin(qJ(5));
t63 = cos(qJ(5));
t5 = t15 * t63 - t25 * t58;
t57 = sin(qJ(6));
t62 = cos(qJ(6));
t110 = -t14 * t62 + t5 * t57;
t109 = t14 * t57 + t5 * t62;
t108 = t15 * t58 + t25 * t63;
t95 = t51 * t58;
t94 = t51 * t63;
t92 = t54 * t59;
t90 = t57 * t63;
t87 = t62 * t63;
t84 = qJ(2) * t52;
t80 = t83 * t56;
t78 = g(1) * t65 + g(2) * t61;
t77 = g(1) * t61 - g(2) * t65;
t76 = t56 * t88 + t86;
t42 = t64 * t80 + (t53 * t55 * t64 - t50 * t60) * t52;
t43 = t60 * t80 + (t50 * t64 + t53 * t91) * t52;
t72 = -t53 * t81 + t55 * t56;
t69 = t72 * t51;
t24 = t43 * t98 + (t42 * t54 + t69) * t59;
t32 = -t42 * t51 + t54 * t72;
t48 = -t56 * t89 + t85;
t68 = -t55 * t76 + t61 * t81;
t37 = -t48 * t60 + t64 * t68;
t38 = t48 * t64 + t60 * t68;
t67 = t61 * t93 + t76 * t83;
t66 = t67 * t51;
t17 = t38 * t98 + (t37 * t54 + t66) * t59;
t27 = -t37 * t51 + t54 * t67;
t6 = -t17 * t58 + t27 * t63;
t75 = g(1) * t6 + g(2) * t108 + g(3) * (-t24 * t58 + t32 * t63);
t16 = -t37 * t82 + t38 * t59 - t66 * t98;
t23 = -t42 * t82 + t43 * t59 - t69 * t98;
t74 = g(1) * t16 - g(2) * t14 + g(3) * t23;
t29 = t42 * t98 - t43 * t92;
t28 = t42 * t59 + t43 * t82;
t22 = t29 * t63 + t43 * t95;
t21 = t37 * t98 - t38 * t92;
t20 = t37 * t59 + t38 * t82;
t19 = -t35 * t98 + t36 * t92;
t18 = -t35 * t59 - t36 * t82;
t11 = t24 * t63 + t32 * t58;
t9 = t21 * t63 + t38 * t95;
t8 = t19 * t63 - t36 * t95;
t7 = t17 * t63 + t27 * t58;
t2 = t16 * t57 + t62 * t7;
t1 = t16 * t62 - t57 * t7;
t3 = [0, t77, t78, g(1) * t47 - g(2) * t48, -g(1) * t46 + g(2) * t76, -t78 * t52, -g(1) * (-pkin(1) * t61 + t65 * t84) - g(2) * (pkin(1) * t65 + t61 * t84) 0, 0, 0, 0, 0, -g(1) * t36 - g(2) * t38, -g(1) * t35 - g(2) * t37, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, g(1) * t14 + g(2) * t16, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, g(1) * t108 - g(2) * t6, 0, 0, 0, 0, 0, -g(1) * t109 - g(2) * t2, g(1) * t110 - g(2) * t1; 0, 0, 0, 0, 0, 0, -g(3) * t56 - t52 * t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t37 + g(2) * t35 - g(3) * t42, g(1) * t38 - g(2) * t36 + g(3) * t43, 0, 0, 0, 0, 0, -g(1) * t21 - g(2) * t19 - g(3) * t29, g(1) * t20 + g(2) * t18 + g(3) * t28, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t8 - g(3) * t22, -g(1) * (-t21 * t58 + t38 * t94) - g(2) * (-t19 * t58 - t36 * t94) - g(3) * (-t29 * t58 + t43 * t94) 0, 0, 0, 0, 0, -g(1) * (t20 * t57 + t62 * t9) - g(2) * (t18 * t57 + t62 * t8) - g(3) * (t22 * t62 + t28 * t57) -g(1) * (t20 * t62 - t57 * t9) - g(2) * (t18 * t62 - t57 * t8) - g(3) * (-t22 * t57 + t28 * t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, g(1) * t17 - g(2) * t15 + g(3) * t24, 0, 0, 0, 0, 0, t74 * t63, -t74 * t58, 0, 0, 0, 0, 0, -g(1) * (-t16 * t87 + t17 * t57) - g(2) * (t14 * t87 - t15 * t57) - g(3) * (-t23 * t87 + t24 * t57) -g(1) * (t16 * t90 + t17 * t62) - g(2) * (-t14 * t90 - t15 * t62) - g(3) * (t23 * t90 + t24 * t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, g(1) * t7 - g(2) * t5 + g(3) * t11, 0, 0, 0, 0, 0, -t75 * t62, t75 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t110 - g(3) * (-t11 * t57 + t23 * t62) g(1) * t2 - g(2) * t109 - g(3) * (-t11 * t62 - t23 * t57);];
taug_reg  = t3;
