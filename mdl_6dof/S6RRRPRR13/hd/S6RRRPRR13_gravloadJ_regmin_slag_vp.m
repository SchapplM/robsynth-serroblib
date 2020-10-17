% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:47:45
% EndTime: 2019-05-07 15:47:49
% DurationCPUTime: 1.10s
% Computational Cost: add. (751->148), mult. (1889->271), div. (0->0), fcn. (2448->16), ass. (0->75)
t53 = sin(qJ(2));
t55 = cos(qJ(2));
t56 = cos(qJ(1));
t75 = cos(pkin(6));
t68 = t56 * t75;
t86 = sin(qJ(1));
t36 = t86 * t53 - t55 * t68;
t37 = t53 * t68 + t86 * t55;
t52 = sin(qJ(3));
t74 = cos(pkin(7));
t87 = cos(qJ(3));
t63 = t74 * t87;
t48 = sin(pkin(7));
t49 = sin(pkin(6));
t77 = t49 * t56;
t72 = t48 * t77;
t12 = t36 * t63 + t37 * t52 + t87 * t72;
t69 = t52 * t74;
t13 = -t36 * t69 + t37 * t87 - t52 * t72;
t70 = t49 * t74;
t26 = t36 * t48 - t56 * t70;
t46 = pkin(13) + qJ(5);
t44 = sin(t46);
t45 = cos(t46);
t4 = t13 * t45 + t26 * t44;
t51 = sin(qJ(6));
t54 = cos(qJ(6));
t94 = -t12 * t54 + t4 * t51;
t93 = t12 * t51 + t4 * t54;
t90 = t13 * t44 - t26 * t45;
t64 = t75 * t86;
t58 = t56 * t53 + t55 * t64;
t71 = t49 * t86;
t89 = -t48 * t71 + t58 * t74;
t88 = pkin(10) * t48;
t85 = t44 * t48;
t84 = t45 * t48;
t83 = t45 * t51;
t82 = t45 * t54;
t47 = sin(pkin(13));
t81 = t47 * t48;
t50 = cos(pkin(13));
t80 = t48 * t50;
t79 = t49 * t53;
t78 = t49 * t55;
t76 = t53 * t48;
t73 = t49 * t76;
t67 = t75 * t48;
t38 = -t53 * t64 + t56 * t55;
t16 = t38 * t52 + t89 * t87;
t65 = -g(1) * t12 + g(2) * t16;
t24 = t52 * t67 + (t87 * t53 + t55 * t69) * t49;
t35 = -t48 * t78 + t75 * t74;
t17 = t38 * t87 - t89 * t52;
t27 = t58 * t48 + t86 * t70;
t6 = -t17 * t44 + t27 * t45;
t62 = g(1) * t6 - g(2) * t90 + g(3) * (-t24 * t44 + t35 * t45);
t23 = t52 * t79 - t63 * t78 - t87 * t67;
t61 = g(1) * t16 + g(2) * t12 + g(3) * t23;
t60 = g(1) * t17 + g(2) * t13 + g(3) * t24;
t19 = -t36 * t52 + t37 * t63;
t21 = t38 * t63 - t58 * t52;
t33 = (t52 * t55 + t53 * t63) * t49;
t59 = g(1) * t21 + g(2) * t19 + g(3) * t33;
t34 = (-t53 * t69 + t87 * t55) * t49;
t22 = -t38 * t69 - t58 * t87;
t20 = -t36 * t87 - t37 * t69;
t18 = t34 * t45 + t44 * t73;
t11 = t24 * t45 + t35 * t44;
t9 = t22 * t45 + t38 * t85;
t8 = t20 * t45 + t37 * t85;
t7 = t17 * t45 + t27 * t44;
t2 = t16 * t51 + t7 * t54;
t1 = t16 * t54 - t7 * t51;
t3 = [0, g(1) * t86 - g(2) * t56, g(1) * t56 + g(2) * t86, 0, 0, 0, 0, 0, g(1) * t37 - g(2) * t38, -g(1) * t36 + g(2) * t58, 0, 0, 0, 0, 0, g(1) * t13 - g(2) * t17, t65, -g(1) * (-t13 * t50 - t26 * t47) - g(2) * (t17 * t50 + t27 * t47) -g(1) * (t13 * t47 - t26 * t50) - g(2) * (-t17 * t47 + t27 * t50) -t65, -g(1) * (-t86 * pkin(1) - t37 * pkin(2) - pkin(3) * t13 + pkin(9) * t77 - qJ(4) * t12) - g(2) * (t56 * pkin(1) + t38 * pkin(2) + t17 * pkin(3) + pkin(9) * t71 + t16 * qJ(4)) + (g(1) * t26 - g(2) * t27) * pkin(10), 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t7, -g(1) * t90 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t93 - g(2) * t2, -g(1) * t94 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t58 + g(2) * t36 - g(3) * t78, g(1) * t38 + g(2) * t37 + g(3) * t79, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t20 - g(3) * t34, t59, -g(1) * (t22 * t50 + t38 * t81) - g(2) * (t20 * t50 + t37 * t81) - g(3) * (t34 * t50 + t47 * t73) -g(1) * (-t22 * t47 + t38 * t80) - g(2) * (-t20 * t47 + t37 * t80) - g(3) * (-t34 * t47 + t50 * t73) -t59, -g(1) * (-t58 * pkin(2) + t22 * pkin(3) + t21 * qJ(4) + t38 * t88) - g(2) * (-t36 * pkin(2) + t20 * pkin(3) + t19 * qJ(4) + t37 * t88) - g(3) * (t34 * pkin(3) + t33 * qJ(4) + (pkin(2) * t55 + pkin(10) * t76) * t49) 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t8 - g(3) * t18, -g(1) * (-t22 * t44 + t38 * t84) - g(2) * (-t20 * t44 + t37 * t84) - g(3) * (-t34 * t44 + t45 * t73) 0, 0, 0, 0, 0, -g(1) * (t21 * t51 + t9 * t54) - g(2) * (t19 * t51 + t8 * t54) - g(3) * (t18 * t54 + t33 * t51) -g(1) * (t21 * t54 - t9 * t51) - g(2) * (t19 * t54 - t8 * t51) - g(3) * (-t18 * t51 + t33 * t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t60, t61 * t50, -t61 * t47, -t60, -g(1) * (-t16 * pkin(3) + t17 * qJ(4)) - g(2) * (-t12 * pkin(3) + t13 * qJ(4)) - g(3) * (-t23 * pkin(3) + t24 * qJ(4)) 0, 0, 0, 0, 0, t61 * t45, -t61 * t44, 0, 0, 0, 0, 0, -g(1) * (-t16 * t82 + t17 * t51) - g(2) * (-t12 * t82 + t13 * t51) - g(3) * (-t23 * t82 + t24 * t51) -g(1) * (t16 * t83 + t17 * t54) - g(2) * (t12 * t83 + t13 * t54) - g(3) * (t23 * t83 + t24 * t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, g(1) * t7 + g(2) * t4 + g(3) * t11, 0, 0, 0, 0, 0, -t62 * t54, t62 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t94 - g(3) * (-t11 * t51 + t23 * t54) g(1) * t2 + g(2) * t93 - g(3) * (-t11 * t54 - t23 * t51);];
taug_reg  = t3;
