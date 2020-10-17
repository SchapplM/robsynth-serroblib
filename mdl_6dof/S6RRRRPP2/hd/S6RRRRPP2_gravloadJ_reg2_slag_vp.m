% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:04:16
% EndTime: 2019-05-07 18:04:18
% DurationCPUTime: 0.48s
% Computational Cost: add. (502->110), mult. (678->136), div. (0->0), fcn. (689->8), ass. (0->70)
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t19 = g(1) * t45 + g(2) * t42;
t39 = qJ(2) + qJ(3);
t36 = sin(t39);
t54 = t19 * t36;
t37 = cos(t39);
t73 = t37 * pkin(3) + t36 * pkin(9);
t8 = -g(3) * t37 + t54;
t90 = -pkin(4) - pkin(5);
t41 = sin(qJ(2));
t89 = pkin(2) * t41;
t46 = -pkin(8) - pkin(7);
t86 = g(2) * t46;
t40 = sin(qJ(4));
t84 = t36 * t40;
t83 = t36 * t42;
t43 = cos(qJ(4));
t82 = t36 * t43;
t81 = t36 * t45;
t80 = t37 * t43;
t79 = t37 * t45;
t78 = t42 * t40;
t77 = t42 * t43;
t76 = t45 * t40;
t75 = t45 * t43;
t74 = t45 * t46;
t72 = qJ(5) * t40;
t71 = qJ(6) * t37;
t70 = t36 * qJ(6);
t69 = t42 * t89;
t68 = t45 * t89;
t44 = cos(qJ(2));
t38 = t44 * pkin(2);
t35 = t38 + pkin(1);
t22 = t45 * t35;
t67 = pkin(3) * t79 + pkin(9) * t81 + t22;
t66 = -pkin(3) - t72;
t15 = t37 * t78 + t75;
t16 = t37 * t77 - t76;
t65 = -t15 * pkin(4) + t16 * qJ(5);
t17 = t37 * t76 - t77;
t18 = t37 * t75 + t78;
t64 = -t17 * pkin(4) + t18 * qJ(5);
t63 = pkin(4) * t80 + t37 * t72 + t73;
t23 = t42 * t37 * pkin(9);
t62 = -t42 * t71 + t23;
t27 = pkin(9) * t79;
t61 = -t45 * t71 + t27;
t60 = t38 + t63;
t59 = -pkin(3) * t36 - t89;
t4 = g(1) * t15 - g(2) * t17;
t58 = g(1) * t42 - g(2) * t45;
t57 = -t35 - t73;
t55 = -t16 * pkin(4) - t15 * qJ(5) - t74;
t53 = t18 * pkin(4) + t17 * qJ(5) + t67;
t2 = g(1) * t17 + g(2) * t15 + g(3) * t84;
t51 = g(1) * t18 + g(2) * t16 + g(3) * t82;
t50 = -g(3) * t44 + t19 * t41;
t49 = (-g(1) * t57 + t86) * t42;
t48 = (pkin(4) * t43 - t66) * t54;
t47 = (g(3) * qJ(6) + t19 * (-t90 * t43 - t66)) * t36;
t24 = pkin(5) * t80;
t20 = qJ(5) * t82;
t14 = g(1) * t83 - g(2) * t81;
t9 = g(3) * t36 + t19 * t37;
t7 = t8 * t43;
t6 = t8 * t40;
t5 = g(1) * t16 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, t58, t19, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t44, -t58 * t41, -t19, -g(1) * (-t42 * pkin(1) + t45 * pkin(7)) - g(2) * (t45 * pkin(1) + t42 * pkin(7)) 0, 0, 0, 0, 0, 0, t58 * t37, -t14, -t19, -g(1) * (-t42 * t35 - t74) - g(2) * (-t42 * t46 + t22) 0, 0, 0, 0, 0, 0, t5, -t4, t14, g(1) * t74 - g(2) * t67 + t49, 0, 0, 0, 0, 0, 0, t5, t14, t4, -g(1) * t55 - g(2) * t53 + t49, 0, 0, 0, 0, 0, 0, t5, t4, -t14, -g(1) * (-t16 * pkin(5) + t55) - g(2) * (t18 * pkin(5) - t45 * t70 + t53) + (-g(1) * (t57 + t70) + t86) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, g(3) * t41 + t19 * t44, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t50 * pkin(2), 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (t45 * t59 + t27) - g(2) * (t42 * t59 + t23) - g(3) * (t38 + t73) 0, 0, 0, 0, 0, 0, t7, -t9, t6, -g(1) * (t27 - t68) - g(2) * (t23 - t69) - g(3) * t60 + t48, 0, 0, 0, 0, 0, 0, t7, t6, t9, -g(1) * (t61 - t68) - g(2) * (t62 - t69) - g(3) * (t24 + t60) + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (-pkin(3) * t81 + t27) - g(2) * (-pkin(3) * t83 + t23) - g(3) * t73, 0, 0, 0, 0, 0, 0, t7, -t9, t6, -g(1) * t27 - g(2) * t23 - g(3) * t63 + t48, 0, 0, 0, 0, 0, 0, t7, t6, t9, -g(1) * t61 - g(2) * t62 - g(3) * (t24 + t63) + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t51, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t51, -g(1) * t64 - g(2) * t65 - g(3) * (-pkin(4) * t84 + t20) 0, 0, 0, 0, 0, 0, t2, -t51, 0, -g(1) * (-t17 * pkin(5) + t64) - g(2) * (-t15 * pkin(5) + t65) - g(3) * (t90 * t84 + t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t1;
