% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPP5
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:28:07
% EndTime: 2019-05-07 18:28:08
% DurationCPUTime: 0.54s
% Computational Cost: add. (511->116), mult. (734->148), div. (0->0), fcn. (775->8), ass. (0->74)
t45 = qJ(3) + qJ(4);
t40 = sin(t45);
t41 = cos(t45);
t51 = cos(qJ(1));
t81 = t51 * t41;
t48 = sin(qJ(1));
t50 = cos(qJ(2));
t83 = t48 * t50;
t15 = t40 * t83 + t81;
t82 = t51 * t40;
t16 = t41 * t83 - t82;
t98 = -t16 * pkin(4) - t15 * qJ(5);
t94 = g(1) * t48;
t61 = -g(2) * t51 + t94;
t92 = g(2) * t48;
t27 = g(1) * t51 + t92;
t47 = sin(qJ(2));
t19 = -g(3) * t50 + t27 * t47;
t97 = -pkin(4) - pkin(5);
t46 = sin(qJ(3));
t96 = pkin(3) * t46;
t95 = pkin(4) * t41;
t90 = g(3) * t47;
t88 = t40 * t47;
t87 = t41 * t47;
t52 = -pkin(9) - pkin(8);
t86 = t47 * t52;
t85 = t48 * t46;
t49 = cos(qJ(3));
t84 = t48 * t49;
t39 = t49 * pkin(3) + pkin(2);
t30 = t50 * t39;
t80 = t51 * t46;
t79 = t51 * t49;
t78 = t51 * pkin(1) + t48 * pkin(7);
t77 = qJ(5) * t40;
t75 = qJ(6) + t52;
t74 = t51 * t86;
t73 = t50 * t80;
t43 = t51 * pkin(7);
t72 = pkin(3) * t80 + t48 * t86 + t43;
t71 = t97 * t40;
t70 = -pkin(1) - t30;
t69 = t75 * t51;
t68 = -t15 * pkin(4) + t16 * qJ(5);
t17 = -t48 * t41 + t50 * t82;
t18 = t48 * t40 + t50 * t81;
t67 = -t17 * pkin(4) + t18 * qJ(5);
t66 = -t39 - t77;
t65 = pkin(3) * t85 + t51 * t30 + t78;
t64 = g(3) * (t30 + (t77 + t95) * t50);
t63 = t50 * pkin(2) + t47 * pkin(8);
t5 = g(1) * t15 - g(2) * t17;
t21 = t46 * t83 + t79;
t58 = t27 * t50;
t56 = t18 * pkin(4) + t17 * qJ(5) + t65;
t2 = g(1) * t17 + g(2) * t15 + g(3) * t88;
t4 = g(1) * t18 + g(2) * t16 + g(3) * t87;
t55 = t70 * t48 + t72;
t36 = pkin(3) * t84;
t54 = -pkin(3) * t73 + t36 + t67;
t20 = t58 + t90;
t53 = -t21 * pkin(3) + t68;
t28 = qJ(5) * t87;
t25 = t61 * t47;
t24 = t50 * t79 + t85;
t23 = -t73 + t84;
t22 = -t49 * t83 + t80;
t12 = t17 * pkin(5);
t9 = t15 * pkin(5);
t8 = t19 * t41;
t7 = t19 * t40;
t6 = g(1) * t16 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, t61, t27, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t50, -t25, -t27, -g(1) * (-t48 * pkin(1) + t43) - g(2) * t78, 0, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t24, -g(1) * t21 - g(2) * t23, t25, -g(1) * t43 - g(2) * (t63 * t51 + t78) - (-pkin(1) - t63) * t94, 0, 0, 0, 0, 0, 0, t6, -t5, t25, -g(1) * t55 - g(2) * (t65 - t74) 0, 0, 0, 0, 0, 0, t6, t25, t5, -g(1) * (t55 + t98) - g(2) * (t56 - t74) 0, 0, 0, 0, 0, 0, t6, t5, -t25, -g(1) * (-t16 * pkin(5) + t72 + t98) - g(2) * (t18 * pkin(5) - t47 * t69 + t56) - (t47 * qJ(6) + t70) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t49, -t19 * t46, -t20, -g(3) * t63 + t27 * (pkin(2) * t47 - pkin(8) * t50) 0, 0, 0, 0, 0, 0, t8, -t7, -t20, -g(3) * (t30 - t86) + t27 * (t39 * t47 + t50 * t52) 0, 0, 0, 0, 0, 0, t8, -t20, t7, -t64 + t52 * t58 + (g(3) * t52 + t27 * (-t66 + t95)) * t47, 0, 0, 0, 0, 0, 0, t8, t7, t20, -t64 + (-g(3) * pkin(5) * t41 + g(1) * t69 + t75 * t92) * t50 + (g(3) * t75 + t27 * (-t97 * t41 - t66)) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t23 + g(2) * t21 + t46 * t90, g(1) * t24 - g(2) * t22 + t49 * t90, 0, 0, 0, 0, 0, 0, 0, 0, t2, t4, 0, -g(1) * t36 + (g(2) * t79 + t20 * t46) * pkin(3), 0, 0, 0, 0, 0, 0, t2, 0, -t4, -g(1) * t54 - g(2) * t53 - g(3) * (t28 + (-pkin(4) * t40 - t96) * t47) 0, 0, 0, 0, 0, 0, t2, -t4, 0, -g(1) * (-t12 + t54) - g(2) * (-t9 + t53) - g(3) * t28 - (t71 - t96) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t4, -g(1) * t67 - g(2) * t68 - g(3) * (-pkin(4) * t88 + t28) 0, 0, 0, 0, 0, 0, t2, -t4, 0, -g(1) * (-t12 + t67) - g(2) * (t68 - t9) - g(3) * (t47 * t71 + t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19;];
taug_reg  = t1;
