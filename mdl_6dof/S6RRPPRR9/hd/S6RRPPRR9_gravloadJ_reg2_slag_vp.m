% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:35:45
% EndTime: 2019-05-06 11:35:47
% DurationCPUTime: 0.52s
% Computational Cost: add. (379->121), mult. (950->171), div. (0->0), fcn. (1134->10), ass. (0->71)
t46 = sin(qJ(2));
t47 = sin(qJ(1));
t50 = cos(qJ(2));
t51 = cos(qJ(1));
t75 = cos(pkin(6));
t70 = t51 * t75;
t26 = t47 * t46 - t50 * t70;
t71 = t47 * t75;
t28 = t51 * t46 + t50 * t71;
t99 = g(1) * t28 + g(2) * t26;
t44 = sin(qJ(6));
t48 = cos(qJ(6));
t27 = t46 * t70 + t47 * t50;
t45 = sin(qJ(5));
t49 = cos(qJ(5));
t43 = sin(pkin(6));
t87 = t43 * t51;
t59 = -t27 * t45 + t49 * t87;
t98 = -t26 * t48 + t44 * t59;
t97 = t26 * t44 + t48 * t59;
t94 = g(3) * t43;
t91 = t43 * t46;
t90 = t43 * t47;
t89 = t43 * t49;
t88 = t43 * t50;
t86 = t44 * t45;
t85 = t44 * t50;
t84 = t45 * t48;
t83 = t48 * t50;
t82 = qJ(3) - pkin(9);
t81 = pkin(2) * t88 + qJ(3) * t91;
t80 = t51 * pkin(1) + pkin(8) * t90;
t79 = t26 * qJ(3);
t78 = t27 * qJ(3);
t77 = t28 * qJ(3);
t29 = -t46 * t71 + t51 * t50;
t76 = t29 * qJ(3);
t74 = t29 * pkin(2) + t80;
t73 = qJ(4) * t88 + t81;
t72 = -t47 * pkin(1) + pkin(8) * t87;
t20 = t26 * pkin(2);
t69 = -t20 + t78;
t22 = t28 * pkin(2);
t68 = -t22 + t76;
t67 = -t26 * qJ(4) - t20;
t66 = -t28 * qJ(4) - t22;
t65 = g(3) * t73;
t64 = -t27 * pkin(2) + t72;
t63 = pkin(5) * t45 - pkin(10) * t49;
t11 = -t29 * t49 + t45 * t90;
t58 = t27 * t49 + t45 * t87;
t62 = g(1) * t58 + g(2) * t11;
t9 = g(1) * t26 - g(2) * t28;
t10 = g(1) * t27 - g(2) * t29;
t61 = g(1) * t51 + g(2) * t47;
t57 = pkin(3) * t90 + t29 * qJ(4) + t74;
t24 = -t75 * t45 + t46 * t89;
t56 = g(1) * t11 - g(2) * t58 - g(3) * t24;
t12 = t29 * t45 + t47 * t89;
t25 = t45 * t91 + t75 * t49;
t55 = g(1) * t12 - g(2) * t59 + g(3) * t25;
t54 = pkin(3) * t87 - t27 * qJ(4) + t64;
t4 = g(3) * t88 - t99;
t7 = g(1) * t29 + g(2) * t27 + g(3) * t91;
t53 = pkin(4) * t90 + t82 * t28 + t57;
t52 = pkin(4) * t87 - t82 * t26 + t54;
t30 = t61 * t43;
t3 = t4 * t49;
t2 = t12 * t48 - t28 * t44;
t1 = -t12 * t44 - t28 * t48;
t5 = [0, 0, 0, 0, 0, 0, g(1) * t47 - g(2) * t51, t61, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, -t30, -g(1) * t72 - g(2) * t80, 0, 0, 0, 0, 0, 0, -t30, -t10, t9, -g(1) * (t64 - t79) - g(2) * (t74 + t77) 0, 0, 0, 0, 0, 0, -t30, t9, t10, -g(1) * (t54 - t79) - g(2) * (t57 + t77) 0, 0, 0, 0, 0, 0, -g(1) * t59 - g(2) * t12, t62, -t9, -g(1) * t52 - g(2) * t53, 0, 0, 0, 0, 0, 0, -g(1) * t97 - g(2) * t2, g(1) * t98 - g(2) * t1, -t62, -g(1) * (pkin(5) * t59 + pkin(10) * t58 + t52) - g(2) * (t12 * pkin(5) + t11 * pkin(10) + t53); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t7, -g(1) * t68 - g(2) * t69 - g(3) * t81, 0, 0, 0, 0, 0, 0, 0, -t7, -t4, -g(1) * (t66 + t76) - g(2) * (t67 + t78) - t65, 0, 0, 0, 0, 0, 0, -t4 * t45, -t3, t7, -g(1) * (t82 * t29 + t66) - g(2) * (t82 * t27 + t67) - g(3) * (-pkin(9) * t91 + t73) 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t84 - t29 * t44) - g(2) * (-t26 * t84 - t27 * t44) - (-t44 * t46 + t45 * t83) * t94, -g(1) * (t28 * t86 - t29 * t48) - g(2) * (t26 * t86 - t27 * t48) - (-t45 * t85 - t46 * t48) * t94, t3, -g(1) * (-t29 * pkin(9) + t68) - g(2) * (-t27 * pkin(9) + t69) - t65 - (-pkin(9) * t46 + t63 * t50) * t94 + t99 * (qJ(4) + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t55, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t48, -t56 * t44, -t55, -g(1) * (-t11 * pkin(5) + t12 * pkin(10)) - g(2) * (pkin(5) * t58 - pkin(10) * t59) - g(3) * (t24 * pkin(5) + t25 * pkin(10)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t98 - g(3) * (-t25 * t44 + t43 * t83) g(1) * t2 - g(2) * t97 - g(3) * (-t25 * t48 - t43 * t85) 0, 0;];
taug_reg  = t5;
