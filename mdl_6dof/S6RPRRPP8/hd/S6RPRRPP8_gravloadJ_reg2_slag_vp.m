% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:52:54
% EndTime: 2019-05-05 21:52:55
% DurationCPUTime: 0.35s
% Computational Cost: add. (204->104), mult. (497->118), div. (0->0), fcn. (522->6), ass. (0->59)
t41 = sin(qJ(3));
t40 = sin(qJ(4));
t45 = cos(qJ(1));
t69 = t45 * t40;
t42 = sin(qJ(1));
t43 = cos(qJ(4));
t73 = t42 * t43;
t15 = t41 * t69 + t73;
t68 = t45 * t43;
t74 = t42 * t40;
t16 = t41 * t68 - t74;
t85 = t16 * pkin(4) + t15 * qJ(5);
t44 = cos(qJ(3));
t70 = t44 * t45;
t78 = g(3) * t41;
t84 = -g(2) * t70 - t78;
t83 = -pkin(1) - pkin(7);
t82 = -pkin(5) - pkin(8);
t81 = pkin(8) * t41;
t80 = g(1) * t42;
t79 = g(2) * t45;
t77 = t40 * t44;
t76 = t41 * t42;
t75 = t41 * t45;
t72 = t42 * t44;
t71 = t43 * t44;
t67 = -pkin(4) - qJ(6);
t66 = pkin(3) * t72 + pkin(8) * t76;
t65 = t45 * pkin(1) + t42 * qJ(2);
t62 = pkin(8) * t72;
t61 = t40 * t72;
t60 = t42 * t71;
t59 = t45 * pkin(7) + t65;
t58 = t83 * t42;
t57 = -qJ(5) * t40 - pkin(3);
t13 = t41 * t74 - t68;
t14 = t41 * t73 + t69;
t56 = -t13 * pkin(4) + t14 * qJ(5);
t55 = t15 * pkin(4) - t16 * qJ(5);
t54 = pkin(4) * t60 + qJ(5) * t61 + t66;
t53 = pkin(3) * t76 + t59;
t52 = g(1) * t15 + g(2) * t13;
t51 = g(1) * t16 + g(2) * t14;
t23 = g(1) * t45 + g(2) * t42;
t22 = -t79 + t80;
t35 = t45 * qJ(2);
t50 = pkin(3) * t75 - pkin(8) * t70 + t35;
t49 = -pkin(4) * t43 + t57;
t48 = t14 * pkin(4) + t13 * qJ(5) + t53;
t2 = g(1) * t13 - g(2) * t15 + g(3) * t77;
t47 = g(1) * t14 - g(2) * t16 + g(3) * t71;
t46 = t58 + t50;
t36 = t44 * pkin(8);
t24 = qJ(5) * t71;
t17 = t23 * t44;
t8 = g(1) * t76 - g(2) * t75 + g(3) * t44;
t7 = g(1) * t60 + t84 * t43;
t6 = g(1) * t61 + t84 * t40;
t1 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, -g(1) * (-t42 * pkin(1) + t35) - g(2) * t65, 0, 0, 0, 0, 0, 0, -t23 * t41, -t17, t22, -g(1) * (t35 + t58) - g(2) * t59, 0, 0, 0, 0, 0, 0, -t51, t52, t17, -g(1) * t46 - g(2) * (t53 - t62) 0, 0, 0, 0, 0, 0, t17, t51, -t52, -g(1) * (t46 + t85) - g(2) * (t48 - t62) 0, 0, 0, 0, 0, 0, t17, -t52, -t51, -g(1) * (-pkin(5) * t70 + t16 * qJ(6) + t50 + t85) - g(2) * (t14 * qJ(6) + t48) + (-g(2) * t82 * t44 - g(1) * t83) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22 * t44 + t78, t8, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, -t8, -g(1) * t66 - g(3) * (-t41 * pkin(3) + t36) - (-pkin(3) * t44 - t81) * t79, 0, 0, 0, 0, 0, 0, -t8, t7, -t6, -g(1) * t54 - g(3) * t36 - t49 * t78 - (t49 * t44 - t81) * t79, 0, 0, 0, 0, 0, 0, -t8, -t6, -t7, -g(1) * (qJ(6) * t60 + t54) - g(3) * (t44 * pkin(5) + t36) + (-pkin(5) * t80 - g(3) * (-qJ(6) * t43 + t49)) * t41 - (t82 * t41 + (t67 * t43 + t57) * t44) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t47, -g(1) * t56 - g(2) * t55 - g(3) * (-pkin(4) * t77 + t24) 0, 0, 0, 0, 0, 0, 0, -t47, t2, -g(1) * (-t13 * qJ(6) + t56) - g(2) * (t15 * qJ(6) + t55) - g(3) * (t67 * t77 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47;];
taug_reg  = t1;
