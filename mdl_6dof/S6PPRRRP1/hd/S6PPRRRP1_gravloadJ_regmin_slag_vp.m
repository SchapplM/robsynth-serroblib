% Calculate minimal parameter regressor of gravitation load for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:11
% EndTime: 2021-01-16 00:50:15
% DurationCPUTime: 0.74s
% Computational Cost: add. (609->96), mult. (1602->167), div. (0->0), fcn. (2097->14), ass. (0->75)
t66 = cos(qJ(5));
t52 = t66 * pkin(5) + pkin(4);
t62 = -qJ(6) - pkin(10);
t64 = sin(qJ(4));
t67 = cos(qJ(4));
t104 = -t52 * t67 + t62 * t64;
t58 = cos(pkin(12));
t60 = cos(pkin(7));
t61 = cos(pkin(6));
t79 = t61 * t60;
t56 = sin(pkin(7));
t57 = sin(pkin(6));
t83 = t57 * t56;
t43 = t58 * t79 - t83;
t55 = sin(pkin(11));
t59 = cos(pkin(11));
t54 = sin(pkin(12));
t87 = t54 * t60;
t30 = t43 * t55 + t59 * t87;
t86 = t54 * t61;
t46 = -t55 * t86 + t59 * t58;
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t14 = t30 * t65 - t46 * t68;
t80 = t61 * t56;
t81 = t60 * t57;
t41 = t58 * t80 + t81;
t88 = t54 * t56;
t28 = t41 * t55 + t59 * t88;
t11 = t14 * t64 + t28 * t67;
t33 = -t43 * t59 + t55 * t87;
t45 = t55 * t58 + t59 * t86;
t20 = t33 * t65 - t45 * t68;
t29 = t41 * t59 - t55 * t88;
t9 = t20 * t64 - t29 * t67;
t44 = t58 * t81 + t80;
t84 = t54 * t68;
t34 = t44 * t65 + t57 * t84;
t42 = t58 * t83 - t79;
t98 = t34 * t64 + t67 * t42;
t70 = g(1) * t11 + g(2) * t9 - g(3) * t98;
t85 = t54 * t65;
t37 = -t44 * t68 + t57 * t85;
t15 = t30 * t68 + t46 * t65;
t21 = t33 * t68 + t45 * t65;
t73 = g(1) * t15 + g(2) * t21;
t100 = g(3) * t37 + t73;
t8 = t14 * t67 - t28 * t64;
t13 = t20 * t67 + t29 * t64;
t23 = t34 * t67 - t64 * t42;
t63 = sin(qJ(5));
t97 = t14 * t63;
t96 = t20 * t63;
t82 = t58 * t60;
t27 = t65 * t80 + (t65 * t82 + t84) * t57;
t95 = t27 * t63;
t77 = t63 * t67;
t75 = t66 * t67;
t50 = pkin(3) * t82 + t54 * pkin(9);
t72 = pkin(3) * t83 - t50 * t61;
t49 = -t54 * pkin(3) + pkin(9) * t82;
t71 = pkin(9) * t83 - t49 * t61;
t69 = g(1) * t8 + g(2) * t13 - g(3) * t23;
t26 = -t68 * t80 + (-t68 * t82 + t85) * t57;
t1 = -g(1) * (t15 * t66 + t63 * t8) - g(2) * (t13 * t63 + t21 * t66) - g(3) * (-t23 * t63 + t26 * t66);
t48 = pkin(3) * t87 - t58 * pkin(9);
t47 = t58 * pkin(3) + pkin(9) * t87;
t38 = -g(3) * t61 + (-g(1) * t55 + g(2) * t59) * t57;
t7 = t100 * t64;
t6 = t70 * t66;
t5 = t70 * t63;
t4 = -g(1) * (-t15 * t75 - t97) - g(2) * (-t21 * t75 - t96) - g(3) * (-t37 * t75 + t95);
t3 = -g(1) * (-t14 * t66 + t15 * t77) - g(2) * (-t20 * t66 + t21 * t77) - g(3) * (t27 * t66 + t37 * t77);
t2 = -g(1) * (-t15 * t63 + t66 * t8) - g(2) * (t13 * t66 - t21 * t63) - g(3) * (-t23 * t66 - t26 * t63);
t10 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, g(3) * t26 + t73, -g(1) * t14 - g(2) * t20 + g(3) * t27, 0, 0, 0, 0, 0, t100 * t67, -t7, 0, 0, 0, 0, 0, t4, t3, t4, t3, t7, -g(1) * (-pkin(5) * t97 - (t59 * t47 - t71 * t55) * t65 + (-t59 * t48 + t72 * t55) * t68 + t104 * t15) - g(2) * (-pkin(5) * t96 - (t55 * t47 + t71 * t59) * t65 + (-t55 * t48 - t72 * t59) * t68 + t104 * t21) - g(3) * (pkin(5) * t95 - (-pkin(9) * t80 - t49 * t57) * t65 + (pkin(3) * t80 + t50 * t57) * t68 + t104 * t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69, 0, 0, 0, 0, 0, -t6, t5, -t6, t5, t69, -g(1) * (t11 * t52 + t8 * t62) - g(2) * (t13 * t62 + t9 * t52) - g(3) * (-t23 * t62 - t52 * t98); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70;];
taug_reg = t10;
