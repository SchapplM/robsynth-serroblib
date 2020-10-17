% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP10
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
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:32:27
% EndTime: 2019-05-08 06:32:30
% DurationCPUTime: 0.68s
% Computational Cost: add. (703->122), mult. (1419->196), div. (0->0), fcn. (1797->12), ass. (0->73)
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t81 = cos(pkin(6));
t96 = cos(qJ(1));
t73 = t81 * t96;
t95 = sin(qJ(1));
t42 = t59 * t73 + t95 * t62;
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t56 = sin(pkin(6));
t78 = t56 * t96;
t27 = t42 * t61 - t58 * t78;
t41 = t95 * t59 - t62 * t73;
t55 = qJ(4) + qJ(5);
t53 = sin(t55);
t54 = cos(t55);
t10 = t27 * t53 - t41 * t54;
t11 = t27 * t54 + t41 * t53;
t57 = sin(qJ(4));
t60 = cos(qJ(4));
t107 = t27 * t57 - t41 * t60;
t106 = t27 * t60 + t41 * t57;
t72 = t81 * t95;
t43 = t96 * t59 + t62 * t72;
t105 = -g(1) * t43 - g(2) * t41;
t44 = -t59 * t72 + t96 * t62;
t77 = t56 * t95;
t31 = t44 * t61 + t58 * t77;
t16 = -t31 * t57 + t43 * t60;
t86 = t56 * t59;
t40 = t81 * t58 + t61 * t86;
t85 = t56 * t62;
t104 = g(2) * t107 - g(3) * (-t40 * t57 - t60 * t85) - g(1) * t16;
t30 = t44 * t58 - t61 * t77;
t76 = -t42 * t58 - t61 * t78;
t68 = -g(2) * t76 - g(3) * (-t58 * t86 + t81 * t61) + g(1) * t30;
t97 = g(3) * t56;
t90 = t42 * t57;
t89 = t44 * t57;
t88 = t53 * t61;
t87 = t54 * t61;
t84 = t57 * t61;
t83 = t60 * t61;
t82 = t61 * t62;
t80 = t53 * t85;
t79 = pkin(4) * t57 + pkin(9);
t14 = t31 * t53 - t43 * t54;
t75 = -g(1) * t10 + g(2) * t14;
t74 = g(1) * t76 + g(2) * t30;
t24 = t40 * t53 + t54 * t85;
t1 = g(1) * t14 + g(2) * t10 + g(3) * t24;
t15 = t31 * t54 + t43 * t53;
t25 = t40 * t54 - t80;
t3 = g(1) * t15 + g(2) * t11 + g(3) * t25;
t19 = -t41 * t88 - t42 * t54;
t21 = -t43 * t88 - t44 * t54;
t32 = -t54 * t86 + t61 * t80;
t69 = g(1) * t21 + g(2) * t19 + g(3) * t32;
t67 = g(1) * t31 + g(2) * t27 + g(3) * t40;
t66 = g(3) * t85 + t105;
t64 = -g(1) * (-t14 * pkin(5) + t15 * qJ(6)) - g(2) * (-t10 * pkin(5) + t11 * qJ(6)) - g(3) * (-t24 * pkin(5) + t25 * qJ(6));
t63 = -pkin(11) - pkin(10);
t52 = t60 * pkin(4) + pkin(3);
t33 = (t53 * t59 + t54 * t82) * t56;
t22 = -t43 * t87 + t44 * t53;
t20 = -t41 * t87 + t42 * t53;
t18 = t66 * t58;
t17 = t31 * t60 + t43 * t57;
t7 = t68 * t54;
t6 = t68 * t53;
t5 = g(1) * t11 - g(2) * t15;
t4 = -g(1) * t22 - g(2) * t20 - g(3) * t33;
t2 = [0, g(1) * t95 - g(2) * t96, g(1) * t96 + g(2) * t95, 0, 0, 0, 0, 0, g(1) * t42 - g(2) * t44, -g(1) * t41 + g(2) * t43, 0, 0, 0, 0, 0, g(1) * t27 - g(2) * t31, t74, 0, 0, 0, 0, 0, g(1) * t106 - g(2) * t17, -g(1) * t107 - g(2) * t16, 0, 0, 0, 0, 0, t5, t75, t5, -t74, -t75, -g(1) * (-t95 * pkin(1) - t42 * pkin(2) - pkin(5) * t11 + pkin(8) * t78 - qJ(6) * t10 - t27 * t52 - t79 * t41 - t63 * t76) - g(2) * (t96 * pkin(1) + t44 * pkin(2) + t15 * pkin(5) + pkin(8) * t77 + t14 * qJ(6) - t30 * t63 + t31 * t52 + t79 * t43); 0, 0, 0, 0, 0, 0, 0, 0, -t66, g(1) * t44 + g(2) * t42 + g(3) * t86, 0, 0, 0, 0, 0, -t66 * t61, t18, 0, 0, 0, 0, 0, -g(1) * (-t43 * t83 + t89) - g(2) * (-t41 * t83 + t90) - (t57 * t59 + t60 * t82) * t97, -g(1) * (t43 * t84 + t44 * t60) - g(2) * (t41 * t84 + t42 * t60) - (-t57 * t82 + t59 * t60) * t97, 0, 0, 0, 0, 0, t4, t69, t4, -t18, -t69, -g(1) * (pkin(4) * t89 + t22 * pkin(5) + t44 * pkin(9) + t21 * qJ(6)) - g(2) * (pkin(4) * t90 + t20 * pkin(5) + t42 * pkin(9) + t19 * qJ(6)) - g(3) * (t33 * pkin(5) + t32 * qJ(6)) - t79 * t59 * t97 + (-t62 * t97 - t105) * (t52 * t61 - t58 * t63 + pkin(2)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t67, 0, 0, 0, 0, 0, t68 * t60, -t68 * t57, 0, 0, 0, 0, 0, t7, -t6, t7, -t67, t6, t67 * t63 + t68 * (pkin(5) * t54 + qJ(6) * t53 + t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, g(1) * t17 + g(2) * t106 - g(3) * (-t40 * t60 + t57 * t85) 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t104 * pkin(4) + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
