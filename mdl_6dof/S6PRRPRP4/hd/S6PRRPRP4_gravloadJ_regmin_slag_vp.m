% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:12:47
% EndTime: 2021-01-16 03:12:51
% DurationCPUTime: 0.76s
% Computational Cost: add. (332->116), mult. (822->196), div. (0->0), fcn. (951->10), ass. (0->86)
t47 = sin(qJ(2));
t48 = cos(qJ(5));
t73 = t48 * pkin(5) + pkin(4) + pkin(8);
t113 = t47 * t73;
t50 = cos(qJ(2));
t45 = sin(qJ(5));
t37 = t45 * pkin(5) + qJ(4);
t40 = qJ(6) + pkin(3) + pkin(9);
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t68 = t37 * t46 + t40 * t49;
t117 = (pkin(2) + t68) * t50 + t113;
t42 = sin(pkin(6));
t88 = t46 * t47;
t44 = cos(pkin(6));
t91 = t44 * t49;
t25 = t42 * t88 - t91;
t96 = t42 * t50;
t116 = g(3) * (t25 * t48 + t45 * t96);
t43 = cos(pkin(10));
t114 = t43 * t47;
t103 = g(3) * t25;
t41 = sin(pkin(10));
t92 = t44 * t47;
t20 = t41 * t92 - t43 * t50;
t22 = t41 * t50 + t43 * t92;
t84 = t49 * t42;
t76 = t43 * t84;
t79 = t41 * t84;
t111 = g(2) * (t22 * t46 + t76) + t103 - g(1) * (t20 * t46 + t79);
t83 = t49 * t47;
t93 = t44 * t46;
t102 = g(3) * (t42 * t83 + t93);
t74 = t44 * t83;
t89 = t46 * t42;
t27 = t74 - t89;
t82 = t49 * t50;
t110 = -t102 - g(2) * (t43 * t27 + t41 * t82);
t109 = t47 * (t37 * t91 - t40 * t93) - t37 * t89 - t40 * t84;
t80 = qJ(4) * t49;
t81 = qJ(4) * t46;
t108 = t47 * (pkin(3) * t93 - t44 * t80) + pkin(3) * t84 + t42 * t81;
t71 = t49 * pkin(3) + t81;
t107 = t50 * (pkin(2) + t71) + pkin(8) * t47;
t101 = g(3) * t42;
t99 = pkin(2) * t114;
t98 = t41 * t44;
t97 = t42 * t47;
t95 = t43 * t44;
t94 = t44 * t45;
t90 = t44 * t50;
t87 = t46 * t48;
t86 = t46 * t50;
t85 = t48 * t50;
t75 = t44 * t87;
t24 = t44 * t88 + t84;
t14 = t43 * t24 + t41 * t86;
t16 = -t41 * t24 + t43 * t86;
t70 = -t46 * pkin(3) + t80;
t69 = t37 * t49 - t40 * t46;
t65 = -t41 * t47 + t43 * t90;
t23 = t41 * t90 + t114;
t64 = t50 * t70;
t63 = t50 * t69;
t62 = t71 * t44;
t28 = t41 * t89;
t56 = g(1) * (t20 * t49 - t28) - g(2) * (t22 * t49 - t43 * t89) - t102;
t33 = t43 * t82;
t55 = -g(1) * (-t41 * t27 + t33) + t110;
t54 = -g(1) * t20 + g(2) * t22 + g(3) * t97;
t53 = g(1) * t23 - g(2) * t65 - g(3) * t96;
t9 = t53 * t49;
t52 = -t47 * t68 + t50 * t73;
t38 = t41 * pkin(2);
t34 = pkin(2) * t95;
t19 = t23 * t46;
t18 = t65 * t46;
t8 = t53 * t46;
t7 = -g(1) * t16 - g(2) * t14 - t103;
t6 = t55 * t48;
t5 = t55 * t45;
t4 = -g(1) * (-t19 * t48 + t20 * t45) - g(2) * (t18 * t48 - t22 * t45) - (-t45 * t47 + t46 * t85) * t101;
t3 = -g(1) * (-t19 * t45 - t20 * t48) - g(2) * (t18 * t45 + t22 * t48) - (t45 * t86 + t47 * t48) * t101;
t2 = -g(1) * (t16 * t48 - t23 * t45) - g(2) * (t14 * t48 + t45 * t65) - t116;
t1 = -g(1) * (-t16 * t45 - t23 * t48) - g(2) * (-t14 * t45 + t48 * t65) - g(3) * (-t25 * t45 + t42 * t85);
t10 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t53, t54, 0, 0, 0, 0, 0, t9, -t8, -t54, -t9, t8, -g(1) * (-t99 + (pkin(8) * t50 - t47 * t71) * t43 - t107 * t98) - g(2) * (-(-pkin(8) * t95 + t71 * t41 + t38) * t47 + (t41 * pkin(8) + t43 * t62 + t34) * t50) - t107 * t101, 0, 0, 0, 0, 0, t3, t4, t3, t4, t9, -g(1) * (-t117 * t98 + t52 * t43 - t99) - g(2) * (t34 * t50 - t38 * t47 + (t50 * t68 + t113) * t95 + t52 * t41) - t117 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t56, 0, -t111, t56, -g(3) * (t70 * t97 + t62) + (-g(1) * t64 + g(2) * t108) * t43 + (-g(1) * t108 - g(2) * t64) * t41, 0, 0, 0, 0, 0, t5, t6, t5, t6, t111, -g(3) * (t44 * t68 + t69 * t97) + (-g(1) * t63 - g(2) * t109) * t43 + (g(1) * t109 - g(2) * t63) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, (-g(1) * ((-t41 * t94 + t43 * t87) * t50 + (-t41 * t75 - t43 * t45) * t47 - t48 * t79) - g(2) * ((t41 * t87 + t43 * t94) * t50 + (-t41 * t45 + t43 * t75) * t47 + t48 * t76) - t116) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t41 * t74 + t28 + t33) + t110;];
taug_reg = t10;
