% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 12:01:40
% EndTime: 2019-05-05 12:01:44
% DurationCPUTime: 1.24s
% Computational Cost: add. (1224->178), mult. (3176->342), div. (0->0), fcn. (3862->14), ass. (0->115)
t73 = sin(qJ(6));
t74 = sin(qJ(5));
t78 = cos(qJ(6));
t79 = cos(qJ(5));
t51 = t73 * t74 - t78 * t79;
t75 = sin(qJ(4));
t43 = t51 * t75;
t125 = 0.2e1 * t43;
t69 = sin(pkin(7));
t76 = sin(qJ(3));
t103 = t69 * t76;
t71 = cos(pkin(7));
t80 = cos(qJ(4));
t46 = t75 * t103 - t80 * t71;
t124 = -0.2e1 * t46;
t123 = 0.2e1 * t46;
t47 = t80 * t103 + t75 * t71;
t122 = -0.2e1 * t47;
t63 = -t79 * pkin(5) - pkin(4);
t121 = 0.2e1 * t63;
t120 = -0.2e1 * t75;
t119 = 0.2e1 * t80;
t118 = pkin(11) + pkin(12);
t117 = pkin(2) * t76;
t81 = cos(qJ(3));
t116 = pkin(2) * t81;
t115 = pkin(4) * t79;
t114 = pkin(10) * t74;
t113 = t46 * pkin(5);
t112 = t73 * pkin(5);
t111 = t78 * pkin(5);
t102 = t69 * t81;
t30 = t79 * t102 + t74 * t47;
t59 = pkin(9) * t103;
t39 = t59 + (-pkin(3) - t116) * t71;
t20 = t46 * pkin(4) - t47 * pkin(11) + t39;
t85 = pkin(9) * t102;
t40 = t85 + (pkin(10) + t117) * t71;
t41 = (-pkin(3) * t81 - pkin(10) * t76 - pkin(2)) * t69;
t24 = t80 * t40 + t75 * t41;
t22 = -pkin(11) * t102 + t24;
t9 = t74 * t20 + t79 * t22;
t7 = -t30 * pkin(12) + t9;
t110 = t78 * t7;
t109 = t80 * pkin(5);
t23 = -t75 * t40 + t80 * t41;
t21 = pkin(4) * t102 - t23;
t108 = t21 * t74;
t107 = t21 * t79;
t31 = -t74 * t102 + t79 * t47;
t106 = t31 * t74;
t105 = t46 * t80;
t64 = t69 ^ 2;
t104 = t64 * t81;
t70 = sin(pkin(6));
t77 = sin(qJ(2));
t101 = t70 * t77;
t82 = cos(qJ(2));
t100 = t70 * t82;
t99 = t71 * t76;
t98 = t71 * t82;
t97 = t74 * t46;
t96 = t74 * t75;
t95 = t74 * t79;
t94 = t74 * t80;
t93 = t75 * t46;
t55 = -t80 * pkin(4) - t75 * pkin(11) - pkin(3);
t89 = t79 * t80;
t86 = pkin(10) * t89;
t32 = t86 + (-pkin(12) * t75 + t55) * t74;
t92 = t78 * t32;
t91 = t79 * t46;
t90 = t79 * t75;
t88 = 0.2e1 * t102;
t87 = t75 * t119;
t84 = t75 * t102;
t83 = t80 * t102;
t8 = t79 * t20 - t74 * t22;
t6 = -t31 * pkin(12) + t113 + t8;
t1 = t78 * t6 - t73 * t7;
t50 = t79 * t55;
t27 = -pkin(12) * t90 + t50 + (-pkin(5) - t114) * t80;
t13 = t78 * t27 - t73 * t32;
t52 = t73 * t79 + t78 * t74;
t72 = cos(pkin(6));
t68 = t80 ^ 2;
t67 = t79 ^ 2;
t66 = t75 ^ 2;
t65 = t74 ^ 2;
t58 = t118 * t79;
t57 = t118 * t74;
t54 = (pkin(5) * t74 + pkin(10)) * t75;
t49 = pkin(2) * t99 + t85;
t48 = t71 * t116 - t59;
t45 = t46 ^ 2;
t44 = -t69 * t100 + t72 * t71;
t42 = t52 * t75;
t37 = t74 * t55 + t86;
t36 = -pkin(10) * t94 + t50;
t34 = -t73 * t57 + t78 * t58;
t33 = -t78 * t57 - t73 * t58;
t29 = t72 * t103 + (t76 * t98 + t77 * t81) * t70;
t28 = -t70 * t81 * t98 + t76 * t101 - t72 * t102;
t19 = t29 * t80 + t44 * t75;
t18 = t29 * t75 - t44 * t80;
t16 = -t73 * t30 + t78 * t31;
t15 = t78 * t30 + t73 * t31;
t14 = t73 * t27 + t92;
t12 = t30 * pkin(5) + t21;
t11 = t19 * t79 + t28 * t74;
t10 = -t19 * t74 + t28 * t79;
t4 = t73 * t10 + t78 * t11;
t3 = t78 * t10 - t73 * t11;
t2 = t73 * t6 + t110;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t100, -t101, 0, 0, 0, 0, 0, -t44 * t102 - t28 * t71, t44 * t103 - t29 * t71, 0, 0, 0, 0, 0, t18 * t102 + t28 * t46, t19 * t102 + t28 * t47, 0, 0, 0, 0, 0, t10 * t46 + t18 * t30, -t11 * t46 + t18 * t31, 0, 0, 0, 0, 0, t18 * t15 + t3 * t46, t18 * t16 - t4 * t46; 0, 1, 0, 0, t64 * t76 ^ 2, 0.2e1 * t76 * t104, 0.2e1 * t69 * t99, t71 * t88, t71 ^ 2, 0.2e1 * pkin(2) * t104 + 0.2e1 * t48 * t71, -0.2e1 * t64 * t117 - 0.2e1 * t49 * t71, t47 ^ 2, t46 * t122, t102 * t122, t46 * t88, t64 * t81 ^ 2, -0.2e1 * t23 * t102 + 0.2e1 * t39 * t46, 0.2e1 * t24 * t102 + 0.2e1 * t39 * t47, t31 ^ 2, -0.2e1 * t31 * t30, t31 * t123, t30 * t124, t45, 0.2e1 * t21 * t30 + 0.2e1 * t8 * t46, 0.2e1 * t21 * t31 - 0.2e1 * t9 * t46, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t123, t15 * t124, t45, 0.2e1 * t1 * t46 + 0.2e1 * t12 * t15, 0.2e1 * t12 * t16 - 0.2e1 * t2 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29, 0, 0, 0, 0, 0, -t28 * t80, t28 * t75, 0, 0, 0, 0, 0, -t10 * t80 + t18 * t96, t11 * t80 + t18 * t90, 0, 0, 0, 0, 0, t18 * t42 - t3 * t80, -t18 * t43 + t4 * t80; 0, 0, 0, 0, 0, 0, t103, t102, t71, t48, -t49, t47 * t75, t47 * t80 - t93, -t84, -t83, 0, -pkin(3) * t46 + pkin(10) * t84 - t39 * t80, -pkin(3) * t47 + pkin(10) * t83 + t39 * t75, t31 * t90 (-t30 * t79 - t106) * t75, -t31 * t80 + t46 * t90, t30 * t80 - t74 * t93, -t105, t36 * t46 - t8 * t80 + (pkin(10) * t30 + t108) * t75, -t37 * t46 + t9 * t80 + (pkin(10) * t31 + t107) * t75, -t16 * t43, t43 * t15 - t16 * t42, -t16 * t80 - t43 * t46, t15 * t80 - t42 * t46, -t105, -t1 * t80 + t12 * t42 + t13 * t46 + t54 * t15, -t12 * t43 - t14 * t46 + t54 * t16 + t2 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t66, t87, 0, 0, 0, pkin(3) * t119, pkin(3) * t120, t67 * t66, -0.2e1 * t66 * t95, t89 * t120, t74 * t87, t68, 0.2e1 * t114 * t66 - 0.2e1 * t36 * t80, 0.2e1 * t66 * pkin(10) * t79 + 0.2e1 * t37 * t80, t43 ^ 2, t42 * t125, t80 * t125, t42 * t119, t68, -0.2e1 * t13 * t80 + 0.2e1 * t54 * t42, 0.2e1 * t14 * t80 - 0.2e1 * t54 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, 0, 0, 0, 0, -t18 * t79, t18 * t74, 0, 0, 0, 0, 0, t18 * t51, t18 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t46, -t102, t23, -t24, t106, -t74 * t30 + t31 * t79, t97, t91, 0, -pkin(4) * t30 - pkin(11) * t97 - t107, -pkin(4) * t31 - pkin(11) * t91 + t108, t16 * t52, -t52 * t15 - t16 * t51, t52 * t46, -t51 * t46, 0, t12 * t51 + t63 * t15 + t33 * t46, t12 * t52 + t63 * t16 - t34 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t80, 0, -t75 * pkin(10), -t80 * pkin(10), t74 * t90 (-t65 + t67) * t75, -t94, -t89, 0, -pkin(10) * t90 + (-pkin(4) * t75 + pkin(11) * t80) * t74, pkin(11) * t89 + (t114 - t115) * t75, -t43 * t52, -t52 * t42 + t43 * t51, -t52 * t80, t51 * t80, 0, -t33 * t80 + t63 * t42 + t54 * t51, t34 * t80 - t63 * t43 + t54 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t65, 0.2e1 * t95, 0, 0, 0, 0.2e1 * t115, -0.2e1 * pkin(4) * t74, t52 ^ 2, -0.2e1 * t52 * t51, 0, 0, 0, t51 * t121, t52 * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, t46, t8, -t9, 0, 0, t16, -t15, t46, t111 * t46 + t1, -t110 + (-t6 - t113) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t96, -t80, t36, -t37, 0, 0, -t43, -t42, -t80, -t109 * t78 + t13, -t92 + (-t27 + t109) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t79, 0, -t74 * pkin(11), -t79 * pkin(11), 0, 0, t52, -t51, 0, t33, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t111, -0.2e1 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t46, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t42, -t80, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t51, 0, t33, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t111, -t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
