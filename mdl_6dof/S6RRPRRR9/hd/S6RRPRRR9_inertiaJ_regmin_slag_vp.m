% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:07:03
% EndTime: 2019-05-06 23:07:08
% DurationCPUTime: 1.33s
% Computational Cost: add. (2351->153), mult. (5421->307), div. (0->0), fcn. (6609->12), ass. (0->110)
t115 = cos(qJ(5));
t80 = sin(pkin(6));
t86 = sin(qJ(2));
t112 = t80 * t86;
t79 = sin(pkin(12));
t81 = cos(pkin(12));
t82 = cos(pkin(6));
t57 = t79 * t112 - t82 * t81;
t58 = t81 * t112 + t79 * t82;
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t38 = t57 * t88 + t85 * t58;
t39 = -t57 * t85 + t58 * t88;
t84 = sin(qJ(5));
t26 = t115 * t38 + t39 * t84;
t123 = -0.2e1 * t26;
t62 = t85 * t79 - t81 * t88;
t71 = -pkin(3) * t81 - pkin(2);
t52 = pkin(4) * t62 + t71;
t122 = 0.2e1 * t52;
t121 = 0.2e1 * t71;
t120 = pkin(1) * t86;
t89 = cos(qJ(2));
t119 = pkin(1) * t89;
t111 = t80 * t89;
t98 = pkin(8) * t111;
t53 = t98 + (qJ(3) + t120) * t82;
t54 = (-pkin(2) * t89 - qJ(3) * t86 - pkin(1)) * t80;
t33 = -t53 * t79 + t81 * t54;
t31 = -pkin(3) * t111 - pkin(9) * t58 + t33;
t34 = t81 * t53 + t79 * t54;
t32 = -pkin(9) * t57 + t34;
t15 = t88 * t31 - t32 * t85;
t99 = pkin(4) * t111;
t11 = -pkin(10) * t39 + t15 - t99;
t16 = t31 * t85 + t32 * t88;
t14 = -pkin(10) * t38 + t16;
t6 = t115 * t11 - t84 * t14;
t4 = pkin(5) * t111 - t6;
t87 = cos(qJ(6));
t118 = t4 * t87;
t117 = t84 * pkin(4);
t96 = t115 * pkin(4);
t73 = -t96 - pkin(5);
t116 = pkin(5) - t73;
t27 = t115 * t39 - t84 * t38;
t83 = sin(qJ(6));
t19 = -t83 * t111 + t27 * t87;
t17 = t19 * t83;
t104 = pkin(9) + qJ(3);
t64 = t104 * t79;
t65 = t104 * t81;
t48 = -t64 * t85 + t65 * t88;
t37 = -pkin(10) * t62 + t48;
t47 = -t64 * t88 - t65 * t85;
t63 = t79 * t88 + t81 * t85;
t91 = -pkin(10) * t63 + t47;
t21 = -t115 * t91 + t37 * t84;
t114 = t21 * t87;
t75 = t80 ^ 2;
t113 = t75 * t89;
t110 = t82 * t86;
t24 = t83 * t26;
t45 = t115 * t62 + t63 * t84;
t42 = t83 * t45;
t46 = t115 * t63 - t84 * t62;
t109 = t83 * t46;
t72 = pkin(11) + t117;
t108 = t83 * t72;
t107 = t83 * t87;
t25 = t87 * t26;
t106 = t87 * t46;
t105 = t87 * t72;
t103 = t79 ^ 2 + t81 ^ 2;
t102 = -0.2e1 * t46 * t45;
t101 = -0.2e1 * t111;
t100 = 0.2e1 * t111;
t97 = qJ(3) * t111;
t95 = t115 * t14;
t94 = -pkin(5) * t46 - pkin(11) * t45;
t93 = -t33 * t79 + t34 * t81;
t92 = -t45 * t72 + t46 * t73;
t66 = pkin(8) * t112;
t56 = t66 + (-pkin(2) - t119) * t82;
t7 = t84 * t11 + t95;
t41 = pkin(3) * t57 + t56;
t28 = pkin(4) * t38 + t41;
t78 = t87 ^ 2;
t77 = t83 ^ 2;
t70 = t75 * t89 ^ 2;
t67 = 0.2e1 * t107;
t61 = pkin(1) * t110 + t98;
t60 = t82 * t119 - t66;
t44 = t46 ^ 2;
t43 = t87 * t45;
t40 = t83 * t106;
t30 = (-t77 + t78) * t46;
t23 = pkin(5) * t45 - pkin(11) * t46 + t52;
t22 = t115 * t37 + t84 * t91;
t20 = t21 * t83;
t18 = t87 * t111 + t27 * t83;
t13 = t22 * t87 + t23 * t83;
t12 = -t22 * t83 + t23 * t87;
t9 = -t18 * t83 + t19 * t87;
t8 = pkin(5) * t26 - pkin(11) * t27 + t28;
t5 = -pkin(11) * t111 + t7;
t3 = t4 * t83;
t2 = t5 * t87 + t8 * t83;
t1 = -t5 * t83 + t8 * t87;
t10 = [1, 0, 0, t75 * t86 ^ 2, 0.2e1 * t86 * t113, 0.2e1 * t80 * t110, t82 * t100, t82 ^ 2, 0.2e1 * pkin(1) * t113 + 0.2e1 * t60 * t82, -0.2e1 * t75 * t120 - 0.2e1 * t61 * t82, -0.2e1 * t33 * t111 + 0.2e1 * t56 * t57, 0.2e1 * t34 * t111 + 0.2e1 * t56 * t58, -0.2e1 * t33 * t58 - 0.2e1 * t34 * t57, t33 ^ 2 + t34 ^ 2 + t56 ^ 2, t39 ^ 2, -0.2e1 * t39 * t38, t39 * t101, t38 * t100, t70, -0.2e1 * t15 * t111 + 0.2e1 * t38 * t41, 0.2e1 * t16 * t111 + 0.2e1 * t39 * t41, t27 ^ 2, t27 * t123, t27 * t101, t26 * t100, t70, -0.2e1 * t111 * t6 + 0.2e1 * t26 * t28, 0.2e1 * t111 * t7 + 0.2e1 * t27 * t28, t19 ^ 2, -0.2e1 * t19 * t18, 0.2e1 * t19 * t26, t18 * t123, t26 ^ 2, 0.2e1 * t1 * t26 + 0.2e1 * t18 * t4, 0.2e1 * t19 * t4 - 0.2e1 * t2 * t26; 0, 0, 0, 0, 0, t112, t111, t82, t60, -t61, -pkin(2) * t57 - t56 * t81 + t79 * t97, -pkin(2) * t58 + t56 * t79 + t81 * t97 (-t57 * t81 + t58 * t79) * qJ(3) + t93, -pkin(2) * t56 + t93 * qJ(3), t39 * t63, -t38 * t63 - t39 * t62, -t63 * t111, t62 * t111, 0, -t47 * t111 + t38 * t71 + t41 * t62, t48 * t111 + t39 * t71 + t41 * t63, t27 * t46, -t26 * t46 - t27 * t45, -t46 * t111, t45 * t111, 0, t21 * t111 + t26 * t52 + t28 * t45, t111 * t22 + t27 * t52 + t28 * t46, t19 * t106 (-t18 * t87 - t17) * t46, t106 * t26 + t19 * t45, -t109 * t26 - t18 * t45, t26 * t45, t1 * t45 + t109 * t4 + t12 * t26 + t18 * t21, t106 * t4 - t13 * t26 + t19 * t21 - t2 * t45; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t81, -0.2e1 * pkin(2) * t79, 0.2e1 * t103 * qJ(3), t103 * qJ(3) ^ 2 + pkin(2) ^ 2, t63 ^ 2, -0.2e1 * t63 * t62, 0, 0, 0, t62 * t121, t63 * t121, t44, t102, 0, 0, 0, t45 * t122, t46 * t122, t78 * t44, -0.2e1 * t44 * t107, 0.2e1 * t45 * t106, t83 * t102, t45 ^ 2, 0.2e1 * t21 * t109 + 0.2e1 * t12 * t45, 0.2e1 * t21 * t106 - 0.2e1 * t13 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t58, 0, t56, 0, 0, 0, 0, 0, t38, t39, 0, 0, 0, 0, 0, t26, t27, 0, 0, 0, 0, 0, t25, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t79, 0, -pkin(2), 0, 0, 0, 0, 0, t62, t63, 0, 0, 0, 0, 0, t45, t46, 0, 0, 0, 0, 0, t43, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t38, -t111, t15, -t16, 0, 0, t27, -t26, -t111, -t111 * t96 + t6, -t95 + (-t11 + t99) * t84, t17, t9, t24, t25, 0, -t108 * t26 + t18 * t73 - t118, -t105 * t26 + t19 * t73 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t62, 0, t47, -t48, 0, 0, t46, -t45, 0, -t21, -t22, t40, t30, t42, t43, 0, t83 * t92 - t114, t87 * t92 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t96, -0.2e1 * t117, t77, t67, 0, 0, 0, -0.2e1 * t73 * t87, 0.2e1 * t73 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, -t111, t6, -t7, t17, t9, t24, t25, 0, -pkin(5) * t18 - pkin(11) * t24 - t118, -pkin(5) * t19 - pkin(11) * t25 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, 0, -t21, -t22, t40, t30, t42, t43, 0, t83 * t94 - t114, t87 * t94 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t96, -t117, t77, t67, 0, 0, 0, t116 * t87, -t116 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t77, t67, 0, 0, 0, 0.2e1 * pkin(5) * t87, -0.2e1 * pkin(5) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t26, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, -t109, t45, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t87, 0, -t108, -t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t87, 0, -t83 * pkin(11), -t87 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
