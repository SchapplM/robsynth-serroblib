% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 14:32:31
% EndTime: 2019-05-07 14:32:35
% DurationCPUTime: 1.45s
% Computational Cost: add. (1272->170), mult. (2929->313), div. (0->0), fcn. (3383->10), ass. (0->113)
t76 = sin(pkin(6));
t81 = sin(qJ(2));
t120 = t76 * t81;
t77 = cos(pkin(6));
t80 = sin(qJ(3));
t84 = cos(qJ(3));
t40 = t80 * t120 - t77 * t84;
t41 = t84 * t120 + t77 * t80;
t79 = sin(qJ(5));
t83 = cos(qJ(5));
t25 = -t40 * t83 + t41 * t79;
t135 = -0.2e1 * t25;
t46 = t79 * t80 + t83 * t84;
t134 = 0.2e1 * t46;
t78 = sin(qJ(6));
t133 = -0.2e1 * t78;
t132 = -0.2e1 * t80;
t82 = cos(qJ(6));
t131 = 0.2e1 * t82;
t130 = 0.2e1 * t84;
t129 = pkin(1) * t81;
t85 = cos(qJ(2));
t62 = t76 * t85;
t99 = pkin(8) * t62;
t35 = t99 + (pkin(9) + t129) * t77;
t36 = (-pkin(2) * t85 - pkin(9) * t81 - pkin(1)) * t76;
t105 = t80 * t35 - t84 * t36;
t58 = pkin(3) * t62;
t19 = t58 + t105;
t10 = pkin(4) * t62 - t41 * pkin(10) + t19;
t23 = t84 * t35 + t80 * t36;
t95 = qJ(4) * t62;
t18 = -t95 + t23;
t16 = t40 * pkin(10) + t18;
t5 = t83 * t10 - t79 * t16;
t3 = -pkin(5) * t62 - t5;
t128 = t3 * t78;
t127 = t3 * t82;
t86 = -pkin(3) - pkin(4);
t51 = t79 * qJ(4) - t83 * t86;
t49 = pkin(5) + t51;
t126 = pkin(5) + t49;
t6 = t79 * t10 + t83 * t16;
t26 = t40 * t79 + t41 * t83;
t21 = t26 * t82 + t78 * t62;
t125 = t21 * t78;
t68 = t84 * pkin(9);
t55 = -t84 * pkin(10) + t68;
t67 = t80 * pkin(9);
t94 = -t80 * pkin(10) + t67;
t28 = t79 * t55 - t83 * t94;
t124 = t28 * t78;
t123 = t28 * t82;
t122 = t41 * t80;
t71 = t76 ^ 2;
t121 = t71 * t85;
t119 = t77 * t81;
t118 = t78 * t25;
t117 = t78 * t46;
t47 = -t79 * t84 + t83 * t80;
t116 = t78 * t47;
t52 = t83 * qJ(4) + t79 * t86;
t50 = -pkin(11) + t52;
t115 = t78 * t50;
t114 = t78 * t79;
t113 = t78 * t82;
t112 = t82 * t25;
t111 = t82 * t46;
t110 = t82 * t47;
t109 = t82 * t50;
t108 = t82 * t79;
t107 = t83 * t78;
t106 = t83 * t82;
t42 = t77 * t85 * pkin(1) - pkin(8) * t120;
t73 = t80 ^ 2;
t104 = t84 ^ 2 + t73;
t103 = -0.2e1 * t47 * t46;
t102 = -0.2e1 * t62;
t101 = 0.2e1 * t62;
t100 = -0.2e1 * t113;
t53 = -t84 * pkin(3) - t80 * qJ(4) - pkin(2);
t98 = t80 * t62;
t97 = t84 * t62;
t96 = t78 * t110;
t34 = -t77 * pkin(2) - t42;
t93 = pkin(9) * t97;
t44 = t84 * pkin(4) - t53;
t17 = t40 * pkin(3) - t41 * qJ(4) + t34;
t92 = -pkin(5) * t47 - pkin(11) * t46;
t91 = -pkin(3) * t80 + t84 * qJ(4);
t90 = t18 * t84 + t19 * t80;
t89 = -t46 * t50 + t47 * t49;
t88 = -t46 * t79 - t47 * t83;
t15 = -t40 * pkin(4) - t17;
t74 = t82 ^ 2;
t72 = t78 ^ 2;
t61 = t71 * t85 ^ 2;
t60 = 0.2e1 * t113;
t54 = pkin(9) * t98;
t45 = t47 ^ 2;
t43 = pkin(1) * t119 + t99;
t29 = t83 * t55 + t79 * t94;
t27 = (t72 - t74) * t47;
t24 = t46 * pkin(5) - t47 * pkin(11) + t44;
t20 = t26 * t78 - t82 * t62;
t12 = t78 * t24 + t82 * t29;
t11 = t82 * t24 - t78 * t29;
t8 = -t78 * t20 + t21 * t82;
t7 = t25 * pkin(5) - t26 * pkin(11) + t15;
t4 = pkin(11) * t62 + t6;
t2 = t82 * t4 + t78 * t7;
t1 = -t78 * t4 + t82 * t7;
t9 = [1, 0, 0, t71 * t81 ^ 2, 0.2e1 * t81 * t121, 0.2e1 * t76 * t119, t77 * t101, t77 ^ 2, 0.2e1 * pkin(1) * t121 + 0.2e1 * t42 * t77, -0.2e1 * t71 * t129 - 0.2e1 * t43 * t77, t41 ^ 2, -0.2e1 * t41 * t40, t41 * t102, t40 * t101, t61, 0.2e1 * t105 * t62 + 0.2e1 * t34 * t40, 0.2e1 * t23 * t62 + 0.2e1 * t34 * t41, 0.2e1 * t17 * t40 + 0.2e1 * t19 * t62, -0.2e1 * t18 * t40 + 0.2e1 * t19 * t41, -0.2e1 * t17 * t41 - 0.2e1 * t18 * t62, t17 ^ 2 + t18 ^ 2 + t19 ^ 2, t26 ^ 2, t26 * t135, t26 * t101, t25 * t102, t61, 0.2e1 * t15 * t25 + 0.2e1 * t5 * t62, 0.2e1 * t15 * t26 - 0.2e1 * t6 * t62, t21 ^ 2, -0.2e1 * t21 * t20, 0.2e1 * t21 * t25, t20 * t135, t25 ^ 2, 0.2e1 * t1 * t25 + 0.2e1 * t3 * t20, -0.2e1 * t2 * t25 + 0.2e1 * t3 * t21; 0, 0, 0, 0, 0, t120, t62, t77, t42, -t43, t122, -t80 * t40 + t41 * t84, -t98, -t97, 0, -pkin(2) * t40 - t34 * t84 + t54, -pkin(2) * t41 + t34 * t80 + t93, -t17 * t84 + t53 * t40 + t54 (-t40 * t84 + t122) * pkin(9) + t90, -t17 * t80 - t53 * t41 - t93, pkin(9) * t90 + t17 * t53, t26 * t47, -t47 * t25 - t26 * t46, t47 * t62, -t46 * t62, 0, t15 * t46 + t44 * t25 - t28 * t62, t15 * t47 + t44 * t26 - t29 * t62, t21 * t110 (-t20 * t82 - t125) * t47, t110 * t25 + t21 * t46, -t25 * t116 - t20 * t46, t25 * t46, t1 * t46 + t11 * t25 + t3 * t116 + t28 * t20, t110 * t3 - t12 * t25 - t2 * t46 + t28 * t21; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t73, t80 * t130, 0, 0, 0, pkin(2) * t130, pkin(2) * t132, -0.2e1 * t53 * t84, 0.2e1 * t104 * pkin(9), t53 * t132, pkin(9) ^ 2 * t104 + t53 ^ 2, t45, t103, 0, 0, 0, t44 * t134, 0.2e1 * t44 * t47, t74 * t45, t45 * t100, t110 * t134, t78 * t103, t46 ^ 2, 0.2e1 * t11 * t46 + 0.2e1 * t28 * t116, 0.2e1 * t110 * t28 - 0.2e1 * t12 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, -t62, -t105, -t23, -0.2e1 * t58 - t105, -t41 * pkin(3) - t40 * qJ(4), -0.2e1 * t95 + t23, -t19 * pkin(3) + t18 * qJ(4), 0, 0, -t26, t25, -t62, -t51 * t62 - t5, -t52 * t62 + t6, -t125, -t8, -t118, -t112, 0, -t25 * t115 + t49 * t20 + t127, -t109 * t25 + t49 * t21 - t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t84, 0, -t67, -t68, -t67, t91, t68, t91 * pkin(9), 0, 0, -t47, t46, 0, t28, t29, -t96, t27, -t117, -t111, 0, t89 * t78 + t123, t82 * t89 - t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t51, 0.2e1 * t52, t72, t60, 0, 0, 0, t49 * t131, t49 * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t41, 0, t19, 0, 0, 0, 0, 0, t83 * t62, -t79 * t62, 0, 0, 0, 0, 0, -t114 * t25 - t83 * t20, -t108 * t25 - t83 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t78, t88 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, -t83, t79, 0, 0, 0, 0, 0, -t106, t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t62, t5, -t6, t125, t8, t118, t112, 0, -pkin(5) * t20 - pkin(11) * t118 - t127, -pkin(5) * t21 - pkin(11) * t112 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t46, 0, -t28, -t29, t96, -t27, t117, t111, 0, t92 * t78 - t123, t82 * t92 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t51, -t52, -t72, t100, 0, 0, 0, -t126 * t82, t126 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t79, 0, 0, 0, 0, 0, t106, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t72, t60, 0, 0, 0, pkin(5) * t131, pkin(5) * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, t25, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, -t116, t46, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t82, 0, -t115, -t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t82, 0, -t78 * pkin(11), -t82 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
