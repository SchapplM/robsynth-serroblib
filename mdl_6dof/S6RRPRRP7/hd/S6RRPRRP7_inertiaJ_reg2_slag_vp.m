% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t77 = sin(qJ(4));
t76 = sin(qJ(5));
t70 = t76 ^ 2;
t79 = cos(qJ(5));
t73 = t79 ^ 2;
t98 = t70 + t73;
t139 = t98 * t77;
t128 = -pkin(2) - pkin(3);
t80 = cos(qJ(4));
t49 = t80 * qJ(3) + t77 * t128;
t44 = -pkin(9) + t49;
t138 = t98 * t44;
t78 = sin(qJ(2));
t72 = t78 ^ 2;
t81 = cos(qJ(2));
t75 = t81 ^ 2;
t137 = t72 + t75;
t37 = t78 * t77 + t81 * t80;
t39 = -t81 * t77 + t78 * t80;
t85 = t77 * t37 + t80 * t39;
t135 = t85 * t79;
t67 = t81 * pkin(7);
t53 = -t81 * pkin(8) + t67;
t94 = (pkin(7) - pkin(8)) * t78;
t16 = t77 * t53 - t80 * t94;
t134 = t16 ^ 2;
t35 = t37 ^ 2;
t133 = 0.2e1 * t37;
t132 = 0.2e1 * t39;
t131 = -0.2e1 * t76;
t130 = -0.2e1 * t78;
t129 = 0.2e1 * t79;
t127 = pkin(9) * t37;
t125 = t37 * pkin(5);
t88 = pkin(5) * t76 - t79 * qJ(6);
t7 = t88 * t39 + t16;
t124 = t7 * t76;
t123 = t7 * t79;
t122 = t76 * pkin(9);
t121 = t78 * pkin(7);
t120 = t79 * pkin(9);
t47 = t77 * qJ(3) - t80 * t128;
t43 = pkin(4) + t47;
t119 = pkin(4) + t43;
t51 = -t81 * pkin(2) - t78 * qJ(3) - pkin(1);
t34 = t81 * pkin(3) - t51;
t10 = t37 * pkin(4) - t39 * pkin(9) + t34;
t18 = t80 * t53 + t77 * t94;
t6 = t76 * t10 + t79 * t18;
t118 = t16 * t76;
t117 = t16 * t79;
t116 = t16 * t80;
t115 = t37 * t44;
t112 = t76 * t39;
t111 = t76 * t44;
t110 = t76 * t77;
t109 = t76 * t79;
t108 = t78 * t81;
t31 = t79 * t39;
t107 = t79 * t44;
t106 = t79 * t77;
t89 = -t79 * pkin(5) - t76 * qJ(6);
t50 = -pkin(4) + t89;
t19 = t47 - t50;
t105 = -t19 + t50;
t104 = t139 * t44;
t103 = pkin(9) * t138;
t102 = t98 * t44 ^ 2;
t101 = t139 * pkin(9);
t100 = t98 * pkin(9) ^ 2;
t99 = t137 * pkin(7) ^ 2;
t97 = t37 * qJ(6);
t96 = t37 * t112;
t36 = t39 ^ 2;
t95 = t36 * t109;
t93 = -t79 * t10 + t76 * t18;
t92 = -pkin(4) * t39 - t127;
t3 = t97 + t6;
t4 = t93 - t125;
t1 = t3 * t79 + t4 * t76;
t2 = t6 * t79 + t76 * t93;
t91 = -t39 * t50 + t127;
t90 = -t78 * pkin(2) + t81 * qJ(3);
t87 = t19 * t39 - t115;
t86 = t39 * t43 - t115;
t84 = t85 * t76;
t74 = t80 ^ 2;
t71 = t77 ^ 2;
t59 = t80 * t79;
t58 = t80 * t76;
t57 = -0.2e1 * t109;
t56 = 0.2e1 * t109;
t46 = 0.2e1 * t137 * pkin(7);
t45 = 0.2e1 * t98 * pkin(9);
t30 = t79 * t37;
t28 = t73 * t36;
t27 = t76 * t37;
t26 = t70 * t36;
t25 = t98 * t71 + t74;
t20 = t76 * t31;
t14 = t31 * t133;
t13 = (t70 - t73) * t39;
t12 = -0.2e1 * t138;
t11 = (-pkin(9) + t44) * t98;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t72, 0.2e1 * t108, 0, t75, 0, 0, 0.2e1 * pkin(1) * t81, pkin(1) * t130, t46, pkin(1) ^ 2 + t99, t72, 0, -0.2e1 * t108, 0, 0, t75, -0.2e1 * t51 * t81, t46, t51 * t130, t51 ^ 2 + t99, t36, -0.2e1 * t39 * t37, 0, t35, 0, 0, t34 * t133, t34 * t132, 0.2e1 * t16 * t39 - 0.2e1 * t18 * t37, t18 ^ 2 + t34 ^ 2 + t134, t28, -0.2e1 * t95, t14, t26, -0.2e1 * t96, t35, 0.2e1 * t16 * t112 - 0.2e1 * t37 * t93, 0.2e1 * t16 * t31 - 0.2e1 * t6 * t37 (-t6 * t76 + t79 * t93) * t132, t6 ^ 2 + t93 ^ 2 + t134, t28, t14, 0.2e1 * t95, t35, 0.2e1 * t96, t26, 0.2e1 * t7 * t112 - 0.2e1 * t4 * t37 (-t3 * t76 + t4 * t79) * t132, 0.2e1 * t3 * t37 - 0.2e1 * t7 * t31, t3 ^ 2 + t4 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, t81, 0, -t121, -t67, 0, 0, 0, t78, 0, 0, -t81, 0, -t121, t90, t67, t90 * pkin(7), 0, 0, -t39, 0, t37, 0, t16, t18, -t49 * t37 + t47 * t39, t16 * t47 + t18 * t49, -t20, t13, -t27, t20, -t30, 0, t86 * t76 + t117, t86 * t79 - t118, -t2, t16 * t43 + t2 * t44, -t20, -t27, -t13, 0, t30, t20, t87 * t76 + t123, -t1, -t87 * t79 + t124, t1 * t44 + t7 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t47, 0.2e1 * t49, 0, t47 ^ 2 + t49 ^ 2, t70, t56, 0, t73, 0, 0, t43 * t129, t43 * t131, t12, t43 ^ 2 + t102, t70, 0, t57, 0, 0, t73, t19 * t129, t12, 0.2e1 * t19 * t76, t19 ^ 2 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, t121, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t18 * t77 - t116, 0, 0, 0, 0, 0, 0, -t84, -t135, 0, t2 * t77 - t116, 0, 0, 0, 0, 0, 0, -t84, 0, t135, t1 * t77 - t7 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t80, t77, 0, -t47 * t80 + t49 * t77, 0, 0, 0, 0, 0, 0, -t59, t58, -t139, -t43 * t80 + t104, 0, 0, 0, 0, 0, 0, -t59, -t139, -t58, -t19 * t80 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 + t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t37, 0, -t16, -t18, 0, 0, t20, -t13, t27, -t20, t30, 0, t92 * t76 - t117, t92 * t79 + t118, t2, -t16 * pkin(4) + t2 * pkin(9), t20, t27, t13, 0, -t30, -t20, -t91 * t76 - t123, t1, t91 * t79 - t124, t1 * pkin(9) + t7 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t47, -t49, 0, 0, -t70, t57, 0, -t73, 0, 0, -t119 * t79, t119 * t76, t11, -t43 * pkin(4) + t103, -t70, 0, t56, 0, 0, -t73, t105 * t79, t11, t105 * t76, t19 * t50 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t77, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t58, t139, t80 * pkin(4) + t101, 0, 0, 0, 0, 0, 0, t59, t139, t58, -t80 * t50 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t70, t56, 0, t73, 0, 0, pkin(4) * t129, pkin(4) * t131, t45, pkin(4) ^ 2 + t100, t70, 0, t57, 0, 0, t73, -0.2e1 * t50 * t79, t45, t50 * t131, t50 ^ 2 + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t112, t37, -t93, -t6, 0, 0, 0, t31, 0, t37, t112, 0, -t93 + 0.2e1 * t125, t89 * t39, 0.2e1 * t97 + t6, -t4 * pkin(5) + t3 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, 0, -t79, 0, -t111, -t107, 0, 0, 0, -t76, 0, 0, t79, 0, -t111, t88, t107, -t88 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t106, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, t106, -t88 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, t79, 0, -t122, -t120, 0, 0, 0, t76, 0, 0, -t79, 0, -t122, -t88, t120, -t88 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t31, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, 0, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;