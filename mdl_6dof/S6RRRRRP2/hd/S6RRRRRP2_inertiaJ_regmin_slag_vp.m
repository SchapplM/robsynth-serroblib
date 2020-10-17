% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:32:49
% EndTime: 2019-05-08 04:32:54
% DurationCPUTime: 1.26s
% Computational Cost: add. (1719->124), mult. (3163->208), div. (0->0), fcn. (3753->8), ass. (0->99)
t104 = cos(qJ(3));
t58 = t104 * pkin(2);
t53 = t58 + pkin(3);
t63 = sin(qJ(4));
t103 = cos(qJ(4));
t64 = sin(qJ(3));
t106 = t64 * pkin(2);
t79 = t103 * t106;
t35 = -t63 * t53 - t79;
t33 = pkin(10) - t35;
t62 = sin(qJ(5));
t60 = t62 ^ 2;
t66 = cos(qJ(5));
t61 = t66 ^ 2;
t84 = t60 + t61;
t89 = t84 * t33;
t107 = t63 * pkin(3);
t51 = pkin(10) + t107;
t120 = t84 * t51;
t67 = cos(qJ(2));
t54 = -t67 * pkin(2) - pkin(1);
t65 = sin(qJ(2));
t69 = -t104 * t67 + t64 * t65;
t30 = t69 * pkin(3) + t54;
t119 = 0.2e1 * t30;
t118 = 0.2e1 * t54;
t117 = -0.2e1 * t62;
t116 = -0.2e1 * t66;
t115 = 0.2e1 * t67;
t114 = pkin(7) + pkin(8);
t38 = t104 * t65 + t64 * t67;
t22 = t103 * t69 + t63 * t38;
t23 = t103 * t38 - t63 * t69;
t11 = t22 * pkin(4) - t23 * pkin(10) + t30;
t42 = t114 * t65;
t43 = t114 * t67;
t26 = -t104 * t43 + t64 * t42;
t16 = -t69 * pkin(9) - t26;
t25 = -t104 * t42 - t64 * t43;
t70 = -t38 * pkin(9) + t25;
t13 = t103 * t16 + t63 * t70;
t5 = t62 * t11 + t66 * t13;
t113 = pkin(4) * t62;
t112 = pkin(10) * t22;
t111 = t22 * pkin(5);
t12 = -t103 * t70 + t63 * t16;
t75 = pkin(5) * t62 - t66 * qJ(6);
t6 = t75 * t23 + t12;
t110 = t6 * t62;
t109 = t6 * t66;
t108 = t62 * pkin(10);
t105 = t66 * pkin(10);
t102 = t12 * t66;
t101 = t22 * t33;
t100 = t22 * t51;
t86 = -t103 * t53 + t63 * t106;
t32 = -pkin(4) + t86;
t99 = t32 * t66;
t57 = t103 * pkin(3);
t52 = -t57 - pkin(4);
t98 = t52 * t66;
t97 = t62 * t23;
t96 = t62 * t33;
t95 = t62 * t51;
t94 = t62 * t66;
t20 = t66 * t23;
t93 = t66 * t33;
t92 = t66 * t51;
t76 = -t66 * pkin(5) - t62 * qJ(6);
t40 = -pkin(4) + t76;
t24 = t40 + t86;
t37 = t40 - t57;
t91 = -t24 - t37;
t90 = -t24 - t40;
t88 = -t37 - t40;
t85 = pkin(10) * t84;
t83 = t22 * qJ(6);
t82 = -0.2e1 * t23 * t22;
t81 = -t66 * t11 + t62 * t13;
t78 = -pkin(4) * t23 - t112;
t2 = t83 + t5;
t3 = t81 - t111;
t1 = t2 * t66 + t3 * t62;
t77 = -t23 * t40 + t112;
t74 = -t23 * t24 + t101;
t73 = t23 * t32 - t101;
t72 = -t23 * t37 + t100;
t71 = t23 * t52 - t100;
t59 = pkin(4) * t66;
t49 = 0.2e1 * t94;
t47 = t52 * t62;
t29 = t32 * t62;
t21 = t23 ^ 2;
t19 = t66 * t22;
t18 = t62 * t22;
t17 = t62 * t20;
t14 = (-t60 + t61) * t23;
t10 = t12 * t62;
t4 = [1, 0, 0, t65 ^ 2, t65 * t115, 0, 0, 0, pkin(1) * t115, -0.2e1 * pkin(1) * t65, t38 ^ 2, -0.2e1 * t38 * t69, 0, 0, 0, t69 * t118, t38 * t118, t21, t82, 0, 0, 0, t22 * t119, t23 * t119, t61 * t21, -0.2e1 * t21 * t94, 0.2e1 * t22 * t20, t62 * t82, t22 ^ 2, 0.2e1 * t12 * t97 - 0.2e1 * t22 * t81, 0.2e1 * t12 * t20 - 0.2e1 * t5 * t22, -0.2e1 * t3 * t22 + 0.2e1 * t6 * t97, 0.2e1 * (-t2 * t62 + t3 * t66) * t23, 0.2e1 * t2 * t22 - 0.2e1 * t6 * t20, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t65, t67, 0, -t65 * pkin(7), -t67 * pkin(7), 0, 0, t38, -t69, 0, t25, t26, 0, 0, t23, -t22, 0, -t12, -t13, t17, t14, t18, t19, 0, t73 * t62 - t102, t73 * t66 + t10, -t74 * t62 - t109, t1, t74 * t66 - t110, t1 * t33 + t6 * t24; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t106, 0, 0, 0, 0, 1, -0.2e1 * t86, 0.2e1 * t35, t60, t49, 0, 0, 0, -0.2e1 * t99, 0.2e1 * t29, t24 * t116, 0.2e1 * t89, t24 * t117, t84 * t33 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t69, 0, t25, t26, 0, 0, t23, -t22, 0, -t12, -t13, t17, t14, t18, t19, 0, t71 * t62 - t102, t71 * t66 + t10, -t72 * t62 - t109, t1, t72 * t66 - t110, t1 * t51 + t6 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t58, -t106, 0, 0, 0, 0, 1, t57 - t86, -t79 + (-pkin(3) - t53) * t63, t60, t49, 0, 0, 0 (-t32 - t52) * t66, t47 + t29, t91 * t66, t120 + t89, t91 * t62, t120 * t33 + t24 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t107, t60, t49, 0, 0, 0, -0.2e1 * t98, 0.2e1 * t47, t37 * t116, 0.2e1 * t120, t37 * t117, t84 * t51 ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, -t12, -t13, t17, t14, t18, t19, 0, t78 * t62 - t102, t78 * t66 + t10, -t77 * t62 - t109, t1, t77 * t66 - t110, t1 * pkin(10) + t6 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t86, t35, t60, t49, 0, 0, 0, t59 - t99, t29 - t113, t90 * t66, t85 + t89, t90 * t62, pkin(10) * t89 + t24 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t57, -t107, t60, t49, 0, 0, 0, t59 - t98, t47 - t113, t88 * t66, t85 + t120, t88 * t62, pkin(10) * t120 + t37 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t60, t49, 0, 0, 0, 0.2e1 * t59, -0.2e1 * t113, t40 * t116, 0.2e1 * t85, t40 * t117, t84 * pkin(10) ^ 2 + t40 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t97, t22, -t81, -t5, -t81 + 0.2e1 * t111, t76 * t23, 0.2e1 * t83 + t5, -t3 * pkin(5) + t2 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t66, 0, -t96, -t93, -t96, -t75, t93, -t75 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t66, 0, -t95, -t92, -t95, -t75, t92, -t75 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t66, 0, -t108, -t105, -t108, -t75, t105, -t75 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t20, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
