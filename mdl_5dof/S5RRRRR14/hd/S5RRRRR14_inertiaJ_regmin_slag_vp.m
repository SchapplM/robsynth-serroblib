% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR14_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:37
% EndTime: 2024-09-27 18:43:38
% DurationCPUTime: 0.62s
% Computational Cost: add. (839->101), mult. (2112->142), div. (0->0), fcn. (2284->10), ass. (0->99)
t86 = cos(pkin(5));
t119 = -0.2e1 * t86;
t118 = 0.2e1 * t86;
t93 = cos(qJ(3));
t117 = pkin(3) * t93;
t85 = sin(pkin(5));
t89 = sin(qJ(3));
t74 = t85 * t89;
t75 = t85 * t93;
t88 = sin(qJ(4));
t92 = cos(qJ(4));
t52 = t88 * t74 - t92 * t75;
t116 = t52 * pkin(4);
t79 = t86 * pkin(3);
t78 = t86 * pkin(4);
t87 = sin(qJ(5));
t115 = t87 * pkin(4);
t114 = t88 * pkin(3);
t113 = sin(qJ(2)) * pkin(1);
t91 = cos(qJ(5));
t80 = t91 * pkin(4);
t106 = t86 * t89;
t65 = t85 * pkin(8) + t113;
t82 = cos(qJ(2)) * pkin(1);
t77 = t82 + pkin(2);
t43 = t77 * t106 + t93 * t65;
t69 = pkin(9) * t75;
t35 = t43 + t69;
t103 = t92 * t35;
t105 = t86 * t93;
t61 = t77 * t105;
t32 = t61 + t79 + (-pkin(9) * t85 - t65) * t89;
t19 = t88 * t32 + t103;
t49 = t52 * pkin(10);
t9 = t19 - t49;
t112 = t91 * t9;
t81 = t92 * pkin(3);
t18 = t92 * t32 - t88 * t35;
t53 = (t88 * t93 + t89 * t92) * t85;
t95 = -t53 * pkin(10) + t78;
t8 = t18 + t95;
t2 = t91 * t8 - t87 * t9;
t27 = t91 * t52 + t87 * t53;
t54 = (-t77 - t117) * t85;
t31 = t54 + t116;
t111 = t2 * t86 + t31 * t27;
t59 = (-pkin(2) - t117) * t85;
t36 = t59 + t116;
t72 = pkin(2) * t105;
t41 = t72 + t79 + (-pkin(8) - pkin(9)) * t74;
t58 = pkin(2) * t106 + pkin(8) * t75;
t46 = t58 + t69;
t22 = t92 * t41 - t88 * t46;
t14 = t22 + t95;
t102 = t92 * t46;
t23 = t88 * t41 + t102;
t20 = t23 - t49;
t5 = t91 * t14 - t87 * t20;
t110 = t36 * t27 + t5 * t86;
t76 = t81 + pkin(4);
t96 = t91 * t114;
t56 = t87 * t76 + t96;
t109 = t56 * t86;
t83 = t85 ^ 2;
t108 = t83 * t89;
t107 = t83 * t93;
t104 = t91 * t20;
t101 = t18 * t86 + t54 * t52;
t100 = t22 * t86 + t59 * t52;
t42 = -t89 * t65 + t61;
t99 = t77 * t107 + t42 * t86;
t57 = -pkin(8) * t74 + t72;
t98 = pkin(2) * t107 + t57 * t86;
t97 = t85 * t118;
t55 = -t87 * t114 + t91 * t76;
t3 = t87 * t8 + t112;
t6 = t87 * t14 + t104;
t84 = t86 ^ 2;
t73 = t83 * t89 ^ 2;
t71 = t86 * t81;
t70 = t86 * t80;
t64 = 0.2e1 * t89 * t107;
t63 = t93 * t97;
t62 = t89 * t97;
t50 = t53 ^ 2;
t48 = t55 * t86;
t45 = t53 * t118;
t44 = t52 * t119;
t38 = t59 * t53;
t34 = t54 * t53;
t30 = -0.2e1 * t53 * t52;
t28 = -t87 * t52 + t91 * t53;
t26 = t28 ^ 2;
t25 = t28 * t118;
t24 = t27 * t119;
t17 = t36 * t28;
t13 = t31 * t28;
t10 = -0.2e1 * t28 * t27;
t1 = [1, 0, 0, 1, 0.2e1 * t82, -0.2e1 * t113, t73, t64, t62, t63, t84, 0.2e1 * t99, -0.2e1 * t77 * t108 - 0.2e1 * t43 * t86, t50, t30, t45, t44, t84, 0.2e1 * t101, -0.2e1 * t19 * t86 + 0.2e1 * t34, t26, t10, t25, t24, t84, 0.2e1 * t111, -0.2e1 * t3 * t86 + 0.2e1 * t13; 0, 0, 0, 1, t82, -t113, t73, t64, t62, t63, t84, t98 + t99, (-t43 - t58) * t86 + (-pkin(2) - t77) * t108, t50, t30, t45, t44, t84, t100 + t101, t34 + t38 + (-t19 - t23) * t86, t26, t10, t25, t24, t84, t110 + t111, t13 + t17 + (-t3 - t6) * t86; 0, 0, 0, 1, 0, 0, t73, t64, t62, t63, t84, 0.2e1 * t98, -0.2e1 * pkin(2) * t108 - 0.2e1 * t58 * t86, t50, t30, t45, t44, t84, 0.2e1 * t100, -0.2e1 * t23 * t86 + 0.2e1 * t38, t26, t10, t25, t24, t84, 0.2e1 * t110, -0.2e1 * t6 * t86 + 0.2e1 * t17; 0, 0, 0, 0, 0, 0, 0, 0, t74, t75, t86, t42, -t43, 0, 0, t53, -t52, t86, t18 + t71, -t103 + (-t32 - t79) * t88, 0, 0, t28, -t27, t86, t2 + t48, -t3 - t109; 0, 0, 0, 0, 0, 0, 0, 0, t74, t75, t86, t57, -t58, 0, 0, t53, -t52, t86, t22 + t71, -t102 + (-t41 - t79) * t88, 0, 0, t28, -t27, t86, t48 + t5, -t6 - t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t81, -0.2e1 * t114, 0, 0, 0, 0, 1, 0.2e1 * t55, -0.2e1 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t52, t86, t18, -t19, 0, 0, t28, -t27, t86, t2 + t70, -t112 + (-t8 - t78) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t52, t86, t22, -t23, 0, 0, t28, -t27, t86, t5 + t70, -t104 + (-t14 - t78) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t81, -t114, 0, 0, 0, 0, 1, t55 + t80, -t96 + (-pkin(4) - t76) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t80, -0.2e1 * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, t86, t2, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, t86, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t55, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t80, -t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t1;
