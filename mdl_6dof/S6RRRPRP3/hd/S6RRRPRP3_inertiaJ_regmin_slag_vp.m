% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:39:45
% EndTime: 2019-05-07 07:39:49
% DurationCPUTime: 1.11s
% Computational Cost: add. (1608->149), mult. (3019->237), div. (0->0), fcn. (3519->8), ass. (0->92)
t75 = sin(pkin(10));
t76 = cos(pkin(10));
t91 = t75 ^ 2 + t76 ^ 2;
t92 = t91 * qJ(4);
t106 = cos(qJ(5));
t77 = sin(qJ(5));
t98 = t77 * t75;
t121 = t106 * t76 - t98;
t120 = -0.2e1 * t121;
t86 = t106 * t75;
t53 = t77 * t76 + t86;
t119 = -0.2e1 * t53;
t107 = cos(qJ(3));
t108 = cos(qJ(2));
t78 = sin(qJ(3));
t79 = sin(qJ(2));
t54 = -t107 * t108 + t78 * t79;
t118 = -0.2e1 * t54;
t112 = t76 * pkin(4);
t87 = t107 * pkin(2);
t68 = -t87 - pkin(3);
t57 = t68 - t112;
t117 = 0.2e1 * t57;
t66 = -pkin(3) - t112;
t116 = 0.2e1 * t66;
t69 = -t108 * pkin(2) - pkin(1);
t115 = 0.2e1 * t69;
t55 = t107 * t79 + t78 * t108;
t35 = t54 * pkin(3) - t55 * qJ(4) + t69;
t61 = (-pkin(8) - pkin(7)) * t79;
t88 = t108 * pkin(7);
t62 = t108 * pkin(8) + t88;
t45 = t107 * t62 + t78 * t61;
t17 = t76 * t35 - t75 * t45;
t99 = t76 * t55;
t14 = t54 * pkin(4) - pkin(9) * t99 + t17;
t100 = t75 * t55;
t18 = t75 * t35 + t76 * t45;
t16 = -pkin(9) * t100 + t18;
t6 = t106 * t16 + t77 * t14;
t90 = t54 * qJ(6);
t3 = t90 + t6;
t113 = t54 * pkin(5);
t83 = -t106 * t14 + t77 * t16;
t4 = t83 - t113;
t114 = t121 * t3 + t4 * t53;
t111 = t78 * pkin(2);
t110 = pkin(3) - t68;
t65 = qJ(4) + t111;
t109 = -pkin(9) - t65;
t72 = t76 * pkin(9);
t49 = t76 * t65 + t72;
t32 = -t109 * t86 + t77 * t49;
t105 = t32 * t54;
t33 = t106 * t49 + t109 * t98;
t104 = t33 * t54;
t58 = t76 * qJ(4) + t72;
t84 = (-pkin(9) - qJ(4)) * t75;
t40 = -t106 * t84 + t77 * t58;
t103 = t40 * t54;
t41 = t106 * t58 + t77 * t84;
t102 = t41 * t54;
t44 = -t107 * t61 + t78 * t62;
t101 = t44 * t76;
t97 = t121 * t33 + t32 * t53;
t96 = t121 * t41 + t40 * t53;
t34 = -pkin(5) * t121 - t53 * qJ(6) + t66;
t27 = -t87 + t34;
t95 = t27 + t34;
t94 = t57 + t66;
t93 = t91 * t65;
t89 = 0.2e1 * t108;
t24 = pkin(4) * t100 + t44;
t82 = -pkin(3) * t55 - qJ(4) * t54;
t7 = -t17 * t75 + t18 * t76;
t81 = -t54 * t65 + t55 * t68;
t50 = t53 ^ 2;
t43 = t53 * t54;
t42 = t121 * t54;
t39 = t121 * t119;
t38 = t44 * t75;
t36 = -pkin(5) * t53 + qJ(6) * t121;
t29 = t121 * t55;
t28 = t53 * t55;
t21 = t29 * t53;
t20 = t24 * t53;
t19 = t24 * t121;
t11 = t121 * t29 - t53 * t28;
t10 = t28 * pkin(5) - t29 * qJ(6) + t24;
t9 = t10 * t53;
t8 = t10 * t121;
t1 = [1, 0, 0, t79 ^ 2, t79 * t89, 0, 0, 0, pkin(1) * t89, -0.2e1 * pkin(1) * t79, t55 ^ 2, t55 * t118, 0, 0, 0, t54 * t115, t55 * t115, 0.2e1 * t44 * t100 + 0.2e1 * t17 * t54, -0.2e1 * t18 * t54 + 0.2e1 * t44 * t99, 0.2e1 * (-t17 * t76 - t18 * t75) * t55, t17 ^ 2 + t18 ^ 2 + t44 ^ 2, t29 ^ 2, -0.2e1 * t29 * t28, 0.2e1 * t29 * t54, t28 * t118, t54 ^ 2, 0.2e1 * t24 * t28 - 0.2e1 * t54 * t83, 0.2e1 * t24 * t29 - 0.2e1 * t6 * t54, 0.2e1 * t10 * t28 - 0.2e1 * t4 * t54, -0.2e1 * t3 * t28 + 0.2e1 * t4 * t29, -0.2e1 * t10 * t29 + 0.2e1 * t3 * t54, t10 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, t79, t108, 0, -t79 * pkin(7), -t88, 0, 0, t55, -t54, 0, -t44, -t45, t81 * t75 - t101, t81 * t76 + t38, t7, t44 * t68 + t7 * t65, t21, t11, t43, t42, 0, t57 * t28 - t105 - t19, t57 * t29 - t104 + t20, t27 * t28 - t105 - t8, -t33 * t28 + t32 * t29 + t114, -t27 * t29 + t104 - t9, t10 * t27 + t3 * t33 + t4 * t32; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t87, -0.2e1 * t111, -0.2e1 * t68 * t76, 0.2e1 * t68 * t75, 0.2e1 * t93, t91 * t65 ^ 2 + t68 ^ 2, t50, -t39, 0, 0, 0, -t121 * t117, t53 * t117, t27 * t120, 0.2e1 * t97, t27 * t119, t27 ^ 2 + t32 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t54, 0, -t44, -t45, t82 * t75 - t101, t82 * t76 + t38, t7, -t44 * pkin(3) + t7 * qJ(4), t21, t11, t43, t42, 0, t66 * t28 - t103 - t19, t66 * t29 - t102 + t20, t34 * t28 - t103 - t8, -t41 * t28 + t40 * t29 + t114, -t34 * t29 + t102 - t9, t10 * t34 + t3 * t41 + t4 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t87, -t111, t110 * t76, -t110 * t75, t92 + t93, -t68 * pkin(3) + t65 * t92, t50, -t39, 0, 0, 0, -t94 * t121, t94 * t53, -t95 * t121, t96 + t97, -t95 * t53, t27 * t34 + t32 * t40 + t33 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t76, -0.2e1 * pkin(3) * t75, 0.2e1 * t92, t91 * qJ(4) ^ 2 + pkin(3) ^ 2, t50, -t39, 0, 0, 0, -t121 * t116, t53 * t116, t34 * t120, 0.2e1 * t96, t34 * t119, t34 ^ 2 + t40 ^ 2 + t41 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t99, 0, t44, 0, 0, 0, 0, 0, t28, t29, t28, 0, -t29, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t75, 0, t68, 0, 0, 0, 0, 0, -t121, t53, -t121, 0, -t53, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t75, 0, -pkin(3), 0, 0, 0, 0, 0, -t121, t53, -t121, 0, -t53, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, t54, -t83, -t6, -t83 + 0.2e1 * t113, -pkin(5) * t29 - t28 * qJ(6), 0.2e1 * t90 + t6, -t4 * pkin(5) + t3 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t121, 0, -t32, -t33, -t32, t36, t33, -t32 * pkin(5) + t33 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t121, 0, -t40, -t41, -t40, t36, t41, -t40 * pkin(5) + t41 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t29, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t1;
