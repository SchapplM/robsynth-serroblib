% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:16:51
% EndTime: 2019-05-05 10:16:55
% DurationCPUTime: 1.26s
% Computational Cost: add. (1312->177), mult. (3297->333), div. (0->0), fcn. (3855->12), ass. (0->107)
t58 = sin(pkin(7));
t68 = cos(qJ(3));
t101 = t58 * t68;
t64 = sin(qJ(3));
t102 = t58 * t64;
t60 = cos(pkin(7));
t63 = sin(qJ(4));
t67 = cos(qJ(4));
t38 = t67 * t102 + t63 * t60;
t62 = sin(qJ(5));
t66 = cos(qJ(5));
t24 = t66 * t101 + t62 * t38;
t123 = -0.2e1 * t24;
t122 = -0.2e1 * t38;
t73 = -t66 * pkin(5) - t62 * qJ(6);
t44 = -pkin(4) + t73;
t121 = -0.2e1 * t44;
t120 = -0.2e1 * t63;
t119 = 0.2e1 * t67;
t118 = pkin(2) * t64;
t117 = pkin(2) * t68;
t116 = pkin(4) * t62;
t115 = pkin(4) * t66;
t114 = pkin(10) * t62;
t113 = pkin(10) * t66;
t37 = t63 * t102 - t67 * t60;
t112 = t37 * pkin(5);
t111 = t62 * pkin(11);
t110 = t66 * pkin(11);
t59 = sin(pkin(6));
t61 = cos(pkin(6));
t65 = sin(qJ(2));
t69 = cos(qJ(2));
t97 = t60 * t69;
t23 = t61 * t102 + (t64 * t97 + t65 * t68) * t59;
t99 = t59 * t69;
t36 = -t58 * t99 + t61 * t60;
t14 = t23 * t63 - t36 * t67;
t109 = t14 * t62;
t108 = t14 * t66;
t82 = pkin(9) * t101;
t33 = t82 + (pkin(10) + t118) * t60;
t34 = (-pkin(3) * t68 - pkin(10) * t64 - pkin(2)) * t58;
t19 = -t63 * t33 + t67 * t34;
t17 = pkin(4) * t101 - t19;
t107 = t17 * t62;
t106 = t17 * t66;
t105 = t24 * t66;
t25 = -t62 * t101 + t66 * t38;
t104 = t25 * t62;
t54 = t58 ^ 2;
t103 = t54 * t68;
t100 = t59 * t65;
t98 = t60 * t64;
t96 = t62 * t37;
t95 = t62 * t63;
t94 = t62 * t66;
t93 = t62 * t67;
t92 = t63 * t37;
t91 = t66 * t37;
t53 = t66 * t63;
t90 = t66 * t67;
t48 = pkin(9) * t102;
t32 = t48 + (-pkin(3) - t117) * t60;
t16 = t37 * pkin(4) - t38 * pkin(11) + t32;
t20 = t67 * t33 + t63 * t34;
t18 = -pkin(11) * t101 + t20;
t5 = t62 * t16 + t66 * t18;
t45 = -t67 * pkin(4) - t63 * pkin(11) - pkin(3);
t30 = pkin(10) * t90 + t62 * t45;
t55 = t62 ^ 2;
t57 = t66 ^ 2;
t89 = t55 + t57;
t88 = t37 * qJ(6);
t87 = t67 * qJ(6);
t86 = 0.2e1 * t101;
t85 = t63 * t119;
t84 = pkin(11) * t96;
t83 = pkin(11) * t91;
t81 = t63 * t101;
t80 = t67 * t101;
t15 = t23 * t67 + t36 * t63;
t22 = -t59 * t68 * t97 + t64 * t100 - t61 * t101;
t7 = t15 * t62 - t22 * t66;
t79 = t14 * t24 - t7 * t37;
t78 = t14 * t95 + t7 * t67;
t77 = -t66 * t16 + t62 * t18;
t2 = t88 + t5;
t3 = t77 - t112;
t76 = t2 * t66 + t3 * t62;
t8 = t15 * t66 + t22 * t62;
t75 = t7 * t62 + t8 * t66;
t74 = t14 * t25 - t8 * t37;
t72 = -pkin(5) * t62 + t66 * qJ(6);
t26 = -t87 + t30;
t42 = t66 * t45;
t27 = -t42 + (pkin(5) + t114) * t67;
t71 = t26 * t66 + t27 * t62;
t56 = t63 ^ 2;
t50 = pkin(11) * t93;
t40 = pkin(2) * t98 + t82;
t39 = t60 * t117 - t48;
t35 = (pkin(10) - t72) * t63;
t29 = -pkin(10) * t93 + t42;
t6 = t24 * pkin(5) - t25 * qJ(6) + t17;
t1 = t14 * t53 + t8 * t67;
t4 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t7 ^ 2 + t8 ^ 2; 0, 0, t99, -t100, 0, 0, 0, 0, 0, -t36 * t101 - t22 * t60, t36 * t102 - t23 * t60, 0, 0, 0, 0, 0, t14 * t101 + t22 * t37, t15 * t101 + t22 * t38, 0, 0, 0, 0, 0, t79, t74, t79, -t8 * t24 + t7 * t25, -t74, t14 * t6 + t8 * t2 + t7 * t3; 0, 1, 0, 0, t54 * t64 ^ 2, 0.2e1 * t64 * t103, 0.2e1 * t58 * t98, t60 * t86, t60 ^ 2, 0.2e1 * pkin(2) * t103 + 0.2e1 * t39 * t60, -0.2e1 * t54 * t118 - 0.2e1 * t40 * t60, t38 ^ 2, t37 * t122, t101 * t122, t37 * t86, t54 * t68 ^ 2, -0.2e1 * t19 * t101 + 0.2e1 * t32 * t37, 0.2e1 * t20 * t101 + 0.2e1 * t32 * t38, t25 ^ 2, t25 * t123, 0.2e1 * t25 * t37, t37 * t123, t37 ^ 2, 0.2e1 * t17 * t24 - 0.2e1 * t37 * t77, 0.2e1 * t17 * t25 - 0.2e1 * t5 * t37, 0.2e1 * t6 * t24 - 0.2e1 * t3 * t37, -0.2e1 * t2 * t24 + 0.2e1 * t3 * t25, 0.2e1 * t2 * t37 - 0.2e1 * t6 * t25, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, 0, 0, 0, 0, -t22 * t67, t22 * t63, 0, 0, 0, 0, 0, t78, t1, t78 (-t62 * t8 + t66 * t7) * t63, -t1, t14 * t35 + t8 * t26 + t7 * t27; 0, 0, 0, 0, 0, 0, t102, t101, t60, t39, -t40, t38 * t63, t38 * t67 - t92, -t81, -t80, 0, -pkin(3) * t37 + pkin(10) * t81 - t32 * t67, -pkin(3) * t38 + pkin(10) * t80 + t32 * t63, t25 * t53 (-t104 - t105) * t63, -t25 * t67 + t37 * t53, t24 * t67 - t62 * t92, -t37 * t67, t29 * t37 + t77 * t67 + (pkin(10) * t24 + t107) * t63, -t30 * t37 + t5 * t67 + (pkin(10) * t25 + t106) * t63, t35 * t24 - t27 * t37 + t3 * t67 + t6 * t95, -t26 * t24 + t27 * t25 + (-t2 * t62 + t3 * t66) * t63, -t2 * t67 - t35 * t25 + t26 * t37 - t6 * t53, t2 * t26 + t3 * t27 + t6 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t56, t85, 0, 0, 0, pkin(3) * t119, pkin(3) * t120, t57 * t56, -0.2e1 * t56 * t94, t90 * t120, t62 * t85, t67 ^ 2, 0.2e1 * t56 * t114 - 0.2e1 * t29 * t67, 0.2e1 * t56 * t113 + 0.2e1 * t30 * t67, 0.2e1 * t27 * t67 + 0.2e1 * t35 * t95, 0.2e1 * (-t26 * t62 + t27 * t66) * t63, -0.2e1 * t26 * t67 - 0.2e1 * t35 * t53, t26 ^ 2 + t27 ^ 2 + t35 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t108, t109, -t108, t75, -t109, pkin(11) * t75 + t14 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, -t101, t19, -t20, t104, -t62 * t24 + t25 * t66, t96, t91, 0, -pkin(4) * t24 - t106 - t84, -pkin(4) * t25 + t107 - t83, t44 * t24 - t6 * t66 - t84 (t104 - t105) * pkin(11) + t76, -t44 * t25 - t6 * t62 + t83, pkin(11) * t76 + t6 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t67, 0, -t63 * pkin(10), -t67 * pkin(10), t62 * t53 (-t55 + t57) * t63, -t93, -t90, 0, t50 + (-t113 - t116) * t63, pkin(11) * t90 + (t114 - t115) * t63, -t35 * t66 + t44 * t95 + t50, t71, -t35 * t62 + (-pkin(11) * t67 - t44 * t63) * t66, pkin(11) * t71 + t35 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t55, 0.2e1 * t94, 0, 0, 0, 0.2e1 * t115, -0.2e1 * t116, t66 * t121, 0.2e1 * t89 * pkin(11), t62 * t121, pkin(11) ^ 2 * t89 + t44 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -t7, 0, t8, -t7 * pkin(5) + t8 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, t37, -t77, -t5, -t77 + 0.2e1 * t112, -pkin(5) * t25 - t24 * qJ(6), 0.2e1 * t88 + t5, -t3 * pkin(5) + t2 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t95, -t67, t29, -t30, t42 + (-0.2e1 * pkin(5) - t114) * t67, t73 * t63, -0.2e1 * t87 + t30, -t27 * pkin(5) + t26 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t66, 0, -t111, -t110, -t111, t72, t110, t72 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t25, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t53, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
