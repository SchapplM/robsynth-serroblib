% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:32:08
% EndTime: 2019-12-05 18:32:11
% DurationCPUTime: 0.77s
% Computational Cost: add. (1143->114), mult. (2556->195), div. (0->0), fcn. (2107->8), ass. (0->92)
t127 = qJD(4) + qJD(5);
t111 = pkin(1) * qJD(2);
t76 = cos(qJ(2));
t106 = t76 * t111;
t74 = sin(qJ(2));
t107 = t74 * t111;
t70 = sin(pkin(9));
t71 = cos(pkin(9));
t49 = t71 * t106 - t70 * t107;
t73 = sin(qJ(4));
t68 = t73 ^ 2;
t75 = cos(qJ(4));
t69 = t75 ^ 2;
t26 = (t69 + t68) * t49;
t126 = t75 * pkin(4);
t125 = cos(qJ(5));
t72 = sin(qJ(5));
t53 = t125 * t73 + t72 * t75;
t35 = t127 * t53;
t121 = t72 * t73;
t99 = t125 * t75;
t52 = -t99 + t121;
t124 = t52 * t35;
t93 = t125 * qJD(5);
t34 = -qJD(4) * t99 + t127 * t121 - t75 * t93;
t123 = t53 * t34;
t122 = t71 * t74;
t120 = t73 * t49;
t119 = t75 * t49;
t109 = t73 * qJD(4);
t105 = pkin(4) * t109;
t48 = (t70 * t76 + t122) * t111;
t37 = t48 + t105;
t65 = t76 * pkin(1) + pkin(2);
t94 = -t70 * t74 * pkin(1) + t71 * t65;
t45 = -pkin(3) - t94;
t40 = t45 - t126;
t117 = t40 * t35 + t37 * t52;
t116 = -t40 * t34 + t37 * t53;
t66 = t75 * qJD(4);
t115 = t45 * t66 + t48 * t73;
t63 = -t71 * pkin(2) - pkin(3);
t54 = t63 - t126;
t114 = t52 * t105 + t54 * t35;
t113 = t53 * t105 - t54 * t34;
t112 = pkin(1) * t122 + t70 * t65;
t110 = qJD(5) * t72;
t108 = pkin(7) + t112;
t104 = pkin(4) * t110;
t103 = t63 * t109;
t102 = t63 * t66;
t101 = t73 * t66;
t100 = t70 * pkin(2) + pkin(7);
t67 = t75 * pkin(8);
t36 = t75 * t108 + t67;
t97 = -pkin(8) - t108;
t89 = t97 * t73;
t81 = t125 * t89;
t16 = -t72 * t36 + t81;
t85 = t72 * t89;
t17 = t125 * t36 + t85;
t87 = qJD(4) * t97;
t77 = t73 * t87 + t119;
t78 = t75 * t87 - t120;
t4 = -qJD(5) * t81 + t36 * t110 - t125 * t77 - t72 * t78;
t5 = -qJD(5) * t85 + t125 * t78 - t36 * t93 - t72 * t77;
t98 = t16 * t34 - t17 * t35 + t4 * t52 - t5 * t53;
t50 = t75 * t100 + t67;
t92 = pkin(8) + t100;
t82 = t92 * t125;
t79 = qJD(4) * t82;
t80 = t73 * t82;
t86 = t72 * t92;
t83 = qJD(4) * t86;
t12 = qJD(5) * t80 + t50 * t110 + t73 * t79 + t75 * t83;
t84 = t73 * t86;
t13 = qJD(5) * t84 - t50 * t93 + t73 * t83 - t75 * t79;
t31 = -t72 * t50 - t80;
t32 = t125 * t50 - t84;
t96 = t12 * t52 - t13 * t53 + t31 * t34 - t32 * t35;
t95 = t45 * t109 - t48 * t75;
t91 = pkin(4) * t93;
t90 = t108 * qJD(4);
t88 = t100 * qJD(4);
t58 = -0.2e1 * t101;
t57 = 0.2e1 * t101;
t51 = 0.2e1 * (-t68 + t69) * qJD(4);
t23 = -0.2e1 * t123;
t22 = 0.2e1 * t124;
t11 = 0.2e1 * t52 * t34 - 0.2e1 * t53 * t35;
t6 = (t125 * t34 - t35 * t72 + (-t125 * t52 + t53 * t72) * qJD(5)) * pkin(4);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t107, -0.2e1 * t106, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t48, -0.2e1 * t49, 0, 0.2e1 * t112 * t49 - 0.2e1 * t94 * t48, t57, t51, 0, t58, 0, 0, 0.2e1 * t95, 0.2e1 * t115, 0.2e1 * t26, 0.2e1 * t108 * t26 + 0.2e1 * t45 * t48, t23, t11, 0, t22, 0, 0, 0.2e1 * t117, 0.2e1 * t116, 0.2e1 * t98, 0.2e1 * t16 * t5 - 0.2e1 * t17 * t4 + 0.2e1 * t40 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t106, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t49, 0, (-t48 * t71 + t49 * t70) * pkin(2), t57, t51, 0, t58, 0, 0, t95 + t103, t102 + t115, t26, t100 * t26 + t48 * t63, t23, t11, 0, t22, 0, 0, t114 + t117, t113 + t116, t96 + t98, t40 * t105 - t17 * t12 + t16 * t13 + t5 * t31 - t4 * t32 + t37 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t51, 0, t58, 0, 0, 0.2e1 * t103, 0.2e1 * t102, 0, 0, t23, t11, 0, t22, 0, 0, 0.2e1 * t114, 0.2e1 * t113, 0.2e1 * t96, 0.2e1 * t54 * t105 - 0.2e1 * t32 * t12 + 0.2e1 * t31 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t35 - t17 * t34 - t4 * t53 - t5 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t53 - t13 * t52 - t31 * t35 - t32 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t123 + 0.2e1 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, -t109, 0, -t75 * t90 - t120, t73 * t90 - t119, 0, 0, 0, 0, -t34, 0, -t35, 0, t5, t4, t6, (t125 * t5 - t4 * t72 + (t125 * t17 - t16 * t72) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, -t109, 0, -t75 * t88, t73 * t88, 0, 0, 0, 0, -t34, 0, -t35, 0, t13, t12, t6, (t125 * t13 - t12 * t72 + (t125 * t32 - t31 * t72) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t66, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t34, 0, (-t125 * t35 - t34 * t72 + (t125 * t53 + t52 * t72) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t104, -0.2e1 * t91, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, -t35, 0, t5, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, -t35, 0, t13, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t91, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
