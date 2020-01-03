% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:13
% EndTime: 2019-12-31 21:51:15
% DurationCPUTime: 0.56s
% Computational Cost: add. (1020->127), mult. (2320->175), div. (0->0), fcn. (1868->6), ass. (0->85)
t122 = qJD(3) + qJD(4);
t79 = 2 * qJD(5);
t121 = pkin(8) + pkin(7);
t78 = cos(qJ(3));
t120 = t78 * pkin(3);
t115 = cos(qJ(4));
t100 = t115 * t78;
t76 = sin(qJ(4));
t77 = sin(qJ(3));
t114 = t76 * t77;
t97 = t115 * qJD(4);
t35 = -qJD(3) * t100 + t122 * t114 - t78 * t97;
t113 = t76 * t78;
t54 = t115 * t77 + t113;
t36 = t122 * t54;
t106 = t77 * qJD(3);
t72 = pkin(3) * t106;
t12 = t36 * pkin(4) + t35 * qJ(5) - t54 * qJD(5) + t72;
t101 = sin(qJ(2)) * pkin(1);
t73 = qJD(2) * t101;
t10 = t12 + t73;
t102 = cos(qJ(2)) * pkin(1);
t53 = -t100 + t114;
t70 = -pkin(2) - t120;
t34 = t53 * pkin(4) - t54 * qJ(5) + t70;
t29 = -t102 + t34;
t119 = t10 * t53 + t29 * t36;
t118 = -t10 * t54 + t29 * t35;
t117 = t12 * t53 + t34 * t36;
t116 = -t12 * t54 + t34 * t35;
t55 = t73 + t72;
t69 = -t102 - pkin(2);
t58 = t69 - t120;
t112 = t58 * t36 + t55 * t53;
t111 = -t58 * t35 + t55 * t54;
t110 = t70 * t36 + t53 * t72;
t109 = -t70 * t35 + t54 * t72;
t74 = t78 * qJD(3);
t108 = t69 * t74 + t77 * t73;
t107 = qJD(4) * t76;
t105 = pkin(2) * t106;
t104 = pkin(2) * t74;
t103 = pkin(3) * t107;
t99 = t121 * qJD(3);
t75 = t78 * pkin(8);
t93 = t101 + pkin(7);
t50 = t78 * t93 + t75;
t90 = -pkin(8) - t93;
t86 = t90 * t77;
t82 = t115 * t86;
t32 = t76 * t50 - t82;
t83 = t76 * t86;
t33 = t115 * t50 + t83;
t84 = t90 * qJD(3);
t95 = qJD(2) * t102;
t91 = t78 * t95;
t80 = t77 * t84 + t91;
t92 = t77 * t95;
t81 = -t78 * t84 + t92;
t8 = -qJD(4) * t82 + t50 * t107 - t115 * t80 + t76 * t81;
t9 = qJD(4) * t83 + t115 * t81 + t50 * t97 + t76 * t80;
t98 = -t32 * t35 - t33 * t36 + t8 * t53 + t9 * t54;
t59 = t78 * pkin(7) + t75;
t94 = t121 * t115;
t88 = qJD(3) * t94;
t89 = t77 * t94;
t22 = qJD(4) * t89 + t59 * t107 + t99 * t113 + t77 * t88;
t23 = t59 * t97 + t78 * t88 + (-qJD(4) * t121 - t99) * t114;
t40 = t76 * t59 + t89;
t41 = -t121 * t114 + t115 * t59;
t96 = t22 * t53 + t23 * t54 - t40 * t35 - t41 * t36;
t87 = t93 * qJD(3);
t85 = t69 * t106 - t78 * t73;
t71 = pkin(3) * t97;
t68 = -t115 * pkin(3) - pkin(4);
t66 = t76 * pkin(3) + qJ(5);
t65 = -0.2e1 * t103;
t62 = t71 + qJD(5);
t61 = 0.2e1 * t77 * t74;
t52 = 0.2e1 * (-t77 ^ 2 + t78 ^ 2) * qJD(3);
t26 = -0.2e1 * t54 * t35;
t13 = pkin(4) * t35 - t36 * qJ(5) - t53 * qJD(5);
t11 = 0.2e1 * t35 * t53 - 0.2e1 * t54 * t36;
t1 = t54 * t103 - t68 * t35 - t66 * t36 - t62 * t53;
t2 = [0, 0, 0, 0, -0.2e1 * t73, -0.2e1 * t95, t61, t52, 0, 0, 0, 0.2e1 * t85, 0.2e1 * t108, t26, t11, 0, 0, 0, 0.2e1 * t112, 0.2e1 * t111, 0.2e1 * t119, 0.2e1 * t98, 0.2e1 * t118, 0.2e1 * t29 * t10 + 0.2e1 * t32 * t9 - 0.2e1 * t33 * t8; 0, 0, 0, 0, -t73, -t95, t61, t52, 0, 0, 0, t85 - t105, -t104 + t108, t26, t11, 0, 0, 0, t110 + t112, t109 + t111, t117 + t119, t96 + t98, t116 + t118, t10 * t34 + t29 * t12 - t33 * t22 + t32 * t23 + t9 * t40 - t8 * t41; 0, 0, 0, 0, 0, 0, t61, t52, 0, 0, 0, -0.2e1 * t105, -0.2e1 * t104, t26, t11, 0, 0, 0, 0.2e1 * t110, 0.2e1 * t109, 0.2e1 * t117, 0.2e1 * t96, 0.2e1 * t116, 0.2e1 * t34 * t12 - 0.2e1 * t41 * t22 + 0.2e1 * t40 * t23; 0, 0, 0, 0, 0, 0, 0, 0, t74, -t106, 0, -t78 * t87 - t92, t77 * t87 - t91, 0, 0, -t35, -t36, 0, -t9, t8, -t9, t1, -t8, t32 * t103 + t33 * t62 - t8 * t66 + t9 * t68; 0, 0, 0, 0, 0, 0, 0, 0, t74, -t106, 0, -pkin(7) * t74, pkin(7) * t106, 0, 0, -t35, -t36, 0, -t23, t22, -t23, t1, -t22, t40 * t103 - t22 * t66 + t23 * t68 + t41 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -0.2e1 * t71, t65, 0, 0.2e1 * t62, 0.2e1 * t68 * t103 + 0.2e1 * t66 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, -t9, t8, -t9, t13, -t8, -t9 * pkin(4) - t8 * qJ(5) + t33 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, -t23, t22, -t23, t13, -t22, -t23 * pkin(4) - t22 * qJ(5) + t41 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t71, -t103, 0, t79 + t71, -pkin(4) * t103 + t62 * qJ(5) + t66 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, qJ(5) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
