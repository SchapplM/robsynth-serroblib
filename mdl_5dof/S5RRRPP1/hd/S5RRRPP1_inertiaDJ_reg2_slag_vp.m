% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:57
% EndTime: 2019-12-31 20:50:00
% DurationCPUTime: 0.79s
% Computational Cost: add. (930->105), mult. (2230->170), div. (0->0), fcn. (1784->6), ass. (0->87)
t115 = cos(pkin(8));
t90 = cos(qJ(3));
t82 = t90 * qJD(4);
t88 = sin(qJ(3));
t117 = pkin(1) * qJD(2);
t91 = cos(qJ(2));
t110 = t91 * t117;
t97 = t90 * t110;
t89 = sin(qJ(2));
t77 = t89 * pkin(1) + pkin(7);
t116 = -qJ(4) - t77;
t98 = qJD(3) * t116;
t43 = t88 * t98 + t82 + t97;
t87 = sin(pkin(8));
t92 = (-qJD(4) - t110) * t88 + t90 * t98;
t12 = -t115 * t92 + t87 * t43;
t13 = t115 * t43 + t87 * t92;
t103 = t115 * t88;
t84 = t90 * qJ(4);
t60 = t90 * t77 + t84;
t36 = -t116 * t103 + t87 * t60;
t126 = t87 * t88;
t37 = t115 * t60 + t116 * t126;
t62 = t87 * t90 + t103;
t58 = t62 * qJD(3);
t102 = t115 * t90;
t114 = t88 * qJD(3);
t59 = qJD(3) * t102 - t87 * t114;
t61 = -t102 + t126;
t137 = t12 * t62 - t13 * t61 + t36 * t59 - t37 * t58;
t142 = 0.2e1 * t137;
t125 = -qJ(4) - pkin(7);
t101 = qJD(3) * t125;
t57 = t88 * t101 + t82;
t93 = -t88 * qJD(4) + t90 * t101;
t33 = -t115 * t93 + t87 * t57;
t34 = t115 * t57 + t87 * t93;
t69 = t90 * pkin(7) + t84;
t44 = -t125 * t103 + t87 * t69;
t45 = t115 * t69 + t125 * t126;
t138 = t33 * t62 - t34 * t61 + t44 * t59 - t45 * t58;
t141 = 0.2e1 * t138;
t140 = t138 + t137;
t139 = 0.2e1 * t62 * t58 + 0.2e1 * t59 * t61;
t136 = 2 * qJD(5);
t135 = t91 * pkin(1);
t80 = pkin(3) * t114;
t18 = t58 * pkin(4) - t59 * qJ(5) - t62 * qJD(5) + t80;
t81 = t89 * t117;
t16 = t18 + t81;
t79 = -t90 * pkin(3) - pkin(2);
t38 = t61 * pkin(4) - t62 * qJ(5) + t79;
t35 = t38 - t135;
t133 = t16 * t61 + t35 * t58;
t132 = -t16 * t62 - t35 * t59;
t131 = t18 * t61 + t38 * t58;
t124 = -t18 * t62 - t38 * t59;
t65 = t81 + t80;
t68 = t79 - t135;
t122 = t68 * t58 + t65 * t61;
t121 = t68 * t59 + t65 * t62;
t120 = t79 * t58 + t61 * t80;
t119 = t79 * t59 + t62 * t80;
t78 = -pkin(2) - t135;
t83 = t90 * qJD(3);
t118 = t78 * t83 + t88 * t81;
t113 = 0.2e1 * t61 * t58;
t112 = pkin(2) * t114;
t111 = pkin(2) * t83;
t109 = t88 * t83;
t108 = t36 * t12 + t37 * t13;
t85 = t88 ^ 2;
t86 = t90 ^ 2;
t107 = (t85 + t86) * t91;
t104 = t44 * t33 + t45 * t34;
t95 = t78 * t114 - t90 * t81;
t94 = t12 * t44 + t13 * t45 + t36 * t33 + t37 * t34;
t76 = -t115 * pkin(3) - pkin(4);
t74 = t87 * pkin(3) + qJ(5);
t72 = -0.2e1 * t109;
t71 = 0.2e1 * t109;
t64 = 0.2e1 * (-t85 + t86) * qJD(3);
t56 = t107 * t117;
t42 = 0.2e1 * t62 * t59;
t32 = (-t115 * t59 - t58 * t87) * pkin(3);
t17 = -qJD(5) * t61 - t74 * t58 + t76 * t59;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t81, -0.2e1 * t110, 0, 0, t71, t64, 0, t72, 0, 0, 0.2e1 * t95, 0.2e1 * t118, 0.2e1 * t56, 0.2e1 * (t77 * t107 + t78 * t89) * t117, t42, -t139, 0, t113, 0, 0, 0.2e1 * t122, 0.2e1 * t121, t142, 0.2e1 * t68 * t65 + 0.2e1 * t108, t42, 0, t139, 0, 0, t113, 0.2e1 * t133, t142, 0.2e1 * t132, 0.2e1 * t35 * t16 + 0.2e1 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t110, 0, 0, t71, t64, 0, t72, 0, 0, t95 - t112, -t111 + t118, t56, (-pkin(2) * t89 + pkin(7) * t107) * t117, t42, -t139, 0, t113, 0, 0, t120 + t122, t119 + t121, t140, t65 * t79 + t68 * t80 + t94, t42, 0, t139, 0, 0, t113, t131 + t133, t140, t124 + t132, t16 * t38 + t35 * t18 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t64, 0, t72, 0, 0, -0.2e1 * t112, -0.2e1 * t111, 0, 0, t42, -t139, 0, t113, 0, 0, 0.2e1 * t120, 0.2e1 * t119, t141, 0.2e1 * t79 * t80 + 0.2e1 * t104, t42, 0, t139, 0, 0, t113, 0.2e1 * t131, t141, 0.2e1 * t124, 0.2e1 * t38 * t18 + 0.2e1 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, -t114, 0, -t88 * t110 - t77 * t83, t77 * t114 - t97, 0, 0, 0, 0, t59, 0, -t58, 0, -t12, -t13, t32, (-t115 * t12 + t13 * t87) * pkin(3), 0, t59, 0, 0, t58, 0, -t12, t17, t13, t37 * qJD(5) + t12 * t76 + t13 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, -t114, 0, -pkin(7) * t83, pkin(7) * t114, 0, 0, 0, 0, t59, 0, -t58, 0, -t33, -t34, t32, (-t115 * t33 + t34 * t87) * pkin(3), 0, t59, 0, 0, t58, 0, -t33, t17, t34, t45 * qJD(5) + t33 * t76 + t34 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t74 * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t59, 0, t65, 0, 0, 0, 0, 0, 0, t58, 0, -t59, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t59, 0, t80, 0, 0, 0, 0, 0, 0, t58, 0, -t59, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
