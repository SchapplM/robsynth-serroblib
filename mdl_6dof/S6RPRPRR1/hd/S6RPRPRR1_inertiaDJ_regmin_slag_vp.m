% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:19
% EndTime: 2019-03-09 03:35:22
% DurationCPUTime: 0.88s
% Computational Cost: add. (1769->135), mult. (3813->238), div. (0->0), fcn. (3851->10), ass. (0->97)
t111 = sin(pkin(11));
t100 = t111 * pkin(3);
t121 = sin(qJ(5));
t122 = cos(qJ(5));
t112 = cos(pkin(11));
t92 = t112 * pkin(3) + pkin(4);
t124 = t121 * t100 - t122 * t92;
t73 = cos(qJ(6));
t69 = t73 ^ 2;
t71 = sin(qJ(6));
t114 = t71 ^ 2 - t69;
t94 = t114 * qJD(6);
t64 = sin(pkin(10)) * pkin(1) + pkin(7);
t113 = qJ(4) + t64;
t72 = sin(qJ(3));
t57 = t113 * t72;
t74 = cos(qJ(3));
t58 = t113 * t74;
t31 = -t111 * t58 - t112 * t57;
t81 = -t111 * t74 - t112 * t72;
t29 = pkin(8) * t81 + t31;
t32 = -t111 * t57 + t112 * t58;
t95 = t111 * t72;
t96 = t112 * t74;
t82 = t95 - t96;
t30 = -t82 * pkin(8) + t32;
t93 = qJD(3) * t113;
t43 = t74 * qJD(4) - t72 * t93;
t44 = -t72 * qJD(4) - t74 * t93;
t27 = -t111 * t43 + t112 * t44;
t63 = qJD(3) * t95;
t84 = qJD(3) * t96 - t63;
t75 = t84 * pkin(8) - t27;
t28 = t111 * t44 + t112 * t43;
t79 = t81 * qJD(3);
t77 = pkin(8) * t79 + t28;
t97 = qJD(5) * t121;
t98 = qJD(5) * t122;
t5 = t121 * t77 + t122 * t75 + t29 * t97 + t30 * t98;
t3 = t5 * t71;
t17 = t121 * t30 - t122 * t29;
t67 = qJD(6) * t73;
t123 = t17 * t67 + t3;
t78 = t122 * t82;
t21 = qJD(5) * t78 - t121 * t79 - t122 * t84 - t81 * t97;
t35 = -t121 * t82 - t122 * t81;
t120 = t35 * t21;
t22 = t35 * qJD(5) + t121 * t84 - t122 * t79;
t119 = t71 * t22;
t118 = t73 * t21;
t117 = t73 * t22;
t34 = -t121 * t81 + t78;
t116 = t35 * t117 - t34 * t118;
t76 = t122 * t100 + t121 * t92;
t46 = t76 * qJD(5);
t53 = -pkin(5) + t124;
t115 = t46 * t71 + t53 * t67;
t110 = qJD(6) * t71;
t109 = t72 * qJD(3);
t108 = t74 * qJD(3);
t107 = t112 ^ 2;
t106 = t71 * t118;
t105 = 0.2e1 * t108;
t104 = pkin(5) * t110;
t103 = pkin(5) * t67;
t66 = pkin(3) * t109;
t102 = t35 * t110;
t101 = t71 * t67;
t65 = -cos(pkin(10)) * pkin(1) - pkin(2);
t99 = t53 * t110 - t46 * t73;
t18 = t121 * t29 + t122 * t30;
t87 = -t74 * pkin(3) + t65;
t38 = t82 * pkin(4) + t87;
t19 = t34 * pkin(5) - t35 * pkin(9) + t38;
t91 = t73 * t18 + t71 * t19;
t90 = t71 * t18 - t73 * t19;
t89 = t21 * t34 - t35 * t22;
t54 = pkin(9) + t76;
t88 = t34 * t54 - t35 * t53;
t85 = -t71 * t21 + t35 * t67;
t10 = t102 + t118;
t11 = t34 * t67 + t119;
t9 = t34 * t110 - t117;
t45 = t124 * qJD(5);
t80 = -t21 * t53 - t22 * t54 + t34 * t45 + t35 * t46;
t39 = -pkin(4) * t79 + t66;
t62 = 0.2e1 * t101;
t60 = -0.2e1 * t94;
t33 = t35 ^ 2;
t13 = t17 * t110;
t8 = -t35 * t94 - t106;
t7 = t22 * pkin(5) + t21 * pkin(9) + t39;
t6 = -0.4e1 * t35 * t101 + t114 * t21;
t4 = t121 * t75 - t122 * t77 - t29 * t98 + t30 * t97;
t2 = -t91 * qJD(6) + t71 * t4 + t73 * t7;
t1 = t90 * qJD(6) + t73 * t4 - t71 * t7;
t12 = [0, 0, 0, 0, t72 * t105, 0.2e1 * (-t72 ^ 2 + t74 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t65 * t109, t65 * t105, 0.2e1 * t27 * t81 - 0.2e1 * t28 * t82 - 0.2e1 * t31 * t84 + 0.2e1 * t32 * t79, 0.2e1 * t31 * t27 + 0.2e1 * t32 * t28 + 0.2e1 * t87 * t66, -0.2e1 * t120, 0.2e1 * t89, 0, 0, 0, 0.2e1 * t38 * t22 + 0.2e1 * t39 * t34, -0.2e1 * t38 * t21 + 0.2e1 * t39 * t35, -0.2e1 * t33 * t101 - 0.2e1 * t69 * t120, 0.4e1 * t35 * t106 + 0.2e1 * t33 * t94, -0.2e1 * t34 * t102 + 0.2e1 * t116, -0.2e1 * t35 * t119 - 0.2e1 * t85 * t34, 0.2e1 * t34 * t22, 0.2e1 * t85 * t17 + 0.2e1 * t2 * t34 - 0.2e1 * t90 * t22 + 0.2e1 * t35 * t3, 0.2e1 * t5 * t73 * t35 + 0.2e1 * t1 * t34 - 0.2e1 * t10 * t17 - 0.2e1 * t91 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27 * t82 - t28 * t81 + t31 * t79 + t32 * t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t73 + t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t63 + (-t82 - t96) * qJD(3)) * t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t108, -t109, 0, -t64 * t108, t64 * t109 (t63 * t112 + (-t107 * t74 + t111 * t81) * qJD(3)) * pkin(3) (t111 * t28 + t112 * t27) * pkin(3), 0, 0, -t21, -t22, 0, -t5, t4, t8, t6, t11, -t9, 0, t13 + (-t88 * qJD(6) - t5) * t73 + t80 * t71, t88 * t110 + t80 * t73 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t108, 0 (-t107 * t109 - t63 * t111) * pkin(3), 0, 0, 0, 0, 0, -t22, t21, 0, 0, 0, 0, 0, t9, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t46, 0.2e1 * t45, t62, t60, 0, 0, 0, 0.2e1 * t99, 0.2e1 * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, t22, -t21, 0, 0, 0, 0, 0, -t9, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, 0, -t5, t4, t8, t6, t11, -t9, 0, t13 + (pkin(5) * t21 - pkin(9) * t22) * t71 + (-t5 + (-pkin(5) * t35 - pkin(9) * t34) * qJD(6)) * t73, t10 * pkin(5) + t9 * pkin(9) + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, 0, 0, 0, 0, 0, t9, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t45, t62, t60, 0, 0, 0, t99 - t104, -t103 + t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t60, 0, 0, 0, -0.2e1 * t104, -0.2e1 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t85, t22, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t110, 0, t71 * t45 - t54 * t67, t54 * t110 + t73 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t110, 0, -pkin(9) * t67, pkin(9) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
