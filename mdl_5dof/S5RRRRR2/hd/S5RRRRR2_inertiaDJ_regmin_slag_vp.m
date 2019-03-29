% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_inertiaDJ_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:31
% EndTime: 2019-03-29 15:26:33
% DurationCPUTime: 0.71s
% Computational Cost: add. (695->85), mult. (2309->184), div. (0->0), fcn. (2039->8), ass. (0->98)
t70 = sin(qJ(4));
t71 = sin(qJ(3));
t120 = t70 * t71;
t72 = sin(qJ(2));
t127 = pkin(1) * t72;
t128 = qJD(3) + qJD(4);
t74 = cos(qJ(4));
t75 = cos(qJ(3));
t76 = cos(qJ(2));
t112 = qJD(2) * t76;
t97 = t75 * t112;
t98 = t71 * t112;
t18 = (t70 * t97 + (t128 * t75 * t72 + t98) * t74) * pkin(1) - t128 * t120 * t127;
t56 = t70 * t75 + t74 * t71;
t44 = t56 * t127;
t73 = cos(qJ(5));
t66 = qJD(5) * t73;
t69 = sin(qJ(5));
t117 = t18 * t69 + t44 * t66;
t125 = t75 * pkin(2);
t32 = t128 * t56;
t105 = t71 * qJD(3);
t101 = pkin(2) * t105;
t55 = -t74 * t75 + t120;
t46 = t55 * t101;
t85 = -t32 * t125 + t46;
t68 = t73 ^ 2;
t114 = t69 ^ 2 - t68;
t88 = t114 * qJD(5);
t126 = pkin(2) * t55;
t111 = qJD(4) * t55;
t31 = -t55 * qJD(3) - t111;
t124 = t56 * t31;
t123 = t56 * t73;
t122 = t69 * t31;
t121 = t69 * t32;
t119 = t73 * t31;
t118 = t73 * t32;
t113 = qJD(2) * t72;
t103 = pkin(1) * t113;
t57 = t101 + t103;
t61 = -t76 * pkin(1) - t125;
t116 = t61 * t32 + t57 * t55;
t115 = t61 * t31 + t57 * t56;
t110 = qJD(4) * t70;
t109 = qJD(4) * t74;
t108 = qJD(5) * t69;
t107 = qJD(5) * t74;
t106 = qJD(5) * t75;
t65 = t75 * qJD(3);
t102 = pkin(1) * t112;
t100 = pkin(2) * t110;
t99 = pkin(2) * t109;
t96 = t73 * t109;
t95 = t70 * t108;
t94 = t69 * t107;
t93 = t69 * t106;
t35 = t44 * t108;
t92 = t69 * t66;
t91 = -0.4e1 * t69 * t123;
t90 = -t18 * t73 + t35;
t87 = t56 * t100;
t86 = t56 * t101 - t31 * t125;
t17 = (t32 * t72 + t70 * t98 - t74 * t97) * pkin(1);
t45 = t55 * t127;
t83 = -t73 * t45 + t69 * t61;
t6 = -t83 * qJD(5) + t69 * t17 + t73 * t57;
t82 = -t69 * t45 - t73 * t61;
t84 = t117 * t56 + t44 * t122 - t82 * t32 + t6 * t55;
t81 = t93 * t126 + t85 * t73;
t24 = -t56 * t66 - t122;
t22 = -t56 * t108 + t119;
t23 = t55 * t66 + t121;
t5 = t82 * qJD(5) + t73 * t17 - t69 * t57;
t80 = t44 * t119 + t18 * t123 - t83 * t32 - t56 * t35 + t5 * t55;
t79 = t23 * t125 - t69 * t46;
t78 = -t96 * t126 + t73 * t87 + (-t70 * t118 - t74 * t119 + t55 * t95 + t56 * t94) * pkin(2);
t77 = t69 * t87 + ((-t55 * t70 - t56 * t74) * t66 + (-t32 * t70 + (-t31 - t111) * t74) * t69) * pkin(2);
t63 = 0.2e1 * t71 * t65;
t62 = 0.2e1 * t92;
t54 = 0.2e1 * (-t71 ^ 2 + t75 ^ 2) * qJD(3);
t53 = -0.2e1 * t88;
t52 = t56 ^ 2;
t51 = (-t76 * t105 - t75 * t113) * pkin(1);
t50 = (t71 * t113 - t76 * t65) * pkin(1);
t49 = (-t73 * t110 - t94) * pkin(2);
t48 = (-t73 * t107 + t69 * t110) * pkin(2);
t26 = 0.2e1 * t124;
t25 = 0.2e1 * t55 * t32;
t21 = -t55 * t108 + t118;
t13 = 0.2e1 * t68 * t124 - 0.2e1 * t52 * t92;
t12 = t69 * t119 - t56 * t88;
t9 = t31 * t91 + 0.2e1 * t52 * t88;
t8 = qJD(5) * t91 - t114 * t31;
t7 = -0.2e1 * t31 * t55 - 0.2e1 * t56 * t32;
t4 = -0.2e1 * t56 * t121 + 0.2e1 * t24 * t55;
t3 = 0.2e1 * t56 * t118 + 0.2e1 * t22 * t55;
t1 = [0, 0, 0, 0, -0.2e1 * t103, -0.2e1 * t102, t63, t54, 0, 0, 0, 0.2e1 * t51, 0.2e1 * t50, t26, t7, 0, 0, 0, 0.2e1 * t116, 0.2e1 * t115, t13, t9, t3, t4, t25, 0.2e1 * t84, 0.2e1 * t80; 0, 0, 0, 0, -t103, -t102, t63, t54, 0, 0, 0, t51, t50, t26, t7, 0, 0, 0, t85 + t116, t86 + t115, t13, t9, t3, t4, t25, t81 + t84, t79 + t80; 0, 0, 0, 0, 0, 0, t63, t54, 0, 0, 0, 0, 0, t26, t7, 0, 0, 0, 0.2e1 * t85, 0.2e1 * t86, t13, t9, t3, t4, t25, 0.2e1 * t81, 0.2e1 * t79; 0, 0, 0, 0, 0, 0, 0, 0, t65, -t105, 0 (-t72 * t65 - t98) * pkin(1) (t72 * t105 - t97) * pkin(1), 0, 0, t31, -t32, 0, -t18, t17, t12, t8, t23, t21, 0, t77 + t90, t78 + t117; 0, 0, 0, 0, 0, 0, 0, 0, t65, -t105, 0, 0, 0, 0, 0, t31, -t32, 0, 0, 0, t12, t8, t23, t21, 0, t77, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t100, -0.2e1 * t99, t62, t53, 0, 0, 0, 0.2e1 * t49, 0.2e1 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t32, 0, -t18, t17, t12, t8, t23, t21, 0, t90, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t32, 0, 0, 0, t12, t8, t23, t21, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t99, t62, t53, 0, 0, 0, t49, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t53, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t24, t32, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t24, t32 (t73 * t105 + t93) * pkin(2) (-t69 * t105 + t73 * t106) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t108, 0 (-t69 * t109 - t70 * t66) * pkin(2) (t95 - t96) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t108, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
