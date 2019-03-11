% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x22]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:34
% EndTime: 2019-03-08 21:26:37
% DurationCPUTime: 1.01s
% Computational Cost: add. (1466->176), mult. (3660->337), div. (0->0), fcn. (3547->10), ass. (0->106)
t58 = sin(pkin(11));
t60 = cos(pkin(11));
t62 = sin(qJ(3));
t65 = cos(qJ(3));
t44 = t58 * t62 - t60 * t65;
t45 = t58 * t65 + t60 * t62;
t54 = -t65 * pkin(3) - pkin(2);
t27 = t44 * pkin(4) - t45 * pkin(9) + t54;
t113 = -qJ(4) - pkin(8);
t49 = t113 * t65;
t89 = t113 * t62;
t32 = -t60 * t49 + t58 * t89;
t64 = cos(qJ(5));
t29 = t64 * t32;
t61 = sin(qJ(5));
t112 = t61 * t27 + t29;
t108 = cos(pkin(6));
t59 = sin(pkin(6));
t63 = sin(qJ(2));
t118 = t59 * t63;
t40 = t108 * t62 + t65 * t118;
t70 = t108 * t65 - t62 * t118;
t23 = t60 * t40 + t58 * t70;
t66 = cos(qJ(2));
t117 = t59 * t66;
t99 = t61 * t117;
t13 = t64 * t23 - t99;
t106 = qJD(2) * t66;
t91 = t59 * t106;
t30 = t70 * qJD(3) + t65 * t91;
t68 = -t40 * qJD(3) - t62 * t91;
t67 = t60 * t30 + t58 * t68;
t74 = t64 * t117 + t61 * t23;
t107 = qJD(2) * t63;
t92 = t59 * t107;
t5 = t74 * qJD(5) - t61 * t92 - t64 * t67;
t101 = qJD(5) * t64;
t6 = qJD(5) * t99 - t23 * t101 - t61 * t67 + t64 * t92;
t124 = qJD(5) * (-t13 * t64 - t61 * t74) + t5 * t61 - t6 * t64;
t110 = qJ(6) * t45;
t38 = t45 * qJD(3);
t104 = qJD(3) * t65;
t105 = qJD(3) * t62;
t39 = t60 * t104 - t58 * t105;
t75 = -qJ(6) * t39 - qJD(6) * t45;
t85 = qJD(3) * t113;
t37 = t65 * qJD(4) + t62 * t85;
t69 = -t62 * qJD(4) + t65 * t85;
t19 = t60 * t37 + t58 * t69;
t55 = pkin(3) * t105;
t20 = t38 * pkin(4) - t39 * pkin(9) + t55;
t87 = -t61 * t19 + t64 * t20;
t1 = t38 * pkin(5) + t75 * t64 + (-t29 + (-t27 + t110) * t61) * qJD(5) + t87;
t94 = t45 * t101;
t97 = t27 * t101 + t64 * t19 + t61 * t20;
t2 = -qJ(6) * t94 + (-qJD(5) * t32 + t75) * t61 + t97;
t86 = t64 * t27 - t61 * t32;
t7 = t44 * pkin(5) - t64 * t110 + t86;
t8 = -t61 * t110 + t112;
t123 = -t1 * t64 - t2 * t61 + (t61 * t7 - t64 * t8) * qJD(5);
t122 = 0.2e1 * qJD(5);
t11 = t58 * t30 - t60 * t68;
t22 = t58 * t40 - t60 * t70;
t121 = t22 * t11;
t120 = t45 * t61;
t119 = t45 * t64;
t116 = t61 * t38;
t115 = t64 * t38;
t114 = t64 * t39;
t56 = t61 ^ 2;
t57 = t64 ^ 2;
t111 = t56 - t57;
t52 = t58 * pkin(3) + pkin(9);
t109 = qJ(6) + t52;
t103 = qJD(3) * t66;
t102 = qJD(5) * t61;
t100 = -0.2e1 * pkin(2) * qJD(3);
t53 = -t60 * pkin(3) - pkin(4);
t98 = t53 * t122;
t96 = pkin(5) * t102;
t95 = t62 * t103;
t93 = t22 * t102;
t90 = t61 * t101;
t88 = -0.4e1 * t61 * t119;
t18 = t58 * t37 - t60 * t69;
t31 = -t58 * t49 - t60 * t89;
t84 = t111 * qJD(5);
t83 = qJD(5) * t109;
t79 = -t13 * t61 + t64 * t74;
t77 = -t38 * t52 + t39 * t53;
t76 = t44 * t52 - t45 * t53;
t73 = t44 * t101 + t116;
t72 = t61 * t39 + t94;
t71 = t45 * t102 - t114;
t47 = -t64 * pkin(5) + t53;
t43 = t45 ^ 2;
t42 = t109 * t64;
t41 = t109 * t61;
t34 = -t61 * qJD(6) - t64 * t83;
t33 = t64 * qJD(6) - t61 * t83;
t26 = -t44 * t102 + t115;
t17 = pkin(5) * t120 + t31;
t9 = t72 * pkin(5) + t18;
t4 = -t112 * qJD(5) + t87;
t3 = t32 * t102 - t97;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t59 ^ 2 * t63 * t106 + 0.2e1 * t23 * t67 + 0.2e1 * t121, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t13 * t5 - 0.2e1 * t6 * t74 + 0.2e1 * t121; 0, 0, -t92, -t91, 0, 0, 0, 0, 0 (-t65 * t107 - t95) * t59 (-t65 * t103 + t62 * t107) * t59, t11 * t45 + t22 * t39 - t23 * t38 - t67 * t44, -t59 * pkin(3) * t95 + t11 * t31 + t22 * t18 + t23 * t19 + t67 * t32 + t54 * t92, 0, 0, 0, 0, 0, t11 * t120 + t72 * t22 - t38 * t74 + t6 * t44, t11 * t119 - t13 * t38 - t22 * t71 + t5 * t44, t124 * t45 + t79 * t39, -t1 * t74 + t11 * t17 + t13 * t2 + t22 * t9 - t5 * t8 + t6 * t7; 0, 0, 0, 0, 0.2e1 * t62 * t104, 0.2e1 * (-t62 ^ 2 + t65 ^ 2) * qJD(3), 0, 0, 0, t62 * t100, t65 * t100, 0.2e1 * t18 * t45 - 0.2e1 * t19 * t44 + 0.2e1 * t31 * t39 - 0.2e1 * t32 * t38, 0.2e1 * t31 * t18 + 0.2e1 * t32 * t19 + 0.2e1 * t54 * t55, 0.2e1 * t57 * t45 * t39 - 0.2e1 * t43 * t90, t111 * t43 * t122 + t39 * t88, 0.2e1 * t45 * t115 - 0.2e1 * t71 * t44, -0.2e1 * t45 * t116 - 0.2e1 * t72 * t44, 0.2e1 * t44 * t38, 0.2e1 * t18 * t120 + 0.2e1 * t72 * t31 + 0.2e1 * t86 * t38 + 0.2e1 * t4 * t44, -0.2e1 * t112 * t38 + 0.2e1 * t18 * t119 + 0.2e1 * t3 * t44 - 0.2e1 * t31 * t71, 0.2e1 * (-t61 * t8 - t64 * t7) * t39 + 0.2e1 * t123 * t45, 0.2e1 * t7 * t1 + 0.2e1 * t17 * t9 + 0.2e1 * t8 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t30, 0 (-t11 * t60 + t67 * t58) * pkin(3), 0, 0, 0, 0, 0, -t11 * t64 + t93, t101 * t22 + t11 * t61, qJD(5) * t79 - t5 * t64 - t6 * t61, pkin(5) * t93 + t11 * t47 + t13 * t33 - t34 * t74 - t6 * t41 - t5 * t42; 0, 0, 0, 0, 0, 0, t104, -t105, 0, -pkin(8) * t104, pkin(8) * t105 (-t38 * t58 - t39 * t60) * pkin(3) (-t18 * t60 + t19 * t58) * pkin(3), t61 * t114 - t45 * t84, qJD(5) * t88 - t111 * t39, t73, t26, 0, -t18 * t64 + t77 * t61 + (t31 * t61 - t64 * t76) * qJD(5), t18 * t61 + t77 * t64 + (t31 * t64 + t61 * t76) * qJD(5) (-t34 * t45 + t39 * t41 + t2 + (-t42 * t45 - t7) * qJD(5)) * t64 + (-t33 * t45 - t39 * t42 - t1 + (-t41 * t45 - t8) * qJD(5)) * t61, -t1 * t41 + t17 * t96 + t2 * t42 + t8 * t33 + t7 * t34 + t9 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t90, -0.2e1 * t84, 0, 0, 0, t61 * t98, t64 * t98, 0.2e1 * t33 * t64 - 0.2e1 * t34 * t61 + 0.2e1 * (t41 * t64 - t42 * t61) * qJD(5), 0.2e1 * t42 * t33 - 0.2e1 * t41 * t34 + 0.2e1 * t47 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, 0, 0, 0, 0, 0, 0, 0, -t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, t26, -t73 (-t56 - t57) * t39, -t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 * t61 + t34 * t64 + (t41 * t61 + t42 * t64) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, t6 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t72, t38, t4, t3, t71 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t102, 0, -t52 * t101, t52 * t102, -pkin(5) * t101, t34 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t101, 0, -t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
