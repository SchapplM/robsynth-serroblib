% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:34
% EndTime: 2019-03-09 05:13:36
% DurationCPUTime: 1.00s
% Computational Cost: add. (2213->154), mult. (4925->252), div. (0->0), fcn. (5070->8), ass. (0->100)
t120 = sin(qJ(3));
t122 = cos(qJ(3));
t65 = sin(pkin(10));
t66 = cos(pkin(10));
t82 = -t120 * t66 - t122 * t65;
t126 = t82 * qJD(2);
t67 = sin(qJ(6));
t63 = t67 ^ 2;
t69 = cos(qJ(6));
t113 = -t69 ^ 2 + t63;
t98 = t113 * qJD(6);
t83 = t120 * t65 - t122 * t66;
t117 = pkin(7) + qJ(2);
t47 = t117 * t65;
t48 = t117 * t66;
t84 = t120 * t48 + t122 * t47;
t125 = t83 * qJD(2) + t84 * qJD(3);
t71 = 2 * qJD(5);
t124 = pkin(4) + pkin(9);
t100 = qJD(3) * t120;
t101 = qJD(3) * t122;
t115 = t66 * t100 + t65 * t101;
t121 = cos(qJ(4));
t68 = sin(qJ(4));
t77 = t83 * qJD(3);
t80 = t68 * t83;
t99 = qJD(4) * t121;
t24 = -qJD(4) * t80 + t121 * t115 - t68 * t77 - t82 * t99;
t112 = qJD(4) * t68;
t28 = -t115 * pkin(8) - t125;
t31 = pkin(8) * t82 - t84;
t85 = t120 * t47 - t122 * t48;
t32 = -t83 * pkin(8) - t85;
t72 = pkin(8) * t77 + t47 * t100 - t48 * t101 + t126;
t8 = t32 * t112 - t121 * t28 - t31 * t99 - t68 * t72;
t5 = -pkin(5) * t24 - t8;
t3 = t5 * t67;
t4 = t5 * t69;
t109 = qJD(6) * t69;
t75 = t121 * t83;
t34 = -t68 * t82 + t75;
t88 = t121 * t32 + t68 * t31;
t16 = -t34 * pkin(5) + t88;
t123 = t16 * t109 + t3;
t119 = t24 * t67;
t118 = t24 * t69;
t97 = pkin(3) * t99;
t50 = t97 + qJD(5);
t54 = pkin(3) * t68 + qJ(5);
t116 = t54 * t109 + t50 * t67;
t107 = qJ(5) * qJD(6);
t114 = qJD(5) * t67 + t69 * t107;
t111 = qJD(6) * t16;
t110 = qJD(6) * t67;
t108 = qJD(6) * t124;
t106 = t67 * t118;
t58 = pkin(3) * t112;
t105 = t67 * t109;
t55 = -t66 * pkin(2) - pkin(1);
t104 = t115 * pkin(3);
t96 = 0.2e1 * (t65 ^ 2 + t66 ^ 2) * qJD(2);
t57 = -t121 * pkin(3) - pkin(4);
t21 = -t121 * t31 + t68 * t32;
t35 = -t121 * t82 - t80;
t36 = t83 * pkin(3) + t55;
t73 = -t35 * qJ(5) + t36;
t14 = t124 * t34 + t73;
t15 = t35 * pkin(5) + t21;
t95 = t14 * t69 + t15 * t67;
t94 = t14 * t67 - t15 * t69;
t93 = -qJ(5) * t24 - qJD(5) * t34;
t23 = -t82 * t112 + t68 * t115 - (-qJD(3) - qJD(4)) * t75;
t92 = t35 * t109 - t23 * t67;
t91 = t35 * t110 + t23 * t69;
t90 = t34 * t109 + t119;
t89 = t34 * t110 - t118;
t53 = -pkin(9) + t57;
t87 = qJD(6) * (t34 * t54 - t35 * t53);
t86 = qJD(6) * (qJ(5) * t34 + t124 * t35);
t81 = t124 * t23 + t93;
t79 = t23 * qJ(5) - t35 * qJD(5) + t104;
t78 = -t24 * t54 - t34 * t50 + t35 * t58;
t76 = -0.2e1 * t77;
t74 = -t23 * t53 + t78;
t9 = t88 * qJD(4) - t121 * t72 + t68 * t28;
t60 = qJD(5) * t69;
t49 = -0.2e1 * t105;
t45 = t50 * t69;
t43 = 0.2e1 * t98;
t33 = t34 ^ 2;
t20 = t34 * pkin(4) + t73;
t19 = -0.2e1 * t35 * t23;
t12 = -t34 * t98 + t106;
t11 = -0.4e1 * t34 * t105 - t113 * t24;
t10 = t24 * pkin(4) + t79;
t7 = t124 * t24 + t79;
t6 = -t23 * pkin(5) + t9;
t2 = -t95 * qJD(6) + t69 * t6 - t67 * t7;
t1 = t94 * qJD(6) - t67 * t6 - t69 * t7;
t13 = [0, 0, 0, 0, 0, t96, qJ(2) * t96, -t82 * t76, 0.2e1 * t83 ^ 2 * qJD(3) + 0.2e1 * t115 * t82, 0, 0, 0, 0.2e1 * t55 * t115, t55 * t76, t19, 0.2e1 * t23 * t34 - 0.2e1 * t24 * t35, 0, 0, 0, 0.2e1 * t34 * t104 + 0.2e1 * t36 * t24, 0.2e1 * t35 * t104 - 0.2e1 * t36 * t23, -0.2e1 * t21 * t23 - 0.2e1 * t24 * t88 + 0.2e1 * t34 * t8 + 0.2e1 * t35 * t9, -0.2e1 * t10 * t34 - 0.2e1 * t20 * t24, -0.2e1 * t10 * t35 + 0.2e1 * t20 * t23, 0.2e1 * t10 * t20 + 0.2e1 * t21 * t9 - 0.2e1 * t8 * t88, 0.2e1 * t24 * t34 * t63 + 0.2e1 * t33 * t105, 0.4e1 * t34 * t106 - 0.2e1 * t33 * t98, 0.2e1 * t35 * t119 + 0.2e1 * t92 * t34, 0.2e1 * t35 * t118 - 0.2e1 * t91 * t34, t19, 0.2e1 * t89 * t16 + 0.2e1 * t2 * t35 + 0.2e1 * t94 * t23 - 0.2e1 * t34 * t4, 0.2e1 * t1 * t35 + 0.2e1 * t90 * t16 + 0.2e1 * t95 * t23 + 0.2e1 * t34 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t77, 0, 0, 0, 0, 0, t24, -t23, 0, -t24, t23, t10, 0, 0, 0, 0, 0, -t92, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t115, 0, t85 * qJD(3) + t126, t125, 0, 0, -t23, -t24, 0, -t9, t8, -t23 * t57 + t78, t9, -t8, t21 * t58 + t50 * t88 - t54 * t8 + t57 * t9, t12, t11, -t91, -t92, 0, t67 * t87 + t74 * t69 + t123, t4 + t69 * t87 + (-t74 - t111) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t58, -0.2e1 * t97, 0, 0.2e1 * t58, 0.2e1 * t50, 0.2e1 * t50 * t54 + 0.2e1 * t57 * t58, t49, t43, 0, 0, 0, 0.2e1 * t116, -0.2e1 * t54 * t110 + 0.2e1 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, -t9, t8, pkin(4) * t23 + t93, t9, -t8, -pkin(4) * t9 - qJ(5) * t8 + qJD(5) * t88, t12, t11, -t91, -t92, 0, t67 * t86 + t81 * t69 + t123, t4 + t69 * t86 + (-t81 - t111) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t97, 0, t58, t71 + t97, -pkin(4) * t58 + qJ(5) * t50 + qJD(5) * t54, t49, t43, 0, 0, 0, t114 + t116, t45 + t60 + (-qJ(5) - t54) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, qJ(5) * t71, t49, t43, 0, 0, 0, 0.2e1 * t114, -0.2e1 * t67 * t107 + 0.2e1 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, t9, 0, 0, 0, 0, 0, -t91, -t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t89, -t23, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t109, 0, -t53 * t110 + t69 * t58, -t53 * t109 - t67 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t109, 0, t67 * t108, t69 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t13;
