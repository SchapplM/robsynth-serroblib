% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:28
% EndTime: 2019-03-08 22:14:31
% DurationCPUTime: 1.12s
% Computational Cost: add. (916->163), mult. (2360->302), div. (0->0), fcn. (2207->10), ass. (0->111)
t70 = cos(qJ(6));
t63 = t70 ^ 2;
t66 = sin(qJ(6));
t111 = t66 ^ 2 - t63;
t129 = qJD(6) * t111;
t132 = 0.2e1 * t129;
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t131 = -t71 * pkin(3) - t68 * qJ(4);
t58 = t71 * qJD(3);
t67 = sin(qJ(5));
t124 = cos(qJ(5));
t94 = qJD(3) * t124;
t130 = -t67 * t58 + t68 * t94;
t41 = t124 * t71 + t67 * t68;
t125 = pkin(8) - pkin(9);
t47 = t125 * t68;
t48 = t125 * t71;
t24 = t124 * t48 + t67 * t47;
t13 = -t41 * qJD(3) * t125 + t24 * qJD(5);
t42 = t124 * t68 - t67 * t71;
t126 = -pkin(3) - pkin(4);
t89 = t124 * t126;
t43 = t67 * qJ(4) + pkin(5) - t89;
t77 = t124 * qJ(4) + t67 * t126;
t44 = -pkin(10) + t77;
t128 = -t13 + (t41 * t44 - t42 * t43) * qJD(6);
t127 = 2 * qJD(4);
t29 = t67 * qJD(4) + t77 * qJD(5);
t123 = t29 * t66;
t122 = t29 * t70;
t105 = t68 * qJD(3);
t19 = t41 * qJD(5) - t67 * t105 - t71 * t94;
t121 = t42 * t19;
t120 = t42 * t66;
t119 = t42 * t70;
t64 = sin(pkin(6));
t69 = sin(qJ(2));
t118 = t64 * t69;
t72 = cos(qJ(2));
t117 = t64 * t72;
t108 = qJD(5) * t67;
t93 = qJD(5) * t124;
t20 = -t71 * t108 + t68 * t93 - t130;
t116 = t66 * t20;
t114 = t70 * t19;
t113 = t70 * t20;
t112 = qJ(4) * t58 + t68 * qJD(4);
t110 = qJD(2) * t69;
t109 = qJD(2) * t72;
t107 = qJD(6) * t66;
t106 = qJD(6) * t70;
t104 = -0.2e1 * pkin(2) * qJD(3);
t103 = -0.2e1 * pkin(5) * qJD(6);
t45 = -pkin(2) + t131;
t102 = pkin(8) * t105;
t101 = pkin(8) * t58;
t100 = t64 * t110;
t99 = t64 * t109;
t98 = t66 * t106;
t95 = 0.4e1 * t66 * t119;
t92 = qJD(6) * (pkin(5) + t43);
t91 = qJD(6) * t124;
t37 = t71 * pkin(4) - t45;
t87 = pkin(5) * t19 - pkin(10) * t20;
t86 = pkin(5) * t42 + pkin(10) * t41;
t16 = t41 * pkin(5) - t42 * pkin(10) + t37;
t85 = t70 * t16 - t66 * t24;
t84 = t66 * t16 + t70 * t24;
t65 = cos(pkin(6));
t31 = t68 * t118 - t65 * t71;
t32 = t71 * t118 + t65 * t68;
t18 = t32 * t124 + t31 * t67;
t82 = t70 * t117 - t18 * t66;
t81 = t66 * t117 + t18 * t70;
t80 = t42 * t106 - t66 * t19;
t79 = -t42 * t107 - t114;
t17 = -t31 * t124 + t32 * t67;
t78 = t124 * t42 + t41 * t67;
t25 = t126 * t105 + t112;
t76 = t131 * qJD(3) + t71 * qJD(4);
t21 = -t31 * qJD(3) + t71 * t99;
t22 = t32 * qJD(3) + t68 * t99;
t75 = t21 * t71 + t22 * t68 + (t31 * t71 - t32 * t68) * qJD(3);
t23 = -t124 * t47 + t67 * t48;
t28 = qJ(4) * t108 - t124 * qJD(4) - qJD(5) * t89;
t74 = -qJD(6) * t23 - t19 * t43 - t20 * t44 + t28 * t41 + t29 * t42;
t73 = t124 * t19 - t20 * t67 + (-t124 * t41 + t42 * t67) * qJD(5);
t49 = 0.2e1 * t98;
t40 = -0.2e1 * t129;
t39 = t42 ^ 2;
t34 = t70 * t108 + t66 * t91;
t33 = t66 * t108 - t70 * t91;
t30 = pkin(3) * t105 - t112;
t27 = (-t72 * t105 - t71 * t110) * t64;
t26 = (t68 * t110 - t72 * t58) * t64;
t15 = t41 * t106 + t116;
t14 = t41 * t107 - t113;
t12 = t48 * t108 + t130 * t125 - t47 * t93;
t11 = t66 * t114 + t129 * t42;
t10 = t20 * pkin(5) + t19 * pkin(10) + t25;
t9 = qJD(6) * t95 - t111 * t19;
t8 = -t17 * qJD(5) + t21 * t124 + t22 * t67;
t7 = t18 * qJD(5) - t22 * t124 + t21 * t67;
t6 = -t17 * t107 + t7 * t70;
t5 = t17 * t106 + t7 * t66;
t4 = t82 * qJD(6) - t66 * t100 + t8 * t70;
t3 = -t81 * qJD(6) - t70 * t100 - t8 * t66;
t2 = -t84 * qJD(6) + t70 * t10 + t66 * t12;
t1 = -t85 * qJD(6) - t66 * t10 + t70 * t12;
t35 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t64 ^ 2 * t69 * t109 + 0.2e1 * t32 * t21 + 0.2e1 * t31 * t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t100, -t99, 0, 0, 0, 0, 0, t27, t26, t27, t75, -t26 (t45 * t110 - t30 * t72) * t64 + t75 * pkin(8), 0, 0, 0, 0, 0 (-t41 * t110 + t20 * t72) * t64 (-t110 * t42 - t19 * t72) * t64, 0, 0, 0, 0, 0, t7 * t120 + t17 * t80 + t20 * t82 + t3 * t41, t7 * t119 + t17 * t79 - t20 * t81 - t4 * t41; 0, 0, 0, 0, 0.2e1 * t68 * t58, 0.2e1 * (-t68 ^ 2 + t71 ^ 2) * qJD(3), 0, 0, 0, t68 * t104, t71 * t104, 0.2e1 * t45 * t105 - 0.2e1 * t30 * t71, 0, -0.2e1 * t30 * t68 - 0.2e1 * t45 * t58, 0.2e1 * t45 * t30, -0.2e1 * t121, 0.2e1 * t19 * t41 - 0.2e1 * t42 * t20, 0, 0, 0, 0.2e1 * t37 * t20 + 0.2e1 * t25 * t41, -0.2e1 * t37 * t19 + 0.2e1 * t25 * t42, -0.2e1 * t63 * t121 - 0.2e1 * t39 * t98, t132 * t39 + t19 * t95, 0.2e1 * t42 * t113 + 0.2e1 * t41 * t79, -0.2e1 * t42 * t116 - 0.2e1 * t41 * t80, 0.2e1 * t41 * t20, 0.2e1 * t13 * t120 + 0.2e1 * t2 * t41 + 0.2e1 * t20 * t85 + 0.2e1 * t23 * t80, 0.2e1 * t1 * t41 + 0.2e1 * t13 * t119 - 0.2e1 * t20 * t84 + 0.2e1 * t23 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t21, -t22, 0, t21, -t22 * pkin(3) + t21 * qJ(4) + t32 * qJD(4), 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, 0, 0, t58, -t105, 0, -t101, t102, -t101, t76, -t102, t76 * pkin(8), 0, 0, t19, t20, 0, t13, -t12, t11, t9, -t15, t14, 0, -t128 * t70 + t74 * t66, t128 * t66 + t74 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, qJ(4) * t127, 0, 0, 0, 0, 0, 0.2e1 * t29, -0.2e1 * t28, t49, t40, 0, 0, 0, -0.2e1 * t107 * t43 + 0.2e1 * t122, -0.2e1 * t106 * t43 - 0.2e1 * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106 * t78 + t66 * t73, t107 * t78 + t70 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t93, 0, 0, 0, 0, 0, t34, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, -t13, t12, -t11, -t9, t15, -t14, 0, -t13 * t70 + t87 * t66 + (t23 * t66 - t70 * t86) * qJD(6), t13 * t66 + t87 * t70 + (t23 * t70 + t66 * t86) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t28, -0.2e1 * t98, t132, 0, 0, 0, t66 * t92 - t122, t70 * t92 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t93, 0, 0, 0, 0, 0, -t34, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t40, 0, 0, 0, t66 * t103, t70 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t80, t20, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, t107, 0, -t106 * t44 + t66 * t28, t107 * t44 + t70 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106 * t67 - t66 * t93, t107 * t67 - t70 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, -t107, 0, -pkin(10) * t106, pkin(10) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t35;
