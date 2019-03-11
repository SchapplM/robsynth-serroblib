% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:35
% EndTime: 2019-03-09 03:49:38
% DurationCPUTime: 1.06s
% Computational Cost: add. (1597->147), mult. (3596->265), div. (0->0), fcn. (3668->8), ass. (0->97)
t74 = cos(qJ(6));
t69 = t74 ^ 2;
t72 = sin(qJ(6));
t109 = t72 ^ 2 - t69;
t126 = qJD(6) * t109;
t129 = 0.2e1 * t126;
t119 = sin(qJ(3));
t121 = cos(qJ(3));
t71 = cos(pkin(10));
t60 = t121 * t71;
t70 = sin(pkin(10));
t93 = t119 * t70 - t60;
t61 = -pkin(2) * t71 - pkin(1);
t128 = -0.2e1 * t61;
t82 = -t119 * t71 - t121 * t70;
t127 = t82 * qJD(3);
t73 = sin(qJ(5));
t120 = cos(qJ(5));
t85 = t120 * t93;
t27 = -t73 * t82 - t85;
t28 = -t120 * t82 + t73 * t93;
t101 = qJD(5) * t120;
t108 = qJD(5) * t73;
t111 = pkin(7) + qJ(2);
t53 = t111 * t70;
t54 = t111 * t71;
t29 = t119 * t54 + t121 * t53;
t23 = pkin(8) * t82 + t29;
t110 = -t119 * t53 + t121 * t54;
t24 = pkin(8) * t93 + t110;
t102 = qJD(3) * t119;
t58 = t70 * t102;
t81 = t82 * qJD(2);
t75 = t58 * pkin(8) + (-pkin(8) * t60 + t110) * qJD(3) - t81;
t103 = qJD(3) * t121;
t21 = qJD(2) * t93 + t54 * t102 + t53 * t103;
t76 = -pkin(8) * t127 - t21;
t4 = t101 * t24 + t108 * t23 - t120 * t75 + t73 * t76;
t123 = -pkin(3) - pkin(4);
t96 = t120 * t123;
t51 = qJ(4) * t73 + pkin(5) - t96;
t84 = qJ(4) * t120 + t123 * t73;
t52 = -pkin(9) + t84;
t125 = -t4 + (t27 * t52 - t28 * t51) * qJD(6);
t124 = 2 * qJD(4);
t122 = t4 * t72;
t89 = -t103 * t71 + t58;
t13 = -qJD(5) * t85 - t108 * t82 + t120 * t89 + t127 * t73;
t118 = t28 * t13;
t117 = t28 * t74;
t35 = t73 * qJD(4) + qJD(5) * t84;
t116 = t35 * t72;
t115 = t35 * t74;
t14 = qJD(5) * t28 + t120 * t127 - t73 * t89;
t114 = t72 * t14;
t113 = t74 * t13;
t112 = t74 * t14;
t63 = qJD(6) * t72;
t64 = qJD(6) * t74;
t107 = -0.2e1 * pkin(5) * qJD(6);
t18 = -pkin(3) * t127 + qJ(4) * t89 + qJD(4) * t82;
t106 = t72 * t64;
t104 = 0.4e1 * t72 * t117;
t100 = qJD(6) * (pkin(5) + t51);
t99 = qJD(6) * t120;
t25 = pkin(3) * t93 + qJ(4) * t82 + t61;
t97 = 0.2e1 * (t70 ^ 2 + t71 ^ 2) * qJD(2);
t95 = pkin(5) * t13 - pkin(9) * t14;
t94 = pkin(5) * t28 + pkin(9) * t27;
t12 = t120 * t24 + t23 * t73;
t19 = -pkin(4) * t93 - t25;
t8 = pkin(5) * t27 - pkin(9) * t28 + t19;
t92 = t12 * t74 + t72 * t8;
t91 = t12 * t72 - t74 * t8;
t88 = -t13 * t72 + t28 * t64;
t87 = -t28 * t63 - t113;
t10 = t27 * t64 + t114;
t86 = t120 * t28 + t27 * t73;
t11 = -t120 * t23 + t24 * t73;
t34 = qJ(4) * t108 - qJD(4) * t120 - qJD(5) * t96;
t78 = -qJD(6) * t11 - t13 * t51 - t14 * t52 + t27 * t34 + t28 * t35;
t15 = pkin(4) * t127 - t18;
t77 = t120 * t13 - t14 * t73 + (-t120 * t27 + t28 * t73) * qJD(5);
t56 = 0.2e1 * t106;
t50 = -0.2e1 * t126;
t41 = t108 * t74 + t72 * t99;
t40 = t108 * t72 - t74 * t99;
t26 = t28 ^ 2;
t22 = qJD(3) * t110 - t81;
t9 = t27 * t63 - t112;
t7 = t113 * t72 + t126 * t28;
t6 = qJD(6) * t104 - t109 * t13;
t5 = t14 * pkin(5) + t13 * pkin(9) + t15;
t3 = -t101 * t23 + t108 * t24 - t120 * t76 - t73 * t75;
t2 = -qJD(6) * t92 + t72 * t3 + t74 * t5;
t1 = qJD(6) * t91 + t74 * t3 - t72 * t5;
t16 = [0, 0, 0, 0, 0, t97, qJ(2) * t97, 0.2e1 * t82 * t89, -0.2e1 * t127 * t82 + 0.2e1 * t89 * t93, 0, 0, 0, t127 * t128, t89 * t128, -0.2e1 * t127 * t25 + 0.2e1 * t18 * t93, 0.2e1 * t110 * t127 + 0.2e1 * t21 * t93 - 0.2e1 * t22 * t82 - 0.2e1 * t29 * t89, 0.2e1 * t18 * t82 + 0.2e1 * t25 * t89, -0.2e1 * t110 * t21 + 0.2e1 * t18 * t25 + 0.2e1 * t22 * t29, -0.2e1 * t118, 0.2e1 * t13 * t27 - 0.2e1 * t14 * t28, 0, 0, 0, 0.2e1 * t14 * t19 + 0.2e1 * t15 * t27, -0.2e1 * t13 * t19 + 0.2e1 * t15 * t28, -0.2e1 * t106 * t26 - 0.2e1 * t118 * t69, t104 * t13 + t129 * t26, 0.2e1 * t112 * t28 + 0.2e1 * t27 * t87, -0.2e1 * t114 * t28 - 0.2e1 * t27 * t88, 0.2e1 * t27 * t14, 0.2e1 * t11 * t88 + 0.2e1 * t122 * t28 - 0.2e1 * t14 * t91 + 0.2e1 * t2 * t27, 0.2e1 * t1 * t27 + 0.2e1 * t11 * t87 + 0.2e1 * t117 * t4 - 0.2e1 * t14 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t89, -t127, 0, t89, t18, 0, 0, 0, 0, 0, -t14, t13, 0, 0, 0, 0, 0, t9, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t127, 0, -t22, t21, -t22, pkin(3) * t89 + qJ(4) * t127 - qJD(4) * t93, -t21, -pkin(3) * t22 - qJ(4) * t21 + qJD(4) * t110, 0, 0, t13, t14, 0, t4, -t3, t7, t6, -t10, t9, 0, -t125 * t74 + t72 * t78, t125 * t72 + t74 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, qJ(4) * t124, 0, 0, 0, 0, 0, 0.2e1 * t35, -0.2e1 * t34, t56, t50, 0, 0, 0, -0.2e1 * t51 * t63 + 0.2e1 * t115, -0.2e1 * t51 * t64 - 0.2e1 * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t86 + t72 * t77, t63 * t86 + t74 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t101, 0, 0, 0, 0, 0, t41, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, -t4, t3, -t7, -t6, t10, -t9, 0, -t4 * t74 + t95 * t72 + (t11 * t72 - t74 * t94) * qJD(6), t122 + t95 * t74 + (t11 * t74 + t72 * t94) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t34, -0.2e1 * t106, t129, 0, 0, 0, t100 * t72 - t115, t100 * t74 + t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t101, 0, 0, 0, 0, 0, -t41, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t50, 0, 0, 0, t72 * t107, t74 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t88, t14, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t63, 0, t34 * t72 - t52 * t64, t34 * t74 + t52 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * t72 - t64 * t73, -t101 * t74 + t63 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t63, 0, -pkin(9) * t64, pkin(9) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t16;
