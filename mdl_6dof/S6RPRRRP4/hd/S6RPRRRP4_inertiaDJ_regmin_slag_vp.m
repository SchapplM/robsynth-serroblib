% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:45
% EndTime: 2019-03-09 06:08:48
% DurationCPUTime: 1.33s
% Computational Cost: add. (3601->199), mult. (7807->328), div. (0->0), fcn. (8119->8), ass. (0->117)
t80 = sin(pkin(10));
t81 = cos(pkin(10));
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t56 = t87 * t80 + t84 * t81;
t129 = pkin(7) + qJ(2);
t61 = t129 * t80;
t62 = t129 * t81;
t99 = -t87 * t61 - t84 * t62;
t37 = -t56 * pkin(8) + t99;
t97 = t84 * t80 - t87 * t81;
t98 = t84 * t61 - t87 * t62;
t38 = -pkin(8) * t97 - t98;
t83 = sin(qJ(4));
t86 = cos(qJ(4));
t28 = t83 * t37 + t86 * t38;
t85 = cos(qJ(5));
t24 = t85 * t28;
t100 = -t83 * t56 - t86 * t97;
t41 = t86 * t56 - t83 * t97;
t68 = -t81 * pkin(2) - pkin(1);
t45 = pkin(3) * t97 + t68;
t29 = -pkin(4) * t100 - t41 * pkin(9) + t45;
t82 = sin(qJ(5));
t127 = t82 * t29 + t24;
t78 = t82 ^ 2;
t79 = t85 ^ 2;
t125 = t78 - t79;
t107 = t125 * qJD(5);
t121 = qJD(4) * t86;
t122 = qJD(4) * t83;
t140 = t97 * qJD(2) - t99 * qJD(3);
t52 = t56 * qJD(3);
t88 = -t52 * pkin(8) - t140;
t51 = t97 * qJD(3);
t90 = -t56 * qJD(2) + t98 * qJD(3);
t89 = t51 * pkin(8) + t90;
t11 = -t37 * t121 + t38 * t122 - t83 * t89 - t86 * t88;
t138 = t52 * pkin(3);
t30 = t100 * qJD(4) - t86 * t51 - t83 * t52;
t31 = t41 * qJD(4) - t83 * t51 + t86 * t52;
t16 = t31 * pkin(4) - t30 * pkin(9) + t138;
t110 = t82 * t11 + t85 * t16;
t124 = qJ(6) * t41;
t96 = -qJ(6) * t30 - qJD(6) * t41;
t1 = t31 * pkin(5) + t96 * t85 + (-t24 + (-t29 + t124) * t82) * qJD(5) + t110;
t10 = -t82 * t124 + t127;
t74 = qJD(5) * t85;
t113 = t41 * t74;
t119 = -t85 * t11 + t82 * t16 + t29 * t74;
t3 = -qJ(6) * t113 + (-qJD(5) * t28 + t96) * t82 + t119;
t109 = -t82 * t28 + t85 * t29;
t75 = t85 * qJ(6);
t7 = -pkin(5) * t100 - t41 * t75 + t109;
t141 = t1 * t85 + t3 * t82 + (t10 * t85 - t7 * t82) * qJD(5);
t139 = 0.2e1 * t68;
t137 = t86 * pkin(3);
t12 = t28 * qJD(4) + t83 * t88 - t86 * t89;
t27 = -t86 * t37 + t83 * t38;
t136 = t12 * t82 + t27 * t74;
t135 = t41 * t30;
t134 = t41 * t82;
t133 = t41 * t85;
t132 = t82 * t31;
t131 = t85 * t30;
t130 = t85 * t31;
t128 = -qJ(6) - pkin(9);
t116 = pkin(3) * t122;
t70 = -pkin(4) - t137;
t126 = t82 * t116 + t70 * t74;
t69 = t83 * pkin(3) + pkin(9);
t123 = -qJ(6) - t69;
t120 = qJD(5) * t82;
t118 = pkin(4) * t120;
t117 = pkin(4) * t74;
t115 = pkin(3) * t121;
t72 = pkin(5) * t120;
t114 = pkin(5) * t74;
t112 = t82 * t74;
t71 = -t85 * pkin(5) - pkin(4);
t111 = -0.4e1 * t82 * t133;
t108 = qJD(5) * t128;
t106 = qJD(5) * t123;
t105 = t85 * t115;
t104 = 0.2e1 * (t80 ^ 2 + t81 ^ 2) * qJD(2);
t101 = -t100 * t69 - t41 * t70;
t95 = -t85 * t116 + t70 * t120;
t94 = t82 * t30 + t113;
t93 = t41 * t120 - t131;
t20 = -t100 * t74 + t132;
t92 = -t100 * t120 - t130;
t91 = t30 * t70 - t31 * t69 + (t100 * t86 + t41 * t83) * qJD(4) * pkin(3);
t73 = t85 * qJD(6);
t67 = 0.2e1 * t112;
t65 = t85 * pkin(9) + t75;
t64 = t128 * t82;
t63 = t71 - t137;
t58 = t72 + t116;
t57 = -0.2e1 * t107;
t54 = t85 * t69 + t75;
t53 = t123 * t82;
t49 = -t82 * qJD(6) + t85 * t108;
t48 = t82 * t108 + t73;
t46 = t48 * t85;
t44 = (-qJD(6) - t115) * t82 + t85 * t106;
t43 = t82 * t106 + t105 + t73;
t42 = t43 * t85;
t39 = t41 ^ 2;
t21 = t27 * t120;
t18 = pkin(5) * t134 + t27;
t17 = -t41 * t107 + t82 * t131;
t13 = qJD(5) * t111 - t125 * t30;
t6 = pkin(5) * t94 + t12;
t5 = -t127 * qJD(5) + t110;
t4 = t28 * t120 - t119;
t2 = t3 * t85;
t8 = [0, 0, 0, 0, 0, t104, qJ(2) * t104, -0.2e1 * t56 * t51, 0.2e1 * t51 * t97 - 0.2e1 * t52 * t56, 0, 0, 0, t52 * t139, -t51 * t139, 0.2e1 * t135, 0.2e1 * t100 * t30 - 0.2e1 * t41 * t31, 0, 0, 0, -0.2e1 * t100 * t138 + 0.2e1 * t45 * t31, 0.2e1 * t41 * t138 + 0.2e1 * t45 * t30, -0.2e1 * t39 * t112 + 0.2e1 * t79 * t135, 0.2e1 * t39 * t107 + t30 * t111, 0.2e1 * t100 * t93 + 0.2e1 * t41 * t130, 0.2e1 * t100 * t94 - 0.2e1 * t41 * t132, -0.2e1 * t100 * t31, -0.2e1 * t100 * t5 + 0.2e1 * t109 * t31 + 0.2e1 * t12 * t134 + 0.2e1 * t94 * t27, -0.2e1 * t100 * t4 + 0.2e1 * t12 * t133 - 0.2e1 * t127 * t31 - 0.2e1 * t93 * t27, 0.2e1 * (-t10 * t82 - t7 * t85) * t30 - 0.2e1 * t141 * t41, 0.2e1 * t7 * t1 + 0.2e1 * t10 * t3 + 0.2e1 * t18 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t51, 0, 0, 0, 0, 0, t31, t30, 0, 0, 0, 0, 0, -t92, -t20 (-t78 - t79) * t30, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t52, 0, t90, t140, 0, 0, t30, -t31, 0, -t12, t11, t17, t13, t20, -t92, 0, t21 + (-t101 * qJD(5) - t12) * t85 + t91 * t82, t101 * t120 + t91 * t85 + t136, t2 + (-t30 * t53 - t41 * t44 + (-t41 * t54 - t7) * qJD(5)) * t85 + (-t30 * t54 - t41 * t43 - t1 + (t41 * t53 - t10) * qJD(5)) * t82, t1 * t53 + t10 * t43 + t18 * t58 + t3 * t54 + t7 * t44 + t6 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * t43 + t85 * t44 + (-t53 * t82 + t54 * t85) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t116, -0.2e1 * t115, t67, t57, 0, 0, 0, 0.2e1 * t95, 0.2e1 * t126, -0.2e1 * t44 * t82 + 0.2e1 * t42 + 0.2e1 * (-t53 * t85 - t54 * t82) * qJD(5), 0.2e1 * t54 * t43 + 0.2e1 * t53 * t44 + 0.2e1 * t63 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t31, 0, -t12, t11, t17, t13, t20, -t92, 0, t21 + (-pkin(4) * t30 - pkin(9) * t31) * t82 + (-t12 + (-pkin(4) * t41 + pkin(9) * t100) * qJD(5)) * t85, pkin(4) * t93 + pkin(9) * t92 + t136, t2 + (-t30 * t64 - t41 * t49 + (-t41 * t65 - t7) * qJD(5)) * t85 + (-t30 * t65 - t41 * t48 - t1 + (t41 * t64 - t10) * qJD(5)) * t82, t1 * t64 + t10 * t48 + t18 * t72 + t3 * t65 + t7 * t49 + t6 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * t48 + t85 * t49 + (-t64 * t82 + t65 * t85) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t115, t67, t57, 0, 0, 0, t95 - t118, -t117 + t126, t42 + t46 + (-t44 - t49) * t82 + ((-t53 - t64) * t85 + (-t54 - t65) * t82) * qJD(5), t43 * t65 + t44 * t64 + t54 * t48 + t53 * t49 + t58 * t71 + t63 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t57, 0, 0, 0, -0.2e1 * t118, -0.2e1 * t117, -0.2e1 * t49 * t82 + 0.2e1 * t46 + 0.2e1 * (-t64 * t85 - t65 * t82) * qJD(5), 0.2e1 * t65 * t48 + 0.2e1 * t64 * t49 + 0.2e1 * t71 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t94, t31, t5, t4, t93 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t74, 0, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t120, 0, -t82 * t115 - t69 * t74, t69 * t120 - t105, -t114, t44 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t120, 0, -pkin(9) * t74, pkin(9) * t120, -t114, t49 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
