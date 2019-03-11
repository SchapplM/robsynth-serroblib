% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:30
% EndTime: 2019-03-08 19:01:33
% DurationCPUTime: 0.98s
% Computational Cost: add. (1295->159), mult. (3765->295), div. (0->0), fcn. (4025->14), ass. (0->117)
t86 = cos(qJ(6));
t75 = t86 ^ 2;
t82 = sin(qJ(6));
t126 = t82 ^ 2 - t75;
t107 = t126 * qJD(6);
t139 = qJD(4) + qJD(5);
t138 = pkin(9) + pkin(10);
t83 = sin(qJ(5));
t84 = sin(qJ(4));
t87 = cos(qJ(5));
t88 = cos(qJ(4));
t60 = t83 * t84 - t87 * t88;
t44 = t139 * t60;
t61 = t83 * t88 + t87 * t84;
t137 = t61 * t44;
t136 = t61 * t82;
t135 = t61 * t86;
t77 = sin(pkin(7));
t85 = sin(qJ(3));
t134 = t77 * t85;
t89 = cos(qJ(3));
t133 = t77 * t89;
t79 = cos(pkin(13));
t80 = cos(pkin(7));
t132 = t79 * t80;
t45 = t139 * t61;
t131 = t82 * t45;
t130 = t86 * t44;
t129 = t86 * t45;
t109 = qJD(4) * t138;
t105 = t88 * t109;
t106 = t84 * t109;
t65 = t138 * t84;
t66 = t138 * t88;
t48 = -t83 * t65 + t87 * t66;
t27 = qJD(5) * t48 + t105 * t87 - t106 * t83;
t47 = t87 * t65 + t83 * t66;
t73 = qJD(6) * t86;
t128 = t27 * t82 + t47 * t73;
t124 = qJD(5) * t83;
t114 = pkin(4) * t124;
t71 = -t87 * pkin(4) - pkin(5);
t127 = t82 * t114 + t71 * t73;
t125 = qJD(3) * t85;
t123 = qJD(5) * t87;
t122 = qJD(6) * t82;
t121 = t84 * qJD(4);
t120 = t88 * qJD(4);
t119 = -0.2e1 * pkin(3) * qJD(4);
t118 = t89 * t132;
t117 = pkin(5) * t122;
t116 = pkin(5) * t73;
t115 = pkin(4) * t121;
t113 = pkin(4) * t123;
t112 = t77 * t125;
t111 = qJD(3) * t133;
t110 = t82 * t73;
t72 = -t88 * pkin(4) - pkin(3);
t108 = -0.4e1 * t82 * t135;
t76 = sin(pkin(13));
t78 = sin(pkin(6));
t81 = cos(pkin(6));
t41 = t81 * t134 + (t85 * t132 + t76 * t89) * t78;
t50 = -t78 * t79 * t77 + t81 * t80;
t28 = -t41 * t84 + t50 * t88;
t29 = t41 * t88 + t50 * t84;
t14 = t83 * t28 + t87 * t29;
t40 = -t81 * t133 + (t76 * t85 - t118) * t78;
t104 = t86 * t14 + t40 * t82;
t103 = t82 * t14 - t40 * t86;
t36 = t60 * pkin(5) - t61 * pkin(11) + t72;
t102 = t86 * t36 - t82 * t48;
t101 = t82 * t36 + t86 * t48;
t52 = t88 * t134 + t84 * t80;
t96 = t84 * t134 - t88 * t80;
t33 = t87 * t52 - t83 * t96;
t70 = t83 * pkin(4) + pkin(11);
t100 = t60 * t70 - t61 * t71;
t99 = -t86 * t114 + t71 * t122;
t98 = t86 * t133 + t82 * t33;
t97 = t82 * t133 - t86 * t33;
t95 = -t82 * t44 + t61 * t73;
t94 = t61 * t122 + t130;
t93 = t60 * t122 - t129;
t34 = -t111 * t81 + (-qJD(3) * t118 + t125 * t76) * t78;
t92 = -qJD(4) * t29 + t34 * t84;
t91 = qJD(4) * t52 + t111 * t84;
t90 = -t44 * t71 - t45 * t70 + (-t60 * t87 + t61 * t83) * qJD(5) * pkin(4);
t68 = 0.2e1 * t110;
t59 = -0.2e1 * t107;
t58 = t61 ^ 2;
t46 = qJD(4) * t96 - t111 * t88;
t42 = t47 * t122;
t35 = t41 * qJD(3);
t32 = t83 * t52 + t87 * t96;
t31 = t60 * t73 + t131;
t26 = t105 * t83 + t106 * t87 + t65 * t123 + t66 * t124;
t23 = t45 * pkin(5) + t44 * pkin(11) + t115;
t22 = -t107 * t61 - t82 * t130;
t21 = qJD(6) * t108 + t126 * t44;
t20 = qJD(4) * t28 - t34 * t88;
t16 = qJD(5) * t33 - t83 * t46 + t87 * t91;
t15 = t123 * t96 + t52 * t124 + t87 * t46 + t83 * t91;
t13 = -t87 * t28 + t83 * t29;
t12 = t32 * t122 - t16 * t86;
t11 = t16 * t82 + t32 * t73;
t10 = t97 * qJD(6) + t112 * t86 + t82 * t15;
t9 = t98 * qJD(6) - t112 * t82 + t86 * t15;
t8 = -t101 * qJD(6) + t86 * t23 + t82 * t26;
t7 = -t102 * qJD(6) - t82 * t23 + t86 * t26;
t6 = qJD(5) * t14 + t83 * t20 - t87 * t92;
t5 = -t28 * t123 + t29 * t124 - t87 * t20 - t83 * t92;
t4 = t13 * t122 - t6 * t86;
t3 = t13 * t73 + t6 * t82;
t2 = -t104 * qJD(6) + t35 * t86 + t82 * t5;
t1 = t103 * qJD(6) - t35 * t82 + t86 * t5;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t35, t34, 0, 0, 0, 0, 0, t40 * t121 - t35 * t88, t40 * t120 + t35 * t84, 0, 0, 0, 0, 0, t35 * t60 + t40 * t45, t35 * t61 - t40 * t44, 0, 0, 0, 0, 0, -t103 * t45 + t13 * t95 + t6 * t136 + t2 * t60, t1 * t60 - t104 * t45 - t13 * t94 + t6 * t135; 0, 0, 0, -t112, -t111, 0, 0, 0, 0, 0 (-t89 * t121 - t88 * t125) * t77 (-t89 * t120 + t84 * t125) * t77, 0, 0, 0, 0, 0 (t60 * t125 - t45 * t89) * t77 (t61 * t125 + t44 * t89) * t77, 0, 0, 0, 0, 0, t10 * t60 + t16 * t136 + t32 * t95 - t45 * t98, t16 * t135 - t32 * t94 + t45 * t97 + t9 * t60; 0, 0, 0, 0, 0, 0.2e1 * t84 * t120, 0.2e1 * (-t84 ^ 2 + t88 ^ 2) * qJD(4), 0, 0, 0, t84 * t119, t88 * t119, -0.2e1 * t137, 0.2e1 * t44 * t60 - 0.2e1 * t61 * t45, 0, 0, 0, 0.2e1 * t60 * t115 + 0.2e1 * t72 * t45, 0.2e1 * t61 * t115 - 0.2e1 * t72 * t44, -0.2e1 * t110 * t58 - 0.2e1 * t75 * t137, 0.2e1 * t58 * t107 - t44 * t108, 0.2e1 * t61 * t129 - 0.2e1 * t60 * t94, -0.2e1 * t61 * t131 - 0.2e1 * t60 * t95, 0.2e1 * t60 * t45, 0.2e1 * t102 * t45 + 0.2e1 * t27 * t136 + 0.2e1 * t47 * t95 + 0.2e1 * t8 * t60, -0.2e1 * t101 * t45 + 0.2e1 * t27 * t135 - 0.2e1 * t47 * t94 + 0.2e1 * t7 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, -t20, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t46, 0, 0, 0, 0, 0, -t16, t15, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, t120, -t121, 0, -pkin(9) * t120, pkin(9) * t121, 0, 0, -t44, -t45, 0, -t27, t26, t22, t21, t31, -t93, 0, t42 + (-t100 * qJD(6) - t27) * t86 + t90 * t82, t100 * t122 + t86 * t90 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t114, -0.2e1 * t113, t68, t59, 0, 0, 0, 0.2e1 * t99, 0.2e1 * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, -t27, t26, t22, t21, t31, -t93, 0, t42 + (pkin(5) * t44 - pkin(11) * t45) * t82 + (-t27 + (-pkin(5) * t61 - pkin(11) * t60) * qJD(6)) * t86, pkin(5) * t94 + pkin(11) * t93 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -t113, t68, t59, 0, 0, 0, t99 - t117, -t116 + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t59, 0, 0, 0, -0.2e1 * t117, -0.2e1 * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t95, t45, t8, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t122, 0, -t113 * t82 - t70 * t73, -t113 * t86 + t122 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t122, 0, -pkin(11) * t73, pkin(11) * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t17;
