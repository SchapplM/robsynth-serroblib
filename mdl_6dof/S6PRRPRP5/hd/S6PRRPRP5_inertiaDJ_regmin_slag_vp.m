% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:16
% EndTime: 2019-03-08 21:49:20
% DurationCPUTime: 1.27s
% Computational Cost: add. (975->186), mult. (2435->323), div. (0->0), fcn. (2057->8), ass. (0->118)
t132 = pkin(4) + pkin(8);
t62 = sin(qJ(3));
t117 = t62 * qJ(4);
t65 = cos(qJ(3));
t67 = -pkin(3) - pkin(9);
t122 = t65 * t67;
t136 = t117 - t122;
t37 = -pkin(2) - t136;
t47 = t132 * t62;
t61 = sin(qJ(5));
t64 = cos(qJ(5));
t135 = t64 * t37 + t61 * t47;
t58 = t65 ^ 2;
t83 = qJD(3) * (t62 ^ 2 - t58);
t55 = t61 ^ 2;
t120 = -t64 ^ 2 + t55;
t82 = t120 * qJD(5);
t118 = qJ(4) * t65;
t106 = t62 * qJD(3);
t81 = pkin(3) * t106 - t62 * qJD(4);
t23 = (pkin(9) * t62 - t118) * qJD(3) + t81;
t52 = t65 * qJD(3);
t50 = pkin(8) * t52;
t41 = pkin(4) * t52 + t50;
t7 = -qJD(5) * t135 - t61 * t23 + t64 * t41;
t134 = 0.2e1 * qJD(4);
t133 = 2 * qJD(6);
t77 = -t61 * pkin(5) + t64 * qJ(6);
t28 = t77 * qJD(5) + t61 * qJD(6);
t78 = pkin(5) * t64 + qJ(6) * t61;
t12 = t28 * t65 + (-t78 - t132) * t106;
t131 = t12 * t61;
t59 = sin(pkin(6));
t63 = sin(qJ(2));
t125 = t59 * t63;
t100 = t62 * t125;
t60 = cos(pkin(6));
t66 = cos(qJ(2));
t115 = qJD(2) * t66;
t90 = t59 * t115;
t17 = -qJD(3) * t100 + (qJD(3) * t60 + t90) * t65;
t130 = t17 * t65;
t26 = t78 * qJD(5) - t64 * qJD(6) + qJD(4);
t129 = t26 * t65;
t31 = t65 * t125 + t60 * t62;
t128 = t31 * t17;
t127 = t31 * t62;
t40 = t132 * t106;
t126 = t40 * t61;
t124 = t59 * t66;
t123 = t62 * t67;
t48 = t132 * t65;
t116 = qJD(2) * t63;
t114 = qJD(3) * t61;
t113 = qJD(3) * t64;
t112 = qJD(4) * t65;
t111 = qJD(5) * t61;
t110 = qJD(5) * t64;
t109 = qJD(5) * t65;
t108 = qJD(5) * t67;
t107 = t48 * qJD(5);
t105 = t62 * qJD(6);
t104 = qJ(4) * qJD(5);
t103 = qJ(6) * qJD(3);
t102 = qJD(3) * qJ(4);
t101 = -0.2e1 * pkin(2) * qJD(3);
t99 = t64 * t124;
t98 = pkin(5) * t52;
t97 = pkin(8) * t106;
t96 = t61 * t109;
t95 = t61 * t108;
t94 = t64 * t109;
t93 = t64 * t108;
t92 = t31 * t111;
t91 = t59 * t116;
t89 = t62 * t52;
t88 = t64 * t106;
t87 = t64 * t52;
t86 = t61 * t110;
t85 = t65 * t103;
t80 = t61 * t88;
t79 = -t65 * pkin(3) - t117;
t13 = t62 * qJ(6) + t135;
t73 = -t61 * t37 + t64 * t47;
t14 = -t62 * pkin(5) - t73;
t76 = t13 * t64 + t14 * t61;
t30 = -t60 * t65 + t100;
t19 = t61 * t124 + t30 * t64;
t20 = t30 * t61 - t99;
t75 = -t19 * t61 + t20 * t64;
t42 = qJ(4) - t77;
t72 = t42 * t65 + t123;
t71 = -t118 - t123;
t10 = t31 * t110 + t17 * t61;
t6 = -t47 * t110 + t37 * t111 - t64 * t23 - t61 * t41;
t70 = t79 * qJD(3) + t112;
t4 = -t6 + t85 + t105;
t5 = -t7 - t98;
t1 = t76 * qJD(5) + t4 * t61 - t5 * t64;
t18 = t31 * qJD(3) + t62 * t90;
t8 = -t18 * t64 - qJD(5) * t99 + (qJD(5) * t30 + t91) * t61;
t9 = t19 * qJD(5) + t18 * t61 + t64 * t91;
t2 = t75 * qJD(5) + t9 * t61 - t8 * t64;
t69 = t130 + t18 * t62 + (t30 * t65 - t127) * qJD(3);
t68 = -t65 * t92 - t8 * t62 + t64 * t130 + (-t64 * t127 + t19 * t65) * qJD(3);
t49 = 0.2e1 * t89;
t46 = t67 * t87;
t43 = -pkin(2) + t79;
t34 = t61 * t106 - t94;
t33 = t62 * t110 + t61 * t52;
t32 = -t62 * t111 + t87;
t29 = -t65 * t102 + t81;
t25 = (t66 * t106 + t65 * t116) * t59;
t24 = (t62 * t116 - t66 * t52) * t59;
t22 = t78 * t65 + t48;
t11 = t17 * t64 - t92;
t3 = (-t31 * t114 + t9) * t62 + (qJD(3) * t20 + t10) * t65;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t59 ^ 2 * t63 * t115 + 0.2e1 * t30 * t18 + 0.2e1 * t128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t19 * t8 + 0.2e1 * t20 * t9 + 0.2e1 * t128; 0, 0, -t91, -t90, 0, 0, 0, 0, 0, -t25, t24, t69, t25, -t24 (t43 * t116 - t29 * t66) * t59 + t69 * pkin(8), 0, 0, 0, 0, 0, t68, -t3, t68, t75 * t106 + (-t61 * t8 - t64 * t9 + (t19 * t64 + t20 * t61) * qJD(5)) * t65, t3, t31 * t12 + t9 * t13 + t8 * t14 + t17 * t22 - t19 * t5 + t20 * t4; 0, 0, 0, 0, t49, -0.2e1 * t83, 0, 0, 0, t62 * t101, t65 * t101, 0, -0.2e1 * t43 * t106 + 0.2e1 * t29 * t65, -0.2e1 * t29 * t62 - 0.2e1 * t43 * t52, 0.2e1 * t43 * t29, -0.2e1 * t55 * t89 + 0.2e1 * t58 * t86, -0.2e1 * t58 * t82 - 0.4e1 * t65 * t80, 0.2e1 * t61 * t83 - 0.2e1 * t62 * t94, 0.2e1 * t62 * t96 + 0.2e1 * t64 * t83, t49, 0.2e1 * (-t48 * t113 + t7) * t62 + 0.2e1 * (t73 * qJD(3) - t61 * t107 - t40 * t64) * t65, 0.2e1 * (t48 * t114 + t6) * t62 + 0.2e1 * (-qJD(3) * t135 - t64 * t107 + t126) * t65, 0.2e1 * (-t22 * t113 - t5) * t62 + 0.2e1 * (-qJD(3) * t14 - t22 * t111 + t12 * t64) * t65, 0.2e1 * t76 * t106 + 0.2e1 * (-t4 * t64 - t5 * t61 + (t13 * t61 - t14 * t64) * qJD(5)) * t65, 0.2e1 * (-t22 * t114 + t4) * t62 + 0.2e1 * (qJD(3) * t13 + t22 * t110 + t131) * t65, 0.2e1 * t22 * t12 + 0.2e1 * t13 * t4 + 0.2e1 * t14 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t17, 0, t18, t17, -t18 * pkin(3) + t17 * qJ(4) + t31 * qJD(4), 0, 0, 0, 0, 0, t10, t11, t10, -t2, -t11, t17 * t42 + t2 * t67 + t31 * t26; 0, 0, 0, 0, 0, 0, t52, -t106, 0, -t50, t97, t70, t50, -t97, t70 * pkin(8), t65 * t82 + t80, -t120 * t106 + 0.4e1 * t65 * t86, t32, -t33, 0, -t126 + t46 + (-t62 * t102 + t112) * t64 + (t48 * t64 + t71 * t61) * qJD(5) (t71 * qJD(5) - t40) * t64 + (t136 * qJD(3) - t107 - t112) * t61, t131 + t46 + (-t42 * t106 + t129) * t64 + (t22 * t64 - t72 * t61) * qJD(5), -t1 (t72 * qJD(5) - t12) * t64 + (qJD(5) * t22 + t129 + (-t42 * t62 + t122) * qJD(3)) * t61, t1 * t67 + t12 * t42 + t22 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, qJ(4) * t134, -0.2e1 * t86, 0.2e1 * t82, 0, 0, 0, 0.2e1 * qJD(4) * t61 + 0.2e1 * t64 * t104, 0.2e1 * qJD(4) * t64 - 0.2e1 * t61 * t104, 0.2e1 * t42 * t110 + 0.2e1 * t26 * t61, 0, 0.2e1 * t42 * t111 - 0.2e1 * t26 * t64, 0.2e1 * t42 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, t50, 0, 0, 0, 0, 0, t32, -t33, t32, 0, t33, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, -t8, 0, t9, -t8 * pkin(5) + t9 * qJ(6) + t20 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t88 + t96, t52, t7, t6, t7 + 0.2e1 * t98 (-pkin(5) * t106 + qJ(6) * t109) * t61 + (t62 * t103 + (pkin(5) * qJD(5) - qJD(6)) * t65) * t64, -t6 + 0.2e1 * t85 + 0.2e1 * t105, -t5 * pkin(5) + t4 * qJ(6) + t13 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t110, 0, -t95, -t93, -t95, -t28, t93, t28 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t110, -t111, 0, t110, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, qJ(6) * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t34, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, 0, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t15;
