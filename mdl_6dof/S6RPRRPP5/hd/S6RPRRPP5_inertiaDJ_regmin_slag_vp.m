% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:42
% EndTime: 2019-03-09 04:44:46
% DurationCPUTime: 1.52s
% Computational Cost: add. (2205->199), mult. (4937->336), div. (0->0), fcn. (4613->6), ass. (0->109)
t79 = cos(qJ(4));
t117 = t79 * qJD(5);
t77 = sin(qJ(4));
t120 = t77 * qJ(5);
t134 = pkin(4) + pkin(5);
t95 = t134 * t79 + t120;
t143 = qJD(4) * t95 - t117;
t131 = cos(qJ(3));
t75 = sin(pkin(9));
t76 = cos(pkin(9));
t78 = sin(qJ(3));
t58 = t131 * t75 + t78 * t76;
t49 = t58 * qJD(3);
t109 = t131 * t76;
t91 = -t78 * t75 + t109;
t141 = t49 * qJ(5) - t91 * qJD(5);
t122 = qJ(5) * t79;
t92 = -t134 * t77 + t122;
t142 = 0.2e1 * t141;
t126 = pkin(7) + qJ(2);
t60 = t126 * t75;
t61 = t126 * t76;
t36 = t131 * t61 - t78 * t60;
t140 = t36 * qJD(3);
t73 = t77 ^ 2;
t74 = t79 ^ 2;
t123 = t73 - t74;
t106 = t123 * qJD(4);
t116 = t79 * qJD(6);
t48 = t91 * qJD(3);
t127 = t79 * t48;
t68 = qJD(4) * t77;
t27 = t58 * t68 - t127;
t44 = t49 * pkin(4);
t110 = t131 * t60;
t17 = qJD(3) * t110 - qJD(2) * t109 + (qJD(2) * t75 + qJD(3) * t61) * t78;
t30 = pkin(3) * t49 - pkin(8) * t48;
t66 = -pkin(2) * t76 - pkin(1);
t31 = -pkin(3) * t91 - pkin(8) * t58 + t66;
t69 = qJD(4) * t79;
t7 = t77 * t17 + t79 * t30 - t31 * t68 - t36 * t69;
t5 = -t44 - t7;
t139 = -t27 * qJ(6) + t58 * t116 - t5;
t124 = t77 * t31 + t79 * t36;
t11 = -qJ(5) * t91 + t124;
t33 = t77 * t36;
t107 = t79 * t31 - t33;
t12 = pkin(4) * t91 - t107;
t6 = t79 * t17 - t77 * t30 - t31 * t69 + t36 * t68;
t3 = -t6 + t141;
t137 = t3 * t77 - t5 * t79 + (t11 * t79 + t12 * t77) * qJD(4);
t99 = pkin(4) * t79 + t120;
t136 = t99 * qJD(4) - t117;
t135 = 0.2e1 * t66;
t81 = 0.2e1 * qJD(5);
t133 = pkin(8) * t49;
t132 = pkin(8) * t91;
t130 = t58 * t48;
t129 = t58 * t77;
t128 = t58 * t79;
t125 = pkin(8) - qJ(6);
t121 = qJ(6) * t77;
t119 = t77 * qJD(5);
t118 = t77 * qJD(6);
t115 = -0.2e1 * pkin(3) * qJD(4);
t114 = pkin(8) * t68;
t113 = pkin(8) * t69;
t112 = t58 * t69;
t111 = t77 * t69;
t63 = t125 * t79;
t108 = -0.4e1 * t77 * t128;
t105 = -pkin(4) * t68 + t119;
t104 = 0.2e1 * (t75 ^ 2 + t76 ^ 2) * qJD(2);
t4 = -t140 + t92 * t48 + (-qJD(2) - t143) * t58;
t54 = pkin(3) + t95;
t103 = -qJD(4) * t54 * t58 + t4;
t102 = -pkin(3) * t48 - t133;
t101 = pkin(3) * t58 - t132;
t35 = t78 * t61 + t110;
t98 = pkin(4) * t77 - t122;
t96 = -t11 * t77 + t12 * t79;
t94 = t48 * t77 + t112;
t28 = t49 * t77 - t69 * t91;
t93 = t49 * t79 + t68 * t91;
t13 = t92 * t58 - t35;
t37 = (-pkin(5) * t77 + t122) * qJD(4) + t105;
t90 = qJD(4) * t13 + t37 * t58 + t48 * t54;
t59 = -pkin(3) - t99;
t18 = t58 * qJD(2) + t140;
t8 = t136 * t58 + t98 * t48 + t18;
t89 = -t8 + (t58 * t59 + t132) * qJD(4);
t87 = qJ(6) * t112 + t58 * t118 + t48 * t121 - t6;
t14 = t98 * t58 + t35;
t46 = qJ(5) * t69 + t105;
t86 = -qJD(4) * t14 + t46 * t58 - t48 * t59 + t133;
t84 = t92 * qJD(4) + t119;
t1 = -pkin(5) * t49 - t139;
t10 = t58 * t121 + t11;
t2 = t87 + t141;
t9 = t33 + (-qJ(6) * t58 - t31) * t79 + t134 * t91;
t83 = -t1 * t79 + t2 * t77 + (t10 * t79 + t77 * t9) * qJD(4);
t82 = t96 * qJD(4) + t3 * t79 + t5 * t77;
t70 = qJ(5) * t81;
t62 = t125 * t77;
t55 = t58 ^ 2;
t47 = qJD(4) * t63 - t118;
t45 = -t125 * t68 - t116;
t21 = (t73 + t74) * t48;
t15 = [0, 0, 0, 0, 0, t104, qJ(2) * t104, 0.2e1 * t130, 0.2e1 * t48 * t91 - 0.2e1 * t49 * t58, 0, 0, 0, t49 * t135, t48 * t135, -0.2e1 * t55 * t111 + 0.2e1 * t74 * t130, 0.2e1 * t55 * t106 + t48 * t108, 0.2e1 * t49 * t128 + 0.2e1 * t27 * t91, -0.2e1 * t49 * t129 + 0.2e1 * t91 * t94, -0.2e1 * t91 * t49, 0.2e1 * t107 * t49 + 0.2e1 * t18 * t129 + 0.2e1 * t94 * t35 - 0.2e1 * t7 * t91, -0.2e1 * t124 * t49 + 0.2e1 * t18 * t128 - 0.2e1 * t27 * t35 - 0.2e1 * t6 * t91, -0.2e1 * t12 * t49 + 0.2e1 * t8 * t129 + 0.2e1 * t94 * t14 + 0.2e1 * t5 * t91, -0.2e1 * t137 * t58 + 0.2e1 * t96 * t48, 0.2e1 * t11 * t49 - 0.2e1 * t8 * t128 + 0.2e1 * t27 * t14 - 0.2e1 * t3 * t91, 0.2e1 * t11 * t3 + 0.2e1 * t12 * t5 + 0.2e1 * t14 * t8, 0.2e1 * t1 * t91 - 0.2e1 * t4 * t129 - 0.2e1 * t94 * t13 - 0.2e1 * t9 * t49, 0.2e1 * t10 * t49 + 0.2e1 * t4 * t128 - 0.2e1 * t27 * t13 - 0.2e1 * t2 * t91, 0.2e1 * (t10 * t77 - t79 * t9) * t48 + 0.2e1 * t83 * t58, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t13 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t48, 0, 0, 0, 0, 0, t93, -t28, t93, -t21, t28, t137, t93, t28, t21, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t49, 0, -t18, t17, -t58 * t106 + t77 * t127, qJD(4) * t108 - t123 * t48, t28, t93, 0, -t18 * t79 + t102 * t77 + (-t101 * t79 + t35 * t77) * qJD(4), t18 * t77 + t102 * t79 + (t101 * t77 + t35 * t79) * qJD(4), -t86 * t77 + t89 * t79, t82, t89 * t77 + t86 * t79, t82 * pkin(8) - t14 * t46 + t59 * t8, t103 * t79 + t47 * t91 - t49 * t62 - t90 * t77, t103 * t77 - t45 * t91 + t49 * t63 + t90 * t79 (-t47 * t58 - t48 * t62 - t2 + (t58 * t63 - t9) * qJD(4)) * t79 + (t45 * t58 + t48 * t63 - t1 + (t58 * t62 + t10) * qJD(4)) * t77, t1 * t62 + t10 * t45 + t13 * t37 + t2 * t63 + t4 * t54 + t47 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t45 - t79 * t47 + (t62 * t77 + t63 * t79) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t111, -0.2e1 * t106, 0, 0, 0, t77 * t115, t79 * t115, 0.2e1 * t46 * t79 + 0.2e1 * t59 * t68, 0, 0.2e1 * t46 * t77 - 0.2e1 * t59 * t69, -0.2e1 * t59 * t46, 0.2e1 * t37 * t79 - 0.2e1 * t54 * t68, 0.2e1 * t37 * t77 + 0.2e1 * t54 * t69, -0.2e1 * t45 * t79 - 0.2e1 * t47 * t77 + 0.2e1 * (-t62 * t79 + t63 * t77) * qJD(4), 0.2e1 * t37 * t54 + 0.2e1 * t45 * t63 + 0.2e1 * t47 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t94, t49, t7, t6, -t5 + t44, -t99 * t48 + (t98 * qJD(4) - t119) * t58, -t6 + t142, -pkin(4) * t5 + qJ(5) * t3 + qJD(5) * t11 (pkin(5) + t134) * t49 + t139, t87 + t142, t48 * t95 + t58 * t84, qJ(5) * t2 + qJD(5) * t10 - t1 * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t69, -t68, 0, t69, t46, -t68, t69, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t68, 0, -t113, t114, -t113, -t136, -t114, -t136 * pkin(8), -t47, t45, t143, qJ(5) * t45 + qJD(5) * t63 - t134 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t70, 0, t81, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t27, 0, t5, -t49, 0, t27, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t113, 0, 0, -t69, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t27, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t69, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t15;
