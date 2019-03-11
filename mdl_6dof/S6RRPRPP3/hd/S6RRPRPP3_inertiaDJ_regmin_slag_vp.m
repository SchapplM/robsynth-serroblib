% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:09
% EndTime: 2019-03-09 09:57:14
% DurationCPUTime: 1.53s
% Computational Cost: add. (2363->237), mult. (5741->413), div. (0->0), fcn. (5019->6), ass. (0->101)
t126 = pkin(4) + qJ(6);
t87 = sin(qJ(2));
t78 = qJD(2) * t87;
t135 = t126 * t78;
t131 = cos(qJ(4));
t107 = qJD(4) * t131;
t86 = sin(qJ(4));
t120 = qJD(4) * t86;
t83 = sin(pkin(9));
t84 = cos(pkin(9));
t134 = t84 * t107 - t83 * t120;
t109 = t131 * t84;
t125 = pkin(8) + qJ(3);
t64 = t125 * t83;
t65 = t125 * t84;
t22 = (qJD(3) * t83 + qJD(4) * t65) * t86 - qJD(3) * t109 + t64 * t107;
t133 = 0.2e1 * t134;
t88 = cos(qJ(2));
t122 = qJD(2) * t88;
t59 = t131 * t83 + t86 * t84;
t28 = t59 * t122 + t134 * t87;
t132 = pkin(5) * t28;
t130 = t83 * t87;
t129 = t83 * t88;
t128 = t84 * t87;
t127 = t84 * t88;
t101 = -pkin(2) * t88 - qJ(3) * t87;
t61 = -pkin(1) + t101;
t57 = t84 * t61;
t37 = -pkin(8) * t128 + t57 + (-pkin(7) * t83 - pkin(3)) * t88;
t72 = pkin(7) * t127;
t46 = t83 * t61 + t72;
t42 = -pkin(8) * t130 + t46;
t124 = t131 * t42 + t86 * t37;
t113 = pkin(7) * t78;
t51 = -t87 * qJD(3) + (pkin(2) * t87 - qJ(3) * t88) * qJD(2);
t40 = t83 * t113 + t84 * t51;
t111 = t83 * t122;
t76 = pkin(7) * t122;
t54 = pkin(3) * t111 + t76;
t60 = pkin(3) * t130 + t87 * pkin(7);
t121 = qJD(3) * t88;
t119 = t88 * qJD(5);
t118 = t88 * qJD(6);
t117 = qJ(5) * qJD(5);
t116 = pkin(7) * t129;
t115 = -0.2e1 * pkin(1) * qJD(2);
t114 = pkin(4) * t78;
t110 = t87 * t122;
t75 = -t84 * pkin(3) - pkin(2);
t108 = qJ(5) * t78;
t24 = (pkin(3) * t87 - pkin(8) * t127) * qJD(2) + t40;
t47 = t83 * t51;
t31 = t47 + (-pkin(7) * t128 - pkin(8) * t129) * qJD(2);
t106 = t42 * t107 + t37 * t120 - t131 * t24 + t86 * t31;
t105 = 0.2e1 * (t83 ^ 2 + t84 ^ 2) * qJD(3);
t12 = qJ(5) * t88 - t124;
t50 = t87 * t109 - t86 * t130;
t103 = -t50 * qJ(5) + t60;
t102 = t131 * t37 - t86 * t42;
t41 = -t84 * t113 + t47;
t100 = -t40 * t83 + t41 * t84;
t49 = t59 * t87;
t99 = -qJ(5) * t28 - qJD(5) * t49;
t98 = -qJ(5) * t134 - t59 * qJD(5);
t53 = t59 * qJD(4);
t58 = t83 * t86 - t109;
t97 = -qJ(5) * t53 - qJD(5) * t58;
t13 = t88 * pkin(4) - t102;
t95 = -t59 * qJ(5) + t75;
t27 = -t109 * t122 + t86 * t111 + t87 * t53;
t94 = -t27 * pkin(5) + t106;
t44 = t131 * t65 - t86 * t64;
t93 = -t22 * t88 - t44 * t78;
t23 = t59 * qJD(3) + t44 * qJD(4);
t43 = t131 * t64 + t86 * t65;
t92 = t23 * t88 - t43 * t78;
t6 = -t37 * t107 + t42 * t120 - t131 * t31 - t86 * t24;
t91 = t27 * qJ(5) - t50 * qJD(5) + t54;
t90 = -t6 + 0.2e1 * t108 - 0.2e1 * t119;
t4 = -t108 + t6 + t119;
t89 = 0.2e1 * qJD(5);
t45 = t57 - t116;
t35 = pkin(4) * t58 + t95;
t26 = -t58 * pkin(5) + t44;
t25 = t59 * pkin(5) + t43;
t19 = t126 * t58 + t95;
t18 = pkin(4) * t53 + t98;
t17 = pkin(4) * t49 + t103;
t16 = pkin(5) * t134 + t23;
t15 = -t53 * pkin(5) - t22;
t14 = t126 * t49 + t103;
t11 = t58 * qJD(6) + t126 * t53 + t98;
t10 = -pkin(5) * t49 - t12;
t9 = t50 * pkin(5) + t88 * qJ(6) + t13;
t8 = pkin(4) * t28 + t91;
t5 = t106 - t114;
t3 = t49 * qJD(6) + t126 * t28 + t91;
t2 = -t4 - t132;
t1 = t118 + t94 - t135;
t7 = [0, 0, 0, 0.2e1 * t110, 0.2e1 * (-t87 ^ 2 + t88 ^ 2) * qJD(2), 0, 0, 0, t87 * t115, t88 * t115, -0.2e1 * t40 * t88 + 0.2e1 * (t45 + 0.2e1 * t116) * t78, 0.2e1 * t41 * t88 + 0.2e1 * (-t46 + 0.2e1 * t72) * t78, 0.2e1 * (-t40 * t84 - t41 * t83) * t87 + 0.2e1 * (-t45 * t84 - t46 * t83) * t122, 0.2e1 * pkin(7) ^ 2 * t110 + 0.2e1 * t40 * t45 + 0.2e1 * t41 * t46, -0.2e1 * t50 * t27, 0.2e1 * t27 * t49 - 0.2e1 * t28 * t50, 0.2e1 * t27 * t88 + 0.2e1 * t50 * t78, 0.2e1 * t28 * t88 - 0.2e1 * t49 * t78, -0.2e1 * t110, 0.2e1 * t102 * t78 + 0.2e1 * t106 * t88 + 0.2e1 * t60 * t28 + 0.2e1 * t54 * t49, -0.2e1 * t124 * t78 - 0.2e1 * t60 * t27 + 0.2e1 * t54 * t50 - 0.2e1 * t6 * t88, 0.2e1 * t12 * t28 - 0.2e1 * t13 * t27 + 0.2e1 * t4 * t49 + 0.2e1 * t5 * t50, 0.2e1 * t13 * t78 - 0.2e1 * t17 * t28 - 0.2e1 * t49 * t8 - 0.2e1 * t5 * t88, -0.2e1 * t12 * t78 + 0.2e1 * t17 * t27 + 0.2e1 * t4 * t88 - 0.2e1 * t50 * t8, 0.2e1 * t12 * t4 + 0.2e1 * t13 * t5 + 0.2e1 * t17 * t8, 0.2e1 * t1 * t50 - 0.2e1 * t10 * t28 - 0.2e1 * t2 * t49 - 0.2e1 * t27 * t9, 0.2e1 * t10 * t78 + 0.2e1 * t14 * t27 - 0.2e1 * t2 * t88 - 0.2e1 * t3 * t50, 0.2e1 * t1 * t88 + 0.2e1 * t14 * t28 + 0.2e1 * t3 * t49 - 0.2e1 * t9 * t78, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t14 * t3; 0, 0, 0, 0, 0, t122, -t78, 0, -t76, t113, t83 * t121 + (t101 * t83 - t72) * qJD(2), t84 * t121 + (t101 * t84 + t116) * qJD(2), t100, -pkin(2) * t76 + (-t45 * t83 + t46 * t84) * qJD(3) + t100 * qJ(3), t134 * t50 - t27 * t59, -t134 * t49 + t27 * t58 - t28 * t59 - t50 * t53, -t134 * t88 + t59 * t78, t53 * t88 - t58 * t78, 0, t28 * t75 + t53 * t60 + t54 * t58 + t92, t134 * t60 - t27 * t75 + t54 * t59 + t93, t12 * t53 + t13 * t134 + t22 * t49 + t23 * t50 - t27 * t43 - t28 * t44 + t4 * t58 + t5 * t59, -t17 * t53 - t18 * t49 - t28 * t35 - t58 * t8 - t92, -t134 * t17 - t18 * t50 + t27 * t35 - t59 * t8 - t93, t12 * t22 + t13 * t23 + t17 * t18 + t35 * t8 - t4 * t44 + t43 * t5, t1 * t59 - t10 * t53 + t134 * t9 - t15 * t49 + t16 * t50 - t2 * t58 - t25 * t27 - t26 * t28, -t11 * t50 - t134 * t14 - t15 * t88 + t19 * t27 + t26 * t78 - t3 * t59, t11 * t49 + t14 * t53 + t16 * t88 + t19 * t28 - t25 * t78 + t3 * t58, t1 * t25 + t10 * t15 + t11 * t14 + t16 * t9 + t19 * t3 + t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, qJ(3) * t105, t59 * t133, -0.2e1 * t134 * t58 - 0.2e1 * t53 * t59, 0, 0, 0, 0.2e1 * t75 * t53, t75 * t133, 0.2e1 * t134 * t43 + 0.2e1 * t22 * t58 + 0.2e1 * t23 * t59 - 0.2e1 * t44 * t53, -0.2e1 * t18 * t58 - 0.2e1 * t35 * t53, -0.2e1 * t134 * t35 - 0.2e1 * t18 * t59, 0.2e1 * t18 * t35 - 0.2e1 * t22 * t44 + 0.2e1 * t23 * t43, 0.2e1 * t134 * t25 - 0.2e1 * t15 * t58 + 0.2e1 * t16 * t59 - 0.2e1 * t26 * t53, -0.2e1 * t11 * t59 - 0.2e1 * t134 * t19, 0.2e1 * t11 * t58 + 0.2e1 * t19 * t53, 0.2e1 * t11 * t19 + 0.2e1 * t15 * t26 + 0.2e1 * t16 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t84 * t122, 0, t76, 0, 0, 0, 0, 0, t28, -t27, 0, -t28, t27, t8, 0, t27, t28, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t134, 0, -t53, -t134, t18, 0, -t134, t53, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t28, t78, -t106, t6, pkin(4) * t27 + t99, t106 - 0.2e1 * t114, t90, -pkin(4) * t5 - qJ(5) * t4 - qJD(5) * t12, -qJD(6) * t50 + t126 * t27 + t99, t90 - t132, -0.2e1 * t118 - t94 + 0.2e1 * t135, qJ(5) * t2 + qJD(5) * t10 - qJD(6) * t9 - t1 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t53, 0, -t23, t22, -pkin(4) * t134 + t97, t23, -t22, -pkin(4) * t23 - qJ(5) * t22 + qJD(5) * t44, -qJD(6) * t59 - t126 * t134 + t97, t15, -t16, qJ(5) * t15 + qJD(5) * t26 - qJD(6) * t25 - t126 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0.2e1 * t117, 0, t89, 0.2e1 * qJD(6), 0.2e1 * qJD(6) * t126 + 0.2e1 * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t78, 0, t5, -t27, 0, -t78, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, 0, 0, t23, t134, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t78, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
