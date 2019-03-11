% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x21]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:25
% EndTime: 2019-03-08 18:54:27
% DurationCPUTime: 1.11s
% Computational Cost: add. (1124->188), mult. (3450->350), div. (0->0), fcn. (3545->12), ass. (0->110)
t58 = sin(qJ(4));
t123 = -0.4e1 * t58;
t61 = cos(qJ(4));
t74 = -t61 * pkin(4) - t58 * pkin(10);
t41 = -pkin(3) + t74;
t60 = cos(qJ(5));
t116 = t60 * t61;
t45 = pkin(9) * t116;
t57 = sin(qJ(5));
t112 = t57 * t41 + t45;
t52 = sin(pkin(7));
t59 = sin(qJ(3));
t120 = t52 * t59;
t55 = cos(pkin(7));
t29 = t58 * t120 - t61 * t55;
t62 = cos(qJ(3));
t119 = t52 * t62;
t86 = qJD(3) * t119;
t122 = -qJD(4) * t29 + t61 * t86;
t49 = t60 ^ 2;
t111 = t57 ^ 2 - t49;
t79 = t111 * qJD(5);
t121 = pkin(9) * t57;
t54 = cos(pkin(12));
t118 = t54 * t55;
t117 = t58 * t60;
t115 = -qJ(6) - pkin(10);
t101 = qJD(5) * t60;
t73 = pkin(4) * t58 - pkin(10) * t61;
t37 = t73 * qJD(4);
t114 = -t41 * t101 - t57 * t37;
t99 = t58 * qJD(4);
t85 = t57 * t99;
t113 = pkin(9) * t85 + t60 * t37;
t48 = t58 ^ 2;
t110 = -t61 ^ 2 + t48;
t109 = qJ(6) * t58;
t108 = qJ(6) * t60;
t107 = qJD(3) * t59;
t106 = qJD(3) * t61;
t104 = qJD(4) * t57;
t103 = qJD(4) * t60;
t102 = qJD(5) * t57;
t100 = qJD(5) * t61;
t98 = t60 * qJD(6);
t97 = t61 * qJD(4);
t96 = -0.2e1 * pkin(3) * qJD(4);
t95 = -0.2e1 * pkin(4) * qJD(5);
t94 = t62 * t118;
t93 = pkin(5) * t102;
t92 = pkin(9) * t97;
t91 = t57 * t100;
t90 = t60 * t100;
t51 = sin(pkin(12));
t53 = sin(pkin(6));
t56 = cos(pkin(6));
t20 = t56 * t120 + (t59 * t118 + t51 * t62) * t53;
t28 = -t53 * t54 * t52 + t56 * t55;
t11 = t20 * t58 - t28 * t61;
t89 = t11 * t102;
t88 = t29 * t102;
t87 = t52 * t107;
t84 = t57 * t97;
t83 = t57 * t101;
t82 = t58 * t97;
t81 = t60 * t97;
t80 = qJD(5) * t115;
t78 = t110 * qJD(4);
t77 = 0.2e1 * t82;
t75 = t57 * t81;
t12 = t20 * t61 + t28 * t58;
t19 = -t56 * t119 + (t51 * t59 - t94) * t53;
t5 = -t12 * t57 + t19 * t60;
t6 = t12 * t60 + t19 * t57;
t72 = -t5 * t60 - t57 * t6;
t30 = t61 * t120 + t58 * t55;
t23 = -t57 * t119 + t60 * t30;
t70 = t60 * t119 + t57 * t30;
t71 = -t23 * t57 + t60 * t70;
t15 = -t56 * t86 + (-qJD(3) * t94 + t107 * t51) * t53;
t3 = t12 * qJD(4) - t15 * t58;
t69 = t11 * t101 + t3 * t57;
t68 = -t3 * t60 + t89;
t21 = t30 * qJD(4) + t58 * t86;
t67 = t29 * t101 + t21 * t57;
t66 = -t21 * t60 + t88;
t65 = -t58 * t102 + t81;
t64 = t60 * t99 + t91;
t63 = t58 * t101 + t84;
t46 = -t60 * pkin(5) - pkin(4);
t43 = t115 * t60;
t42 = t115 * t57;
t38 = (pkin(5) * t57 + pkin(9)) * t58;
t36 = t60 * t41;
t27 = -t57 * qJD(6) + t60 * t80;
t26 = t57 * t80 + t98;
t25 = t63 * pkin(5) + t92;
t24 = -t57 * t109 + t112;
t18 = -t58 * t108 + t36 + (-pkin(5) - t121) * t61;
t16 = t20 * qJD(3);
t14 = -t112 * qJD(5) + t113;
t13 = t64 * pkin(9) + t114;
t10 = (-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t117 + (-qJD(6) * t58 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t61) * t57 - t114;
t9 = -t55 * t84 - t30 * t101 + (t60 * t107 + (t59 * t99 + (qJD(5) - t106) * t62) * t57) * t52;
t8 = t70 * qJD(5) - t60 * t122 - t57 * t87;
t7 = -t58 * t98 + (pkin(5) * t58 - t61 * t108) * qJD(4) + (-t45 + (-t41 + t109) * t57) * qJD(5) + t113;
t4 = -qJD(4) * t11 - t15 * t61;
t2 = t5 * qJD(5) + t16 * t57 + t4 * t60;
t1 = -t6 * qJD(5) + t16 * t60 - t4 * t57;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t5 * t1 + 0.2e1 * t11 * t3 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t70 + t11 * t21 + t2 * t23 + t3 * t29 + t5 * t9 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t29 * t21 - 0.2e1 * t23 * t8 - 0.2e1 * t70 * t9; 0, 0, 0, -t16, t15, 0, 0, 0, 0, 0, -t16 * t61 + t19 * t99, t16 * t58 + t19 * t97, 0, 0, 0, 0, 0 (t11 * t104 - t1) * t61 + (qJD(4) * t5 + t69) * t58 (t11 * t103 + t2) * t61 + (-qJD(4) * t6 - t68) * t58, t72 * t97 + (-t1 * t60 - t2 * t57 + (t5 * t57 - t6 * t60) * qJD(5)) * t58, t1 * t18 + t6 * t10 + t11 * t25 + t2 * t24 + t3 * t38 + t5 * t7; 0, 0, 0, -t87, -t86, 0, 0, 0, 0, 0 (-t59 * t106 - t62 * t99) * t52 (t58 * t107 - t62 * t97) * t52, 0, 0, 0, 0, 0 (t29 * t104 - t9) * t61 + (-qJD(4) * t70 + t67) * t58 (t29 * t103 - t8) * t61 + (-qJD(4) * t23 - t66) * t58, t71 * t97 + (t57 * t8 - t60 * t9 + (-t23 * t60 - t57 * t70) * qJD(5)) * t58, t23 * t10 + t9 * t18 + t21 * t38 - t8 * t24 + t29 * t25 - t7 * t70; 0, 0, 0, 0, 0, t77, -0.2e1 * t78, 0, 0, 0, t58 * t96, t61 * t96, -0.2e1 * t48 * t83 + 0.2e1 * t49 * t82, t75 * t123 + 0.2e1 * t48 * t79, 0.2e1 * t110 * t103 + 0.2e1 * t58 * t91, -0.2e1 * t57 * t78 + 0.2e1 * t58 * t90, -0.2e1 * t82, 0.2e1 * t36 * t99 - 0.2e1 * t14 * t61 + 0.2e1 * (t48 * t101 + t57 * t82) * pkin(9), -0.2e1 * t13 * t61 - 0.2e1 * t112 * t99 + 0.2e1 * (-t48 * t102 + t60 * t77) * pkin(9), 0.2e1 * (-t18 * t60 - t24 * t57) * t97 + 0.2e1 * (-t10 * t57 - t60 * t7 + (t18 * t57 - t24 * t60) * qJD(5)) * t58, 0.2e1 * t24 * t10 + 0.2e1 * t18 * t7 + 0.2e1 * t38 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, t68, t69, t72 * qJD(5) - t1 * t57 + t2 * t60, pkin(5) * t89 + t1 * t42 - t2 * t43 + t6 * t26 + t5 * t27 + t3 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t122, 0, 0, 0, 0, 0, t66, t67, t71 * qJD(5) - t9 * t57 - t8 * t60, pkin(5) * t88 + t21 * t46 + t23 * t26 - t27 * t70 + t9 * t42 + t8 * t43; 0, 0, 0, 0, 0, 0, 0, t97, -t99, 0, -t92, pkin(9) * t99, -t58 * t79 + t75, -t111 * t97 + t83 * t123, t85 - t90, t64, 0 (pkin(10) * t116 + (-pkin(4) * t60 + t121) * t58) * qJD(5) + (t74 * t57 - t45) * qJD(4) (pkin(9) * t117 + t73 * t57) * qJD(5) + (t61 * t121 + t74 * t60) * qJD(4) (-t42 * t97 - t27 * t58 + t10 + (t43 * t58 - t18) * qJD(5)) * t60 + (t43 * t97 - t26 * t58 - t7 + (t42 * t58 - t24) * qJD(5)) * t57, -t10 * t43 + t18 * t27 + t24 * t26 + t25 * t46 + t38 * t93 + t7 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t83, -0.2e1 * t79, 0, 0, 0, t57 * t95, t60 * t95, 0.2e1 * t26 * t60 - 0.2e1 * t27 * t57 + 0.2e1 * (-t42 * t60 + t43 * t57) * qJD(5), -0.2e1 * t43 * t26 + 0.2e1 * t42 * t27 + 0.2e1 * t46 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t8, 0, t9 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t63, t99, t14, t13, -t65 * pkin(5), t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t102, 0, -pkin(10) * t101, pkin(10) * t102, -pkin(5) * t101, t27 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t17;
