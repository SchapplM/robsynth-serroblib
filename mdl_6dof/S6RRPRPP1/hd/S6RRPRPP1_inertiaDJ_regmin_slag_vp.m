% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:47:51
% EndTime: 2019-03-09 09:47:55
% DurationCPUTime: 1.47s
% Computational Cost: add. (3708->209), mult. (8240->382), div. (0->0), fcn. (7771->8), ass. (0->114)
t104 = sin(qJ(4));
t106 = cos(qJ(4));
t140 = qJD(4) * t106;
t102 = sin(pkin(9));
t103 = cos(pkin(9));
t105 = sin(qJ(2));
t107 = cos(qJ(2));
t84 = t102 * t107 + t103 * t105;
t134 = t84 * t140;
t142 = qJD(2) * t107;
t143 = qJD(2) * t105;
t78 = -t102 * t143 + t103 * t142;
t159 = t104 * t78 + t134;
t101 = sin(pkin(10));
t149 = cos(pkin(10));
t126 = t149 * t106;
t141 = qJD(4) * t104;
t77 = qJD(4) * t126 - t101 * t141;
t148 = t101 * t104;
t158 = t126 - t148;
t127 = t149 * t104;
t83 = t101 * t106 + t127;
t100 = t106 ^ 2;
t150 = -t104 ^ 2 + t100;
t123 = t150 * qJD(4);
t157 = 2 * qJD(6);
t117 = -qJ(5) * t78 - qJD(5) * t84;
t156 = -qJ(3) - pkin(7);
t128 = qJD(2) * t156;
t73 = t107 * qJD(3) + t105 * t128;
t74 = -t105 * qJD(3) + t107 * t128;
t42 = t102 * t74 + t103 * t73;
t76 = t84 * qJD(2);
t98 = pkin(2) * t143;
t43 = t76 * pkin(3) - t78 * pkin(8) + t98;
t133 = -t107 * pkin(2) - pkin(1);
t82 = t102 * t105 - t103 * t107;
t53 = t82 * pkin(3) - t84 * pkin(8) + t133;
t136 = t104 * t43 + t106 * t42 + t53 * t140;
t87 = t156 * t105;
t88 = t156 * t107;
t59 = t102 * t87 - t103 * t88;
t10 = -qJ(5) * t134 + (-qJD(4) * t59 + t117) * t104 + t136;
t129 = -t104 * t42 + t106 * t43;
t154 = qJ(5) * t84;
t54 = t106 * t59;
t8 = t76 * pkin(4) + t117 * t106 + (-t54 + (-t53 + t154) * t104) * qJD(4) + t129;
t4 = t149 * t10 + t101 * t8;
t49 = t106 * t53;
t20 = t82 * pkin(4) - t104 * t59 - t106 * t154 + t49;
t152 = t104 * t84;
t155 = t104 * t53 + t54;
t25 = -qJ(5) * t152 + t155;
t14 = t101 * t20 + t149 * t25;
t94 = t102 * pkin(2) + pkin(8);
t151 = qJ(5) + t94;
t146 = t104 * t106;
t145 = t82 * qJD(6);
t144 = t83 * qJD(6);
t139 = t76 * qJ(6) + t4;
t138 = -0.2e1 * pkin(1) * qJD(2);
t96 = -t103 * pkin(2) - pkin(3);
t137 = 0.2e1 * qJD(4) * t96;
t97 = pkin(4) * t141;
t135 = t84 * t141;
t131 = t104 * t140;
t3 = -t101 * t10 + t149 * t8;
t124 = qJD(4) * t151;
t109 = -t104 * qJD(5) - t106 * t124;
t64 = t106 * qJD(5) - t104 * t124;
t35 = t101 * t64 - t149 * t109;
t36 = t101 * t109 + t149 * t64;
t79 = t151 * t106;
t51 = t101 * t79 + t151 * t127;
t52 = -t151 * t148 + t149 * t79;
t130 = t51 * t35 + t52 * t36;
t41 = t102 * t73 - t103 * t74;
t58 = -t102 * t88 - t103 * t87;
t125 = -0.4e1 * t84 * t146;
t40 = pkin(4) * t152 + t58;
t121 = t41 * t84 + t58 * t78;
t120 = t76 * t84 + t78 * t82;
t119 = -t76 * t94 + t78 * t96;
t118 = t82 * t94 - t84 * t96;
t28 = t159 * pkin(4) + t41;
t75 = t83 * qJD(4);
t116 = -0.2e1 * t158 * t75 + 0.2e1 * t83 * t77;
t86 = -t106 * pkin(4) + t96;
t115 = t104 * t76 + t82 * t140;
t13 = -t101 * t25 + t149 * t20;
t26 = -t77 * t84 - t83 * t78;
t27 = -t158 * t78 + t84 * t75;
t46 = t83 * t84;
t47 = t158 * t84;
t114 = t52 * t26 - t51 * t27 + t35 * t47 - t36 * t46;
t113 = t158 * t27 + t83 * t26 - t77 * t46 + t75 * t47;
t112 = -t158 * t35 + t36 * t83 + t51 * t75 + t52 * t77;
t110 = -t59 * t76 + t121;
t108 = 0.2e1 * t158 * t36 + 0.2e1 * t35 * t83 + 0.2e1 * t51 * t77 - 0.2e1 * t52 * t75;
t95 = -t149 * pkin(4) - pkin(5);
t91 = t101 * pkin(4) + qJ(6);
t80 = t84 ^ 2;
t50 = t106 * t76 - t82 * t141;
t44 = -pkin(5) * t158 - t83 * qJ(6) + t86;
t31 = t75 * pkin(5) - t77 * qJ(6) - t144 + t97;
t17 = t46 * pkin(5) - t47 * qJ(6) + t40;
t16 = -t155 * qJD(4) + t129;
t15 = t59 * t141 - t136;
t12 = -t82 * pkin(5) - t13;
t11 = t82 * qJ(6) + t14;
t5 = -t26 * pkin(5) + t27 * qJ(6) - t47 * qJD(6) + t28;
t2 = -t76 * pkin(5) - t3;
t1 = t139 + t145;
t6 = [0, 0, 0, 0.2e1 * t105 * t142, 0.2e1 * (-t105 ^ 2 + t107 ^ 2) * qJD(2), 0, 0, 0, t105 * t138, t107 * t138, -0.2e1 * t42 * t82 + 0.2e1 * t110, 0.2e1 * t133 * t98 + 0.2e1 * t58 * t41 + 0.2e1 * t59 * t42, 0.2e1 * t100 * t84 * t78 - 0.2e1 * t80 * t131, -0.2e1 * t80 * t123 + t78 * t125, 0.2e1 * t120 * t106 - 0.2e1 * t82 * t135, -0.2e1 * t120 * t104 - 0.2e1 * t82 * t134, 0.2e1 * t82 * t76, 0.2e1 * t110 * t104 + 0.2e1 * t58 * t134 + 0.2e1 * t16 * t82 + 0.2e1 * t49 * t76, 0.2e1 * t121 * t106 - 0.2e1 * t58 * t135 + 0.2e1 * t15 * t82 - 0.2e1 * t155 * t76, 0.2e1 * t13 * t27 + 0.2e1 * t14 * t26 - 0.2e1 * t3 * t47 - 0.2e1 * t4 * t46, 0.2e1 * t13 * t3 + 0.2e1 * t14 * t4 + 0.2e1 * t40 * t28, -0.2e1 * t12 * t76 - 0.2e1 * t17 * t26 - 0.2e1 * t2 * t82 + 0.2e1 * t5 * t46, -0.2e1 * t1 * t46 + 0.2e1 * t11 * t26 - 0.2e1 * t12 * t27 + 0.2e1 * t2 * t47, 0.2e1 * t1 * t82 + 0.2e1 * t11 * t76 + 0.2e1 * t17 * t27 - 0.2e1 * t5 * t47, 0.2e1 * t11 * t1 + 0.2e1 * t12 * t2 + 0.2e1 * t17 * t5; 0, 0, 0, 0, 0, t142, -t143, 0, -pkin(7) * t142, pkin(7) * t143 (-t102 * t76 - t103 * t78) * pkin(2) (t102 * t42 - t103 * t41) * pkin(2), t84 * t123 + t78 * t146, qJD(4) * t125 + t150 * t78, t115, t50, 0, -t41 * t106 + t119 * t104 + (t104 * t58 - t118 * t106) * qJD(4), t41 * t104 + t119 * t106 + (t118 * t104 + t106 * t58) * qJD(4), -t13 * t77 - t14 * t75 + t158 * t4 - t3 * t83 + t114, -t13 * t35 + t14 * t36 + t28 * t86 - t3 * t51 + t4 * t52 + t40 * t97, -t158 * t5 + t17 * t75 - t44 * t26 + t31 * t46 - t35 * t82 - t51 * t76, t1 * t158 - t11 * t75 + t12 * t77 + t2 * t83 + t114, -t17 * t77 + t44 * t27 - t31 * t47 + t36 * t82 - t5 * t83 + t52 * t76, t1 * t52 + t11 * t36 + t12 * t35 + t17 * t31 + t2 * t51 + t5 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t131, 0.2e1 * t123, 0, 0, 0, t104 * t137, t106 * t137, t108, 0.2e1 * t86 * t97 + 0.2e1 * t130, -0.2e1 * t158 * t31 + 0.2e1 * t44 * t75, t108, -0.2e1 * t31 * t83 - 0.2e1 * t44 * t77, 0.2e1 * t44 * t31 + 0.2e1 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, 0, 0, 0, 0, t50, -t115, t113, -t13 * t75 + t14 * t77 + t158 * t3 + t4 * t83, t158 * t76 - t75 * t82, t113, t83 * t76 + t77 * t82, t1 * t83 + t11 * t77 + t12 * t75 - t158 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, 0, 0, 0, t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, 0, 0, 0, t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106 * t78 - t135, -t159, t76, t16, t15 (t101 * t26 + t149 * t27) * pkin(4) (t101 * t4 + t149 * t3) * pkin(4) (pkin(5) - t95) * t76 + t3, -qJD(6) * t46 + t91 * t26 - t95 * t27, t91 * t76 + t139 + 0.2e1 * t145, t11 * qJD(6) + t1 * t91 + t2 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t141, 0, -t94 * t140, t94 * t141 (-t101 * t75 - t149 * t77) * pkin(4) (t101 * t36 - t149 * t35) * pkin(4), -t35, qJD(6) * t158 - t91 * t75 + t95 * t77, t36, t52 * qJD(6) + t35 * t95 + t36 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t140, 0 (t101 * t77 - t149 * t75) * pkin(4), -t75, 0, t77, t75 * t95 + t77 * t91 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t91 * t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26, 0, t27, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t75, 0, -t77, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t27, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
