% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:17
% EndTime: 2019-03-09 04:52:21
% DurationCPUTime: 1.81s
% Computational Cost: add. (1194->230), mult. (2509->357), div. (0->0), fcn. (1772->4), ass. (0->129)
t62 = cos(qJ(4));
t118 = t62 * qJD(5);
t60 = sin(qJ(4));
t128 = t60 * qJ(5);
t143 = pkin(4) + pkin(5);
t77 = t143 * t62 + t128;
t153 = t77 * qJD(4) - t118;
t114 = qJ(5) * qJD(3);
t63 = cos(qJ(3));
t50 = t63 * t114;
t61 = sin(qJ(3));
t150 = t61 * qJD(5) + t50;
t152 = t143 * qJD(3);
t151 = 0.2e1 * t150;
t65 = -pkin(1) - pkin(7);
t130 = qJ(5) * t62;
t74 = -t143 * t60 + t130;
t149 = t65 + t74;
t56 = t60 ^ 2;
t58 = t62 ^ 2;
t134 = t56 - t58;
t94 = t134 * qJD(4);
t121 = t60 * qJD(5);
t82 = pkin(4) * t60 - t130;
t18 = t82 * qJD(4) - t121;
t138 = t63 * t18;
t83 = t62 * pkin(4) + t128;
t35 = -pkin(3) - t83;
t140 = t35 * t61;
t141 = t63 * pkin(8);
t75 = -t65 + t82;
t16 = t75 * t63;
t148 = (t140 + t141) * qJD(3) - qJD(4) * t16 - t138;
t117 = t62 * qJD(6);
t125 = qJD(4) * t60;
t136 = pkin(8) - qJ(6);
t17 = -t136 * t125 - t117;
t120 = t60 * qJD(6);
t40 = t136 * t62;
t19 = qJD(4) * t40 - t120;
t39 = t136 * t60;
t147 = -(t39 * t62 - t40 * t60) * qJD(4) - t17 * t62 - t19 * t60;
t145 = t83 * qJD(4) - t118;
t144 = 2 * qJD(2);
t66 = 0.2e1 * qJD(5);
t142 = pkin(8) * t61;
t139 = t61 * t65;
t137 = t63 * t65;
t135 = t61 * t118 + t62 * t50;
t133 = t56 + t58;
t57 = t61 ^ 2;
t59 = t63 ^ 2;
t132 = t57 - t59;
t131 = t57 + t59;
t129 = qJ(6) * t63;
t127 = qJD(3) * t60;
t126 = qJD(3) * t62;
t124 = qJD(4) * t61;
t53 = qJD(4) * t62;
t123 = qJD(4) * t63;
t122 = qJD(4) * t65;
t119 = t61 * qJD(3);
t116 = t63 * qJD(3);
t115 = qJ(2) * qJD(3);
t113 = -0.2e1 * pkin(3) * qJD(4);
t43 = t60 * t139;
t44 = t62 * t139;
t88 = pkin(3) * t63 + t142;
t33 = t88 * qJD(3) + qJD(2);
t87 = t61 * pkin(3) - t141;
t34 = qJ(2) + t87;
t99 = t65 * t116;
t112 = t60 * t33 + t34 * t53 + t62 * t99;
t32 = t60 * t34;
t13 = t61 * qJ(5) + t32 + t44;
t111 = pkin(4) * t116;
t110 = pkin(8) * t125;
t109 = pkin(8) * t53;
t108 = t60 * t124;
t107 = t60 * t123;
t106 = t60 * t122;
t105 = t62 * t123;
t104 = t62 * t122;
t103 = t60 * t116;
t102 = t60 * t53;
t101 = t62 * t119;
t100 = t61 * t116;
t98 = t61 * t114;
t97 = qJ(5) * t123;
t96 = t63 * t133;
t95 = -t62 * t34 + t43;
t93 = t132 * qJD(3);
t92 = 0.2e1 * t100;
t91 = t61 * t104 + t34 * t125 - t62 * t33 + t60 * t99;
t90 = t60 * t101;
t3 = -t149 * t119 - t153 * t63;
t31 = pkin(3) + t77;
t89 = -t31 * t123 + t3;
t10 = t60 * t129 + t13;
t9 = t43 + (-t34 - t129) * t62 - t143 * t61;
t85 = t10 * t62 + t60 * t9;
t84 = -t10 * t60 + t62 * t9;
t14 = -t61 * pkin(4) + t95;
t81 = t13 * t62 + t14 * t60;
t80 = t13 * t60 - t14 * t62;
t73 = -qJ(6) * t125 + t117;
t27 = t61 * t53 + t103;
t6 = t61 * t106 - t112;
t8 = -t75 * t119 + t145 * t63;
t72 = -t8 + (t35 * t63 - t142) * qJD(4);
t71 = -qJ(6) * t101 - t91;
t12 = t149 * t63;
t15 = t74 * qJD(4) + t121;
t70 = -qJD(4) * t12 + t31 * t119 - t63 * t15;
t4 = -t6 + t150;
t5 = t91 - t111;
t68 = -t80 * qJD(4) + t4 * t62 + t5 * t60;
t67 = qJ(6) * t105 + t63 * t120 + (-qJ(6) * qJD(3) - t122) * t61 * t60 + t112;
t55 = qJ(5) * t66;
t28 = -t60 * t119 + t105;
t26 = t131 * t53;
t25 = t101 + t107;
t24 = -t62 * t116 + t108;
t23 = t131 * t125;
t22 = qJD(3) * t96;
t11 = (-0.1e1 + t133) * t92;
t2 = t67 + t150;
t1 = (-t73 - t152) * t63 - t71;
t7 = [0, 0, 0, 0, t144, qJ(2) * t144, -0.2e1 * t100, 0.2e1 * t93, 0, 0, 0, 0.2e1 * qJD(2) * t61 + 0.2e1 * t63 * t115, 0.2e1 * qJD(2) * t63 - 0.2e1 * t61 * t115, -0.2e1 * t58 * t100 - 0.2e1 * t59 * t102, 0.2e1 * t59 * t94 + 0.4e1 * t63 * t90, -0.2e1 * t61 * t107 - 0.2e1 * t132 * t126, -0.2e1 * t61 * t105 + 0.2e1 * t60 * t93, t92, -0.2e1 * t59 * t104 - 0.2e1 * t91 * t61 + 0.2e1 * (-t95 + 0.2e1 * t43) * t116, 0.2e1 * t59 * t106 + 0.2e1 * t6 * t61 + 0.2e1 * (-t32 + t44) * t116, 0.2e1 * (-t16 * t127 - t5) * t61 + 0.2e1 * (-qJD(3) * t14 + t16 * t53 + t8 * t60) * t63, 0.2e1 * t80 * t119 + 0.2e1 * (-t81 * qJD(4) - t4 * t60 + t5 * t62) * t63, 0.2e1 * (t16 * t126 + t4) * t61 + 0.2e1 * (qJD(3) * t13 + t16 * t125 - t8 * t62) * t63, 0.2e1 * t13 * t4 + 0.2e1 * t14 * t5 + 0.2e1 * t16 * t8, 0.2e1 * (t12 * t127 - t1) * t61 + 0.2e1 * (-qJD(3) * t9 - t12 * t53 - t3 * t60) * t63, 0.2e1 * (-t12 * t126 + t2) * t61 + 0.2e1 * (qJD(3) * t10 - t12 * t125 + t3 * t62) * t63, 0.2e1 * t84 * t119 + 0.2e1 * (t85 * qJD(4) - t1 * t62 + t2 * t60) * t63, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t12 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t23, -t26, 0, -t23 (t81 * qJD(3) - t8) * t63 + (qJD(3) * t16 + t68) * t61, -t26, -t23, 0 (qJD(3) * t85 + t3) * t63 + (-qJD(3) * t12 + qJD(4) * t84 + t1 * t60 + t2 * t62) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t116, 0, -t65 * t119, -t99, -t63 * t94 - t90, -0.4e1 * t63 * t102 + t134 * t119, t27, -t24, 0 (-t60 * t137 - t88 * t62) * qJD(4) + (t87 * t60 - t44) * qJD(3) (-t62 * t137 + t88 * t60) * qJD(4) + (-t62 * t141 + (pkin(3) * t62 + t60 * t65) * t61) * qJD(3), -t148 * t60 + t72 * t62, t68, t148 * t62 + t72 * t60, t68 * pkin(8) + t16 * t18 + t8 * t35, -t39 * t116 - t19 * t61 + t70 * t60 + t89 * t62, t40 * t116 + t17 * t61 + t89 * t60 - t70 * t62 (t39 * t119 - t19 * t63 - t2 + (t40 * t63 - t9) * qJD(4)) * t62 + (-t40 * t119 + t17 * t63 - t1 + (t39 * t63 + t10) * qJD(4)) * t60, t1 * t39 + t10 * t17 + t12 * t15 + t9 * t19 + t2 * t40 + t3 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t116, 0, 0, 0, 0, 0, -t25, -t28, -t25, t22, t28, -t138 + (pkin(8) * t96 + t140) * qJD(3), -t25, t28, -t22 (t15 + (t39 * t60 + t40 * t62) * qJD(3)) * t63 + (-qJD(3) * t31 - t147) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t102, -0.2e1 * t94, 0, 0, 0, t60 * t113, t62 * t113, 0.2e1 * t35 * t125 - 0.2e1 * t18 * t62, 0, -0.2e1 * t18 * t60 - 0.2e1 * t35 * t53, 0.2e1 * t35 * t18, -0.2e1 * t31 * t125 + 0.2e1 * t15 * t62, 0.2e1 * t15 * t60 + 0.2e1 * t31 * t53, 0.2e1 * t147, 0.2e1 * t31 * t15 + 0.2e1 * t40 * t17 + 0.2e1 * t39 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t28, t116, -t91, t6, -t91 + 0.2e1 * t111 (pkin(4) * t119 - t97) * t62 + (t98 + (pkin(4) * qJD(4) - qJD(5)) * t63) * t60, -t6 + t151, -t5 * pkin(4) + t4 * qJ(5) + t13 * qJD(5) (t73 + 0.2e1 * t152) * t63 + t71, t67 + t151 (-t119 * t143 + t97) * t62 + (-t98 + (-qJD(4) * t143 + qJD(5)) * t63) * t60, t2 * qJ(5) + t10 * qJD(5) - t1 * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t24, -t27, 0, -t24, -t27 * pkin(4) - qJ(5) * t108 + t135, -t27, -t24, 0, -t103 * t143 - t77 * t124 + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t125, 0, -t109, t110, -t109, -t145, -t110, -t145 * pkin(8), -t19, t17, t153, t17 * qJ(5) + t40 * qJD(5) - t143 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t55, 0, t66, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t25, 0, t5, -t116, 0, t25, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t109, 0, 0, -t53, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t25, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, t53, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
