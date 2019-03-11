% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRPP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:29
% EndTime: 2019-03-08 22:58:36
% DurationCPUTime: 2.06s
% Computational Cost: add. (1357->257), mult. (3637->419), div. (0->0), fcn. (3058->8), ass. (0->139)
t81 = cos(qJ(4));
t70 = qJD(4) * t81;
t79 = sin(qJ(3));
t129 = t79 * t70;
t82 = cos(qJ(3));
t142 = t82 * qJD(3);
t78 = sin(qJ(4));
t177 = t78 * t142 + t129;
t145 = qJD(4) * t82;
t128 = t81 * t145;
t68 = t79 * qJD(3);
t176 = t78 * t68 - t128;
t175 = -0.4e1 * t79;
t77 = -pkin(4) - qJ(6);
t174 = t77 * t79;
t149 = sin(pkin(6));
t117 = sin(qJ(2)) * t149;
t108 = qJD(2) * t117;
t150 = cos(pkin(6));
t89 = t82 * t117 + t150 * t79;
t170 = -qJD(4) * t89 + t108;
t106 = t149 * cos(qJ(2));
t35 = t79 * t117 - t150 * t82;
t97 = qJD(2) * t106;
t23 = -qJD(3) * t35 + t82 * t97;
t172 = -qJD(4) * t106 + t23;
t7 = -t170 * t81 + t172 * t78;
t8 = t170 * t78 + t172 * t81;
t173 = t7 * t78 + t8 * t81;
t151 = t78 * qJ(5);
t105 = -t81 * pkin(4) - t151;
t69 = qJD(5) * t81;
t171 = t105 * qJD(4) + t69;
t103 = t77 * t81 - t151;
t75 = t81 ^ 2;
t154 = t78 ^ 2 - t75;
t116 = t154 * qJD(4);
t167 = pkin(9) * t82;
t109 = pkin(3) * t79 - t167;
t46 = t109 * qJD(3);
t166 = t79 * pkin(9);
t110 = -t82 * pkin(3) - t166;
t51 = -pkin(2) + t110;
t157 = t78 * t46 + t51 * t70;
t130 = t78 * t145;
t92 = t81 * t68 + t130;
t16 = t92 * pkin(8) - t157;
t138 = qJ(5) * qJD(4);
t120 = t79 * t138;
t64 = pkin(8) * t142;
t113 = t177 * pkin(4) + t78 * t120 + t64;
t139 = qJ(5) * qJD(3);
t121 = t82 * t139;
t15 = (-qJD(5) * t79 - t121) * t81 + t113;
t47 = -pkin(3) + t105;
t169 = (t47 * t79 + t167) * qJD(4) - t15;
t168 = pkin(5) + pkin(9);
t152 = qJ(5) * t81;
t102 = qJ(6) * t78 - t152;
t144 = qJD(6) * t78;
t10 = t102 * t142 + (t144 + (qJ(6) * qJD(4) - qJD(5)) * t81) * t79 + t113;
t164 = t10 * t78;
t163 = t10 * t81;
t161 = t78 * t79;
t160 = t78 * t82;
t159 = t79 * t81;
t158 = t81 * t82;
t156 = pkin(4) * t161 + t79 * pkin(8);
t62 = pkin(8) * t158;
t155 = t78 * t51 + t62;
t74 = t79 ^ 2;
t153 = -t82 ^ 2 + t74;
t148 = qJD(3) * t78;
t147 = qJD(3) * t81;
t146 = qJD(4) * t78;
t143 = t78 * qJD(5);
t141 = t82 * qJD(5);
t140 = t82 * qJD(6);
t137 = qJ(5) * qJD(5);
t136 = pkin(5) * t158;
t61 = pkin(8) * t160;
t135 = -0.2e1 * pkin(2) * qJD(3);
t134 = -0.2e1 * pkin(3) * qJD(4);
t133 = pkin(4) * t68;
t132 = pkin(9) * t146;
t131 = t79 * t146;
t40 = -pkin(3) + t103;
t127 = t40 * t70;
t124 = t78 * t70;
t123 = t79 * t142;
t122 = t81 * t142;
t119 = t79 * t139;
t118 = t81 * t51 - t61;
t115 = t153 * qJD(3);
t114 = 0.2e1 * t123;
t112 = -t176 * pkin(8) + t51 * t146 - t81 * t46;
t111 = t78 * t122;
t29 = t82 * qJ(5) - t155;
t72 = t82 * pkin(4);
t30 = -t118 + t72;
t104 = t29 * t78 + t30 * t81;
t26 = -t78 * t106 + t89 * t81;
t101 = t8 * qJ(5) + t26 * qJD(5);
t100 = -t81 * qJD(6) - t143;
t99 = qJD(3) * t106;
t24 = t89 * qJD(3) + t79 * t97;
t12 = t24 * t78 + t35 * t70;
t13 = -t35 * t146 + t24 * t81;
t65 = pkin(4) * t146;
t27 = t102 * qJD(4) + t100 + t65;
t96 = t40 * t146 - t27 * t81;
t25 = t81 * t106 + t89 * t78;
t95 = 0.2e1 * t35 * t24 + 0.2e1 * t25 * t7 + 0.2e1 * t26 * t8;
t28 = t102 * t79 + t156;
t94 = -qJD(4) * t28 - t40 * t142;
t91 = -t26 * t146 + t25 * t70 + t173;
t90 = -pkin(5) * t131 + t112;
t11 = -t119 + t16 + t141;
t14 = t112 - t133;
t87 = t104 * qJD(4) - t11 * t81 + t14 * t78;
t33 = -t79 * t152 + t156;
t34 = -t81 * t138 - t143 + t65;
t86 = -qJD(4) * t33 - t34 * t79 + (-t47 * t82 + t166) * qJD(3);
t85 = 0.2e1 * t119 - 0.2e1 * t141 - t16;
t1 = (t35 * t148 + t7) * t82 + (-qJD(3) * t25 + t12) * t79;
t2 = (-t35 * t147 - t8) * t82 + (qJD(3) * t26 - t13) * t79;
t84 = -t26 * t129 + t25 * t122 + t7 * t159 + (-t26 * t142 + (-qJD(4) * t25 - t8) * t79) * t78;
t83 = 0.2e1 * qJD(5);
t66 = pkin(9) * t70;
t53 = t168 * t81;
t52 = t168 * t78;
t45 = pkin(5) * t70 + t66;
t44 = t168 * t146;
t36 = t122 - t131;
t22 = -pkin(5) * t161 - t29;
t19 = t82 * qJ(6) + t61 + t72 + (pkin(5) * t79 - t51) * t81;
t9 = -t141 + (-pkin(5) * t159 - t61) * qJD(4) + (-pkin(5) * t160 + (-pkin(8) * t81 + qJ(5)) * t79) * qJD(3) + t157;
t5 = t140 + (t136 + t174) * qJD(3) + t90;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, t95; 0, 0, -t108, -t97, 0, 0, 0, 0, 0, -t82 * t108 - t79 * t99, t79 * t108 - t82 * t99, 0, 0, 0, 0, 0, t1, -t2, t84, -t1, t2, -t26 * t11 + t25 * t14 + t35 * t15 + t24 * t33 - t8 * t29 + t7 * t30, t84, t2, t1, t35 * t10 + t7 * t19 + t8 * t22 + t24 * t28 + t25 * t5 + t26 * t9; 0, 0, 0, 0, t114, -0.2e1 * t115, 0, 0, 0, t79 * t135, t82 * t135, 0.2e1 * t75 * t123 - 0.2e1 * t74 * t124, t111 * t175 + 0.2e1 * t74 * t116, 0.2e1 * t79 * t130 + 0.2e1 * t153 * t147, -0.2e1 * t78 * t115 + 0.2e1 * t79 * t128, -0.2e1 * t123, 0.2e1 * t112 * t82 + 0.2e1 * t118 * t68 + 0.2e1 * (t78 * t114 + t74 * t70) * pkin(8), -0.2e1 * t16 * t82 - 0.2e1 * t155 * t68 + 0.2e1 * (t81 * t114 - t74 * t146) * pkin(8), 0.2e1 * t104 * t142 + 0.2e1 * (t11 * t78 + t14 * t81 + (t29 * t81 - t30 * t78) * qJD(4)) * t79, 0.2e1 * (-t33 * t148 - t14) * t82 + 0.2e1 * (qJD(3) * t30 - t15 * t78 - t33 * t70) * t79, 0.2e1 * (-t33 * t147 + t11) * t82 + 0.2e1 * (-qJD(3) * t29 + t33 * t146 - t15 * t81) * t79, 0.2e1 * t29 * t11 + 0.2e1 * t30 * t14 + 0.2e1 * t33 * t15, 0.2e1 * (t19 * t81 - t22 * t78) * t142 + 0.2e1 * (t5 * t81 - t78 * t9 + (-t19 * t78 - t22 * t81) * qJD(4)) * t79, 0.2e1 * (-t28 * t147 - t9) * t82 + 0.2e1 * (qJD(3) * t22 + t28 * t146 - t163) * t79, 0.2e1 * (t28 * t148 + t5) * t82 + 0.2e1 * (-qJD(3) * t19 + t28 * t70 + t164) * t79, 0.2e1 * t28 * t10 + 0.2e1 * t19 * t5 + 0.2e1 * t22 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, 0, 0, 0, 0, 0, -t13, t12, t91, t13, -t12, t24 * t47 + t35 * t34 + ((t25 * t81 - t26 * t78) * qJD(4) + t173) * pkin(9), t91, -t12, -t13, t24 * t40 + t25 * t45 - t26 * t44 + t35 * t27 + t7 * t52 + t8 * t53; 0, 0, 0, 0, 0, 0, t142, -t68, 0, -t64, pkin(8) * t68, -t79 * t116 + t111, t124 * t175 - t154 * t142, t176, t92, 0 (pkin(9) * t158 + (-pkin(3) * t81 + pkin(8) * t78) * t79) * qJD(4) + (t110 * t78 - t62) * qJD(3) (pkin(8) * t159 + t109 * t78) * qJD(4) + (t110 * t81 + t61) * qJD(3), t87, -t169 * t81 + t86 * t78, t169 * t78 + t86 * t81, t87 * pkin(9) + t15 * t47 + t33 * t34 (t52 * t142 + t45 * t79 + t9 + (-t53 * t79 + t19) * qJD(4)) * t81 + (-t53 * t142 + t44 * t79 + t5 + (-t52 * t79 - t22) * qJD(4)) * t78, -t164 + t44 * t82 + t94 * t81 + (qJD(3) * t53 + t96) * t79, -t163 + t45 * t82 + (-qJD(3) * t52 + t127) * t79 + (t27 * t79 - t94) * t78, t10 * t40 + t19 * t45 - t22 * t44 + t28 * t27 + t5 * t52 + t9 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t124, -0.2e1 * t116, 0, 0, 0, t78 * t134, t81 * t134, 0, -0.2e1 * t47 * t146 + 0.2e1 * t34 * t81, -0.2e1 * t34 * t78 - 0.2e1 * t47 * t70, 0.2e1 * t47 * t34, -0.2e1 * t44 * t81 + 0.2e1 * t45 * t78 + 0.2e1 * (t52 * t81 - t53 * t78) * qJD(4), -0.2e1 * t27 * t78 - 0.2e1 * t127, 0.2e1 * t96, 0.2e1 * t40 * t27 - 0.2e1 * t53 * t44 + 0.2e1 * t52 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, t7, t8, -t7 * pkin(4) + t101, 0, t8, -t7, -t25 * qJD(6) + t7 * t77 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t177, t68, -t112, t16 (-pkin(4) * t142 - t120) * t81 + (-t121 + (pkin(4) * qJD(4) - qJD(5)) * t79) * t78, t112 - 0.2e1 * t133, t85, -t14 * pkin(4) - t11 * qJ(5) - t29 * qJD(5), t103 * t142 + ((-t77 * t78 - t152) * qJD(4) + t100) * t79, -pkin(5) * t177 + t85, -0.2e1 * t140 + (-t136 - 0.2e1 * t174) * qJD(3) - t90, t9 * qJ(5) + t22 * qJD(5) - t19 * qJD(6) + t5 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t146, 0, -t66, t132, t171, t66, -t132, t171 * pkin(9), qJD(4) * t103 - t144 + t69, -t44, -t45, -t44 * qJ(5) + t53 * qJD(5) - t52 * qJD(6) + t45 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0.2e1 * t137, 0, t83, 0.2e1 * qJD(6), -0.2e1 * t77 * qJD(6) + 0.2e1 * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t68, 0, t14, t36, 0, -t68, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, t66, t70, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, t68, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, 0, 0, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
