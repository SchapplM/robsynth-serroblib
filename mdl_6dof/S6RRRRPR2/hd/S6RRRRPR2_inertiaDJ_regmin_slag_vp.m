% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:59:14
% EndTime: 2019-03-09 21:59:20
% DurationCPUTime: 2.49s
% Computational Cost: add. (6233->240), mult. (13963->409), div. (0->0), fcn. (13732->10), ass. (0->156)
t120 = sin(pkin(11));
t121 = cos(pkin(11));
t172 = t120 ^ 2 + t121 ^ 2;
t202 = t172 * qJD(5);
t203 = 0.2e1 * t202;
t123 = sin(qJ(4));
t124 = sin(qJ(3));
t128 = cos(qJ(3));
t127 = cos(qJ(4));
t168 = qJD(4) * t127;
t176 = t124 * t127;
t201 = ((t123 * t128 + t176) * qJD(3) + t124 * t168) * pkin(2);
t112 = t128 * pkin(2) + pkin(3);
t184 = pkin(2) * qJD(3);
t165 = t128 * t184;
t107 = t123 * t124 * pkin(2);
t166 = t124 * t184;
t174 = qJD(4) * t107 + t123 * t166;
t67 = -(qJD(4) * t112 + t165) * t127 + t174;
t66 = qJD(5) - t67;
t200 = t172 * t66;
t163 = pkin(3) * t168;
t106 = qJD(5) + t163;
t199 = t172 * t106;
t125 = sin(qJ(2));
t196 = pkin(7) + pkin(8);
t162 = qJD(2) * t196;
t94 = t125 * t162;
t129 = cos(qJ(2));
t95 = t129 * t162;
t146 = -t124 * t95 - t128 * t94;
t101 = t196 * t125;
t102 = t196 * t129;
t143 = -t128 * t101 - t124 * t102;
t93 = t124 * t129 + t128 * t125;
t60 = -t93 * pkin(9) + t143;
t198 = t60 * qJD(3) + t146;
t194 = t121 * pkin(5);
t117 = t121 * pkin(10);
t193 = t127 * pkin(3);
t137 = t93 * qJD(2);
t130 = -pkin(9) * t137 + t198;
t169 = qJD(4) * t123;
t131 = -t93 * qJD(3) - t137;
t140 = t124 * t125 - t128 * t129;
t70 = t123 * t93 + t127 * t140;
t73 = (-qJD(2) - qJD(3)) * t140;
t33 = -t70 * qJD(4) + t123 * t131 + t127 * t73;
t183 = t120 * t33;
t142 = t124 * t101 - t128 * t102;
t50 = t142 * qJD(3) + t124 * t94 - t128 * t95;
t132 = -t73 * pkin(9) + t50;
t37 = t127 * t132;
t61 = -t140 * pkin(9) - t142;
t12 = pkin(5) * t183 + t123 * t130 + t61 * t168 + t60 * t169 - t37;
t71 = -t123 * t140 + t127 * t93;
t182 = t120 * t71;
t42 = t123 * t61 - t127 * t60;
t32 = pkin(5) * t182 + t42;
t122 = sin(qJ(6));
t126 = cos(qJ(6));
t175 = t126 * t120;
t92 = t122 * t121 + t175;
t87 = t92 * qJD(6);
t177 = t122 * t120;
t91 = -t126 * t121 + t177;
t192 = t12 * t91 + t32 * t87;
t86 = t91 * qJD(6);
t191 = t12 * t92 - t32 * t86;
t34 = t71 * qJD(4) + t123 * t73 - t127 * t131;
t171 = qJD(2) * t125;
t114 = pkin(2) * t171;
t65 = -t131 * pkin(3) + t114;
t19 = t34 * pkin(4) - t33 * qJ(5) - t71 * qJD(5) + t65;
t20 = -t123 * t132 - t127 * t130 - t60 * t168 + t61 * t169;
t8 = t120 * t19 - t121 * t20;
t113 = -t129 * pkin(2) - pkin(1);
t81 = t140 * pkin(3) + t113;
t41 = t70 * pkin(4) - t71 * qJ(5) + t81;
t43 = t123 * t60 + t127 * t61;
t25 = t120 * t41 + t121 * t43;
t68 = t112 * t169 + t201;
t85 = -t127 * t112 - pkin(4) + t107;
t80 = t85 - t194;
t190 = t68 * t91 + t80 * t87;
t189 = t68 * t92 - t80 * t86;
t164 = pkin(3) * t169;
t110 = -pkin(4) - t194;
t98 = t110 - t193;
t187 = t91 * t164 + t98 * t87;
t186 = t92 * t164 - t98 * t86;
t181 = t121 * t33;
t170 = qJD(2) * t129;
t21 = t43 * qJD(4) - t37 + ((-t124 * t170 - t128 * t171) * pkin(9) + t198) * t123;
t180 = t21 * t121;
t179 = t68 * t121;
t178 = qJD(6) * t71;
t167 = -0.2e1 * pkin(1) * qJD(2);
t7 = t120 * t20 + t121 * t19;
t24 = -t120 * t43 + t121 * t41;
t160 = (-pkin(3) - t112) * qJD(4);
t109 = t123 * pkin(3) + qJ(5);
t157 = t172 * t109;
t155 = t121 * t164;
t3 = -t7 * t120 + t8 * t121;
t154 = t21 * t71 + t33 * t42;
t153 = -t120 * t24 + t121 * t25;
t22 = t70 * pkin(5) - t71 * t117 + t24;
t23 = -pkin(10) * t182 + t25;
t152 = t122 * t23 - t126 * t22;
t151 = t122 * t22 + t126 * t23;
t84 = pkin(2) * t176 + t123 * t112 + qJ(5);
t76 = (-pkin(10) - t84) * t120;
t77 = t121 * t84 + t117;
t150 = t122 * t77 - t126 * t76;
t149 = t122 * t76 + t126 * t77;
t89 = (-pkin(10) - t109) * t120;
t90 = t121 * t109 + t117;
t148 = t122 * t90 - t126 * t89;
t147 = t122 * t89 + t126 * t90;
t100 = t121 * qJ(5) + t117;
t99 = (-pkin(10) - qJ(5)) * t120;
t145 = t126 * t100 + t122 * t99;
t144 = t122 * t100 - t126 * t99;
t138 = t113 * t93;
t136 = -pkin(4) * t33 - qJ(5) * t34 - qJD(5) * t70;
t135 = t33 * t85 - t34 * t84 - t66 * t70 + t68 * t71;
t111 = -pkin(4) - t193;
t134 = -t106 * t70 - t109 * t34 + t111 * t33 + t71 * t164;
t103 = t120 * t164;
t79 = t110 * t87;
t78 = t110 * t86;
t72 = -0.2e1 * t92 * t86;
t62 = t68 * t120;
t54 = -t92 * qJD(5) - t145 * qJD(6);
t53 = t91 * qJD(5) + t144 * qJD(6);
t49 = -t143 * qJD(3) - t146;
t48 = 0.2e1 * t86 * t91 - 0.2e1 * t92 * t87;
t47 = -t147 * qJD(6) - t92 * t106;
t46 = t148 * qJD(6) + t91 * t106;
t45 = t91 * t71;
t44 = t92 * t71;
t29 = -t149 * qJD(6) - t92 * t66;
t28 = t150 * qJD(6) + t91 * t66;
t27 = -t91 * t34 - t87 * t70;
t26 = t92 * t34 - t86 * t70;
t18 = t21 * t120;
t14 = t33 * t175 - t177 * t178 + (t122 * t33 + t126 * t178) * t121;
t13 = -t92 * t178 - t91 * t33;
t9 = t13 * t92 + t45 * t86;
t6 = -pkin(10) * t183 + t8;
t5 = t34 * pkin(5) - pkin(10) * t181 + t7;
t4 = -t13 * t91 - t92 * t14 + t86 * t44 + t45 * t87;
t2 = -t151 * qJD(6) - t122 * t6 + t126 * t5;
t1 = t152 * qJD(6) - t122 * t5 - t126 * t6;
t10 = [0, 0, 0, 0.2e1 * t125 * t170, 0.2e1 * (-t125 ^ 2 + t129 ^ 2) * qJD(2), 0, 0, 0, t125 * t167, t129 * t167, 0.2e1 * t93 * t73, 0.2e1 * t93 * t131 - 0.2e1 * t73 * t140, 0, 0, 0, 0.2e1 * qJD(3) * t138 + 0.2e1 * (t125 * pkin(2) * t140 + t138) * qJD(2), 0.2e1 * t113 * t73 + 0.2e1 * t93 * t114, 0.2e1 * t71 * t33, -0.2e1 * t33 * t70 - 0.2e1 * t34 * t71, 0, 0, 0, 0.2e1 * t81 * t34 + 0.2e1 * t65 * t70, 0.2e1 * t81 * t33 + 0.2e1 * t65 * t71, 0.2e1 * t120 * t154 + 0.2e1 * t24 * t34 + 0.2e1 * t7 * t70, 0.2e1 * t121 * t154 - 0.2e1 * t25 * t34 - 0.2e1 * t8 * t70, 0.2e1 * (-t24 * t33 - t7 * t71) * t121 + 0.2e1 * (-t25 * t33 - t71 * t8) * t120, 0.2e1 * t21 * t42 + 0.2e1 * t24 * t7 + 0.2e1 * t25 * t8, -0.2e1 * t45 * t13, -0.2e1 * t13 * t44 + 0.2e1 * t14 * t45, 0.2e1 * t13 * t70 - 0.2e1 * t45 * t34, -0.2e1 * t14 * t70 - 0.2e1 * t44 * t34, 0.2e1 * t70 * t34, 0.2e1 * t12 * t44 + 0.2e1 * t32 * t14 - 0.2e1 * t152 * t34 + 0.2e1 * t2 * t70, 0.2e1 * t1 * t70 - 0.2e1 * t12 * t45 + 0.2e1 * t32 * t13 - 0.2e1 * t151 * t34; 0, 0, 0, 0, 0, t170, -t171, 0, -pkin(7) * t170, pkin(7) * t171, 0, 0, t73, t131, 0, t50, t49, 0, 0, t33, -t34, 0, -t21, t20, t120 * t135 - t180, t121 * t135 + t18, t3, t21 * t85 + t42 * t68 + (t25 * t66 + t8 * t84) * t121 + (-t24 * t66 - t7 * t84) * t120, t9, t4, t26, t27, 0, t80 * t14 - t150 * t34 + t29 * t70 + t68 * t44 + t192, t80 * t13 - t149 * t34 + t28 * t70 - t68 * t45 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t166, -0.2e1 * t165, 0, 0, 0, 0, 0, -0.2e1 * t68, 0.2e1 * t67, -0.2e1 * t179, 0.2e1 * t62, 0.2e1 * t200, 0.2e1 * t200 * t84 + 0.2e1 * t85 * t68, t72, t48, 0, 0, 0, 0.2e1 * t190, 0.2e1 * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t131, 0, t50, t49, 0, 0, t33, -t34, 0, -t21, t20, t120 * t134 - t180, t121 * t134 + t18, t3, t106 * t153 + t109 * t3 + t21 * t111 + t164 * t42, t9, t4, t26, t27, 0, t98 * t14 - t148 * t34 + t164 * t44 + t47 * t70 + t192, t98 * t13 - t147 * t34 - t164 * t45 + t46 * t70 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, -t165, 0, 0, 0, 0, 0, t123 * t160 - t201 (t160 - t165) * t127 + t174 (-t68 - t164) * t121, t103 + t62, t199 + t200, t68 * t111 + t157 * t66 + t164 * t85 + t199 * t84, t72, t48, 0, 0, 0, t187 + t190, t186 + t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t164, -0.2e1 * t163, -0.2e1 * t155, 0.2e1 * t103, 0.2e1 * t199, 0.2e1 * t106 * t157 + 0.2e1 * t111 * t164, t72, t48, 0, 0, 0, 0.2e1 * t187, 0.2e1 * t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t34, 0, -t21, t20, t120 * t136 - t180, t121 * t136 + t18, t3, -t21 * pkin(4) + qJ(5) * t3 + qJD(5) * t153, t9, t4, t26, t27, 0, t110 * t14 - t144 * t34 + t54 * t70 + t192, t110 * t13 - t145 * t34 + t53 * t70 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t67, -t179, t62, t202 + t200, -t68 * pkin(4) + qJ(5) * t200 + t202 * t84, t72, t48, 0, 0, 0, t79 + t190, -t78 + t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, -t163, -t155, t103, t202 + t199, -pkin(4) * t164 + qJ(5) * t199 + t109 * t202, t72, t48, 0, 0, 0, t79 + t187, -t78 + t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, qJ(5) * t203, t72, t48, 0, 0, 0, 0.2e1 * t79, -0.2e1 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, t181, 0, t21, 0, 0, 0, 0, 0, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, 0, 0, 0, t87, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, 0, 0, 0, 0, 0, t87, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, t34, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t87, 0, t29, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t87, 0, t47, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t87, 0, t54, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
