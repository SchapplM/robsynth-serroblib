% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:17
% EndTime: 2021-01-16 02:21:27
% DurationCPUTime: 2.13s
% Computational Cost: add. (2287->291), mult. (6128->407), div. (0->0), fcn. (4647->10), ass. (0->168)
t108 = sin(pkin(11));
t110 = cos(pkin(11));
t113 = sin(qJ(3));
t116 = cos(qJ(3));
t86 = t108 * t116 + t110 * t113;
t214 = t86 * qJD(2);
t218 = qJD(6) + t214;
t112 = sin(qJ(6));
t115 = cos(qJ(6));
t182 = t110 * t116;
t158 = qJD(2) * t182;
t171 = qJD(2) * t113;
t78 = t108 * t171 - t158;
t59 = qJD(3) * t112 - t115 * t78;
t220 = t218 * t59;
t152 = t112 * t218;
t165 = qJD(2) * qJD(3);
t155 = t116 * t165;
t156 = t113 * t165;
t93 = t108 * t156;
t70 = t110 * t155 - t93;
t139 = t115 * t70 - t152 * t218;
t114 = sin(qJ(2));
t109 = sin(pkin(6));
t173 = qJD(1) * t109;
t162 = t114 * t173;
t169 = qJD(3) * t113;
t219 = pkin(3) * t169 - t162;
t217 = -0.2e1 * qJD(3);
t117 = cos(qJ(2));
t161 = t117 * t173;
t194 = qJ(4) + pkin(8);
t154 = qJD(3) * t194;
t74 = t116 * qJD(4) - t113 * t154;
t75 = -t113 * qJD(4) - t116 * t154;
t191 = t108 * t74 - t110 * t75 - t86 * t161;
t211 = -t108 * t113 + t182;
t190 = t108 * t75 + t110 * t74 - t211 * t161;
t77 = t214 ^ 2;
t215 = -t78 ^ 2 - t77;
t89 = qJD(2) * pkin(8) + t162;
t151 = qJ(4) * qJD(2) + t89;
t111 = cos(pkin(6));
t172 = qJD(1) * t111;
t157 = t113 * t172;
t54 = t151 * t116 + t157;
t47 = t108 * t54;
t97 = t116 * t172;
t53 = -t151 * t113 + t97;
t24 = t110 * t53 - t47;
t177 = -qJD(5) + t24;
t213 = -qJD(6) + t218;
t168 = qJD(3) * t116;
t83 = -t108 * t169 + t110 * t168;
t212 = qJ(5) * t83 + qJD(5) * t86 - t219;
t188 = t110 * t54;
t51 = qJD(3) * pkin(3) + t53;
t20 = t108 * t51 + t188;
t17 = -qJD(3) * qJ(5) - t20;
t204 = pkin(5) * t78;
t11 = -t17 - t204;
t22 = t108 * t53 + t188;
t102 = -pkin(3) * t110 - pkin(4);
t98 = -pkin(9) + t102;
t210 = t98 * t70 + (t11 - t22 + t204) * t218;
t170 = qJD(2) * t114;
t185 = t109 * t114;
t84 = t111 * t113 + t116 * t185;
t125 = t84 * qJD(3);
t184 = t109 * t117;
t159 = qJD(2) * t184;
t150 = t113 * t159;
t121 = -t125 - t150;
t131 = t111 * t116 - t113 * t185;
t149 = t116 * t159;
t52 = qJD(3) * t131 + t149;
t23 = t108 * t121 + t110 * t52;
t209 = t109 * (t117 * t70 - t170 * t214) + t23 * qJD(3);
t21 = t108 * t52 - t110 * t121;
t80 = t86 * qJD(3);
t69 = qJD(2) * t80;
t208 = (-t117 * t69 + t78 * t170) * t109 - t21 * qJD(3);
t206 = pkin(4) + pkin(9);
t205 = pkin(4) * t69;
t203 = pkin(5) * t214;
t43 = t108 * t84 - t110 * t131;
t141 = qJD(4) + t161;
t33 = (-t113 * t89 + t97) * qJD(3) + (-qJ(4) * t169 + t141 * t116) * qJD(2);
t34 = (-t116 * t89 - t157) * qJD(3) + (-qJ(4) * t168 - t141 * t113) * qJD(2);
t7 = t108 * t33 - t110 * t34;
t202 = t43 * t7;
t91 = t194 * t113;
t92 = t194 * t116;
t57 = t108 * t92 + t110 * t91;
t201 = t57 * t7;
t200 = pkin(3) * t113;
t199 = t11 * t211;
t103 = -pkin(3) * t116 - pkin(2);
t71 = t103 * qJD(2) + qJD(4) - t161;
t122 = -qJ(5) * t214 + t71;
t32 = pkin(4) * t78 + t122;
t198 = t32 * t214;
t140 = -qJ(5) * t86 + t103;
t35 = -t206 * t211 + t140;
t197 = t35 * t70;
t196 = t59 * t78;
t61 = qJD(3) * t115 + t112 * t78;
t195 = t61 * t78;
t193 = pkin(5) * t83 + t191;
t192 = -pkin(5) * t80 + t190;
t8 = t108 * t34 + t110 * t33;
t160 = t109 * t170;
t76 = pkin(3) * t156 + qJD(1) * t160;
t189 = qJD(2) * pkin(2);
t166 = qJD(6) * t115;
t167 = qJD(6) * t112;
t36 = -qJD(3) * t167 + t112 * t69 + t78 * t166;
t187 = t36 * t115;
t119 = qJD(2) ^ 2;
t183 = t109 * t119;
t118 = qJD(3) ^ 2;
t181 = t118 * t113;
t180 = t118 * t116;
t176 = t203 - t177;
t174 = t113 ^ 2 - t116 ^ 2;
t104 = pkin(3) * t171;
t164 = t211 * t166;
t163 = t114 * t183;
t19 = t110 * t51 - t47;
t153 = t78 * qJ(5) + t104;
t148 = -qJ(5) * t70 + t76;
t147 = qJD(5) - t19;
t146 = -t206 * t80 + t212;
t145 = -pkin(4) * t80 + t212;
t142 = -t211 * t70 + t218 * t80;
t6 = -qJD(3) * qJD(5) - t8;
t10 = -t206 * qJD(3) + t147 + t203;
t18 = t206 * t78 + t122;
t1 = t10 * t115 - t112 * t18;
t2 = t10 * t112 + t115 * t18;
t58 = -t108 * t91 + t110 * t92;
t138 = -qJD(3) * t22 + t7;
t135 = -t112 * t43 + t115 * t184;
t134 = t112 * t184 + t115 * t43;
t3 = -pkin(5) * t69 - t6;
t132 = t3 + (-qJD(6) * t98 + t206 * t214 + t153) * t218;
t130 = t189 * qJD(2);
t41 = pkin(5) * t86 + t57;
t129 = -t11 * t80 + t211 * t3 + t41 * t70;
t128 = -qJD(5) * t214 + t148;
t126 = -t115 * t218 ^ 2 - t112 * t70;
t124 = t189 * t217;
t44 = t108 * t131 + t110 * t84;
t123 = t21 * t214 - t23 * t78 + t43 * t70 - t44 * t69;
t120 = -t190 * t78 + t191 * t214 + t57 * t70 - t58 * t69 + t7 * t86;
t100 = pkin(3) * t108 + qJ(5);
t72 = qJD(3) * t78;
t63 = t115 * t69;
t45 = -pkin(4) * t211 + t140;
t42 = pkin(5) * t211 + t58;
t38 = pkin(4) * t214 + t153;
t37 = t61 * qJD(6) - t63;
t16 = -qJD(3) * pkin(4) + t147;
t14 = t128 + t205;
t9 = t206 * t69 + t128;
t5 = pkin(5) * t70 + t7;
t4 = t115 * t5;
t12 = [0, 0, -t163, -t117 * t183, 0, 0, 0, 0, 0, -t116 * t163 + (-t125 - 0.2e1 * t150) * qJD(3), t113 * t163 + (-t52 - t149) * qJD(3), t208, -t209, t123, -t19 * t21 + t20 * t23 + t202 + t44 * t8 + (-t117 * t76 + t71 * t170) * t109, t123, -t208, t209, t16 * t21 - t17 * t23 + t202 - t44 * t6 + (-t117 * t14 + t32 * t170) * t109, 0, 0, 0, 0, 0, (qJD(6) * t135 - t112 * t160 + t115 * t21) * t218 + t134 * t70 + t23 * t59 + t44 * t37, -(qJD(6) * t134 + t112 * t21 + t115 * t160) * t218 + t135 * t70 + t23 * t61 + t44 * t36; 0, 0, 0, 0, 0.2e1 * t113 * t155, -0.2e1 * t174 * t165, t180, -t181, 0, -pkin(8) * t180 + t113 * t124, pkin(8) * t181 + t116 * t124, -t78 * t162 + t103 * t69 + t71 * t80 - t76 * t211 + (t78 * t200 - t191) * qJD(3), -t214 * t162 + t103 * t70 + t71 * t83 + t76 * t86 + (t200 * t214 - t190) * qJD(3), -t19 * t83 - t20 * t80 + t211 * t8 + t120, t103 * t76 - t191 * t19 + t190 * t20 + t219 * t71 + t58 * t8 + t201, t16 * t83 + t17 * t80 - t211 * t6 + t120, t191 * qJD(3) + t14 * t211 + t145 * t78 - t32 * t80 - t45 * t69, t190 * qJD(3) - t14 * t86 + t145 * t214 - t32 * t83 - t45 * t70, t14 * t45 - t145 * t32 + t191 * t16 - t190 * t17 - t58 * t6 + t201, -t61 * t164 + (-t211 * t36 + t61 * t80) * t112, (-t112 * t59 + t115 * t61) * t80 - (-t112 * t37 + t187 + (-t112 * t61 - t115 * t59) * qJD(6)) * t211, t112 * t142 - t164 * t218 + t36 * t86 + t61 * t83, t167 * t211 * t218 + t115 * t142 - t37 * t86 - t59 * t83, t218 * t83 + t70 * t86, t1 * t83 + t42 * t37 + t4 * t86 + t192 * t59 + (t146 * t218 - t9 * t86 - t197) * t112 + (t193 * t218 + t129) * t115 + ((-t112 * t41 - t115 * t35) * t218 - t2 * t86 - t112 * t199) * qJD(6), -t2 * t83 + t42 * t36 + t192 * t61 + (-t197 - (qJD(6) * t10 + t9) * t86 - qJD(6) * t199 + (-qJD(6) * t41 + t146) * t218) * t115 + (-(-qJD(6) * t18 + t5) * t86 + (qJD(6) * t35 - t193) * t218 - t129) * t112; 0, 0, 0, 0, -t113 * t119 * t116, t174 * t119, 0, 0, 0, t113 * t130, t116 * t130, -t78 * t104 - t214 * t71 - t138, qJD(3) * t24 - t104 * t214 + t71 * t78 - t8, (t20 - t22) * t214 + (-t19 + t24) * t78 + (-t108 * t69 - t110 * t70) * pkin(3), t19 * t22 - t20 * t24 + (t108 * t8 - t110 * t7 - t71 * t171) * pkin(3), -t100 * t69 + t102 * t70 + (-t17 - t22) * t214 + (t16 + t177) * t78, t38 * t78 + t138 + t198, -t32 * t78 + t38 * t214 + (0.2e1 * qJD(5) - t24) * qJD(3) + t8, -t100 * t6 + t102 * t7 - t16 * t22 + t177 * t17 - t32 * t38, -t152 * t61 + t187, (-t218 * t61 - t37) * t115 + (-t36 + t220) * t112, t139 + t195, t126 - t196, t218 * t78, t1 * t78 + t100 * t37 + t132 * t112 + t115 * t210 + t176 * t59, t100 * t36 - t112 * t210 + t132 * t115 + t176 * t61 - t2 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t214 * qJD(3), -t93 + (-t78 + t158) * qJD(3), t215, t19 * t214 + t20 * t78 + t76, t215, t214 * t217, -t70 + t72, t205 - t17 * t78 + (-qJD(5) - t16) * t214 + t148, 0, 0, 0, 0, 0, t126 + t196, t195 - t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 + t72, -t214 * t78, -t77 - t118, qJD(3) * t17 + t198 + t7, 0, 0, 0, 0, 0, -qJD(3) * t59 + t139, -qJD(3) * t61 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t59, -t59 ^ 2 + t61 ^ 2, t36 + t220, t213 * t61 + t63, t70, -t11 * t61 - t112 * t9 + t2 * t213 + t4, t1 * t213 + t11 * t59 - t112 * t5 - t115 * t9;];
tauc_reg = t12;
