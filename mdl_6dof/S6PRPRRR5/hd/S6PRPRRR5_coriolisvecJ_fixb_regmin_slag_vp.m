% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:35
% EndTime: 2019-03-08 20:43:41
% DurationCPUTime: 1.89s
% Computational Cost: add. (2166->255), mult. (4906->367), div. (0->0), fcn. (3655->10), ass. (0->163)
t105 = sin(qJ(5));
t106 = sin(qJ(4));
t109 = cos(qJ(5));
t110 = cos(qJ(4));
t73 = t105 * t110 + t106 * t109;
t66 = t73 * qJD(2);
t207 = qJD(6) + t66;
t104 = sin(qJ(6));
t108 = cos(qJ(6));
t170 = qJD(2) * t110;
t147 = t109 * t170;
t171 = qJD(2) * t106;
t65 = t105 * t171 - t147;
t99 = qJD(4) + qJD(5);
t55 = -t104 * t65 - t108 * t99;
t209 = t207 * t55;
t208 = qJD(6) - t207;
t166 = qJD(5) * t105;
t168 = qJD(4) * t106;
t206 = -t105 * t168 - t106 * t166;
t111 = cos(qJ(2));
t102 = sin(pkin(6));
t174 = qJD(1) * t102;
t150 = t111 * t174;
t125 = qJD(3) - t150;
t205 = t66 * t99;
t167 = qJD(4) * t110;
t129 = -pkin(4) * t167 - t125;
t107 = sin(qJ(2));
t181 = t102 * t107;
t153 = qJD(2) * t181;
t133 = qJD(1) * t153;
t112 = -pkin(2) - pkin(8);
t64 = t112 * qJD(2) + t125;
t144 = pkin(9) * qJD(2) - t64;
t103 = cos(pkin(6));
t173 = qJD(1) * t103;
t148 = t106 * t173;
t31 = t106 * t133 + (-t144 * t110 - t148) * qJD(4);
t45 = -pkin(9) * t170 + t110 * t64 - t148;
t39 = qJD(4) * pkin(4) + t45;
t204 = (qJD(5) * t39 + t31) * t109;
t146 = t110 * t173;
t46 = -pkin(9) * t171 + t106 * t64 + t146;
t186 = t105 * t46;
t18 = t109 * t39 - t186;
t16 = -pkin(5) * t99 - t18;
t154 = t107 * t174;
t194 = pkin(9) - t112;
t76 = t194 * t106;
t77 = t194 * t110;
t51 = -t105 * t76 + t109 * t77;
t68 = t194 * t168;
t69 = qJD(4) * t77;
t193 = t51 * qJD(5) - t105 * t68 + t109 * t69 + t73 * t154;
t75 = t110 * t133;
t32 = t75 + (t144 * t106 - t146) * qJD(4);
t141 = t105 * t32 - t46 * t166;
t3 = t141 + t204;
t172 = qJD(2) * qJ(3);
t74 = t154 + t172;
t60 = pkin(4) * t171 + t74;
t30 = pkin(5) * t66 + pkin(10) * t65 + t60;
t142 = t105 * t31 - t109 * t32;
t184 = t109 * t46;
t19 = t105 * t39 + t184;
t4 = t19 * qJD(5) + t142;
t178 = t109 * t110;
t126 = t99 * t178;
t191 = t206 * qJD(2);
t41 = qJD(2) * t126 + t191;
t72 = t105 * t106 - t178;
t91 = t106 * pkin(4) + qJ(3);
t43 = pkin(5) * t73 + pkin(10) * t72 + t91;
t165 = qJD(5) * t109;
t47 = -t105 * t167 - t106 * t165 - t109 * t168 - t110 * t166;
t52 = -t105 * t77 - t109 * t76;
t202 = -(qJD(6) * t30 + t3) * t73 + t16 * t47 + (-qJD(6) * t43 + t193) * t207 - t4 * t72 - t52 * t41;
t201 = t16 * t66;
t200 = t16 * t72;
t199 = t41 * t72;
t198 = t43 * t41;
t48 = t126 + t206;
t197 = t48 * t99;
t196 = t207 * t65;
t195 = t65 * t66;
t192 = t52 * qJD(5) - t105 * t69 - t109 * t68 - t72 * t154;
t190 = qJD(2) * pkin(2);
t189 = t104 * t205;
t188 = t104 * t41;
t187 = t104 * t207;
t185 = t108 * t41;
t183 = t111 * t74;
t163 = qJD(6) * t108;
t164 = qJD(6) * t104;
t23 = -t108 * t205 + t99 * t163 + t65 * t164;
t182 = t23 * t104;
t180 = t102 * t111;
t114 = qJD(2) ^ 2;
t179 = t102 * t114;
t176 = t106 ^ 2 - t110 ^ 2;
t113 = qJD(4) ^ 2;
t175 = -t113 - t114;
t169 = qJD(2) * t111;
t162 = qJD(2) * qJD(4);
t161 = pkin(4) * t170;
t160 = t72 * t164;
t159 = t207 * t163;
t158 = t107 * t179;
t157 = t111 * t179;
t156 = -pkin(4) * t99 - t39;
t17 = pkin(10) * t99 + t19;
t124 = t104 * t17 - t108 * t30;
t155 = -t124 * t65 + t16 * t164;
t152 = t102 * t169;
t145 = t110 * t162;
t42 = -pkin(5) * t65 + pkin(10) * t66;
t93 = pkin(4) * t105 + pkin(10);
t139 = qJD(6) * t93 + t161 + t42;
t138 = t108 * t207;
t137 = qJD(6) * t73 + qJD(2);
t7 = t104 * t30 + t108 * t17;
t136 = t4 * t104 + t16 * t163 - t7 * t65;
t135 = t110 * t153;
t134 = t106 * t153;
t20 = t105 * t45 + t184;
t132 = pkin(4) * t166 - t20;
t131 = -t74 + t154;
t130 = pkin(5) * t48 - pkin(10) * t47 - t129;
t57 = t104 * t99 - t108 * t65;
t128 = -t23 * t72 + t47 * t57;
t127 = -t207 * t47 + t199;
t119 = -t103 * t110 + t106 * t180;
t62 = -t103 * t106 - t110 * t180;
t35 = -t105 * t119 - t109 * t62;
t36 = t105 * t62 - t109 * t119;
t123 = t60 * t65 - t142;
t122 = t60 * t66 - t141;
t121 = -t104 * t36 + t108 * t181;
t120 = t104 * t181 + t108 * t36;
t70 = (qJD(3) + t150) * qJD(2);
t118 = t131 - t172;
t58 = pkin(4) * t145 + t70;
t117 = t125 * qJD(2) - t112 * t113 + t70;
t21 = t109 * t45 - t186;
t116 = t201 - t93 * t41 + (-pkin(4) * t165 + t21) * t207;
t94 = -pkin(4) * t109 - pkin(5);
t71 = t125 - t190;
t50 = t119 * qJD(4) + t135;
t49 = t62 * qJD(4) + t134;
t44 = t47 * t99;
t33 = t65 ^ 2 - t66 ^ 2;
t28 = -t191 + (-t147 - t65) * t99;
t24 = qJD(6) * t57 - t189;
t13 = t41 * pkin(5) + pkin(10) * t205 + t58;
t12 = t108 * t13;
t11 = t36 * qJD(5) + t105 * t49 - t109 * t50;
t10 = -t35 * qJD(5) + t105 * t50 + t109 * t49;
t9 = t138 * t207 + t57 * t65 + t188;
t8 = -t187 * t207 - t55 * t65 + t185;
t5 = t57 * t138 + t182;
t1 = (t23 - t209) * t108 + (-t207 * t57 - t24) * t104;
t2 = [0, 0, -t158, -t157, t158, t157 (t107 * t70 + (t183 + (t71 - t150) * t107) * qJD(2)) * t102, 0, 0, 0, 0, 0, t106 * t157 + (t50 + t135) * qJD(4), t110 * t157 + (-t49 - t134) * qJD(4), 0, 0, 0, 0, 0, -t11 * t99 + (t107 * t41 + t66 * t169) * t102, -t10 * t99 + (-t107 * t205 - t65 * t169) * t102, 0, 0, 0, 0, 0 (-qJD(6) * t120 - t104 * t10 + t108 * t152) * t207 + t121 * t41 + t11 * t55 + t35 * t24 -(qJD(6) * t121 + t108 * t10 + t104 * t152) * t207 - t120 * t41 + t11 * t57 + t35 * t23; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t70 * qJ(3) + t74 * qJD(3) + (-t183 + (-t71 - t190) * t107) * t174, -0.2e1 * t106 * t145, 0.2e1 * t176 * t162, -t113 * t106, -t113 * t110, 0, t106 * t117 - t118 * t167, t110 * t117 + t118 * t168, t205 * t72 - t47 * t65, t205 * t73 - t47 * t66 + t48 * t65 + t199, t44, -t197, 0, -t129 * t66 - t192 * t99 + t41 * t91 + t48 * t60 + t58 * t73, t129 * t65 + t193 * t99 - t205 * t91 + t47 * t60 - t58 * t72, t108 * t128 + t57 * t160 (-t104 * t57 - t108 * t55) * t47 + (t182 + t108 * t24 + (-t104 * t55 + t108 * t57) * qJD(6)) * t72, -t108 * t127 + t160 * t207 + t23 * t73 + t57 * t48, t104 * t127 + t159 * t72 - t24 * t73 - t55 * t48, t207 * t48 + t41 * t73, t12 * t73 + t51 * t24 - t124 * t48 + t192 * t55 + (t198 + t130 * t207 + (-t17 * t73 - t207 * t52 - t200) * qJD(6)) * t108 + t202 * t104, t51 * t23 - t7 * t48 + t192 * t57 + (-t198 - (-qJD(6) * t17 + t13) * t73 + qJD(6) * t200 + (qJD(6) * t52 - t130) * t207) * t104 + t202 * t108; 0, 0, 0, 0, 0, -t114, t131 * qJD(2), 0, 0, 0, 0, 0, t175 * t106, t175 * t110, 0, 0, 0, 0, 0, -qJD(2) * t66 + t44, qJD(2) * t65 - t197, 0, 0, 0, 0, 0, -t73 * t188 + t24 * t72 - t47 * t55 + (-t104 * t48 - t108 * t137) * t207, -t73 * t185 + (t104 * t137 - t108 * t48) * t207 - t128; 0, 0, 0, 0, 0, 0, 0, t110 * t114 * t106, -t176 * t114, 0, 0, 0, -t74 * t170 + t75, -t131 * t171, -t195, t33, 0, t28, 0, -t66 * t161 + t20 * t99 + (t156 * t105 - t184) * qJD(5) + t123, t65 * t161 + t21 * t99 + (t156 * qJD(5) - t31) * t109 + t122, t5, t1, t9, t8, t196, t94 * t24 + t132 * t55 + (-t139 * t207 - t4) * t108 + t116 * t104 + t155, t108 * t116 + t132 * t57 + t139 * t187 + t94 * t23 + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t195, t33, 0, t28, 0, t123 + (-qJD(5) + t99) * t19, t18 * t99 + t122 - t204, t5, t1, t9, t8, t196, -pkin(5) * t24 - t4 * t108 - (-t104 * t18 + t108 * t42) * t207 - t19 * t55 + t104 * t201 + (-t159 - t188) * pkin(10) + t155, -pkin(5) * t23 + (t104 * t42 + t108 * t18) * t207 - t19 * t57 + t108 * t201 + (t164 * t207 - t185) * pkin(10) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t55, -t55 ^ 2 + t57 ^ 2, t23 + t209, -t208 * t57 + t189, t41, -t104 * t3 - t16 * t57 - t208 * t7 + t12, -t104 * t13 - t108 * t3 + t124 * t208 + t16 * t55;];
tauc_reg  = t2;
