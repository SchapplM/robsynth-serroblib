% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:22
% EndTime: 2019-03-09 02:21:28
% DurationCPUTime: 2.38s
% Computational Cost: add. (3409->289), mult. (8401->398), div. (0->0), fcn. (6544->10), ass. (0->159)
t149 = cos(qJ(4));
t142 = cos(pkin(11));
t189 = qJD(1) * t142;
t129 = t149 * t189;
t140 = sin(pkin(11));
t146 = sin(qJ(4));
t195 = t146 * t140;
t176 = qJD(1) * t195;
t109 = t129 - t176;
t224 = qJD(5) + qJD(6);
t236 = t109 - t224;
t130 = sin(pkin(10)) * pkin(1) + qJ(3);
t123 = t130 * qJD(1);
t136 = t142 * qJD(2);
t91 = t136 + (-pkin(7) * qJD(1) - t123) * t140;
t101 = t140 * qJD(2) + t142 * t123;
t92 = pkin(7) * t189 + t101;
t39 = t146 * t91 + t149 * t92;
t235 = t39 * qJD(4);
t144 = sin(qJ(6));
t145 = sin(qJ(5));
t147 = cos(qJ(6));
t148 = cos(qJ(5));
t119 = t144 * t148 + t147 * t145;
t213 = t236 * t119;
t117 = t149 * t140 + t146 * t142;
t110 = t117 * qJD(1);
t183 = t148 * qJD(4);
t93 = t145 * t110 - t183;
t95 = t145 * qJD(4) + t148 * t110;
t164 = t144 * t93 - t147 * t95;
t41 = t144 * t95 + t147 * t93;
t234 = t164 * t41;
t187 = qJD(5) * t145;
t197 = t145 * t109;
t233 = t187 - t197;
t232 = t164 ^ 2 - t41 ^ 2;
t107 = qJD(5) - t109;
t102 = qJD(6) + t107;
t184 = qJD(6) * t147;
t185 = qJD(6) * t144;
t128 = qJD(4) * t129;
t103 = -qJD(4) * t176 + t128;
t54 = qJD(5) * t183 + t148 * t103 - t110 * t187;
t55 = qJD(5) * t95 + t145 * t103;
t9 = -t144 * t55 + t147 * t54 - t93 * t184 - t185 * t95;
t231 = t41 * t102 + t9;
t35 = qJD(4) * pkin(8) + t39;
t122 = -cos(pkin(10)) * pkin(1) - t142 * pkin(3) - pkin(2);
t108 = qJD(1) * t122 + qJD(3);
t53 = -t109 * pkin(4) - t110 * pkin(8) + t108;
t22 = t145 * t53 + t148 * t35;
t14 = -t93 * pkin(9) + t22;
t12 = t14 * t185;
t227 = -t146 * t92 + t149 * t91;
t34 = -qJD(4) * pkin(4) - t227;
t26 = t93 * pkin(5) + t34;
t230 = t26 * t41 + t12;
t112 = t117 * qJD(4);
t104 = qJD(1) * t112;
t116 = -t149 * t142 + t195;
t155 = t116 * qJD(3);
t28 = -qJD(1) * t155 + qJD(4) * t227;
t64 = t104 * pkin(4) - t103 * pkin(8);
t59 = t148 * t64;
t153 = -qJD(5) * t22 - t145 * t28 + t59;
t2 = t104 * pkin(5) - t54 * pkin(9) + t153;
t186 = qJD(5) * t148;
t160 = t145 * t64 + t148 * t28 + t53 * t186 - t187 * t35;
t5 = -t55 * pkin(9) + t160;
t180 = -t144 * t5 + t147 * t2;
t21 = -t145 * t35 + t148 * t53;
t13 = -t95 * pkin(9) + t21;
t11 = t107 * pkin(5) + t13;
t210 = t147 * t14;
t4 = t144 * t11 + t210;
t229 = -qJD(6) * t4 + t26 * t164 + t180;
t152 = qJD(6) * t164 - t144 * t54 - t147 * t55;
t228 = -t102 * t164 + t152;
t72 = t119 * t117;
t219 = pkin(7) + t130;
t113 = t219 * t140;
t114 = t219 * t142;
t226 = -t149 * t113 - t146 * t114;
t111 = t116 * qJD(4);
t193 = t148 * t111;
t225 = -t117 * t187 - t193;
t118 = t144 * t145 - t147 * t148;
t214 = t236 * t118;
t223 = -t102 * t214 - t119 * t104;
t97 = t148 * t104;
t222 = t107 * t225 + t117 * t97;
t221 = pkin(8) + pkin(9);
t220 = -t112 * t164 + t9 * t116;
t196 = t145 * t111;
t201 = t117 * t148;
t202 = t117 * t145;
t19 = -t185 * t202 + (t224 * t201 - t196) * t147 + t225 * t144;
t218 = -t19 * t102 - t72 * t104;
t217 = t95 * t112 + t54 * t116;
t77 = t110 * pkin(4) - t109 * pkin(8);
t216 = t145 * t77 + t148 * t227;
t75 = -t146 * t113 + t149 * t114;
t66 = t148 * t75;
t67 = t116 * pkin(4) - t117 * pkin(8) + t122;
t215 = t145 * t67 + t66;
t212 = t110 * t41;
t211 = t110 * t93;
t182 = qJD(1) * qJD(3);
t29 = t117 * t182 + t235;
t208 = t29 * t148;
t207 = t164 * t110;
t206 = t54 * t145;
t205 = t93 * t107;
t204 = t95 * t107;
t203 = t95 * t110;
t198 = t145 * t104;
t190 = t140 ^ 2 + t142 ^ 2;
t188 = qJD(4) * t111;
t179 = qJD(5) * t221;
t177 = t117 * t186;
t175 = qJD(6) * t11 + t5;
t173 = qJD(1) * t190;
t172 = t107 * t148;
t171 = t213 * t102 - t118 * t104;
t170 = t233 * pkin(5) - t39;
t124 = t221 * t145;
t169 = -pkin(9) * t197 + qJD(6) * t124 + t145 * t179 + t216;
t125 = t221 * t148;
t70 = t148 * t77;
t168 = t110 * pkin(5) + qJD(6) * t125 - t145 * t227 + t70 + (-pkin(9) * t109 + t179) * t148;
t167 = -t112 * t41 + t116 * t152;
t18 = t118 * t111 - t224 * t72;
t73 = t118 * t117;
t166 = -t18 * t102 + t73 * t104;
t165 = -t112 * t93 - t116 * t55;
t162 = (-t140 * t123 + t136) * t140 - t101 * t142;
t161 = -t233 * t107 + t97;
t46 = t226 * qJD(4) - t155;
t78 = t112 * pkin(4) + t111 * pkin(8);
t159 = t145 * t78 + t148 * t46 + t67 * t186 - t187 * t75;
t158 = t177 - t196;
t156 = -pkin(8) * t104 + t107 * t34;
t151 = -t107 * t158 - t117 * t198;
t47 = qJD(3) * t117 + qJD(4) * t75;
t134 = -t148 * pkin(5) - pkin(4);
t106 = t112 * qJD(4);
t79 = t104 * t116;
t71 = t148 * t78;
t63 = t148 * t67;
t48 = pkin(5) * t202 - t226;
t25 = pkin(5) * t158 + t47;
t24 = -pkin(9) * t202 + t215;
t23 = t116 * pkin(5) - pkin(9) * t201 - t145 * t75 + t63;
t17 = t55 * pkin(5) + t29;
t7 = -pkin(9) * t158 + t159;
t6 = pkin(9) * t193 + t112 * pkin(5) - t145 * t46 + t71 + (-t66 + (pkin(9) * t117 - t67) * t145) * qJD(5);
t3 = t147 * t11 - t144 * t14;
t1 = [0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t173 (t130 * t173 - t162) * qJD(3), t103 * t117 - t110 * t111, -t116 * t103 - t104 * t117 - t109 * t111 - t112 * t110, -t188, -t106, 0, -t47 * qJD(4) + t122 * t104 + t108 * t112, -t46 * qJD(4) + t122 * t103 - t108 * t111, -t95 * t193 + (t54 * t148 - t187 * t95) * t117 -(-t145 * t95 - t148 * t93) * t111 + (-t206 - t148 * t55 + (t145 * t93 - t148 * t95) * qJD(5)) * t117, t217 + t222, t151 + t165, t107 * t112 + t79 (-t186 * t75 + t71) * t107 + t63 * t104 + (-t186 * t35 + t59) * t116 + t21 * t112 + t47 * t93 - t226 * t55 + t34 * t177 + ((-qJD(5) * t67 - t46) * t107 - t75 * t104 + (-qJD(5) * t53 - t28) * t116 + t29 * t117 - t34 * t111) * t145, -t159 * t107 - t215 * t104 - t160 * t116 - t22 * t112 + t47 * t95 - t226 * t54 - t34 * t193 + (-t187 * t34 + t208) * t117, -t164 * t18 - t9 * t73, -t152 * t73 + t164 * t19 - t18 * t41 - t9 * t72, -t166 + t220, t167 + t218, t102 * t112 + t79 (-t144 * t7 + t147 * t6) * t102 + (-t144 * t24 + t147 * t23) * t104 + t180 * t116 + t3 * t112 + t25 * t41 - t48 * t152 + t17 * t72 + t26 * t19 + ((-t144 * t23 - t147 * t24) * t102 - t4 * t116) * qJD(6), -t4 * t112 + t12 * t116 - t17 * t73 + t26 * t18 - t25 * t164 + t48 * t9 + (-(-qJD(6) * t24 + t6) * t102 - t23 * t104 - t2 * t116) * t144 + (-(qJD(6) * t23 + t7) * t102 - t24 * t104 - t175 * t116) * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, t188, 0, 0, 0, 0, 0, t151 - t165, t217 - t222, 0, 0, 0, 0, 0, -t167 + t218, t166 + t220; 0, 0, 0, 0, 0, 0, -t190 * qJD(1) ^ 2, t162 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t110 * qJD(4), t128 + (t109 - t176) * qJD(4), 0, 0, 0, 0, 0, t161 - t211, -t107 ^ 2 * t148 - t198 - t203, 0, 0, 0, 0, 0, t171 - t212, t207 + t223; 0, 0, 0, 0, 0, 0, 0, 0, -t109 * t110, -t109 ^ 2 + t110 ^ 2, t128 + (-t109 - t176) * qJD(4), 0, 0, -t108 * t110 + t235 - t29, -t108 * t109 + t116 * t182, t172 * t95 + t206 (t54 - t205) * t148 + (-t55 - t204) * t145, t107 * t172 + t198 - t203, t161 + t211, -t107 * t110, -pkin(4) * t55 - t21 * t110 - t208 - t39 * t93 + (-pkin(8) * t186 - t70) * t107 + (t107 * t227 + t156) * t145, -pkin(4) * t54 + t22 * t110 + t29 * t145 - t39 * t95 + (pkin(8) * t187 + t216) * t107 + t156 * t148, t9 * t119 - t164 * t214, -t9 * t118 + t119 * t152 - t164 * t213 - t214 * t41, t207 - t223, t171 + t212, -t102 * t110 (-t147 * t124 - t144 * t125) * t104 - t134 * t152 + t17 * t118 - t3 * t110 + t170 * t41 - t213 * t26 + (t144 * t169 - t147 * t168) * t102 -(-t144 * t124 + t147 * t125) * t104 + t134 * t9 + t17 * t119 + t4 * t110 - t170 * t164 + t214 * t26 + (t144 * t168 + t147 * t169) * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t93, -t93 ^ 2 + t95 ^ 2, t54 + t205, t204 - t55, t104, t22 * t107 - t34 * t95 + t153, t21 * t107 + t34 * t93 - t160, -t234, t232, t231, t228, t104 -(-t144 * t13 - t210) * t102 + (-t102 * t185 + t147 * t104 - t95 * t41) * pkin(5) + t229 (-t14 * t102 - t2) * t144 + (t13 * t102 - t175) * t147 + (-t102 * t184 - t144 * t104 + t164 * t95) * pkin(5) + t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, t232, t231, t228, t104, t4 * t102 + t229, t3 * t102 - t144 * t2 - t147 * t175 + t230;];
tauc_reg  = t1;
