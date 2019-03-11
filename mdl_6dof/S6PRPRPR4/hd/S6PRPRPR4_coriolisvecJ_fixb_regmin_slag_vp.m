% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:22
% EndTime: 2019-03-08 19:41:29
% DurationCPUTime: 2.70s
% Computational Cost: add. (3252->307), mult. (8727->455), div. (0->0), fcn. (7201->12), ass. (0->162)
t146 = cos(pkin(11));
t224 = cos(qJ(4));
t181 = t224 * t146;
t134 = qJD(2) * t181;
t143 = sin(pkin(11));
t149 = sin(qJ(4));
t196 = t149 * t143;
t177 = qJD(2) * t196;
t110 = -t134 + t177;
t106 = qJD(6) + t110;
t160 = t181 - t196;
t144 = sin(pkin(6));
t152 = cos(qJ(2));
t200 = t144 * t152;
t155 = t160 * t200;
t222 = pkin(8) + qJ(3);
t127 = t222 * t143;
t129 = t222 * t146;
t226 = -t224 * t127 - t149 * t129;
t233 = -qJD(1) * t155 + t160 * qJD(3) + t226 * qJD(4);
t115 = t160 * qJD(4);
t123 = t224 * t143 + t149 * t146;
t116 = t123 * qJD(4);
t150 = sin(qJ(2));
t192 = qJD(1) * t144;
t179 = t150 * t192;
t232 = t116 * pkin(4) - t115 * qJ(5) - t123 * qJD(5) - t179;
t125 = qJD(2) * qJ(3) + t179;
t147 = cos(pkin(6));
t191 = qJD(1) * t147;
t133 = t146 * t191;
t213 = pkin(8) * qJD(2);
t86 = t133 + (-t125 - t213) * t143;
t99 = t146 * t125 + t143 * t191;
t87 = t146 * t213 + t99;
t46 = t149 * t86 + t224 * t87;
t231 = t46 * qJD(4);
t230 = qJD(6) - t106;
t151 = cos(qJ(6));
t112 = qJD(2) * t123;
t142 = sin(pkin(12));
t145 = cos(pkin(12));
t93 = -t145 * qJD(4) + t142 * t112;
t229 = t151 * t93;
t148 = sin(qJ(6));
t122 = t151 * t142 + t148 * t145;
t208 = t106 * t122;
t95 = t142 * qJD(4) + t145 * t112;
t164 = t148 * t93 - t151 * t95;
t228 = t106 * t164;
t218 = -t233 * t142 + t232 * t145;
t217 = t232 * t142 + t233 * t145;
t156 = t123 * t200;
t90 = -t149 * t127 + t224 * t129;
t215 = -qJD(1) * t156 + qJD(3) * t123 + qJD(4) * t90;
t45 = -t149 * t87 + t224 * t86;
t178 = t152 * t192;
t119 = (qJD(3) + t178) * qJD(2);
t227 = t160 * t119;
t171 = qJD(3) - t178;
t104 = qJD(2) * t116;
t194 = t151 * t145;
t197 = t148 * t142;
t120 = -t194 + t197;
t209 = t106 * t120;
t225 = -t122 * t104 + t209 * t106;
t131 = qJD(4) * t134;
t103 = -qJD(4) * t177 + t131;
t22 = -t164 * qJD(6) + t122 * t103;
t107 = t110 ^ 2;
t223 = t145 * pkin(9);
t221 = pkin(9) + qJ(5);
t220 = t116 * pkin(5) - t115 * t223 + t218;
t206 = t115 * t142;
t219 = pkin(9) * t206 - t217;
t18 = t227 + (qJD(5) + t45) * qJD(4);
t130 = qJD(2) * t179;
t44 = t104 * pkin(4) - t103 * qJ(5) - t112 * qJD(5) + t130;
t7 = t142 * t44 + t145 * t18;
t41 = qJD(4) * qJ(5) + t46;
t138 = -t146 * pkin(3) - pkin(2);
t105 = t138 * qJD(2) + t171;
t55 = t110 * pkin(4) - t112 * qJ(5) + t105;
t13 = t142 * t55 + t145 * t41;
t75 = t112 * pkin(4) + t110 * qJ(5);
t20 = t142 * t75 + t145 * t45;
t216 = pkin(5) * t206 + t215;
t76 = -pkin(4) * t160 - t123 * qJ(5) + t138;
t39 = t142 * t76 + t145 * t90;
t188 = qJD(6) * t151;
t214 = t103 * t194 - t93 * t188;
t212 = qJD(2) * pkin(2);
t48 = t148 * t95 + t229;
t211 = t112 * t48;
t210 = t164 * t112;
t207 = t110 * t142;
t204 = t123 * t142;
t203 = t123 * t145;
t202 = t142 * t103;
t201 = t144 * t150;
t153 = qJD(2) ^ 2;
t199 = t144 * t153;
t198 = t145 * t103;
t193 = t143 ^ 2 + t146 ^ 2;
t190 = qJD(2) * t150;
t189 = qJD(6) * t123;
t186 = t150 * t199;
t185 = t152 * t199;
t6 = -t142 * t18 + t145 * t44;
t4 = t104 * pkin(5) - pkin(9) * t198 + t6;
t5 = -pkin(9) * t202 + t7;
t184 = -t148 * t5 + t151 * t4;
t180 = t144 * t190;
t12 = -t142 * t41 + t145 * t55;
t19 = -t142 * t45 + t145 * t75;
t38 = -t142 * t90 + t145 * t76;
t176 = t193 * t119;
t23 = t123 * t119 + t231;
t174 = -t120 * t104 - t208 * t106;
t173 = t148 * t4 + t151 * t5;
t8 = t110 * pkin(5) - t95 * pkin(9) + t12;
t9 = -t93 * pkin(9) + t13;
t172 = t148 * t9 - t151 * t8;
t2 = t148 * t8 + t151 * t9;
t170 = -t12 * t142 + t13 * t145;
t169 = t143 * (-t143 * t125 + t133) - t146 * t99;
t26 = -pkin(5) * t160 - pkin(9) * t203 + t38;
t28 = -pkin(9) * t204 + t39;
t168 = -t148 * t28 + t151 * t26;
t167 = t148 * t26 + t151 * t28;
t108 = -t143 * t201 + t147 * t146;
t109 = t147 * t143 + t146 * t201;
t68 = t149 * t108 + t224 * t109;
t56 = -t142 * t68 - t145 * t200;
t57 = -t142 * t200 + t145 * t68;
t166 = -t148 * t57 + t151 * t56;
t165 = t148 * t56 + t151 * t57;
t163 = -qJD(6) * t95 - t202;
t161 = t224 * t108 - t149 * t109;
t128 = t221 * t145;
t159 = t112 * pkin(5) + qJD(5) * t142 + qJD(6) * t128 + t110 * t223 + t19;
t126 = t221 * t142;
t158 = pkin(9) * t207 - qJD(5) * t145 + qJD(6) * t126 + t20;
t37 = -qJD(4) * pkin(4) + qJD(5) - t45;
t157 = -t103 * t226 + t115 * t37 + t123 * t23;
t21 = t163 * t148 + t214;
t154 = -pkin(4) * t103 - qJ(5) * t104 + (-qJD(5) + t37) * t110;
t137 = -t145 * pkin(5) - pkin(4);
t124 = t171 - t212;
t72 = t120 * t123;
t71 = t122 * t123;
t62 = pkin(5) * t204 - t226;
t43 = qJD(2) * t156 + qJD(4) * t68;
t42 = qJD(2) * t155 + qJD(4) * t161;
t33 = t142 * t180 + t145 * t42;
t32 = -t142 * t42 + t145 * t180;
t31 = t122 * t115 + t188 * t203 - t189 * t197;
t30 = -t120 * t115 - t122 * t189;
t29 = -pkin(5) * t207 + t46;
t27 = t93 * pkin(5) + t37;
t14 = pkin(5) * t202 + t23;
t1 = [0, 0, -t186, -t185, -t146 * t186, t143 * t186, t193 * t185 (-t108 * t143 + t109 * t146) * t119 + (t124 * t150 + (-t169 - t179) * t152) * t144 * qJD(2), 0, 0, 0, 0, 0, -t43 * qJD(4) + (-t104 * t152 + t110 * t190) * t144, -t42 * qJD(4) + (-t103 * t152 + t112 * t190) * t144, t56 * t104 + t32 * t110 - t161 * t202 + t43 * t93, -t57 * t104 - t33 * t110 - t161 * t198 + t43 * t95, -t32 * t95 - t33 * t93 + (-t142 * t57 - t145 * t56) * t103, t12 * t32 + t13 * t33 - t161 * t23 + t37 * t43 + t6 * t56 + t7 * t57, 0, 0, 0, 0, 0 (-t165 * qJD(6) - t148 * t33 + t151 * t32) * t106 + t166 * t104 + t43 * t48 - t161 * t22 -(t166 * qJD(6) + t148 * t32 + t151 * t33) * t106 - t165 * t104 - t43 * t164 - t161 * t21; 0, 0, 0, 0, 0, 0, t171 * qJD(2) * t193 + t176, -t169 * qJD(3) + qJ(3) * t176 + (t169 * t152 + (-t124 - t212) * t150) * t192, t103 * t123 + t112 * t115, t103 * t160 - t123 * t104 - t115 * t110 - t112 * t116, t115 * qJD(4), -t116 * qJD(4), 0, t138 * t104 + t105 * t116 - t215 * qJD(4) + (-qJD(2) * t160 - t110) * t179, -t233 * qJD(4) + t138 * t103 + t105 * t115, t38 * t104 + t218 * t110 + t12 * t116 + t157 * t142 - t160 * t6 + t215 * t93, -t39 * t104 - t217 * t110 - t13 * t116 + t157 * t145 + t160 * t7 + t215 * t95, -t218 * t95 - t217 * t93 + (-t103 * t38 - t115 * t12 - t123 * t6) * t145 + (-t103 * t39 - t115 * t13 - t123 * t7) * t142, t218 * t12 + t217 * t13 + t215 * t37 - t226 * t23 + t6 * t38 + t7 * t39, -t164 * t30 - t21 * t72, t164 * t31 - t21 * t71 + t72 * t22 - t30 * t48, -t72 * t104 + t30 * t106 - t116 * t164 - t160 * t21, -t71 * t104 - t31 * t106 - t48 * t116 + t160 * t22, -t104 * t160 + t106 * t116, t168 * t104 - t184 * t160 - t172 * t116 + t62 * t22 + t14 * t71 + t27 * t31 + t216 * t48 + (t219 * t148 + t220 * t151) * t106 + (-t167 * t106 + t160 * t2) * qJD(6), -t167 * t104 + t173 * t160 - t2 * t116 + t62 * t21 - t14 * t72 + t27 * t30 - t216 * t164 + (-t220 * t148 + t219 * t151) * t106 + (-t168 * t106 - t160 * t172) * qJD(6); 0, 0, 0, 0, 0, 0, -t193 * t153, t169 * qJD(2) + t130, 0, 0, 0, 0, 0, 0.2e1 * t112 * qJD(4), t131 + (-t110 - t177) * qJD(4), t145 * t104 - t142 * t107 - t112 * t93, -t142 * t104 - t145 * t107 - t112 * t95 (t142 * t95 - t145 * t93) * t110 + (-t142 ^ 2 - t145 ^ 2) * t103, t170 * t110 - t37 * t112 + t7 * t142 + t6 * t145, 0, 0, 0, 0, 0, t174 - t211, t210 + t225; 0, 0, 0, 0, 0, 0, 0, 0, t112 * t110, t112 ^ 2 - t107, t131 + (t110 - t177) * qJD(4), 0, 0, -t105 * t112 - t23 + t231, t105 * t110 - t227, -t19 * t110 - t12 * t112 + t142 * t154 - t23 * t145 - t46 * t93, t20 * t110 + t13 * t112 + t23 * t142 + t145 * t154 - t46 * t95, t19 * t95 + t20 * t93 + (-qJD(5) * t93 - t110 * t12 + t7) * t145 + (qJD(5) * t95 - t110 * t13 - t6) * t142, -t23 * pkin(4) - t12 * t19 - t13 * t20 - t37 * t46 + t170 * qJD(5) + (-t6 * t142 + t7 * t145) * qJ(5), t21 * t122 + t164 * t209, -t21 * t120 - t122 * t22 + t164 * t208 + t209 * t48, t210 - t225, t174 + t211, -t106 * t112 (-t151 * t126 - t148 * t128) * t104 + t137 * t22 + t14 * t120 + t172 * t112 - t29 * t48 + t208 * t27 + (t148 * t158 - t151 * t159) * t106 -(-t148 * t126 + t151 * t128) * t104 + t137 * t21 + t14 * t122 + t2 * t112 + t29 * t164 - t209 * t27 + (t148 * t159 + t151 * t158) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110 * t95 + t202, -t110 * t93 + t198, -t93 ^ 2 - t95 ^ 2, t12 * t95 + t13 * t93 + t23, 0, 0, 0, 0, 0, t22 - t228, -t106 * t229 + (-t106 * t95 + t163) * t148 + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164 * t48, t164 ^ 2 - t48 ^ 2, t48 * t106 + t21, -t22 - t228, t104, t27 * t164 - t230 * t2 + t184, t230 * t172 + t27 * t48 - t173;];
tauc_reg  = t1;
