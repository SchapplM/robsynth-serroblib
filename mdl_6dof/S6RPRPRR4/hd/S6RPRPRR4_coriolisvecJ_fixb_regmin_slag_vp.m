% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:27
% EndTime: 2019-03-09 03:46:34
% DurationCPUTime: 2.48s
% Computational Cost: add. (2201->304), mult. (4952->430), div. (0->0), fcn. (3164->8), ass. (0->186)
t127 = sin(pkin(10)) * pkin(1) + pkin(7);
t248 = pkin(4) + t127;
t146 = sin(qJ(6));
t147 = sin(qJ(5));
t149 = cos(qJ(6));
t150 = cos(qJ(5));
t106 = t146 * t150 + t147 * t149;
t148 = sin(qJ(3));
t161 = t106 * t148;
t252 = qJD(5) + qJD(6);
t245 = -qJD(1) * t161 - t252 * t106;
t217 = qJD(3) * t147;
t151 = cos(qJ(3));
t218 = qJD(1) * t151;
t103 = t150 * t218 + t217;
t194 = t147 * t218;
t215 = qJD(3) * t150;
t105 = -t194 + t215;
t175 = t103 * t146 - t149 * t105;
t32 = t149 * t103 + t105 * t146;
t261 = t175 * t32;
t112 = t127 * qJD(1);
t78 = -t151 * qJD(2) + t148 * t112;
t260 = qJD(4) + t78;
t210 = qJD(5) * t151;
t196 = t147 * t210;
t259 = t148 * t215 + t196;
t258 = t175 ^ 2 - t32 ^ 2;
t152 = -pkin(3) - pkin(8);
t220 = qJD(1) * t148;
t222 = pkin(4) * t220 + t260;
t41 = t152 * qJD(3) + t222;
t128 = -cos(pkin(10)) * pkin(1) - pkin(2);
t168 = -qJ(4) * t148 + t128;
t75 = t152 * t151 + t168;
t53 = t75 * qJD(1);
t17 = t147 * t41 + t150 * t53;
t14 = -pkin(9) * t103 + t17;
t209 = qJD(6) * t146;
t12 = t14 * t209;
t141 = qJD(3) * qJ(4);
t79 = t148 * qJD(2) + t151 * t112;
t60 = pkin(4) * t218 + t79;
t50 = t141 + t60;
t25 = pkin(5) * t103 + t50;
t257 = t25 * t32 + t12;
t126 = qJD(5) + t220;
t120 = qJD(6) + t126;
t208 = qJD(6) * t149;
t207 = qJD(1) * qJD(3);
t191 = t148 * t207;
t57 = -qJD(5) * t103 + t147 * t191;
t211 = qJD(5) * t150;
t58 = qJD(3) * t211 - qJD(5) * t194 - t150 * t191;
t9 = -t103 * t208 - t105 * t209 - t146 * t58 + t149 * t57;
t256 = t120 * t32 + t9;
t16 = -t147 * t53 + t150 * t41;
t13 = -pkin(9) * t105 + t16;
t11 = pkin(5) * t126 + t13;
t239 = t14 * t149;
t2 = t11 * t146 + t239;
t130 = t151 * t207;
t125 = pkin(3) * t191;
t178 = pkin(8) * t148 - qJ(4) * t151;
t213 = qJD(4) * t148;
t159 = qJD(3) * t178 - t213;
t44 = qJD(1) * t159 + t125;
t206 = qJD(2) * qJD(3);
t214 = qJD(3) * t151;
t69 = t112 * t214 + t148 * t206;
t45 = pkin(4) * t130 + t69;
t188 = -t147 * t44 + t150 * t45;
t157 = -qJD(5) * t17 + t188;
t4 = pkin(5) * t130 - pkin(9) * t57 + t157;
t205 = -t147 * t45 - t150 * t44 - t41 * t211;
t212 = qJD(5) * t147;
t164 = -t212 * t53 - t205;
t5 = -pkin(9) * t58 + t164;
t201 = -t146 * t5 + t149 * t4;
t255 = -qJD(6) * t2 + t25 * t175 + t201;
t156 = qJD(6) * t175 - t146 * t57 - t149 * t58;
t254 = -t120 * t175 + t156;
t65 = -t141 - t79;
t62 = -qJD(3) * pkin(3) + t260;
t143 = t151 ^ 2;
t195 = t150 * t210;
t229 = t126 * t148;
t251 = -(qJD(1) * t143 - t229) * t217 - t126 * t195;
t250 = t252 * t151;
t249 = t9 * t148 - t175 * t214;
t247 = pkin(9) - t152;
t225 = t149 * t150;
t227 = t146 * t147;
t174 = -t225 + t227;
t216 = qJD(3) * t148;
t20 = t106 * t250 - t174 * t216;
t81 = t174 * t151;
t246 = t20 * t120 + t81 * t130;
t219 = qJD(1) * t150;
t199 = t148 * t219;
t244 = -t146 * t212 - t147 * t209 + t149 * t199 - t220 * t227 + t252 * t225;
t243 = t105 * t214 + t57 * t148;
t93 = t248 * t148;
t83 = t147 * t93;
t242 = t150 * t75 + t83;
t135 = pkin(3) * t220;
t80 = qJD(1) * t178 + t135;
t241 = t147 * t60 + t150 * t80;
t240 = t259 * t126;
t238 = t148 * t58;
t140 = qJD(3) * qJD(4);
t233 = t112 * t216 - t151 * t206;
t56 = -t140 + t233;
t37 = -pkin(4) * t191 - t56;
t237 = t37 * t147;
t236 = t37 * t150;
t235 = t57 * t150;
t202 = -pkin(5) * t150 - pkin(4);
t234 = pkin(5) * t211 - t202 * t220 + t260;
t92 = -pkin(3) * t151 + t168;
t66 = qJD(1) * t92;
t232 = t103 * t126;
t231 = t103 * t151;
t230 = t105 * t126;
t228 = t126 * t152;
t226 = t147 * t148;
t224 = t150 * t151;
t153 = qJD(3) ^ 2;
t223 = t153 * t148;
t138 = t153 * t151;
t94 = t248 * t151;
t142 = t148 ^ 2;
t221 = t142 - t143;
t113 = qJD(1) * t128;
t204 = t150 * t229;
t154 = qJD(1) ^ 2;
t203 = t151 * t154 * t148;
t200 = t143 * t219;
t197 = t126 * t211;
t193 = pkin(9) * t151 - t75;
t110 = t247 * t150;
t192 = t244 * t120;
t190 = qJD(6) * t11 + t5;
t187 = -t147 * t80 + t150 * t60;
t134 = pkin(3) * t216;
t63 = t134 + t159;
t88 = t248 * t214;
t186 = -t147 * t63 + t150 * t88;
t183 = t245 * t120 - t174 * t130;
t182 = qJD(3) * t79 - t69;
t181 = qJD(3) * t78 - t233;
t109 = t247 * t147;
t169 = pkin(5) * t151 - pkin(9) * t226;
t180 = qJD(1) * t169 - qJD(6) * t109 - t247 * t212 + t187;
t179 = pkin(9) * t199 + t252 * t110 + t241;
t84 = t150 * t93;
t23 = pkin(5) * t148 + t147 * t193 + t84;
t24 = -pkin(9) * t224 + t242;
t177 = t146 * t23 + t149 * t24;
t176 = -0.2e1 * qJD(3) * t66;
t172 = 0.2e1 * qJD(3) * t113;
t171 = t126 * t147;
t167 = t148 * t156 - t214 * t32;
t162 = -qJ(4) * t214 - t213;
t67 = qJD(1) * t162 + t125;
t89 = t134 + t162;
t166 = qJD(1) * t89 + t127 * t153 + t67;
t165 = t148 * t50 + t152 * t214;
t163 = t147 * t88 + t150 * t63 + t93 * t211 - t212 * t75;
t19 = qJD(3) * t161 + t174 * t250;
t82 = t106 * t151;
t160 = -t120 * t19 + t130 * t82;
t155 = t148 * t69 - t151 * t56 + (t148 * t65 + t151 * t62) * qJD(3);
t129 = pkin(5) * t147 + qJ(4);
t117 = t150 * t130;
t116 = t148 * t130;
t108 = -qJ(4) * t218 + t135;
t87 = t248 * t216;
t68 = pkin(5) * t224 + t94;
t54 = t66 * t220;
t30 = -pkin(5) * t196 + (-t127 + t202) * t216;
t21 = pkin(5) * t58 + t37;
t7 = t259 * pkin(9) + t163;
t6 = t169 * qJD(3) + (t150 * t193 - t83) * qJD(5) + t186;
t1 = t11 * t149 - t14 * t146;
t3 = [0, 0, 0, 0, 0.2e1 * t116, -0.2e1 * t221 * t207, t138, -t223, 0, -t127 * t138 + t148 * t172, t127 * t223 + t151 * t172, t155, t148 * t176 + t151 * t166, -t148 * t166 + t151 * t176, t127 * t155 + t66 * t89 + t67 * t92, -t57 * t147 * t151 + (t147 * t216 - t195) * t105 (-t103 * t147 + t105 * t150) * t216 + (t147 * t58 - t235 + (t103 * t150 + t105 * t147) * qJD(5)) * t151, t243 + t251, -t238 + (-t200 - t231) * qJD(3) + t240, t126 * t214 + t116, t186 * t126 - t87 * t103 + t94 * t58 + (-t215 * t50 + t188) * t148 + (-t126 * t242 - t148 * t17) * qJD(5) + (-t50 * t212 + t236 + ((-t147 * t75 + t84) * qJD(1) + t16) * qJD(3)) * t151, -t163 * t126 - t87 * t105 + t94 * t57 + ((qJD(3) * t50 + qJD(5) * t53) * t147 + t205) * t148 + (-t50 * t211 - t237 + (-qJD(1) * t242 - t17) * qJD(3)) * t151, -t175 * t19 - t82 * t9, -t156 * t82 - t175 * t20 - t19 * t32 + t81 * t9, -t160 + t249, t167 + t246, t120 * t214 + t116 (-t146 * t7 + t149 * t6) * t120 + t201 * t148 + t30 * t32 - t68 * t156 - t21 * t81 - t25 * t20 + (-t120 * t177 - t148 * t2) * qJD(6) + ((-t146 * t24 + t149 * t23) * qJD(1) + t1) * t214, t12 * t148 + t25 * t19 - t21 * t82 - t30 * t175 + t68 * t9 + (-(-qJD(6) * t24 + t6) * t120 - t4 * t148) * t146 + (-(qJD(6) * t23 + t7) * t120 - t190 * t148) * t149 + (-qJD(1) * t177 - t2) * t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223, -t138, 0, t223, t138, -t148 * t56 - t151 * t69 + (t148 * t62 - t151 * t65) * qJD(3), 0, 0, 0, 0, 0, t238 + (-t200 + t231) * qJD(3) + t240, t243 - t251, 0, 0, 0, 0, 0, -t167 + t246, t160 + t249; 0, 0, 0, 0, -t203, t221 * t154, 0, 0, 0, -t113 * t220 + t182, -t113 * t218 - t181, 0, -t108 * t218 - t182 + t54, 0.2e1 * t140 + (t108 * t148 + t151 * t66) * qJD(1) + t181, -pkin(3) * t69 - qJ(4) * t56 - t108 * t66 - t260 * t65 - t62 * t79, -t105 * t171 + t235 (-t58 - t230) * t150 + (-t57 + t232) * t147, -t126 * t212 + t117 + (-t105 * t151 - t126 * t226) * qJD(1), -t197 + (-t204 + (t103 - t217) * t151) * qJD(1), -t126 * t218, qJ(4) * t58 + t237 - t187 * t126 + t222 * t103 + (-t147 * t228 + t150 * t50) * qJD(5) + (t150 * t165 - t16 * t151) * qJD(1), qJ(4) * t57 + t236 + t241 * t126 + t222 * t105 + (-t147 * t50 - t150 * t228) * qJD(5) + (-t147 * t165 + t17 * t151) * qJD(1), -t174 * t9 - t175 * t245, -t106 * t9 - t156 * t174 + t175 * t244 - t245 * t32, t175 * t218 + t183, -t192 + (-qJD(3) * t106 + t32) * t218, -t120 * t218, -t129 * t156 + t21 * t106 + t234 * t32 + t244 * t25 + (t146 * t179 - t149 * t180) * t120 + ((t109 * t146 - t110 * t149) * qJD(3) - t1) * t218, -t21 * t174 + t129 * t9 - t234 * t175 + t245 * t25 + (t146 * t180 + t149 * t179) * t120 + (-(-t109 * t149 - t110 * t146) * qJD(3) + t2) * t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, -t142 * t154 - t153, qJD(3) * t65 + t54 + t69, 0, 0, 0, 0, 0, -qJD(3) * t103 - t126 * t171 + t117, -t197 - qJD(3) * t105 + (-t147 * t214 - t204) * qJD(1), 0, 0, 0, 0, 0, -qJD(3) * t32 + t183, -t192 + (-t106 * t218 + t175) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t103, -t103 ^ 2 + t105 ^ 2, t57 + t232, -t58 + t230, t130, -t105 * t50 + t126 * t17 + t157, t103 * t50 + t126 * t16 - t164, -t261, t258, t256, t254, t130 -(-t13 * t146 - t239) * t120 + (-t105 * t32 - t120 * t209 + t130 * t149) * pkin(5) + t255 (-t120 * t14 - t4) * t146 + (t120 * t13 - t190) * t149 + (t105 * t175 - t120 * t208 - t130 * t146) * pkin(5) + t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261, t258, t256, t254, t130, t120 * t2 + t255, t1 * t120 - t146 * t4 - t149 * t190 + t257;];
tauc_reg  = t3;
