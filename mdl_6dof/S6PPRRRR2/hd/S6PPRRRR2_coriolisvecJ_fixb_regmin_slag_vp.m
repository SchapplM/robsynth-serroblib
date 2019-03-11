% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:29
% EndTime: 2019-03-08 19:05:36
% DurationCPUTime: 3.09s
% Computational Cost: add. (2922->339), mult. (7794->520), div. (0->0), fcn. (6694->14), ass. (0->189)
t142 = sin(qJ(4));
t146 = cos(qJ(4));
t176 = pkin(4) * t142 - pkin(10) * t146;
t113 = t176 * qJD(4);
t139 = cos(pkin(6));
t124 = qJD(1) * t139 + qJD(2);
t134 = sin(pkin(13));
t136 = sin(pkin(6));
t143 = sin(qJ(3));
t147 = cos(qJ(3));
t137 = cos(pkin(13));
t138 = cos(pkin(7));
t224 = t137 * t138;
t154 = (t134 * t147 + t143 * t224) * t136;
t135 = sin(pkin(7));
t227 = t135 * t143;
t58 = qJD(1) * t154 + t124 * t227;
t264 = t113 - t58;
t141 = sin(qJ(5));
t145 = cos(qJ(5));
t206 = t145 * qJD(4);
t216 = qJD(3) * t142;
t104 = t141 * t216 - t206;
t56 = qJD(3) * pkin(9) + t58;
t217 = qJD(1) * t136;
t194 = t137 * t217;
t84 = t124 * t138 - t135 * t194;
t253 = -t142 * t56 + t146 * t84;
t33 = -qJD(4) * pkin(4) - t253;
t26 = t104 * pkin(5) + t33;
t207 = t141 * qJD(4);
t106 = t145 * t216 + t207;
t140 = sin(qJ(6));
t144 = cos(qJ(6));
t62 = t144 * t104 + t140 * t106;
t263 = t26 * t62;
t167 = t140 * t104 - t144 * t106;
t262 = t167 * t62;
t211 = qJD(5) * t141;
t189 = t142 * t211;
t250 = t146 * t206 - t189;
t226 = t135 * t147;
t251 = (-t134 * t143 + t147 * t224) * t136;
t261 = t139 * t226 + t251;
t260 = t167 ^ 2 - t62 ^ 2;
t205 = t146 * qJD(3);
t125 = -qJD(5) + t205;
t122 = -qJD(6) + t125;
t208 = qJD(6) * t144;
t209 = qJD(6) * t140;
t203 = qJD(4) * qJD(5);
t81 = qJD(3) * t250 + t145 * t203;
t187 = t146 * t207;
t210 = qJD(5) * t145;
t188 = t142 * t210;
t157 = t187 + t188;
t82 = qJD(3) * t157 + t141 * t203;
t24 = -t104 * t208 - t106 * t209 - t140 * t82 + t144 * t81;
t259 = -t122 * t62 + t24;
t204 = qJD(3) * qJD(4);
t127 = t142 * t204;
t116 = -t146 * pkin(4) - t142 * pkin(10) - pkin(3);
t180 = t138 * t194;
t195 = t134 * t217;
t57 = -t143 * t195 + t147 * (t124 * t135 + t180);
t48 = qJD(3) * t116 - t57;
t239 = t141 * t48;
t36 = t142 * t84 + t146 * t56;
t34 = qJD(4) * pkin(10) + t36;
t14 = t145 * t34 + t239;
t256 = qJD(1) * t251 + t124 * t226;
t53 = t256 * qJD(3);
t16 = qJD(4) * t253 + t146 * t53;
t215 = qJD(3) * t143;
t193 = t135 * t215;
t214 = qJD(3) * t147;
t54 = t124 * t193 + t180 * t215 + t195 * t214;
t44 = qJD(3) * t113 + t54;
t182 = t141 * t16 - t145 * t44;
t151 = -qJD(5) * t14 - t182;
t2 = pkin(5) * t127 - t81 * pkin(11) + t151;
t160 = t141 * t44 + t145 * t16 + t48 * t210 - t34 * t211;
t3 = -pkin(11) * t82 + t160;
t198 = -t140 * t3 + t144 * t2;
t11 = -t104 * pkin(11) + t14;
t240 = t11 * t144;
t13 = -t141 * t34 + t145 * t48;
t10 = -pkin(11) * t106 + t13;
t8 = -pkin(5) * t125 + t10;
t5 = t140 * t8 + t240;
t258 = -qJD(6) * t5 + t26 * t167 + t198;
t150 = qJD(6) * t167 - t140 * t81 - t144 * t82;
t257 = t122 * t167 + t150;
t222 = t141 * t146;
t255 = -t142 * pkin(9) * t207 - t264 * t145 - t57 * t222;
t220 = t145 * t146;
t254 = t116 * t210 + t264 * t141 - t57 * t220;
t252 = qJD(3) * t58 - t54;
t108 = t140 * t145 + t144 * t141;
t87 = t108 * t142;
t249 = qJD(5) + qJD(6);
t185 = qJD(6) * t8 + t3;
t9 = t11 * t209;
t248 = t140 * t2 + t144 * t185 - t9;
t247 = pkin(10) + pkin(11);
t126 = pkin(9) * t220;
t165 = pkin(5) * t142 - pkin(11) * t220;
t246 = -t165 * qJD(4) - (-t126 + (pkin(11) * t142 - t116) * t141) * qJD(5) + t255;
t245 = -t157 * pkin(11) + (-t142 * t206 - t146 * t211) * pkin(9) + t254;
t110 = t176 * qJD(3);
t244 = t141 * t110 + t145 * t253;
t107 = t140 * t141 - t144 * t145;
t159 = t107 * t146;
t243 = qJD(3) * t159 - t249 * t107;
t242 = (-t205 + t249) * t108;
t241 = qJD(3) * pkin(3);
t212 = qJD(4) * t146;
t213 = qJD(4) * t142;
t17 = t142 * t53 + t56 * t212 + t84 * t213;
t237 = t17 * t141;
t236 = t17 * t145;
t235 = t33 * t141;
t234 = t81 * t141;
t230 = t104 * t125;
t229 = t106 * t125;
t228 = t125 * t145;
t149 = qJD(3) ^ 2;
t225 = t135 * t149;
t223 = t141 * t142;
t221 = t142 * t145;
t219 = t141 * t116 + t126;
t132 = t142 ^ 2;
t218 = -t146 ^ 2 + t132;
t200 = t142 * t227;
t199 = t146 * t227;
t196 = qJD(5) * t247;
t192 = t135 * t214;
t191 = t141 * t205;
t190 = t125 * t211;
t181 = t145 * t110 - t141 * t253;
t179 = t146 * t192;
t178 = t142 * t192;
t177 = -t36 + (-t191 + t211) * pkin(5);
t119 = t247 * t141;
t175 = pkin(11) * t191 - qJD(6) * t119 - t141 * t196 - t244;
t120 = t247 * t145;
t174 = t165 * qJD(3) + qJD(6) * t120 + t145 * t196 + t181;
t68 = t139 * t227 + t154;
t90 = -t135 * t136 * t137 + t138 * t139;
t47 = t90 * t142 + t68 * t146;
t22 = -t47 * t141 - t145 * t261;
t23 = -t141 * t261 + t145 * t47;
t173 = -t140 * t23 + t144 * t22;
t172 = t140 * t22 + t144 * t23;
t103 = t145 * t116;
t66 = -pkin(11) * t221 + t103 + (-pkin(9) * t141 - pkin(5)) * t146;
t75 = -pkin(11) * t223 + t219;
t171 = t140 * t66 + t144 * t75;
t94 = t142 * t138 + t199;
t162 = t141 * t226 - t145 * t94;
t73 = -t141 * t94 - t145 * t226;
t170 = t140 * t162 + t144 * t73;
t169 = t140 * t73 - t144 * t162;
t46 = t68 * t142 - t90 * t146;
t166 = qJD(3) * t132 - t125 * t146;
t148 = qJD(4) ^ 2;
t164 = pkin(9) * t148 - t252;
t55 = -t57 - t241;
t163 = qJD(4) * (t55 + t57 - t241);
t93 = -t146 * t138 + t200;
t130 = -pkin(5) * t145 - pkin(4);
t114 = (pkin(5) * t141 + pkin(9)) * t142;
t88 = t107 * t142;
t83 = pkin(5) * t157 + pkin(9) * t212;
t72 = qJD(4) * t94 + t178;
t71 = -qJD(4) * t93 + t179;
t60 = t68 * qJD(3);
t59 = t261 * qJD(3);
t40 = -t209 * t223 + (t249 * t221 + t187) * t144 + t250 * t140;
t39 = -qJD(4) * t159 - t249 * t87;
t32 = qJD(5) * t162 - t141 * t71 + t145 * t193;
t31 = qJD(5) * t73 + t141 * t193 + t145 * t71;
t21 = -qJD(4) * t46 + t59 * t146;
t20 = qJD(4) * t47 + t59 * t142;
t12 = pkin(5) * t82 + t17;
t7 = qJD(5) * t22 + t60 * t141 + t145 * t21;
t6 = -qJD(5) * t23 - t21 * t141 + t145 * t60;
t4 = -t140 * t11 + t144 * t8;
t1 = [0, 0, 0, -t60 * qJD(3), -t59 * qJD(3), 0, 0, 0, 0, 0, -t20 * qJD(4) + (-t146 * t60 - t213 * t261) * qJD(3), -t21 * qJD(4) + (t142 * t60 - t212 * t261) * qJD(3), 0, 0, 0, 0, 0, t104 * t20 - t125 * t6 + t127 * t22 + t46 * t82, t106 * t20 + t125 * t7 - t127 * t23 + t46 * t81, 0, 0, 0, 0, 0 -(-qJD(6) * t172 - t140 * t7 + t144 * t6) * t122 + t173 * t127 + t20 * t62 - t46 * t150 (qJD(6) * t173 + t140 * t6 + t144 * t7) * t122 - t172 * t127 - t20 * t167 + t46 * t24; 0, 0, 0, -t143 * t225, -t147 * t225, 0, 0, 0, 0, 0, -t149 * t199 + (-t72 - t178) * qJD(4), t149 * t200 + (-t71 - t179) * qJD(4), 0, 0, 0, 0, 0, t104 * t72 - t125 * t32 + t127 * t73 + t82 * t93, t106 * t72 + t125 * t31 + t127 * t162 + t81 * t93, 0, 0, 0, 0, 0 -(-qJD(6) * t169 - t140 * t31 + t144 * t32) * t122 + t170 * t127 + t72 * t62 - t93 * t150 (qJD(6) * t170 + t140 * t32 + t144 * t31) * t122 - t169 * t127 - t72 * t167 + t93 * t24; 0, 0, 0, t252 (-t256 + t57) * qJD(3), 0.2e1 * t146 * t127, -0.2e1 * t218 * t204, t148 * t146, -t148 * t142, 0, t142 * t163 - t146 * t164, t142 * t164 + t146 * t163, t250 * t106 + t81 * t221 (-t104 * t145 - t106 * t141) * t212 + (-t234 - t145 * t82 + (t104 * t141 - t106 * t145) * qJD(5)) * t142, t125 * t189 - t81 * t146 + (t106 * t142 + t145 * t166) * qJD(4), t125 * t188 + t82 * t146 + (-t104 * t142 - t141 * t166) * qJD(4) (-t125 - t205) * t213 (t116 * t211 + t255) * t125 + ((pkin(9) * t104 + t235) * qJD(4) + (t239 + (pkin(9) * t125 + t34) * t145) * qJD(5) + t182) * t146 + (t33 * t210 + pkin(9) * t82 - t57 * t104 + t237 + ((-pkin(9) * t222 + t103) * qJD(3) + t13) * qJD(4)) * t142, t254 * t125 + (t33 * t206 + (qJD(4) * t106 - t190) * pkin(9) + t160) * t146 + (-t33 * t211 + pkin(9) * t81 - t57 * t106 + t236 + (-pkin(9) * t228 - t219 * qJD(3) - t14) * qJD(4)) * t142, -t167 * t39 - t24 * t88, -t150 * t88 + t167 * t40 - t24 * t87 - t39 * t62, -t39 * t122 - t24 * t146 + (-qJD(3) * t88 - t167) * t213, t40 * t122 - t150 * t146 + (-qJD(3) * t87 - t62) * t213 (-t122 - t205) * t213, -t198 * t146 + t83 * t62 - t114 * t150 + t12 * t87 + t26 * t40 + (t245 * t140 + t246 * t144) * t122 + (t122 * t171 + t146 * t5) * qJD(6) + (-t57 * t62 + ((-t140 * t75 + t144 * t66) * qJD(3) + t4) * qJD(4)) * t142, t248 * t146 - t83 * t167 + t114 * t24 - t12 * t88 + t26 * t39 + ((qJD(6) * t66 + t245) * t144 + (-qJD(6) * t75 - t246) * t140) * t122 + (t57 * t167 + (-qJD(3) * t171 - t5) * qJD(4)) * t142; 0, 0, 0, 0, 0, -t146 * t149 * t142, t218 * t149, 0, 0, 0, qJD(4) * t36 - t55 * t216 - t17 (-qJD(3) * t55 - t53) * t146, -t106 * t228 + t234 (t81 + t230) * t145 + (-t82 + t229) * t141, -t125 * t210 + (t125 * t220 + (-t106 + t207) * t142) * qJD(3), t190 + (-t125 * t222 + (t104 + t206) * t142) * qJD(3), t125 * t216, -pkin(4) * t82 - t236 + t181 * t125 - t36 * t104 + (pkin(10) * t228 + t235) * qJD(5) + (-t13 * t142 + (-pkin(10) * t213 - t146 * t33) * t141) * qJD(3), -pkin(4) * t81 + t237 - t244 * t125 - t36 * t106 + (-t141 * pkin(10) * t125 + t33 * t145) * qJD(5) + (-t33 * t220 + (-pkin(10) * t206 + t14) * t142) * qJD(3), t108 * t24 - t167 * t243, -t107 * t24 + t108 * t150 + t167 * t242 - t243 * t62, -t243 * t122 + (qJD(4) * t108 + t167) * t216, t242 * t122 + (-qJD(4) * t107 + t62) * t216, t122 * t216, t12 * t107 - t130 * t150 + t177 * t62 + t242 * t26 + (t140 * t175 + t144 * t174) * t122 + ((-t119 * t144 - t120 * t140) * qJD(4) - t4) * t216, t12 * t108 + t130 * t24 - t177 * t167 + t243 * t26 + (-t140 * t174 + t144 * t175) * t122 + (-(-t119 * t140 + t120 * t144) * qJD(4) + t5) * t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106 * t104, -t104 ^ 2 + t106 ^ 2, t81 - t230, -t229 - t82, t127, -t33 * t106 - t14 * t125 + t151, t104 * t33 - t125 * t13 - t160, -t262, t260, t259, t257, t127 (-t10 * t140 - t240) * t122 + (-t106 * t62 + t122 * t209 + t127 * t144) * pkin(5) + t258, t263 + t9 + (t11 * t122 - t2) * t140 + (-t10 * t122 - t185) * t144 + (t106 * t167 + t122 * t208 - t127 * t140) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, t260, t259, t257, t127, -t122 * t5 + t258, -t122 * t4 - t248 + t263;];
tauc_reg  = t1;
