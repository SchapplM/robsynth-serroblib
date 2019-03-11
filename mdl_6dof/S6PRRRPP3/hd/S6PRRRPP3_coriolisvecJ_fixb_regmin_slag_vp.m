% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRPP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:29
% EndTime: 2019-03-08 22:58:37
% DurationCPUTime: 2.90s
% Computational Cost: add. (3180->408), mult. (7792->518), div. (0->0), fcn. (5333->8), ass. (0->203)
t137 = sin(qJ(3));
t140 = cos(qJ(3));
t103 = -t140 * pkin(3) - t137 * pkin(9) - pkin(2);
t136 = sin(qJ(4));
t138 = sin(qJ(2));
t139 = cos(qJ(4));
t212 = qJD(4) * t139;
t133 = sin(pkin(6));
t219 = qJD(1) * t133;
t141 = cos(qJ(2));
t227 = t140 * t141;
t167 = pkin(3) * t137 - pkin(9) * t140;
t97 = t167 * qJD(3);
t278 = -(t136 * t138 + t139 * t227) * t219 + t103 * t212 + t136 * t97;
t208 = qJD(2) * qJD(3);
t124 = t137 * t208;
t260 = pkin(4) + qJ(6);
t277 = t260 * t124;
t191 = t137 * t212;
t211 = t136 * qJD(3);
t276 = t140 * t211 + t191;
t134 = cos(pkin(6));
t232 = t134 * t137;
t116 = qJD(1) * t232;
t190 = t138 * t219;
t99 = qJD(2) * pkin(8) + t190;
t73 = t140 * t99 + t116;
t62 = qJD(3) * pkin(9) + t73;
t189 = t141 * t219;
t74 = qJD(2) * t103 - t189;
t23 = t136 * t62 - t139 * t74;
t216 = qJD(2) * t137;
t92 = t139 * t216 + t211;
t162 = t92 * pkin(5) + t23;
t223 = qJD(5) + t162;
t209 = t140 * qJD(2);
t119 = -qJD(4) + t209;
t24 = t136 * t74 + t139 * t62;
t19 = t119 * qJ(5) - t24;
t175 = pkin(4) * t124;
t213 = qJD(4) * t136;
t217 = qJD(2) * t133;
t195 = t141 * t217;
t215 = qJD(3) * t137;
t218 = qJD(1) * t140;
t41 = -t99 * t215 + (qJD(3) * t134 + t195) * t218;
t71 = (t97 + t190) * qJD(2);
t180 = -t136 * t41 + t139 * t71 - t62 * t212 - t74 * t213;
t5 = -t175 - t180;
t275 = -t19 * t119 + t5;
t228 = t139 * t140;
t123 = pkin(8) * t228;
t274 = qJD(4) * t123 + t103 * t213 - t139 * t97 - (t136 * t227 - t138 * t139) * t219;
t273 = t140 * qJD(5) - t278;
t72 = t134 * t218 - t137 * t99;
t117 = t119 ^ 2;
t88 = t92 ^ 2;
t272 = -t88 - t117;
t108 = qJD(5) * t119;
t118 = qJ(5) * t124;
t271 = t118 - t108;
t270 = t136 * qJD(5) + t73 + (t136 * t209 - t213) * pkin(4);
t196 = t138 * t217;
t234 = t133 * t141;
t199 = t139 * t234;
t172 = t140 * t195;
t235 = t133 * t138;
t79 = -t134 * t140 + t137 * t235;
t49 = -qJD(3) * t79 + t172;
t80 = t140 * t235 + t232;
t12 = -qJD(4) * t199 + t136 * t196 + t49 * t139 - t80 * t213;
t173 = t137 * t195;
t50 = qJD(3) * t80 + t173;
t52 = -t136 * t234 + t80 * t139;
t187 = t136 * t216;
t210 = t139 * qJD(3);
t194 = t140 * t210;
t207 = qJD(3) * qJD(4);
t67 = -qJD(2) * t194 + qJD(4) * t187 - t139 * t207;
t147 = -t12 * t119 + t52 * t124 - t50 * t92 + t79 * t67;
t11 = qJD(4) * t52 + t49 * t136 - t139 * t196;
t51 = t80 * t136 + t199;
t68 = qJD(2) * t276 + t136 * t207;
t90 = t187 - t210;
t148 = t11 * t119 - t51 * t124 + t50 * t90 + t79 * t68;
t269 = t90 ^ 2;
t268 = pkin(5) + pkin(9);
t267 = t90 * pkin(5);
t61 = -qJD(3) * pkin(3) - t72;
t146 = -t92 * qJ(5) + t61;
t17 = t260 * t90 + t146;
t266 = t17 * t92;
t25 = t90 * pkin(4) + t146;
t265 = t25 * t92;
t171 = t137 * t189;
t214 = qJD(3) * t140;
t42 = qJD(2) * t171 + qJD(3) * t116 + t99 * t214;
t151 = t67 * qJ(5) - t92 * qJD(5) + t42;
t6 = t68 * pkin(4) + t151;
t264 = t6 * t136;
t263 = t6 * t139;
t7 = t260 * t119 + t223;
t262 = t7 * t119;
t261 = t92 * t90;
t192 = t137 * t213;
t197 = -pkin(8) * t136 - pkin(4);
t205 = pkin(5) * t228;
t259 = -pkin(5) * t192 + t140 * qJD(6) + (t205 + (-qJ(6) + t197) * t137) * qJD(3) + t274;
t230 = t136 * t140;
t122 = pkin(8) * t230;
t206 = pkin(5) * t230;
t229 = t137 * t139;
t258 = (-pkin(5) * t229 - t122) * qJD(4) + (-t206 + (-pkin(8) * t139 + qJ(5)) * t137) * qJD(3) - t273;
t105 = t268 * t139;
t94 = t167 * qJD(2);
t182 = -t136 * t72 + t139 * t94;
t257 = -(-t260 * t137 + t205) * qJD(2) + t182 + qJD(4) * t105;
t256 = -qJ(5) * t215 + (t137 * t210 + t140 * t213) * pkin(8) + t273;
t251 = t136 * t94 + t139 * t72;
t255 = (qJ(5) * t137 - t206) * qJD(2) + t251 + t268 * t213;
t254 = t197 * t215 + t274;
t238 = qJ(5) * t139;
t164 = qJ(6) * t136 - t238;
t155 = t164 * t140;
t253 = qJD(2) * t155 - t164 * qJD(4) + t139 * qJD(6) + t270;
t252 = qJ(5) * t212 - t209 * t238 + t270;
t249 = qJ(5) * t68;
t248 = qJ(5) * t90;
t247 = qJD(2) * pkin(2);
t246 = t119 * t90;
t245 = t119 * t92;
t243 = t42 * t136;
t242 = t42 * t139;
t241 = t61 * t136;
t240 = t67 * t136;
t239 = t136 * t103 + t123;
t237 = t119 * t136;
t236 = t119 * t139;
t143 = qJD(2) ^ 2;
t233 = t133 * t143;
t231 = t136 * t137;
t142 = qJD(3) ^ 2;
t226 = t142 * t137;
t225 = t142 * t140;
t224 = -qJD(5) - t23;
t16 = t24 - t267;
t222 = -qJD(6) - t16;
t221 = pkin(4) * t231 + t137 * pkin(8);
t131 = t137 ^ 2;
t220 = -t140 ^ 2 + t131;
t204 = pkin(9) * t237;
t203 = pkin(9) * t236;
t202 = pkin(9) * t215;
t201 = pkin(9) * t210;
t198 = t138 * t233;
t193 = t119 * t213;
t184 = -t136 * qJ(5) - pkin(3);
t183 = -t124 + t261;
t181 = -t136 * t71 - t139 * t41 - t74 * t212 + t62 * t213;
t179 = t139 * t103 - t122;
t177 = t90 * t189;
t176 = t92 * t189;
t174 = t276 * pkin(4) + pkin(8) * t214 + qJ(5) * t192;
t169 = -t108 - t181;
t65 = t140 * qJ(5) - t239;
t100 = -t189 - t247;
t166 = -t100 - t189;
t18 = t119 * pkin(4) - t224;
t165 = t136 * t19 + t139 * t18;
t163 = qJD(2) * t131 - t119 * t140;
t161 = -t67 * pkin(5) - t180;
t160 = -t68 * pkin(5) - t181;
t3 = t90 * qJD(6) + t260 * t68 + t151;
t159 = t3 * t136 + t17 * t212;
t158 = -t3 * t139 + t17 * t213;
t157 = -t24 * t119 + t180;
t156 = t11 * t92 - t12 * t90 - t51 * t67 - t52 * t68;
t150 = -t17 * t90 + t160;
t149 = qJD(3) * (-t166 - t247);
t32 = -t67 - t246;
t145 = t161 - t277;
t144 = -t245 - t68;
t130 = t140 * pkin(4);
t115 = 0.2e1 * t118;
t104 = t268 * t136;
t101 = -t139 * pkin(4) + t184;
t83 = -t260 * t139 + t184;
t76 = -qJ(5) * t229 + t221;
t66 = t130 - t179;
t60 = t164 * t137 + t221;
t47 = pkin(4) * t92 + t248;
t44 = -pkin(5) * t231 - t65;
t43 = t140 * qJ(6) + t122 + t130 + (pkin(5) * t137 - t103) * t139;
t33 = t260 * t92 + t248;
t31 = (-qJ(5) * t214 - qJD(5) * t137) * t139 + t174;
t30 = -pkin(4) * t216 - t182;
t29 = -qJ(5) * t216 - t251;
t14 = qJD(3) * t155 + (qJD(6) * t136 + (qJ(6) * qJD(4) - qJD(5)) * t139) * t137 + t174;
t9 = qJD(6) - t19 - t267;
t4 = -t118 - t169;
t2 = t160 + t271;
t1 = qJD(6) * t119 + t145;
t8 = [0, 0, -t198, -t141 * t233, 0, 0, 0, 0, 0, -t140 * t198 + (-t50 - t173) * qJD(3), t137 * t198 + (-t49 - t172) * qJD(3), 0, 0, 0, 0, 0, t148, -t147, t156, -t148, t147, t11 * t18 - t12 * t19 + t25 * t50 - t4 * t52 + t5 * t51 + t6 * t79, t156, t147, t148, t1 * t51 + t11 * t7 + t12 * t9 + t17 * t50 + t2 * t52 + t3 * t79; 0, 0, 0, 0, 0.2e1 * t140 * t124, -0.2e1 * t220 * t208, t225, -t226, 0, -pkin(8) * t225 + t137 * t149, pkin(8) * t226 + t140 * t149, t92 * t194 + (-t67 * t139 - t92 * t213) * t137 (-t136 * t92 - t139 * t90) * t214 + (t240 - t139 * t68 + (t136 * t90 - t139 * t92) * qJD(4)) * t137, t119 * t192 + t67 * t140 + (t137 * t92 + t139 * t163) * qJD(3), t119 * t191 + t68 * t140 + (-t136 * t163 - t137 * t90) * qJD(3) (-t119 - t209) * t215, t274 * t119 + ((pkin(8) * t90 + t241) * qJD(3) - t180) * t140 + (-t177 + t61 * t212 + pkin(8) * t68 + t243 + (-pkin(8) * t237 + qJD(2) * t179 - t23) * qJD(3)) * t137, t278 * t119 + (t61 * t210 + (qJD(3) * t92 - t193) * pkin(8) - t181) * t140 + (-t176 - t61 * t213 - pkin(8) * t67 + t242 + (-pkin(8) * t236 - t239 * qJD(2) - t24) * qJD(3)) * t137, t65 * t68 - t66 * t67 + t254 * t92 + t256 * t90 + t165 * t214 + (t136 * t4 + t139 * t5 + (-t136 * t18 + t139 * t19) * qJD(4)) * t137, -t31 * t90 - t76 * t68 + (-t211 * t25 - t5) * t140 - t254 * t119 + (t177 - t25 * t212 - t264 + (qJD(2) * t66 + t18) * qJD(3)) * t137, -t31 * t92 + t76 * t67 + (-t210 * t25 + t4) * t140 + t256 * t119 + (t176 + t25 * t213 - t263 + (-qJD(2) * t65 - t19) * qJD(3)) * t137, t4 * t65 + t5 * t66 + t6 * t76 + (t31 - t171) * t25 + t256 * t19 + t254 * t18, -t43 * t67 - t44 * t68 + t259 * t92 - t258 * t90 + (-t136 * t9 + t139 * t7) * t214 + (t1 * t139 - t136 * t2 + (-t136 * t7 - t139 * t9) * qJD(4)) * t137, -t14 * t92 + t60 * t67 + (-t17 * t210 - t2) * t140 - t258 * t119 + (t176 + (qJD(2) * t44 + t9) * qJD(3) + t158) * t137, t14 * t90 + t60 * t68 + (t17 * t211 + t1) * t140 + t259 * t119 + (-t177 + (-qJD(2) * t43 - t7) * qJD(3) + t159) * t137, t1 * t43 + t2 * t44 + t3 * t60 + t258 * t9 + t259 * t7 + (t14 - t171) * t17; 0, 0, 0, 0, -t137 * t143 * t140, t220 * t143, 0, 0, 0, t73 * qJD(3) - t100 * t216 - t42, t166 * t209, -t92 * t236 - t240 (-t67 + t246) * t139 + (-t68 + t245) * t136, -t119 * t212 + (t119 * t228 + (-t92 + t211) * t137) * qJD(2), t193 + (-t119 * t230 + (t90 + t210) * t137) * qJD(2), t119 * t216, -pkin(3) * t68 - t242 + t182 * t119 - t73 * t90 + (t203 + t241) * qJD(4) + (t23 * t137 + (-t140 * t61 - t202) * t136) * qJD(2), pkin(3) * t67 + t243 - t251 * t119 - t73 * t92 + (t61 * t139 - t204) * qJD(4) + (-t61 * t228 + (t24 - t201) * t137) * qJD(2), -t29 * t90 - t30 * t92 + (-t4 - t119 * t18 + (qJD(4) * t92 - t68) * pkin(9)) * t139 + ((qJD(4) * t90 - t67) * pkin(9) + t275) * t136, -t101 * t68 + t30 * t119 + t263 + t252 * t90 + (-t136 * t25 - t203) * qJD(4) + (-t137 * t18 + (t140 * t25 + t202) * t136) * qJD(2), t101 * t67 - t29 * t119 - t264 + t252 * t92 + (-t139 * t25 + t204) * qJD(4) + (t25 * t228 + (t19 + t201) * t137) * qJD(2), t6 * t101 - t18 * t30 - t19 * t29 - t252 * t25 + (qJD(4) * t165 + t5 * t136 - t4 * t139) * pkin(9), -t104 * t67 - t105 * t68 + t257 * t92 + t255 * t90 + (t2 - t262) * t139 + (t9 * t119 + t1) * t136, t83 * t67 + t253 * t92 + t255 * t119 + (t17 * t228 + (qJD(3) * t105 - t9) * t137) * qJD(2) - t159, t83 * t68 - t253 * t90 + t257 * t119 + (-t17 * t230 + (-qJD(3) * t104 + t7) * t137) * qJD(2) + t158, t1 * t104 + t2 * t105 - t253 * t17 - t255 * t9 + t257 * t7 + t3 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, t88 - t269, t32, t144, t124, -t61 * t92 + t157, t23 * t119 + t61 * t90 + t181, pkin(4) * t67 - t249 + (-t19 - t24) * t92 + (t18 + t224) * t90, t47 * t90 - t157 - 0.2e1 * t175 + t265, t119 * t224 - t25 * t90 + t47 * t92 + t115 + t169, -pkin(4) * t5 - qJ(5) * t4 - t18 * t24 + t19 * t224 - t25 * t47, -t249 + t260 * t67 + (t9 + t222) * t92 + (t7 - t223) * t90, -t119 * t162 + t33 * t92 - 0.2e1 * t108 + t115 + t150, -t266 - t33 * t90 + (-0.2e1 * qJD(6) - t16) * t119 + 0.2e1 * t277 - t161, t2 * qJ(5) - t1 * t260 - t17 * t33 + t222 * t7 + t223 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t183, t272, t265 + t275, t32, t272, t183, t266 + (qJD(6) + t9) * t119 + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t124 + t261, -t117 - t269, t150 - t262 + t271;];
tauc_reg  = t8;
