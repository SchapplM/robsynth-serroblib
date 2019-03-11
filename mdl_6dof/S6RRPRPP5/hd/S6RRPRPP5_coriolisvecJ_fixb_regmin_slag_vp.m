% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:17
% EndTime: 2019-03-09 10:06:25
% DurationCPUTime: 3.34s
% Computational Cost: add. (3150->389), mult. (6899->490), div. (0->0), fcn. (3806->4), ass. (0->205)
t149 = sin(qJ(2));
t223 = qJD(1) * t149;
t273 = -t223 - qJD(4);
t148 = sin(qJ(4));
t150 = cos(qJ(4));
t221 = qJD(2) * t148;
t151 = cos(qJ(2));
t222 = qJD(1) * t151;
t93 = t150 * t222 + t221;
t178 = t273 * t93;
t210 = qJD(1) * qJD(2);
t198 = t149 * t210;
t51 = qJD(4) * t93 - t148 * t198;
t266 = t51 + t178;
t199 = t148 * t222;
t219 = qJD(2) * t150;
t95 = -t199 + t219;
t244 = t273 * t95;
t215 = qJD(4) * t150;
t52 = qJD(2) * t215 - qJD(4) * t199 - t150 * t198;
t277 = t52 + t244;
t281 = t148 * t277 - t266 * t150;
t280 = t51 - t178;
t203 = t150 * t223;
t279 = t203 + t215;
t260 = pkin(3) + pkin(7);
t278 = -t52 + t244;
t197 = t151 * t210;
t259 = pkin(4) + pkin(5);
t276 = t259 * t197;
t274 = qJ(5) * t197 - qJD(5) * t273;
t134 = pkin(7) * t223;
t229 = pkin(3) * t223 + qJD(3) + t134;
t272 = -0.2e1 * t210;
t271 = t52 * qJ(6) + t93 * qJD(6);
t122 = t273 * qJ(5);
t153 = -pkin(2) - pkin(8);
t196 = -qJ(3) * t149 - pkin(1);
t87 = t153 * t151 + t196;
t65 = t87 * qJD(1);
t67 = t153 * qJD(2) + t229;
t28 = t148 * t67 + t150 * t65;
t24 = -t122 + t28;
t186 = pkin(4) * t197;
t216 = qJD(4) * t148;
t130 = pkin(2) * t198;
t183 = pkin(8) * t149 - qJ(3) * t151;
t217 = qJD(3) * t149;
t158 = qJD(2) * t183 - t217;
t44 = qJD(1) * t158 + t130;
t128 = pkin(7) * t197;
t86 = pkin(3) * t197 + t128;
t195 = -t148 * t44 + t150 * t86 - t65 * t215 - t67 * t216;
t8 = -t186 - t195;
t270 = -t24 * t273 - t8;
t136 = pkin(7) * t222;
t101 = pkin(3) * t222 + t136;
t145 = qJD(2) * qJ(3);
t79 = t145 + t101;
t175 = qJ(5) * t95 - t79;
t21 = -t259 * t93 + qJD(6) + t175;
t269 = (qJD(6) + t21) * t95;
t268 = t279 * t273;
t124 = t273 ^ 2;
t92 = t95 ^ 2;
t267 = -t124 - t92;
t27 = -t148 * t65 + t150 * t67;
t227 = qJD(5) - t27;
t264 = 0.2e1 * t274;
t263 = qJD(5) * t150 - t229;
t238 = qJ(5) * t150;
t262 = t259 * t148 - t238;
t239 = qJ(5) * t148;
t166 = t259 * t150 + t239;
t261 = t93 ^ 2;
t31 = pkin(4) * t93 - t175;
t258 = t31 * t95;
t257 = t95 * t93;
t230 = qJ(6) + t153;
t188 = qJD(4) * t230;
t139 = pkin(2) * t223;
t75 = qJD(1) * t183 + t139;
t194 = t101 * t150 - t148 * t75;
t235 = t148 * t149;
t208 = qJ(6) * t235;
t211 = qJD(6) * t150;
t256 = (-t259 * t151 - t208) * qJD(1) - t194 - t148 * t188 + t211;
t212 = qJD(6) * t148;
t252 = t148 * t101 + t150 * t75;
t33 = qJ(5) * t222 + t252;
t255 = qJ(6) * t203 + t150 * t188 + t212 - t33;
t254 = -t273 * t166 - t263;
t185 = pkin(4) * t150 + t239;
t253 = t185 * t273 + t263;
t114 = t260 * t149;
t251 = t148 * t114 + t150 * t87;
t250 = qJ(5) * t52;
t249 = qJ(5) * t93;
t248 = qJD(2) * pkin(2);
t144 = qJD(2) * qJD(3);
t105 = pkin(7) * t198 - t144;
t72 = -pkin(3) * t198 - t105;
t11 = t52 * pkin(4) + t51 * qJ(5) - t95 * qJD(5) + t72;
t247 = t11 * t148;
t246 = t11 * t150;
t243 = t150 * t51;
t242 = t151 * t95;
t241 = t72 * t148;
t240 = t72 * t150;
t237 = qJ(6) * t151;
t236 = t273 * t153;
t234 = t149 * t150;
t155 = qJD(1) ^ 2;
t233 = t151 * t155;
t154 = qJD(2) ^ 2;
t232 = t154 * t149;
t231 = t154 * t151;
t19 = qJ(6) * t95 + t27;
t228 = qJD(5) - t19;
t115 = t260 * t151;
t146 = t149 ^ 2;
t147 = t151 ^ 2;
t224 = t146 - t147;
t111 = -pkin(2) * t151 + t196;
t80 = qJD(1) * t111;
t220 = qJD(2) * t149;
t218 = qJD(2) * t151;
t214 = qJD(4) * t151;
t213 = qJD(5) * t148;
t209 = -t148 * t86 - t150 * t44 - t67 * t215;
t37 = t149 * qJ(5) + t251;
t207 = t148 * t236;
t206 = t149 * t233;
t205 = t150 * t236;
t204 = t148 * t223;
t202 = t153 * t218;
t103 = t273 * t216;
t201 = t148 * t214;
t200 = t150 * t214;
t193 = t114 * t150 - t148 * t87;
t192 = pkin(1) * t272;
t191 = qJD(3) - t248;
t117 = t148 * t197;
t190 = -qJD(2) * t95 - t117;
t119 = t150 * t197;
t189 = -qJD(2) * t93 + t119;
t20 = qJ(6) * t93 + t28;
t184 = -pkin(4) * t148 + t238;
t22 = pkin(4) * t273 + t227;
t182 = t148 * t22 + t150 * t24;
t102 = t260 * t218;
t138 = pkin(2) * t220;
t59 = t138 + t158;
t181 = t102 * t150 - t114 * t216 - t148 * t59 - t87 * t215;
t180 = -0.2e1 * qJD(2) * t80;
t177 = -qJD(1) * t147 - t149 * t273;
t176 = t273 * t148;
t174 = t204 * t273 + t103 + t189;
t173 = -t190 - t268;
t161 = -qJ(3) * t218 - t217;
t62 = qJD(1) * t161 + t130;
t78 = t138 + t161;
t172 = pkin(7) * t154 + qJD(1) * t78 + t62;
t5 = -pkin(5) * t52 - t11;
t171 = -t148 * t5 - t21 * t215;
t170 = -t150 * t5 + t21 * t216;
t169 = qJ(6) * t51 - t195;
t168 = -t273 * t28 + t195;
t165 = -t197 + t257;
t164 = t216 * t65 + t209;
t163 = t148 * t102 + t114 * t215 + t150 * t59 - t216 * t87;
t14 = t259 * t273 + t228;
t17 = -t122 + t20;
t4 = -t164 + t274;
t2 = t4 + t271;
t157 = t169 - t276;
t3 = -qJD(6) * t95 + t157;
t162 = t2 * t148 - t150 * t3 + t279 * t17 + (t204 + t216) * t14;
t10 = qJ(5) * t218 + t149 * qJD(5) + t163;
t160 = -t27 * t273 + t164;
t109 = t134 + t191;
t112 = -t136 - t145;
t156 = -t105 * t151 + (t109 * t151 + (t112 + t136) * t149) * qJD(2);
t110 = qJ(3) - t184;
t108 = t230 * t150;
t107 = t230 * t148;
t106 = t153 * t119;
t100 = t260 * t220;
t98 = -qJ(3) * t222 + t139;
t83 = -qJ(3) - t262;
t66 = t80 * t223;
t58 = t151 * t185 + t115;
t42 = pkin(4) * t95 + t249;
t41 = -t151 * t166 - t115;
t38 = -pkin(4) * t149 - t193;
t35 = -pkin(4) * t222 - t194;
t34 = t150 * t237 + t37;
t32 = -t259 * t95 - t249;
t29 = t148 * t237 - t259 * t149 - t193;
t25 = (qJD(4) * t184 + t213) * t151 + (-t185 - t260) * t220;
t18 = (t262 * qJD(4) - t213) * t151 + (t166 + t260) * t220;
t13 = -pkin(4) * t218 - t181;
t7 = t151 * t211 + (-t149 * t219 - t201) * qJ(6) + t10;
t6 = -qJD(2) * t208 + (qJ(6) * t215 - t259 * qJD(2) + t212) * t151 - t181;
t1 = [0, 0, 0, 0.2e1 * t149 * t197, t224 * t272, t231, -t232, 0, -pkin(7) * t231 + t149 * t192, pkin(7) * t232 + t151 * t192, t156, t149 * t180 + t151 * t172, -t149 * t172 + t151 * t180, pkin(7) * t156 + t111 * t62 + t78 * t80, -t95 * t200 + (t151 * t51 + t220 * t95) * t148 (-t148 * t93 + t150 * t95) * t220 + (t148 * t52 + t243 + (t148 * t95 + t150 * t93) * qJD(4)) * t151, t273 * t200 - t149 * t51 + (t148 * t177 + t242) * qJD(2), -t273 * t201 - t149 * t52 + (t150 * t177 - t151 * t93) * qJD(2) (-t273 + t223) * t218, -t181 * t273 - t100 * t93 + t115 * t52 + (-t219 * t79 + t195) * t149 + (-t79 * t216 + t240 + (qJD(1) * t193 + t27) * qJD(2)) * t151, t163 * t273 - t100 * t95 - t115 * t51 + ((qJD(2) * t79 + qJD(4) * t65) * t148 + t209) * t149 + (-t79 * t215 - t241 + (-qJD(1) * t251 - t28) * qJD(2)) * t151, t13 * t273 + t25 * t93 + t52 * t58 + (-t219 * t31 - t8) * t149 + (-t31 * t216 + t246 + (-qJD(1) * t38 - t22) * qJD(2)) * t151, -t10 * t93 + t13 * t95 - t37 * t52 - t38 * t51 + t182 * t220 + (-t148 * t8 - t150 * t4 + (t148 * t24 - t150 * t22) * qJD(4)) * t151, -t10 * t273 - t25 * t95 + t51 * t58 + (-t221 * t31 + t4) * t149 + (t31 * t215 + t247 + (qJD(1) * t37 + t24) * qJD(2)) * t151, t10 * t24 + t11 * t58 + t13 * t22 + t25 * t31 + t37 * t4 + t38 * t8, t273 * t6 - t18 * t93 - t41 * t52 + (t21 * t219 - t3) * t149 + ((-qJD(1) * t29 - t14) * qJD(2) + t170) * t151, -t273 * t7 + t18 * t95 - t41 * t51 + (t21 * t221 + t2) * t149 + ((qJD(1) * t34 + t17) * qJD(2) + t171) * t151, t29 * t51 + t34 * t52 - t6 * t95 + t7 * t93 + (-t14 * t148 - t150 * t17) * t220 + (t148 * t3 + t150 * t2 + (t14 * t150 - t148 * t17) * qJD(4)) * t151, t14 * t6 + t17 * t7 + t18 * t21 + t2 * t34 + t29 * t3 + t41 * t5; 0, 0, 0, -t206, t224 * t155, 0, 0, 0, t155 * pkin(1) * t149, pkin(1) * t233 ((-t112 - t145) * t149 + (-t109 + t191) * t151) * qJD(1), -t222 * t98 + t66, 0.2e1 * t144 + (t149 * t98 + t151 * t80) * qJD(1), -qJ(3) * t105 - qJD(3) * t112 - t80 * t98 + (-t112 * t149 + (-t109 - t248) * t151) * qJD(1) * pkin(7), t176 * t95 - t243, t280 * t148 + t278 * t150, t103 + t119 + (t235 * t273 - t242) * qJD(1), t222 * t93 - t117 + t268, t273 * t222, t106 + qJ(3) * t52 + t241 + t194 * t273 + t229 * t93 + (t150 * t79 + t207) * qJD(4) + (-t27 * t151 + t234 * t79) * qJD(1), -qJ(3) * t51 + t240 - t252 * t273 + t229 * t95 + (-t79 * t148 + t205) * qJD(4) + (t28 * t151 + (-t149 * t79 - t202) * t148) * qJD(1), t247 + t110 * t52 - t273 * t35 + t106 - t253 * t93 + (t150 * t31 + t207) * qJD(4) + (t151 * t22 + t234 * t31) * qJD(1), t33 * t93 - t35 * t95 + (-t24 * t223 + t153 * t51 + t8 + (-t153 * t93 - t24) * qJD(4)) * t150 + (-t22 * t223 - t153 * t52 - t4 + (t153 * t95 - t22) * qJD(4)) * t148, -t246 + t110 * t51 + t273 * t33 + t253 * t95 + (t148 * t31 - t205) * qJD(4) + (-t151 * t24 + (t149 * t31 + t202) * t148) * qJD(1), t11 * t110 - t22 * t35 - t24 * t33 - t253 * t31 + (qJD(4) * t182 + t148 * t4 - t150 * t8) * t153, -t52 * t83 + t254 * t93 - t256 * t273 + (-t21 * t234 + (qJD(2) * t108 + t14) * t151) * qJD(1) + t171, -t51 * t83 - t254 * t95 - t255 * t273 + (-t21 * t235 + (qJD(2) * t107 - t17) * t151) * qJD(1) - t170, t107 * t52 - t108 * t51 + t255 * t93 + t256 * t95 + t162, t107 * t2 - t108 * t3 - t256 * t14 + t255 * t17 - t254 * t21 + t5 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, -t146 * t155 - t154, qJD(2) * t112 + t128 + t66, 0, 0, 0, 0, 0, -t176 * t273 + t189, -t124 * t150 + t190, t174, -t281, t173, -qJD(2) * t31 + t270 * t150 + (-t22 * t273 + t4) * t148, t174, t173, t281, qJD(2) * t21 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257, t92 - t261, -t266, -t277, t197, -t79 * t95 + t168, t79 * t93 + t160, -t42 * t93 + t168 + 0.2e1 * t186 - t258, pkin(4) * t51 - t250 + (t24 - t28) * t95 + (t22 - t227) * t93, -t31 * t93 + t42 * t95 - t160 + t264, -pkin(4) * t8 + qJ(5) * t4 - t22 * t28 + t227 * t24 - t31 * t42, -t20 * t273 + t32 * t93 - t169 + t269 + 0.2e1 * t276, t19 * t273 + t21 * t93 - t32 * t95 - t164 + t264 + t271, t250 - t259 * t51 + (-t17 + t20) * t95 + (-t14 + t228) * t93, qJ(5) * t2 - t14 * t20 + t17 * t228 - t21 * t32 - t259 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, -t266, t267, t258 - t270, t165, t267, t266, t17 * t273 + t157 - t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, -t280, -t92 - t261, t14 * t95 - t17 * t93 + t5;];
tauc_reg  = t1;
