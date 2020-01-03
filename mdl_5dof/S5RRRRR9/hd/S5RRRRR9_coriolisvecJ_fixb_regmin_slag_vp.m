% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:52
% EndTime: 2019-12-31 22:30:08
% DurationCPUTime: 5.46s
% Computational Cost: add. (4615->386), mult. (11589->564), div. (0->0), fcn. (8425->8), ass. (0->203)
t185 = sin(qJ(4));
t191 = cos(qJ(2));
t256 = qJD(1) * t191;
t169 = -qJD(3) + t256;
t190 = cos(qJ(3));
t186 = sin(qJ(3));
t255 = qJD(2) * t186;
t187 = sin(qJ(2));
t257 = qJD(1) * t187;
t141 = t190 * t257 + t255;
t152 = -pkin(2) * t191 - pkin(7) * t187 - pkin(1);
t132 = t152 * qJD(1);
t178 = pkin(6) * t256;
t158 = qJD(2) * pkin(7) + t178;
t96 = t190 * t132 - t158 * t186;
t66 = -pkin(8) * t141 + t96;
t58 = -pkin(3) * t169 + t66;
t189 = cos(qJ(4));
t239 = t186 * t257;
t246 = t190 * qJD(2);
t139 = t239 - t246;
t274 = t186 * t132;
t97 = t158 * t190 + t274;
t67 = -pkin(8) * t139 + t97;
t63 = t189 * t67;
t20 = t185 * t58 + t63;
t88 = t189 * t139 + t141 * t185;
t312 = pkin(9) * t88;
t15 = t20 - t312;
t184 = sin(qJ(5));
t248 = qJD(5) * t184;
t13 = t15 * t248;
t188 = cos(qJ(5));
t209 = t139 * t185 - t189 * t141;
t45 = t184 * t209 - t188 * t88;
t157 = -qJD(2) * pkin(2) + pkin(6) * t257;
t109 = pkin(3) * t139 + t157;
t56 = pkin(4) * t88 + t109;
t316 = -t56 * t45 + t13;
t245 = qJD(1) * qJD(2);
t172 = t187 * t245;
t231 = t191 * t245;
t252 = qJD(3) * t186;
t236 = t187 * t252;
t244 = qJD(2) * qJD(3);
t106 = -qJD(1) * t236 + (t231 + t244) * t190;
t213 = pkin(2) * t187 - pkin(7) * t191;
t150 = t213 * qJD(2);
t133 = qJD(1) * t150;
t217 = pkin(6) * t172;
t264 = -t190 * t133 - t186 * t217;
t199 = -t97 * qJD(3) - t264;
t27 = pkin(3) * t172 - pkin(8) * t106 + t199;
t253 = qJD(2) * t191;
t238 = t186 * t253;
t251 = qJD(3) * t190;
t200 = t187 * t251 + t238;
t107 = t200 * qJD(1) + t186 * t244;
t206 = t132 * t251 + t186 * t133 - t158 * t252;
t194 = -t190 * t217 + t206;
t37 = -pkin(8) * t107 + t194;
t227 = -t185 * t37 + t189 * t27;
t197 = -t20 * qJD(4) + t227;
t249 = qJD(4) * t189;
t250 = qJD(4) * t185;
t32 = t189 * t106 - t185 * t107 - t139 * t249 - t141 * t250;
t2 = pkin(4) * t172 - pkin(9) * t32 + t197;
t315 = -t184 * t2 + t316;
t211 = t184 * t88 + t188 * t209;
t308 = t211 * t45;
t305 = t211 ^ 2 - t45 ^ 2;
t162 = -qJD(4) + t169;
t155 = -qJD(5) + t162;
t195 = t209 * qJD(4) - t185 * t106 - t189 * t107;
t247 = qJD(5) * t188;
t8 = t184 * t195 + t188 * t32 + t209 * t248 - t88 * t247;
t304 = t155 * t45 + t8;
t224 = -t185 * t27 - t189 * t37 - t58 * t249 + t67 * t250;
t3 = pkin(9) * t195 - t224;
t241 = -t184 * t3 + t188 * t2;
t314 = t56 * t211 + t241;
t196 = t211 * qJD(5) - t184 * t32 + t188 * t195;
t299 = t155 * t211 + t196;
t147 = t213 * qJD(1);
t126 = t186 * t147;
t292 = pkin(7) + pkin(8);
t240 = qJD(3) * t292;
t271 = t187 * t190;
t272 = t186 * t191;
t307 = -t186 * t240 - t126 - (-pkin(6) * t271 - pkin(8) * t272) * qJD(1);
t270 = t190 * t191;
t207 = pkin(3) * t187 - pkin(8) * t270;
t260 = pkin(6) * t239 + t190 * t147;
t313 = t207 * qJD(1) + t190 * t240 + t260;
t311 = pkin(9) * t209;
t309 = t209 * t88;
t142 = t185 * t186 - t189 * t190;
t202 = t142 * t191;
t293 = qJD(3) + qJD(4);
t266 = qJD(1) * t202 - t293 * t142;
t143 = t185 * t190 + t186 * t189;
t265 = (-t256 + t293) * t143;
t306 = t209 ^ 2 - t88 ^ 2;
t303 = -t162 * t88 + t32;
t302 = t109 * t88 + t224;
t61 = t185 * t67;
t19 = t189 * t58 - t61;
t14 = t19 + t311;
t12 = -pkin(4) * t162 + t14;
t283 = t188 * t15;
t5 = t184 * t12 + t283;
t301 = -t5 * qJD(5) + t314;
t300 = t109 * t209 + t197;
t298 = t162 * t209 + t195;
t297 = -0.2e1 * t245;
t296 = t313 * t189;
t119 = t143 * t187;
t214 = -t178 + (-t186 * t256 + t252) * pkin(3);
t159 = t292 * t186;
t160 = t292 * t190;
t261 = -t185 * t159 + t189 * t160;
t295 = -t159 * t249 - t160 * t250 - t313 * t185 + t307 * t189;
t294 = t191 * t246 - t236;
t291 = pkin(6) * t186;
t92 = t188 * t142 + t143 * t184;
t290 = -t92 * qJD(5) - t265 * t184 + t266 * t188;
t93 = -t142 * t184 + t143 * t188;
t289 = t93 * qJD(5) + t266 * t184 + t265 * t188;
t288 = t189 * t66 - t61;
t287 = t265 * pkin(4) + t214;
t171 = pkin(6) * t270;
t259 = t186 * t152 + t171;
t273 = t186 * t187;
t102 = -pkin(8) * t273 + t259;
t138 = t190 * t152;
t95 = -pkin(8) * t271 + t138 + (-pkin(3) - t291) * t191;
t285 = t189 * t102 + t185 * t95;
t284 = t188 * t12;
t282 = t106 * t186;
t281 = t139 * t169;
t280 = t141 * t169;
t279 = t157 * t186;
t278 = t157 * t190;
t277 = t169 * t190;
t276 = t184 * t185;
t275 = t185 * t188;
t193 = qJD(1) ^ 2;
t269 = t191 * t193;
t192 = qJD(2) ^ 2;
t268 = t192 * t187;
t267 = t192 * t191;
t263 = t186 * t150 + t152 * t251;
t254 = qJD(2) * t187;
t262 = t190 * t150 + t254 * t291;
t151 = pkin(3) * t273 + t187 * pkin(6);
t182 = t187 ^ 2;
t258 = -t191 ^ 2 + t182;
t242 = pkin(3) * qJD(4) * t155;
t110 = t200 * pkin(3) + pkin(6) * t253;
t176 = -pkin(3) * t190 - pkin(2);
t235 = t191 * t252;
t233 = t169 * t251;
t85 = pkin(3) * t107 + pkin(6) * t231;
t229 = qJD(5) * t12 + t3;
t49 = t207 * qJD(2) + (-t171 + (pkin(8) * t187 - t152) * t186) * qJD(3) + t262;
t52 = -t200 * pkin(8) + (-t187 * t246 - t235) * pkin(6) + t263;
t226 = -t185 * t52 + t189 * t49;
t225 = -t185 * t66 - t63;
t223 = -t102 * t185 + t189 * t95;
t222 = pkin(1) * t297;
t220 = -t189 * t159 - t160 * t185;
t219 = t139 + t246;
t218 = -t141 + t255;
t73 = -pkin(9) * t142 + t261;
t216 = pkin(4) * t257 + t266 * pkin(9) + t261 * qJD(4) + qJD(5) * t73 + t307 * t185 + t296;
t72 = -pkin(9) * t143 + t220;
t215 = -t265 * pkin(9) + qJD(5) * t72 + t295;
t120 = t142 * t187;
t36 = -pkin(4) * t191 + pkin(9) * t120 + t223;
t38 = -pkin(9) * t119 + t285;
t212 = t184 * t36 + t188 * t38;
t68 = t188 * t119 - t120 * t184;
t69 = -t119 * t184 - t120 * t188;
t208 = qJD(1) * t182 - t169 * t191;
t175 = pkin(3) * t189 + pkin(4);
t205 = pkin(3) * t275 + t175 * t184;
t204 = -pkin(3) * t276 + t175 * t188;
t203 = -t102 * t250 + t185 * t49 + t189 * t52 + t95 * t249;
t114 = pkin(4) * t142 + t176;
t98 = pkin(4) * t119 + t151;
t60 = pkin(3) * t141 - pkin(4) * t209;
t54 = -t250 * t273 + (t293 * t271 + t238) * t189 + t294 * t185;
t53 = -qJD(2) * t202 - t293 * t119;
t39 = pkin(4) * t54 + t110;
t18 = -pkin(4) * t195 + t85;
t17 = t288 + t311;
t16 = t225 + t312;
t11 = t69 * qJD(5) + t184 * t53 + t188 * t54;
t10 = -t68 * qJD(5) - t184 * t54 + t188 * t53;
t7 = -pkin(9) * t54 + t203;
t6 = pkin(4) * t254 - pkin(9) * t53 - qJD(4) * t285 + t226;
t4 = -t15 * t184 + t284;
t1 = [0, 0, 0, 0.2e1 * t191 * t172, t258 * t297, t267, -t268, 0, -pkin(6) * t267 + t187 * t222, pkin(6) * t268 + t191 * t222, t106 * t271 + t294 * t141, (-t139 * t190 - t141 * t186) * t253 + (-t282 - t107 * t190 + (t139 * t186 - t141 * t190) * qJD(3)) * t187, t169 * t236 - t106 * t191 + (t141 * t187 + t190 * t208) * qJD(2), t187 * t233 + t107 * t191 + (-t139 * t187 - t186 * t208) * qJD(2), (-t169 - t256) * t254, -(-t152 * t252 + t262) * t169 + (t157 * t251 + pkin(6) * t107 + (qJD(1) * t138 + t96) * qJD(2)) * t187 + ((pkin(6) * t139 + t279) * qJD(2) + (t274 + (pkin(6) * t169 + t158) * t190) * qJD(3) + t264) * t191, (-pkin(6) * t235 + t263) * t169 + t206 * t191 + (pkin(6) * t106 - t157 * t252) * t187 + ((pkin(6) * t141 + t278) * t191 + (-pkin(6) * t277 - qJD(1) * t259 - t97) * t187) * qJD(2), -t120 * t32 - t209 * t53, -t119 * t32 - t120 * t195 + t209 * t54 - t53 * t88, -t162 * t53 - t191 * t32 + (-qJD(1) * t120 - t209) * t254, t162 * t54 - t191 * t195 + (-qJD(1) * t119 - t88) * t254, (-t162 - t256) * t254, -t226 * t162 - t227 * t191 + t110 * t88 - t151 * t195 + t85 * t119 + t109 * t54 + (t162 * t285 + t191 * t20) * qJD(4) + (qJD(1) * t223 + t19) * t254, t203 * t162 - t224 * t191 - t110 * t209 + t151 * t32 - t85 * t120 + t109 * t53 + (-t285 * qJD(1) - t20) * t254, -t10 * t211 + t69 * t8, t10 * t45 + t11 * t211 + t196 * t69 - t68 * t8, -t10 * t155 - t8 * t191 + (qJD(1) * t69 - t211) * t254, t11 * t155 - t191 * t196 + (-qJD(1) * t68 + t45) * t254, (-t155 - t256) * t254, -(-t184 * t7 + t188 * t6) * t155 - t241 * t191 - t39 * t45 - t98 * t196 + t18 * t68 + t56 * t11 + (t155 * t212 + t191 * t5) * qJD(5) + ((-t184 * t38 + t188 * t36) * qJD(1) + t4) * t254, t56 * t10 - t13 * t191 + t18 * t69 - t39 * t211 + t98 * t8 + ((-qJD(5) * t38 + t6) * t155 + t2 * t191) * t184 + ((qJD(5) * t36 + t7) * t155 + t229 * t191) * t188 + (-qJD(1) * t212 - t5) * t254; 0, 0, 0, -t187 * t269, t258 * t193, 0, 0, 0, t193 * pkin(1) * t187, pkin(1) * t269, -t141 * t277 + t282, (t106 + t281) * t190 + (-t107 + t280) * t186, -t233 + (t169 * t270 + t187 * t218) * qJD(1), t169 * t252 + (-t169 * t272 + t187 * t219) * qJD(1), t169 * t257, -pkin(2) * t107 + t260 * t169 + (pkin(7) * t277 + t279) * qJD(3) + ((-pkin(7) * t255 - t96) * t187 + (-pkin(6) * t219 - t279) * t191) * qJD(1), -pkin(2) * t106 - t126 * t169 + (-pkin(7) * t169 * t186 + t278) * qJD(3) + (-t157 * t270 + (-pkin(7) * t246 + t97) * t187 + (t169 * t271 + t191 * t218) * pkin(6)) * qJD(1), t32 * t143 - t209 * t266, -t142 * t32 + t143 * t195 + t209 * t265 - t266 * t88, -t266 * t162 + (qJD(2) * t143 + t209) * t257, t265 * t162 + (-qJD(2) * t142 + t88) * t257, t162 * t257, t85 * t142 - t176 * t195 + t214 * t88 + (t160 * t249 + (-qJD(4) * t159 + t307) * t185 + t296) * t162 + t265 * t109 + (qJD(2) * t220 - t19) * t257, t85 * t143 + t176 * t32 - t214 * t209 + t295 * t162 + t266 * t109 + (-qJD(2) * t261 + t20) * t257, -t211 * t290 + t8 * t93, t196 * t93 + t211 * t289 + t290 * t45 - t8 * t92, -t290 * t155 + (qJD(2) * t93 + t211) * t257, t289 * t155 + (-qJD(2) * t92 - t45) * t257, t155 * t257, -t114 * t196 + t18 * t92 + t289 * t56 - t287 * t45 + (t184 * t215 + t188 * t216) * t155 + ((-t184 * t73 + t188 * t72) * qJD(2) - t4) * t257, t114 * t8 + t18 * t93 + t290 * t56 - t287 * t211 + (-t184 * t216 + t188 * t215) * t155 + (-(t184 * t72 + t188 * t73) * qJD(2) + t5) * t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141 * t139, -t139 ^ 2 + t141 ^ 2, t106 - t281, -t107 - t280, t172, -t141 * t157 - t169 * t97 + t199, t139 * t157 - t169 * t96 - t194, -t309, t306, t303, t298, t172, t225 * t162 + (-t141 * t88 + t162 * t250 + t172 * t189) * pkin(3) + t300, -t288 * t162 + (t141 * t209 + t162 * t249 - t172 * t185) * pkin(3) + t302, t308, t305, t304, t299, t172, t204 * t172 + (t16 * t188 - t17 * t184) * t155 + t60 * t45 - (-t184 * t189 - t275) * t242 + (t155 * t205 - t5) * qJD(5) + t314, -t205 * t172 - t188 * t3 - (t16 * t184 + t17 * t188) * t155 + t60 * t211 + (t188 * t189 - t276) * t242 + (t155 * t204 - t284) * qJD(5) + t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t309, t306, t303, t298, t172, -t162 * t20 + t300, -t162 * t19 + t302, t308, t305, t304, t299, t172, (-t14 * t184 - t283) * t155 + (t155 * t248 + t172 * t188 - t209 * t45) * pkin(4) + t301, (t15 * t155 - t2) * t184 + (-t14 * t155 - t229) * t188 + (t155 * t247 - t172 * t184 - t209 * t211) * pkin(4) + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t308, t305, t304, t299, t172, -t155 * t5 + t301, -t155 * t4 - t188 * t229 + t315;];
tauc_reg = t1;
