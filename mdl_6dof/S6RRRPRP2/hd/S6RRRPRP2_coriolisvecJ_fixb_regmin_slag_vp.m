% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:37:18
% EndTime: 2019-03-09 16:37:30
% DurationCPUTime: 4.29s
% Computational Cost: add. (10606->414), mult. (27279->537), div. (0->0), fcn. (20273->8), ass. (0->231)
t184 = cos(qJ(5));
t251 = qJD(5) * t184;
t182 = sin(qJ(3));
t183 = sin(qJ(2));
t185 = cos(qJ(2));
t305 = cos(qJ(3));
t206 = -t182 * t183 + t305 * t185;
t139 = t206 * qJD(1);
t154 = t182 * t185 + t305 * t183;
t141 = t154 * qJD(1);
t180 = sin(pkin(10));
t273 = cos(pkin(10));
t235 = t273 * t139 - t180 * t141;
t268 = t235 * t184;
t322 = t251 - t268;
t181 = sin(qJ(5));
t252 = qJD(5) * t181;
t321 = -t235 * t181 + t252;
t309 = qJD(5) - t235;
t320 = t309 ^ 2;
t204 = t180 * t139 + t273 * t141;
t248 = qJD(2) + qJD(3);
t100 = t181 * t248 + t184 * t204;
t233 = t100 * t309;
t194 = t141 * t248;
t232 = t309 * t181;
t121 = t248 * t206;
t240 = t273 * t121;
t244 = qJD(3) * t305;
t253 = qJD(3) * t182;
t310 = qJD(1) * (t240 + t180 * (-t154 * qJD(2) - t183 * t244 - t185 * t253));
t318 = qJD(5) * t248 + t310;
t307 = pkin(7) + pkin(8);
t162 = t307 * t185;
t157 = qJD(1) * t162;
t146 = t305 * t157;
t161 = t307 * t183;
t155 = qJD(1) * t161;
t294 = qJD(2) * pkin(2);
t149 = -t155 + t294;
t208 = -t182 * t149 - t146;
t266 = t139 * qJ(4);
t97 = -t208 + t266;
t243 = t273 * t97;
t135 = t141 * qJ(4);
t142 = t182 * t157;
t234 = t305 * t149 - t142;
t96 = -t135 + t234;
t88 = pkin(3) * t248 + t96;
t56 = t180 * t88 + t243;
t54 = pkin(9) * t248 + t56;
t174 = -t185 * pkin(2) - pkin(1);
t160 = qJD(1) * t174;
t123 = -t139 * pkin(3) + qJD(4) + t160;
t63 = -pkin(4) * t235 - pkin(9) * t204 + t123;
t28 = t181 * t63 + t184 * t54;
t16 = qJ(6) * t309 + t28;
t193 = t121 * qJD(1);
t247 = qJD(2) * t307;
t229 = qJD(1) * t247;
t151 = t185 * t229;
t130 = t305 * t151;
t213 = -t157 * t244 - t130;
t150 = t183 * t229;
t262 = t182 * t150;
t311 = -qJ(4) * t193 - t141 * qJD(4) - t149 * t253 + t213 + t262;
t198 = t149 * t244 - t305 * t150 - t182 * t151 - t157 * t253;
t49 = -qJ(4) * t194 + t139 * qJD(4) + t198;
t22 = t311 * t180 + t273 * t49;
t249 = qJD(1) * qJD(2);
t242 = t183 * t249;
t168 = pkin(2) * t242;
t191 = pkin(3) * t194 + t168;
t84 = t180 * t193 + t273 * t194;
t39 = t84 * pkin(4) - pkin(9) * t310 + t191;
t238 = t181 * t22 - t184 * t39 + t54 * t251 + t63 * t252;
t306 = pkin(5) * t84;
t4 = t238 - t306;
t317 = t16 * t309 - t4;
t316 = -0.2e1 * t249;
t210 = t181 * t39 + t184 * t22 + t63 * t251 - t54 * t252;
t278 = t84 * qJ(6);
t2 = qJD(6) * t309 + t210 + t278;
t315 = t4 * t181 + t2 * t184;
t102 = t182 * t155 - t146 - t266;
t256 = -t305 * t155 - t142;
t103 = -t135 + t256;
t239 = t273 * t182;
t295 = pkin(2) * qJD(3);
t276 = t273 * t102 - t180 * t103 + (t305 * t180 + t239) * t295;
t263 = t180 * t182;
t132 = (t273 * t305 - t263) * t295;
t69 = t180 * t102 + t273 * t103;
t275 = t132 - t69;
t91 = t180 * t97;
t55 = t273 * t88 - t91;
t53 = -pkin(4) * t248 - t55;
t98 = t181 * t204 - t184 * t248;
t31 = t98 * pkin(5) - t100 * qJ(6) + t53;
t21 = t180 * t49 - t273 * t311;
t46 = -t318 * t184 + t204 * t252;
t47 = t318 * t181 + t204 * t251;
t7 = t47 * pkin(5) + t46 * qJ(6) - t100 * qJD(6) + t21;
t314 = -t7 * t184 + t31 * t252;
t313 = t21 * t184 - t53 * t252;
t312 = t321 * pkin(5) - t322 * qJ(6) - t181 * qJD(6);
t308 = t100 ^ 2;
t304 = t204 * pkin(5);
t303 = t141 * pkin(3);
t302 = t180 * pkin(3);
t301 = t7 * t181;
t299 = -t181 * t47 - t98 * t251;
t61 = t273 * t96 - t91;
t73 = pkin(4) * t204 - pkin(9) * t235 + t303;
t298 = t181 * t73 + t184 * t61;
t254 = qJD(1) * t183;
t175 = pkin(2) * t254;
t70 = t175 + t73;
t297 = t181 * t70 + t184 * t69;
t119 = t180 * t154 - t206 * t273;
t120 = t273 * t154 + t180 * t206;
t219 = -pkin(3) * t206 + t174;
t77 = t119 * pkin(4) - t120 * pkin(9) + t219;
t110 = -t154 * qJ(4) - t305 * t161 - t182 * t162;
t207 = t182 * t161 - t305 * t162;
t111 = qJ(4) * t206 - t207;
t79 = t180 * t110 + t273 * t111;
t296 = t181 * t77 + t184 * t79;
t293 = t100 * t31;
t292 = t100 * t98;
t291 = t204 * t98;
t290 = t235 * t31;
t289 = t235 * t98;
t173 = t305 * pkin(2) + pkin(3);
t134 = pkin(2) * t239 + t180 * t173;
t128 = pkin(9) + t134;
t288 = t128 * t84;
t170 = pkin(9) + t302;
t287 = t170 * t84;
t286 = t181 * t61;
t80 = t181 * t84;
t122 = t248 * t154;
t86 = -t180 * t122 + t240;
t285 = t181 * t86;
t284 = t181 * t98;
t283 = t184 * t47;
t81 = t184 * t84;
t282 = t184 * t86;
t280 = t46 * t181;
t279 = t53 * t235;
t277 = t312 + t276;
t60 = t180 * t96 + t243;
t274 = -t60 + t312;
t272 = t100 * t204;
t271 = t100 * t184;
t270 = t309 * t204;
t267 = t120 * t184;
t265 = t141 * t139;
t264 = t160 * t141;
t187 = qJD(1) ^ 2;
t260 = t185 * t187;
t186 = qJD(2) ^ 2;
t259 = t186 * t183;
t258 = t186 * t185;
t27 = -t181 * t54 + t184 * t63;
t257 = qJD(6) - t27;
t255 = t183 ^ 2 - t185 ^ 2;
t177 = t183 * t294;
t245 = t273 * pkin(3);
t241 = t122 * pkin(3) + t177;
t156 = t183 * t247;
t158 = t185 * t247;
t203 = -t305 * t156 - t182 * t158 - t161 * t244 - t162 * t253;
t66 = -t122 * qJ(4) + qJD(4) * t206 + t203;
t200 = qJD(3) * t207 + t182 * t156 - t305 * t158;
t67 = -t121 * qJ(4) - t154 * qJD(4) + t200;
t35 = t180 * t66 - t273 * t67;
t26 = t181 * t69 - t184 * t70 - t304;
t237 = t181 * t132 - t26;
t236 = pkin(1) * t316;
t78 = -t273 * t110 + t180 * t111;
t230 = t21 * t181 + t204 * t28 + t53 * t251;
t171 = -t245 - pkin(4);
t227 = t21 * t120 - t79 * t84;
t226 = t184 * pkin(5) + t181 * qJ(6);
t225 = pkin(5) * t181 - qJ(6) * t184;
t224 = t204 * t56 + t235 * t55;
t223 = -t279 - t288;
t222 = -t279 - t287;
t15 = -pkin(5) * t309 + t257;
t221 = t15 * t184 - t16 * t181;
t220 = t15 * t181 + t16 * t184;
t133 = -pkin(2) * t263 + t273 * t173;
t218 = -t309 * t321 + t81;
t217 = t322 * t309 + t80;
t216 = t15 * t204 + t314;
t215 = -t16 * t204 + t31 * t268 - t301;
t127 = -pkin(4) - t133;
t214 = -t204 * t27 - t313;
t212 = t28 * t309 - t238;
t211 = -t120 * t252 + t282;
t36 = t180 * t67 + t273 * t66;
t85 = t180 * t121 + t273 * t122;
t43 = t85 * pkin(4) - t86 * pkin(9) + t241;
t209 = t181 * t43 + t184 * t36 + t77 * t251 - t79 * t252;
t205 = -t128 * t252 + t184 * t132;
t202 = t15 * t251 - t16 * t321 + t315;
t201 = qJD(5) * t221 + t315;
t199 = -t280 - t283 + (t271 + t284) * qJD(5);
t197 = -t160 * t139 - t198;
t148 = -t226 + t171;
t124 = t127 - t226;
t108 = t204 * qJ(6);
t104 = -t139 ^ 2 + t141 ^ 2;
t89 = -t139 * t248 + t193;
t71 = pkin(5) * t100 + qJ(6) * t98;
t44 = t120 * t225 + t78;
t34 = -t119 * pkin(5) + t181 * t79 - t184 * t77;
t33 = qJ(6) * t119 + t296;
t30 = t309 * t98 - t46;
t25 = t108 + t297;
t24 = -t184 * t73 + t286 - t304;
t23 = t108 + t298;
t13 = t184 * t233 - t280;
t12 = t217 - t272;
t11 = t218 + t291;
t9 = t225 * t86 + (qJD(5) * t226 - qJD(6) * t184) * t120 + t35;
t8 = (-t46 + t289) * t184 - t100 * t232 + t299;
t6 = -t85 * pkin(5) + t296 * qJD(5) + t181 * t36 - t184 * t43;
t5 = t85 * qJ(6) + t119 * qJD(6) + t209;
t1 = [0, 0, 0, 0.2e1 * t185 * t242, t255 * t316, t258, -t259, 0, -pkin(7) * t258 + t183 * t236, pkin(7) * t259 + t185 * t236, t141 * t121 + t154 * t193, t121 * t139 - t141 * t122 + t248 * qJD(1) * (-t154 ^ 2 + t206 ^ 2) t121 * t248, -t122 * t248, 0, 0.2e1 * t160 * t122 - t139 * t177 - t168 * t206 + t200 * t248, t160 * t121 + t141 * t177 + t154 * t168 + t174 * t193 - t203 * t248, -t22 * t119 + t204 * t35 + t235 * t36 + t310 * t78 - t55 * t86 - t56 * t85 + t227, t123 * t241 + t191 * t219 + t21 * t78 + t22 * t79 - t55 * t35 + t56 * t36, t100 * t211 - t46 * t267 (-t100 * t181 - t184 * t98) * t86 + (t280 - t283 + (-t271 + t284) * qJD(5)) * t120, t100 * t85 - t46 * t119 + t211 * t309 + t84 * t267, -t120 * t80 - t47 * t119 - t98 * t85 + (-t120 * t251 - t285) * t309, t119 * t84 + t309 * t85, -t238 * t119 + t27 * t85 + t35 * t98 + t78 * t47 + ((-qJD(5) * t79 + t43) * t309 + t77 * t84 + t53 * qJD(5) * t120) * t184 + ((-qJD(5) * t77 - t36) * t309 + t53 * t86 + t227) * t181, t35 * t100 - t210 * t119 + t313 * t120 - t209 * t309 - t28 * t85 + t53 * t282 - t296 * t84 - t78 * t46, t31 * t285 - t6 * t309 - t4 * t119 - t15 * t85 - t34 * t84 + t44 * t47 + t9 * t98 + (t31 * t251 + t301) * t120, t6 * t100 - t33 * t47 - t34 * t46 - t5 * t98 + t221 * t86 + (-qJD(5) * t220 - t2 * t181 + t4 * t184) * t120, -t9 * t100 + t2 * t119 + t314 * t120 + t16 * t85 - t31 * t282 + t309 * t5 + t33 * t84 + t44 * t46, t15 * t6 + t16 * t5 + t2 * t33 + t31 * t9 + t34 * t4 + t44 * t7; 0, 0, 0, -t183 * t260, t255 * t187, 0, 0, 0, t187 * pkin(1) * t183, pkin(1) * t260, -t265, t104, t89, 0, 0, t139 * t175 - t264 + (-qJD(3) * t149 + t150) * t182 + t213 + (t146 + (-t155 - t295) * t182) * t248, t256 * t248 + (-t141 * t254 - t244 * t248) * pkin(2) + t197, -t133 * t310 - t134 * t84 + t276 * t204 + t275 * t235 + t224, t22 * t134 - t21 * t133 - t123 * (t175 + t303) + t275 * t56 - t276 * t55, t13, t8, t12, t11, -t270, t127 * t47 + t276 * t98 + t223 * t181 + ((-qJD(5) * t128 - t70) * t184 - t275 * t181) * t309 + t214, -t127 * t46 + t223 * t184 + t276 * t100 + (-t205 + t297) * t309 + t230, t124 * t47 + t277 * t98 + (-t288 - t290) * t181 + (-t128 * t251 - t237) * t309 + t216, t25 * t98 + (-t132 * t98 - t15 * t235) * t184 + t237 * t100 + t199 * t128 + t202, t124 * t46 + (-qJD(5) * t31 + t288) * t184 - t277 * t100 + (t205 - t25) * t309 + t215, t7 * t124 + t128 * t201 + t132 * t220 - t15 * t26 - t16 * t25 + t277 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t265, t104, t89, 0, 0, -qJD(2) * t208 - t130 + t262 - t264, t234 * t248 + t197, -t204 * t60 - t235 * t61 - t245 * t310 - t84 * t302 + t224, t55 * t60 - t56 * t61 + (-t123 * t141 + t180 * t22 - t273 * t21) * pkin(3), t13, t8, t12, t11, -t270, t171 * t47 - t60 * t98 + t222 * t181 + (t286 + (-qJD(5) * t170 - t73) * t184) * t309 + t214, -t60 * t100 - t171 * t46 + t222 * t184 + (t170 * t252 + t298) * t309 + t230, t148 * t47 + t274 * t98 + (-t287 - t290) * t181 + (-t170 * t251 + t24) * t309 + t216, -t24 * t100 - t15 * t268 + t170 * t199 + t23 * t98 + t202, t170 * t81 - t23 * t309 + t148 * t46 - t274 * t100 + (-t170 * t232 - t184 * t31) * qJD(5) + t215, t7 * t148 - t15 * t24 - t16 * t23 + t170 * t201 + t274 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204 ^ 2 - t235 ^ 2, t204 * t55 - t235 * t56 + t191, 0, 0, 0, 0, 0, t218 - t291, -t184 * t320 - t272 - t80, -t232 * t309 - t291 + t81 (t46 + t289) * t184 + t181 * t233 + t299, t217 + t272, -t31 * t204 + t317 * t184 + (t15 * t309 + t2) * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292, -t98 ^ 2 + t308, t30, -t47 + t233, t84, -t100 * t53 + t212, t27 * t309 + t53 * t98 - t210, -t71 * t98 + t212 - t293 + 0.2e1 * t306, pkin(5) * t46 - qJ(6) * t47 + (t16 - t28) * t100 + (t15 - t257) * t98, 0.2e1 * t278 + t71 * t100 - t31 * t98 + (0.2e1 * qJD(6) - t27) * t309 + t210, -pkin(5) * t4 + qJ(6) * t2 - t15 * t28 + t16 * t257 - t31 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84 + t292, t30, -t308 - t320, t293 - t317;];
tauc_reg  = t1;
