% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:11
% EndTime: 2019-03-08 22:34:21
% DurationCPUTime: 3.31s
% Computational Cost: add. (2348->357), mult. (5712->516), div. (0->0), fcn. (4052->10), ass. (0->211)
t276 = pkin(4) + pkin(8);
t144 = qJD(3) * qJ(4);
t155 = cos(qJ(3));
t238 = qJD(2) * t155;
t152 = sin(qJ(2));
t147 = sin(pkin(6));
t244 = qJD(1) * t147;
t217 = t152 * t244;
t110 = qJD(2) * pkin(8) + t217;
t151 = sin(qJ(3));
t148 = cos(pkin(6));
t243 = qJD(1) * t148;
t70 = t155 * t110 + t151 * t243;
t55 = pkin(4) * t238 + t70;
t44 = t144 + t55;
t154 = cos(qJ(5));
t150 = sin(qJ(5));
t236 = qJD(3) * t150;
t98 = t154 * t238 + t236;
t26 = pkin(5) * t98 + t44;
t209 = t150 * t238;
t234 = qJD(3) * t154;
t100 = -t209 + t234;
t149 = sin(qJ(6));
t153 = cos(qJ(6));
t40 = t100 * t149 + t153 * t98;
t294 = t26 * t40;
t101 = t149 * t154 + t150 * t153;
t167 = t101 * t151;
t284 = qJD(5) + qJD(6);
t272 = -qJD(2) * t167 - t284 * t101;
t182 = t153 * t100 - t149 * t98;
t293 = t182 * t40;
t69 = t151 * t110 - t155 * t243;
t292 = qJD(4) + t69;
t156 = cos(qJ(2));
t241 = qJD(2) * t147;
t203 = qJD(1) * t241;
t190 = t156 * t203;
t237 = qJD(3) * t148;
t291 = qJD(1) * t237 + t190;
t290 = t182 ^ 2 - t40 ^ 2;
t240 = qJD(2) * t151;
t135 = qJD(5) + t240;
t127 = qJD(6) + t135;
t227 = qJD(6) * t153;
t228 = qJD(6) * t149;
t226 = qJD(2) * qJD(3);
t201 = t151 * t226;
t65 = -qJD(5) * t98 + t150 * t201;
t230 = qJD(5) * t154;
t66 = qJD(3) * t230 - qJD(5) * t209 - t154 * t201;
t9 = -t100 * t228 - t149 * t66 + t153 * t65 - t98 * t227;
t289 = t127 * t40 + t9;
t63 = -t144 - t70;
t233 = qJD(3) * t155;
t107 = t276 * t233;
t117 = t276 * t151;
t231 = qJD(5) * t150;
t253 = t150 * t156;
t235 = qJD(3) * t151;
t139 = pkin(3) * t235;
t186 = pkin(9) * t151 - qJ(4) * t155;
t232 = qJD(4) * t151;
t164 = t186 * qJD(3) - t232;
t72 = t139 + t164;
t157 = -pkin(3) - pkin(9);
t200 = -qJ(4) * t151 - pkin(2);
t94 = t157 * t155 + t200;
t287 = -t150 * t107 - t117 * t230 - t154 * t72 + t94 * t231 + (t151 * t253 + t152 * t154) * t244;
t62 = -qJD(3) * pkin(3) + t292;
t250 = t154 * t156;
t286 = t154 * t107 - t150 * t72 - (-t150 * t152 + t151 * t250) * t244;
t103 = t150 * t117;
t262 = t154 * t94 + t103;
t247 = pkin(4) * t240 + t292;
t137 = t155 * t226;
t38 = t157 * qJD(3) + t247;
t216 = t156 * t244;
t56 = qJD(2) * t94 - t216;
t18 = t150 * t38 + t154 * t56;
t35 = t110 * t233 + t291 * t151;
t27 = pkin(4) * t137 + t35;
t246 = pkin(3) * t201 + t152 * t203;
t36 = qJD(2) * t164 + t246;
t198 = -t150 * t36 + t154 * t27;
t162 = -t18 * qJD(5) + t198;
t2 = pkin(5) * t137 - pkin(10) * t65 + t162;
t225 = -t150 * t27 - t154 * t36 - t38 * t230;
t170 = -t56 * t231 - t225;
t3 = -pkin(10) * t66 + t170;
t218 = -t149 * t3 + t153 * t2;
t13 = -pkin(10) * t98 + t18;
t268 = t13 * t153;
t17 = -t150 * t56 + t154 * t38;
t12 = -pkin(10) * t100 + t17;
t7 = pkin(5) * t135 + t12;
t5 = t149 * t7 + t268;
t283 = -t5 * qJD(6) - t26 * t182 + t218;
t10 = t182 * qJD(6) + t149 * t65 + t153 * t66;
t282 = t127 * t182 - t10;
t11 = t13 * t228;
t206 = qJD(6) * t7 + t3;
t281 = t149 * t2 + t153 * t206 - t11;
t158 = qJD(3) ^ 2;
t168 = -qJ(4) * t233 - t232;
t47 = qJD(2) * t168 + t246;
t82 = t139 + t168;
t280 = qJD(2) * (-t82 + t217) - pkin(8) * t158 - t47;
t214 = t156 * t241;
t159 = qJD(2) ^ 2;
t256 = t147 * t159;
t221 = t152 * t256;
t257 = t147 * t152;
t222 = t151 * t257;
t50 = -qJD(3) * t222 + (t214 + t237) * t155;
t279 = qJD(3) * (t155 * t214 + t50) - t151 * t221;
t192 = t151 * t214;
t89 = t148 * t151 + t155 * t257;
t51 = qJD(3) * t89 + t192;
t278 = (t51 + t192) * qJD(3) + t155 * t221;
t277 = t284 * t155;
t254 = t150 * t151;
t175 = pkin(5) * t155 - pkin(10) * t254;
t207 = pkin(10) * t155 - t94;
t275 = t175 * qJD(3) + (t207 * t154 - t103) * qJD(5) + t286;
t229 = qJD(5) * t155;
t212 = t150 * t229;
t274 = -(t151 * t234 + t212) * pkin(10) + t287;
t273 = pkin(10) - t157;
t213 = t154 * t240;
t252 = t153 * t154;
t255 = t149 * t150;
t271 = -t149 * t231 - t150 * t228 + t153 * t213 - t240 * t255 + t252 * t284;
t140 = pkin(3) * t240;
t77 = t186 * qJD(2) + t140;
t270 = t150 * t55 + t154 * t77;
t269 = qJD(2) * pkin(2);
t267 = t135 * t98;
t266 = t154 * t65;
t143 = qJD(3) * qJD(4);
t224 = t110 * t235 - t291 * t155;
t29 = -t143 + t224;
t23 = -pkin(4) * t201 - t29;
t265 = t23 * t150;
t264 = t23 * t154;
t219 = -pkin(5) * t154 - pkin(4);
t263 = pkin(5) * t230 - t219 * t240 + t292;
t261 = t100 * t135;
t260 = t100 * t155;
t259 = t135 * t151;
t258 = t135 * t157;
t251 = t154 * t155;
t249 = t158 * t151;
t248 = t158 * t155;
t118 = t276 * t155;
t145 = t151 ^ 2;
t146 = t155 ^ 2;
t245 = t145 - t146;
t114 = -pkin(3) * t155 + t200;
t242 = qJD(2) * t114;
t239 = qJD(2) * t152;
t223 = t154 * t259;
t220 = t151 * t159 * t155;
t215 = t147 * t239;
t211 = t135 * t230;
t210 = t154 * t229;
t113 = t273 * t154;
t204 = t271 * t127;
t197 = -t150 * t77 + t154 * t55;
t181 = -t252 + t255;
t191 = t272 * t127 - t181 * t137;
t112 = t273 * t150;
t189 = t175 * qJD(2) - qJD(6) * t112 - t273 * t231 + t197;
t188 = pkin(10) * t213 + t284 * t113 + t270;
t104 = t154 * t117;
t28 = pkin(5) * t151 + t207 * t150 + t104;
t34 = -pkin(10) * t251 + t262;
t185 = t149 * t28 + t153 * t34;
t88 = -t148 * t155 + t222;
t172 = t147 * t250 - t150 * t88;
t52 = t147 * t253 + t154 * t88;
t184 = t149 * t172 + t153 * t52;
t183 = t149 * t52 - t153 * t172;
t180 = -qJD(2) * t146 + t259;
t178 = t135 * t150;
t177 = qJD(3) * t70 - t35;
t176 = qJD(3) * t69 - t224;
t171 = t151 * t44 + t157 * t233;
t111 = -t216 - t269;
t166 = qJD(3) * (t111 + t216 - t269);
t71 = -t216 + t242;
t165 = qJD(3) * (-t216 - t71 - t242);
t160 = t151 * t35 - t155 * t29 + (t151 * t63 + t155 * t62) * qJD(3);
t136 = pkin(5) * t150 + qJ(4);
t124 = t154 * t137;
t123 = t151 * t137;
t106 = t276 * t235;
t105 = -qJ(4) * t238 + t140;
t87 = pkin(5) * t251 + t118;
t79 = t101 * t155;
t78 = t181 * t155;
t58 = t71 * t240;
t57 = -pkin(5) * t212 + (-pkin(8) + t219) * t235;
t22 = t101 * t277 - t181 * t235;
t21 = qJD(3) * t167 + t181 * t277;
t16 = pkin(5) * t66 + t23;
t15 = qJD(5) * t52 + t150 * t51 + t154 * t215;
t14 = qJD(5) * t172 - t150 * t215 + t154 * t51;
t4 = -t13 * t149 + t153 * t7;
t1 = [0, 0, -t221, -t156 * t256, 0, 0, 0, 0, 0, -t278, -t279 (t151 * t51 + t155 * t50 + (-t151 * t89 + t155 * t88) * qJD(3)) * qJD(2), t278, t279, -t29 * t89 + t35 * t88 - t50 * t63 + t51 * t62 + (-t156 * t47 + t71 * t239) * t147, 0, 0, 0, 0, 0, t135 * t14 + t137 * t52 + t50 * t98 + t66 * t89, t100 * t50 - t135 * t15 + t137 * t172 + t65 * t89, 0, 0, 0, 0, 0 (-qJD(6) * t183 + t14 * t153 - t149 * t15) * t127 + t184 * t137 + t50 * t40 + t89 * t10 -(qJD(6) * t184 + t14 * t149 + t15 * t153) * t127 - t183 * t137 + t50 * t182 + t89 * t9; 0, 0, 0, 0, 0.2e1 * t123, -0.2e1 * t245 * t226, t248, -t249, 0, -pkin(8) * t248 + t151 * t166, pkin(8) * t249 + t155 * t166 (-t145 - t146) * t190 + t160, t151 * t165 - t155 * t280, t151 * t280 + t155 * t165, t114 * t47 + t71 * t82 + (-t152 * t71 + (-t151 * t62 + t155 * t63) * t156) * t244 + t160 * pkin(8), -t150 * t155 * t65 + (t150 * t235 - t210) * t100 (t100 * t154 - t150 * t98) * t235 + (t150 * t66 - t266 + (t100 * t150 + t154 * t98) * qJD(5)) * t155, -t135 * t210 + t151 * t65 + (t150 * t180 + t260) * qJD(3), t135 * t212 - t151 * t66 + (t154 * t180 - t155 * t98) * qJD(3), t135 * t233 + t123, -t106 * t98 + t118 * t66 + (-t44 * t234 + t198) * t151 + t286 * t135 + (-t135 * t262 - t18 * t151) * qJD(5) + (-t98 * t216 - t44 * t231 + t264 + ((-t150 * t94 + t104) * qJD(2) + t17) * qJD(3)) * t155, -t106 * t100 + t118 * t65 + ((qJD(3) * t44 + qJD(5) * t56) * t150 + t225) * t151 + t287 * t135 + (-t100 * t216 - t44 * t230 - t265 + (-t262 * qJD(2) - t18) * qJD(3)) * t155, t182 * t21 - t79 * t9, t10 * t79 + t182 * t22 - t21 * t40 + t78 * t9, t127 * t21 + t151 * t9 + (-qJD(2) * t79 + t182) * t233, -t10 * t151 + t127 * t22 + (qJD(2) * t78 - t40) * t233, t127 * t233 + t123, t218 * t151 + t57 * t40 + t87 * t10 - t16 * t78 - t26 * t22 + (t274 * t149 + t275 * t153) * t127 + (-t127 * t185 - t151 * t5) * qJD(6) + (-t40 * t216 + ((-t149 * t34 + t153 * t28) * qJD(2) + t4) * qJD(3)) * t155, -t281 * t151 + t57 * t182 + t87 * t9 - t16 * t79 + t26 * t21 + ((-qJD(6) * t28 + t274) * t153 + (qJD(6) * t34 - t275) * t149) * t127 + (-t182 * t216 + (-qJD(2) * t185 - t5) * qJD(3)) * t155; 0, 0, 0, 0, -t220, t245 * t159, 0, 0, 0, -t111 * t240 + t177, -t111 * t238 - t176, 0, -t105 * t238 - t177 + t58, 0.2e1 * t143 + (t105 * t151 + t155 * t71) * qJD(2) + t176, -pkin(3) * t35 - qJ(4) * t29 - t105 * t71 - t292 * t63 - t62 * t70, -t100 * t178 + t266 (-t66 - t261) * t154 + (-t65 + t267) * t150, -t135 * t231 + t124 + (-t135 * t254 - t260) * qJD(2), -t211 + (-t223 + (t98 - t236) * t155) * qJD(2), -t135 * t238, qJ(4) * t66 + t265 - t197 * t135 + t247 * t98 + (-t150 * t258 + t44 * t154) * qJD(5) + (t154 * t171 - t17 * t155) * qJD(2), qJ(4) * t65 + t264 + t270 * t135 + t247 * t100 + (-t44 * t150 - t154 * t258) * qJD(5) + (-t150 * t171 + t18 * t155) * qJD(2), -t181 * t9 + t182 * t272, t10 * t181 - t101 * t9 - t182 * t271 - t272 * t40, -t182 * t238 + t191, -t204 + (-qJD(3) * t101 + t40) * t238, -t127 * t238, t136 * t10 + t16 * t101 + t263 * t40 + t271 * t26 + (t149 * t188 - t153 * t189) * t127 + ((t112 * t149 - t113 * t153) * qJD(3) - t4) * t238, -t16 * t181 + t136 * t9 + t263 * t182 + t272 * t26 + (t149 * t189 + t153 * t188) * t127 + (-(-t112 * t153 - t113 * t149) * qJD(3) + t5) * t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, -t145 * t159 - t158, qJD(3) * t63 + t35 + t58, 0, 0, 0, 0, 0, -qJD(3) * t98 - t135 * t178 + t124, -t211 - qJD(3) * t100 + (-t150 * t233 - t223) * qJD(2), 0, 0, 0, 0, 0, -qJD(3) * t40 + t191, -t204 + (-t101 * t238 - t182) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100 * t98, t100 ^ 2 - t98 ^ 2, t65 + t267, -t66 + t261, t137, -t100 * t44 + t135 * t18 + t162, t135 * t17 + t44 * t98 - t170, t293, t290, t289, t282, t137 -(-t12 * t149 - t268) * t127 + (-t100 * t40 - t127 * t228 + t137 * t153) * pkin(5) + t283, t294 + t11 + (-t127 * t13 - t2) * t149 + (t12 * t127 - t206) * t153 + (-t100 * t182 - t127 * t227 - t137 * t149) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, t290, t289, t282, t137, t127 * t5 + t283, t127 * t4 - t281 + t294;];
tauc_reg  = t1;
