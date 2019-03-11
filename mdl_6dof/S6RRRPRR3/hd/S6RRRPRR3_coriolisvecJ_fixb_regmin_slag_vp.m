% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:22
% EndTime: 2019-03-09 18:13:34
% DurationCPUTime: 4.12s
% Computational Cost: add. (5919->365), mult. (14352->463), div. (0->0), fcn. (10606->8), ass. (0->203)
t191 = cos(qJ(6));
t183 = qJD(2) + qJD(3);
t245 = qJD(5) - t183;
t246 = qJD(6) * t191;
t187 = sin(qJ(6));
t247 = qJD(6) * t187;
t193 = cos(qJ(2));
t290 = cos(qJ(3));
t238 = t290 * t193;
t219 = qJD(1) * t238;
t189 = sin(qJ(3));
t190 = sin(qJ(2));
t251 = qJD(1) * t190;
t237 = t189 * t251;
t134 = -t219 + t237;
t188 = sin(qJ(5));
t192 = cos(qJ(5));
t152 = t189 * t193 + t290 * t190;
t252 = qJD(1) * t152;
t294 = t134 * t188 + t192 * t252;
t248 = qJD(5) * t192;
t249 = qJD(5) * t188;
t261 = t189 * t190;
t217 = t183 * t261;
t254 = t183 * t219;
t93 = qJD(1) * t217 - t254;
t107 = t183 * t152;
t94 = t107 * qJD(1);
t31 = t134 * t248 + t188 * t94 - t192 * t93 - t249 * t252;
t15 = t191 * t31 + t245 * t246 - t247 * t294;
t215 = -t187 * t245 - t191 * t294;
t16 = -qJD(6) * t215 + t187 * t31;
t295 = -t192 * t134 + t188 * t252;
t302 = qJD(6) + t295;
t284 = t215 * t302;
t67 = t187 * t294 - t191 * t245;
t285 = t67 * t302;
t313 = (t15 - t285) * t191 + (-t16 + t284) * t187;
t129 = t252 * pkin(9);
t291 = -pkin(8) - pkin(7);
t164 = t291 * t190;
t158 = qJD(1) * t164;
t279 = qJD(2) * pkin(2);
t144 = t158 + t279;
t165 = t291 * t193;
t160 = qJD(1) * t165;
t262 = t189 * t160;
t97 = t290 * t144 + t262;
t236 = qJD(4) - t97;
t303 = -t129 + t236;
t277 = t15 * t187;
t308 = t191 * t302;
t9 = t215 * t308 - t277;
t32 = qJD(5) * t294 - t188 * t93 - t192 * t94;
t271 = t187 * t32;
t310 = t215 * t294 + t302 * t308 + t271;
t143 = t290 * t160;
t104 = t189 * t158 - t143;
t250 = qJD(3) * t189;
t218 = pkin(2) * t250 - t104;
t105 = t290 * t158 + t262;
t235 = qJD(3) * t290;
t255 = pkin(2) * t235 + qJD(4) - t105;
t307 = qJD(6) - t302;
t306 = pkin(5) * t294;
t177 = -t193 * pkin(2) - pkin(1);
t163 = t177 * qJD(1);
t74 = t134 * pkin(3) - qJ(4) * t252 + t163;
t57 = -pkin(4) * t134 - t74;
t20 = pkin(5) * t295 - pkin(10) * t294 + t57;
t194 = -pkin(3) - pkin(4);
t55 = t194 * t183 + t303;
t180 = t183 * qJ(4);
t287 = t134 * pkin(9);
t98 = t189 * t144 - t143;
t66 = t98 + t287;
t61 = t180 + t66;
t27 = t188 * t55 + t192 * t61;
t24 = pkin(10) * t245 + t27;
t12 = t187 * t20 + t191 * t24;
t288 = t12 * t294;
t283 = t302 * t294;
t305 = -t129 + t255;
t304 = t287 - t218;
t282 = t294 * t295;
t301 = t294 ^ 2 - t295 ^ 2;
t179 = t183 * qJD(4);
t239 = qJD(2) * t291;
t220 = qJD(1) * t239;
t145 = t190 * t220;
t146 = t193 * t220;
t221 = -t144 * t235 - t290 * t145 - t189 * t146 - t160 * t250;
t51 = t179 - t221;
t35 = pkin(9) * t94 + t51;
t52 = t144 * t250 + t189 * t145 - t290 * t146 - t160 * t235;
t37 = pkin(9) * t93 + t52;
t6 = t188 * t35 - t192 * t37 + t61 * t248 + t55 * t249;
t201 = t294 * t57 + t6;
t5 = t188 * t37 + t192 * t35 + t55 * t248 - t61 * t249;
t199 = -t295 * t57 + t5;
t216 = t187 * t24 - t191 * t20;
t240 = t6 * t191 - t216 * t294;
t300 = t245 * t294 - t32;
t299 = t245 * t295 + t31;
t269 = t191 * t32;
t298 = t294 * t67 + t269;
t244 = qJD(1) * qJD(2);
t296 = -0.2e1 * t244;
t109 = t189 * t164 - t290 * t165;
t151 = -t238 + t261;
t101 = t151 * t188 + t152 * t192;
t159 = t190 * t239;
t161 = t193 * t239;
t58 = t290 * t159 + t189 * t161 + t164 * t235 + t165 * t250;
t45 = pkin(9) * t107 + t58;
t106 = -qJD(2) * t238 - t193 * t235 + t217;
t59 = t109 * qJD(3) + t189 * t159 - t290 * t161;
t46 = t106 * pkin(9) + t59;
t108 = -t290 * t164 - t189 * t165;
t79 = -t152 * pkin(9) + t108;
t80 = pkin(9) * t151 + t109;
t47 = t188 * t80 - t192 * t79;
t13 = -qJD(5) * t47 + t188 * t46 + t192 * t45;
t211 = t192 * t151 - t152 * t188;
t26 = -t188 * t61 + t192 * t55;
t23 = -pkin(5) * t245 - t26;
t96 = t151 * pkin(3) - t152 * qJ(4) + t177;
t70 = -pkin(4) * t151 - t96;
t33 = -pkin(5) * t211 - pkin(10) * t101 + t70;
t42 = qJD(5) * t211 - t106 * t192 + t107 * t188;
t48 = t188 * t79 + t192 * t80;
t293 = (qJD(6) * t20 + t5) * t211 + t6 * t101 - (qJD(6) * t33 + t13) * t302 + t23 * t42 - t48 * t32;
t292 = t252 ^ 2;
t289 = pkin(4) * t252;
t286 = t33 * t32;
t176 = -t290 * pkin(2) - pkin(3);
t173 = -pkin(4) + t176;
t174 = pkin(2) * t189 + qJ(4);
t210 = t173 * t192 - t174 * t188;
t281 = -qJD(5) * t210 + t188 * t304 - t192 * t305;
t209 = t173 * t188 + t174 * t192;
t280 = qJD(5) * t209 + t188 * t305 + t192 * t304;
t278 = t101 * t23;
t276 = t183 * t58;
t275 = t183 * t59;
t274 = t183 * t97;
t273 = t183 * t98;
t270 = t187 * t302;
t268 = t191 * t215;
t213 = -qJ(4) * t188 + t192 * t194;
t267 = qJD(5) * t213 - t188 * t66 + t192 * t303;
t214 = qJ(4) * t192 + t188 * t194;
t266 = qJD(5) * t214 + t188 * t303 + t192 * t66;
t264 = t252 * t134;
t196 = qJD(1) ^ 2;
t259 = t193 * t196;
t195 = qJD(2) ^ 2;
t258 = t195 * t190;
t257 = t195 * t193;
t95 = pkin(3) * t252 + t134 * qJ(4);
t253 = t190 ^ 2 - t193 ^ 2;
t243 = pkin(2) * t251;
t242 = t190 * t279;
t234 = t190 * t244;
t233 = pkin(10) * t302 + t306;
t230 = t302 * t23;
t226 = t245 * t302;
t225 = pkin(1) * t296;
t111 = -pkin(10) + t209;
t64 = -t95 - t289;
t25 = -pkin(10) * t295 - t306 + t64;
t224 = qJD(6) * t111 - t243 + t25;
t157 = -pkin(10) + t214;
t223 = qJD(6) * t157 + t25;
t222 = t245 ^ 2;
t78 = t243 + t95;
t208 = -t246 * t302 - t271;
t207 = -t247 * t302 + t269;
t49 = t107 * pkin(3) + t106 * qJ(4) - t152 * qJD(4) + t242;
t206 = -t252 * t74 - t52;
t205 = -t134 * t74 - t221;
t204 = -t163 * t252 - t52;
t203 = t163 * t134 + t221;
t44 = pkin(2) * t234 + t94 * pkin(3) + t93 * qJ(4) - qJD(4) * t252;
t200 = -pkin(10) * t32 + (t23 + t26) * t302;
t34 = -pkin(4) * t107 - t49;
t21 = -pkin(4) * t94 - t44;
t198 = -t111 * t32 + t281 * t302 - t230;
t197 = -t157 * t32 - t267 * t302 - t230;
t156 = pkin(5) - t213;
t110 = pkin(5) - t210;
t82 = t180 + t98;
t77 = -pkin(3) * t183 + t236;
t71 = -t134 ^ 2 + t292;
t62 = t254 + (t134 - t237) * t183;
t60 = -t78 - t289;
t43 = qJD(5) * t101 - t106 * t188 - t192 * t107;
t14 = qJD(5) * t48 + t188 * t45 - t192 * t46;
t10 = pkin(5) * t43 - pkin(10) * t42 + t34;
t7 = t270 * t302 - t298;
t3 = pkin(5) * t32 - pkin(10) * t31 + t21;
t2 = t191 * t3;
t1 = [0, 0, 0, 0.2e1 * t193 * t234, t253 * t296, t257, -t258, 0, -pkin(7) * t257 + t190 * t225, pkin(7) * t258 + t193 * t225, -t106 * t252 - t152 * t93, t106 * t134 - t107 * t252 + t151 * t93 - t152 * t94, -t106 * t183, -t107 * t183, 0, t107 * t163 + t177 * t94 - t275 + (qJD(1) * t151 + t134) * t242, -t106 * t163 - t177 * t93 + 0.2e1 * t252 * t242 - t276, t107 * t74 + t134 * t49 + t151 * t44 + t94 * t96 - t275, -t106 * t77 - t107 * t82 - t108 * t93 - t109 * t94 - t134 * t58 - t151 * t51 + t152 * t52 + t252 * t59, t106 * t74 - t152 * t44 - t252 * t49 + t93 * t96 + t276, t108 * t52 + t109 * t51 + t44 * t96 + t49 * t74 + t58 * t82 + t59 * t77, t101 * t31 + t294 * t42, -t101 * t32 + t211 * t31 - t294 * t43 - t295 * t42, t42 * t245, -t43 * t245, 0, -t14 * t245 - t21 * t211 + t295 * t34 + t32 * t70 + t43 * t57, t101 * t21 - t13 * t245 + t294 * t34 + t31 * t70 + t42 * t57, -t42 * t268 + (t15 * t191 + t215 * t247) * t101 (t187 * t215 - t191 * t67) * t42 + (-t277 - t16 * t191 + (t187 * t67 + t268) * qJD(6)) * t101, t101 * t207 - t15 * t211 - t215 * t43 + t308 * t42, t101 * t208 + t16 * t211 - t42 * t270 - t43 * t67, -t211 * t32 + t302 * t43, -t2 * t211 - t216 * t43 + t14 * t67 + t47 * t16 + (t10 * t302 + t286 + (t211 * t24 - t302 * t48 + t278) * qJD(6)) * t191 + t293 * t187, -t12 * t43 - t14 * t215 + t47 * t15 + (-(-qJD(6) * t48 + t10) * t302 - t286 + (-qJD(6) * t24 + t3) * t211 - qJD(6) * t278) * t187 + t293 * t191; 0, 0, 0, -t190 * t259, t253 * t196, 0, 0, 0, t196 * pkin(1) * t190, pkin(1) * t259, t264, t71, t62, 0, 0, t104 * t183 + (-t134 * t251 - t183 * t250) * pkin(2) + t204, t105 * t183 + (-t183 * t235 - t251 * t252) * pkin(2) + t203, -t134 * t78 - t183 * t218 + t206, -t174 * t94 - t176 * t93 + (t218 + t82) * t252 + (t77 - t255) * t134, t255 * t183 + t252 * t78 + t179 + t205, t174 * t51 + t176 * t52 + t218 * t77 + t255 * t82 - t74 * t78, -t282, -t301, -t299, -t300, 0, -t245 * t280 - t295 * t60 + t201, t245 * t281 - t294 * t60 + t199, t9, -t313, -t310, t7, t283, t110 * t16 + t198 * t187 - t224 * t308 + t280 * t67 + t240, t110 * t15 - t288 - t280 * t215 + (t224 * t302 - t6) * t187 + t198 * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, t71, t62, 0, 0, t204 + t273, t203 + t274, -t134 * t95 + t206 + t273, pkin(3) * t93 - qJ(4) * t94 + (t82 - t98) * t252 + (t77 - t236) * t134, t252 * t95 + 0.2e1 * t179 + t205 - t274, -pkin(3) * t52 + qJ(4) * t51 + t236 * t82 - t74 * t95 - t77 * t98, -t282, -t301, -t299, -t300, 0, -t245 * t266 - t295 * t64 + t201, -t245 * t267 - t294 * t64 + t199, t9, -t313, -t310, t7, t283, t156 * t16 + t187 * t197 - t223 * t308 + t266 * t67 + t240, -t288 + t156 * t15 - t266 * t215 + (t223 * t302 - t6) * t187 + t197 * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, t62, -t183 ^ 2 - t292, -t183 * t82 - t206, 0, 0, 0, 0, 0, -t188 * t222 - t252 * t295, -t192 * t222 - t252 * t294, 0, 0, 0, 0, 0, -t252 * t308 + (-t187 * t226 - t16) * t192 + (t245 * t67 + t208) * t188, t252 * t270 + (-t191 * t226 - t15) * t192 + (-t215 * t245 - t207) * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, t301, t299, t300, 0, t245 * t27 - t201, t245 * t26 - t199, -t9, t313, t310, -t187 * t302 ^ 2 + t298, -t283, -pkin(5) * t16 + t187 * t200 - t233 * t308 - t27 * t67 - t240, -pkin(5) * t15 + t288 + t27 * t215 + (t233 * t302 + t6) * t187 + t200 * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215 * t67, t215 ^ 2 - t67 ^ 2, t15 + t285, -t16 - t284, t32, -t307 * t12 - t187 * t5 + t23 * t215 + t2, -t187 * t3 - t191 * t5 + t307 * t216 + t23 * t67;];
tauc_reg  = t1;
