% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:45:56
% EndTime: 2019-03-09 16:46:08
% DurationCPUTime: 4.05s
% Computational Cost: add. (6222->414), mult. (14759->508), div. (0->0), fcn. (10162->6), ass. (0->215)
t172 = sin(qJ(3));
t173 = sin(qJ(2));
t175 = cos(qJ(2));
t278 = cos(qJ(3));
t136 = t172 * t175 + t278 * t173;
t235 = qJD(1) * t136;
t292 = qJD(5) + t235;
t223 = t278 * t175;
t206 = qJD(1) * t223;
t234 = qJD(1) * t173;
t222 = t172 * t234;
t122 = -t206 + t222;
t167 = qJD(2) + qJD(3);
t171 = sin(qJ(5));
t174 = cos(qJ(5));
t99 = -t174 * t122 + t167 * t171;
t288 = t292 * t99;
t98 = t167 * t136;
t180 = qJD(1) * t98;
t231 = qJD(5) * t174;
t232 = qJD(5) * t171;
t51 = -t122 * t231 + t167 * t232 - t171 * t180;
t295 = t51 - t288;
t212 = t292 ^ 2;
t244 = t172 * t173;
t203 = t167 * t244;
t237 = t167 * t206;
t90 = qJD(1) * t203 - t237;
t88 = t174 * t90;
t294 = -t171 * t212 - t88;
t218 = qJD(3) * t278;
t280 = -pkin(8) - pkin(7);
t147 = t280 * t175;
t142 = qJD(1) * t147;
t126 = t172 * t142;
t146 = t280 * t173;
t140 = qJD(1) * t146;
t96 = t278 * t140 + t126;
t293 = pkin(2) * t218 - t96;
t179 = t167 * t235;
t228 = qJD(1) * qJD(2);
t291 = -0.2e1 * t228;
t101 = t122 * t171 + t167 * t174;
t165 = t167 * qJ(4);
t275 = t122 * pkin(4);
t129 = t278 * t142;
t268 = qJD(2) * pkin(2);
t131 = t140 + t268;
t94 = t172 * t131 - t129;
t73 = t94 - t275;
t64 = t165 + t73;
t35 = pkin(5) * t99 - qJ(6) * t101 + t64;
t224 = qJD(2) * t280;
t207 = qJD(1) * t224;
t132 = t173 * t207;
t133 = t175 * t207;
t233 = qJD(3) * t172;
t208 = -t131 * t218 - t278 * t132 - t172 * t133 - t142 * t233;
t47 = -t167 * qJD(4) + t208;
t32 = -pkin(4) * t180 - t47;
t89 = t174 * t180;
t52 = qJD(5) * t101 - t89;
t8 = t52 * pkin(5) + t51 * qJ(6) - t101 * qJD(6) + t32;
t290 = -t8 * t171 - t35 * t231;
t289 = t32 * t171 + t64 * t231;
t251 = qJD(4) + t293;
t287 = -t174 * t8 + t35 * t232;
t286 = -t167 * t99 + t294;
t227 = pkin(2) * t233;
t95 = t172 * t140 - t129;
t205 = -t95 + t227;
t210 = t292 * t101;
t93 = -t278 * t131 - t126;
t239 = qJD(4) + t93;
t196 = -pkin(5) * t174 - qJ(6) * t171 - pkin(4);
t285 = pkin(5) * t231 + qJ(6) * t232 - t174 * qJD(6) - t196 * t235 + qJD(4);
t104 = -t172 * t146 + t278 * t147;
t284 = t171 * pkin(5) - qJ(6) * t174;
t283 = t101 ^ 2;
t282 = t235 ^ 2;
t281 = pkin(3) + pkin(9);
t279 = pkin(5) * t90;
t277 = pkin(2) * t172;
t276 = pkin(5) * t122;
t274 = t235 * pkin(4);
t217 = t173 * t228;
t154 = pkin(2) * t217;
t40 = pkin(3) * t179 + t90 * qJ(4) - qJD(4) * t235 + t154;
t17 = pkin(9) * t180 + t40;
t50 = t131 * t233 + t172 * t132 - t278 * t133 - t142 * t218;
t37 = -pkin(4) * t90 + t50;
t240 = t274 + t239;
t56 = -t281 * t167 + t240;
t160 = -pkin(2) * t175 - pkin(1);
t145 = t160 * qJD(1);
t182 = -qJ(4) * t235 + t145;
t59 = t281 * t122 + t182;
t215 = t171 * t17 - t174 * t37 + t59 * t231 + t56 * t232;
t3 = t215 + t279;
t2 = t3 * t174;
t118 = t235 * pkin(9);
t91 = pkin(3) * t235 + qJ(4) * t122;
t83 = pkin(2) * t234 + t91;
t62 = t118 + t83;
t77 = t95 - t275;
t272 = t171 * t77 + t174 * t62;
t69 = t118 + t91;
t271 = t171 * t73 + t174 * t69;
t135 = -t223 + t244;
t197 = -qJ(4) * t136 + t160;
t75 = t281 * t135 + t197;
t103 = -t278 * t146 - t172 * t147;
t84 = t136 * pkin(4) + t103;
t270 = t171 * t84 + t174 * t75;
t269 = qJ(6) * t90;
t267 = t101 * t35;
t266 = t101 * t99;
t21 = t171 * t56 + t174 * t59;
t11 = qJ(6) * t292 + t21;
t265 = t11 * t235;
t141 = t173 * t224;
t143 = t175 * t224;
t60 = -t278 * t141 - t172 * t143 - t146 * t218 - t147 * t233;
t264 = t167 * t60;
t263 = t167 * t94;
t262 = t171 * t90;
t261 = t171 * t98;
t260 = t174 * t51;
t259 = t174 * t73;
t258 = t174 * t98;
t257 = t281 * t90;
t61 = -t104 * qJD(3) + t172 * t141 - t278 * t143;
t256 = t61 * t167;
t255 = t90 * t135;
t254 = t285 + t293;
t253 = t93 + t285;
t252 = t274 + t251;
t249 = t101 * t167;
t248 = t292 * t122;
t247 = t235 * t122;
t246 = t235 * t174;
t245 = t135 * t171;
t178 = qJD(1) ^ 2;
t243 = t175 * t178;
t177 = qJD(2) ^ 2;
t242 = t177 * t173;
t241 = t177 * t175;
t20 = -t171 * t59 + t174 * t56;
t238 = qJD(6) - t20;
t236 = t173 ^ 2 - t175 ^ 2;
t230 = qJD(5) * t281;
t163 = t173 * t268;
t159 = -t278 * pkin(2) - pkin(3);
t155 = -pkin(9) + t159;
t226 = t155 * t88;
t225 = t174 * t257;
t221 = t155 * t231;
t220 = t174 * t230;
t192 = t174 * t17 + t171 * t37 + t56 * t231 - t59 * t232;
t1 = qJD(6) * t292 + t192 - t269;
t10 = -pkin(5) * t292 + t238;
t219 = -t10 * t235 - t1;
t216 = -t21 * t122 + t32 * t174;
t214 = t292 * t64;
t213 = pkin(1) * t291;
t211 = t174 * t292;
t144 = qJ(4) + t284;
t202 = t10 * t174 - t11 * t171;
t201 = t10 * t171 + t11 * t174;
t82 = -pkin(3) * t167 + t239;
t86 = -t165 - t94;
t200 = t82 * t122 - t235 * t86;
t199 = -t10 * t122 + t35 * t246 - t290;
t198 = t20 * t122 + t64 * t246 + t289;
t195 = t21 * t292 - t215;
t194 = t135 * t231 + t261;
t97 = -qJD(2) * t223 - t175 * t218 + t203;
t193 = qJ(4) * t97 - qJD(4) * t136 + t163;
t31 = t281 * t98 + t193;
t44 = -t97 * pkin(4) + t61;
t191 = t171 * t44 + t174 * t31 + t84 * t231 - t75 * t232;
t190 = t171 * t235 * t35 + t11 * t122 + t287;
t80 = pkin(3) * t122 + t182;
t189 = t235 * t80 + t50;
t188 = -t145 * t235 - t50;
t187 = t145 * t122 + t208;
t185 = -t122 * t80 - t47;
t184 = -t174 * t212 + t262;
t65 = t237 + (t122 - t222) * t167;
t181 = qJD(5) * t201 + t1 * t171 - t2;
t156 = qJ(4) + t277;
t130 = t144 + t277;
t117 = t122 * qJ(6);
t92 = pkin(3) * t135 + t197;
t85 = -pkin(4) * t135 - t104;
t76 = -t122 ^ 2 + t282;
t74 = t90 * t136;
t63 = pkin(5) * t101 + qJ(6) * t99;
t49 = t135 * t196 - t104;
t46 = pkin(3) * t98 + t193;
t43 = -pkin(4) * t98 - t60;
t39 = -pkin(5) * t136 + t171 * t75 - t174 * t84;
t38 = qJ(6) * t136 + t270;
t25 = t171 * t69 - t259 + t276;
t24 = -t117 + t271;
t23 = t171 * t62 - t174 * t77 + t276;
t22 = -t117 + t272;
t16 = -t122 * t99 + t184;
t15 = t101 * t122 + t294;
t12 = -t171 * t210 - t260;
t9 = t196 * t98 + (t284 * qJD(5) - qJD(6) * t171) * t135 - t60;
t6 = (-t52 - t210) * t174 + (t51 + t288) * t171;
t5 = pkin(5) * t97 + t270 * qJD(5) + t171 * t31 - t174 * t44;
t4 = -qJ(6) * t97 + qJD(6) * t136 + t191;
t7 = [0, 0, 0, 0.2e1 * t175 * t217, t236 * t291, t241, -t242, 0, -pkin(7) * t241 + t173 * t213, pkin(7) * t242 + t175 * t213, -t235 * t97 - t74, t97 * t122 - t136 * t179 - t235 * t98 + t255, -t97 * t167, -t98 * t167, 0, t122 * t163 + t135 * t154 + t145 * t98 + t160 * t180 - t256, -t145 * t97 - t160 * t90 + 0.2e1 * t235 * t163 + t264, -t103 * t90 + t104 * t180 + t60 * t122 + t47 * t135 + t50 * t136 + t235 * t61 - t82 * t97 + t86 * t98, -t46 * t122 - t40 * t135 - t179 * t92 - t80 * t98 + t256, -t136 * t40 - t235 * t46 + t80 * t97 + t90 * t92 - t264, t103 * t50 + t104 * t47 + t40 * t92 + t46 * t80 + t60 * t86 + t61 * t82, t101 * t194 - t51 * t245 (t101 * t174 - t171 * t99) * t98 + (-t171 * t52 - t260 + (-t101 * t171 - t174 * t99) * qJD(5)) * t135, -t101 * t97 - t136 * t51 + t194 * t292 - t90 * t245, -t174 * t255 - t136 * t52 + t97 * t99 + (-t135 * t232 + t258) * t292, -t292 * t97 - t74, -t215 * t136 - t20 * t97 + t43 * t99 + t85 * t52 + ((-qJD(5) * t84 - t31) * t292 + t75 * t90 + t64 * qJD(5) * t135) * t171 + ((-qJD(5) * t75 + t44) * t292 - t84 * t90 - t32 * t135 - t64 * t98) * t174, t43 * t101 + t289 * t135 - t192 * t136 - t191 * t292 + t21 * t97 + t64 * t261 + t270 * t90 - t85 * t51, t10 * t97 + t287 * t135 - t136 * t3 - t35 * t258 - t292 * t5 + t39 * t90 + t49 * t52 + t9 * t99, t101 * t5 - t38 * t52 - t39 * t51 - t4 * t99 + t201 * t98 + (qJD(5) * t202 + t1 * t174 + t171 * t3) * t135, t1 * t136 - t101 * t9 - t11 * t97 + t290 * t135 - t35 * t261 + t292 * t4 - t38 * t90 + t49 * t51, t1 * t38 + t10 * t5 + t11 * t4 + t3 * t39 + t35 * t9 + t49 * t8; 0, 0, 0, -t173 * t243, t236 * t178, 0, 0, 0, t178 * pkin(1) * t173, pkin(1) * t243, t247, t76, t65, 0, 0, t167 * t95 + (-t122 * t234 - t167 * t233) * pkin(2) + t188, t96 * t167 + (-t167 * t218 - t234 * t235) * pkin(2) + t187, -t251 * t122 - t156 * t180 - t159 * t90 + t205 * t235 + t200, t122 * t83 + t167 * t205 + t189, t251 * t167 + t235 * t83 + t185, -t156 * t47 + t159 * t50 + t205 * t82 - t251 * t86 - t80 * t83, t12, t6, t15, t16, t248, -t226 + t156 * t52 + t252 * t99 + ((-t77 + t227) * t174 + (-qJD(5) * t155 + t62) * t171) * t292 + t198, -t156 * t51 + (-t221 + t272) * t292 + t252 * t101 + (t155 * t90 - t227 * t292 - t214) * t171 + t216, -t226 + t130 * t52 + t254 * t99 + (-t155 * t232 + t174 * t227 + t23) * t292 + t199, -t101 * t23 + t22 * t99 + t2 + (-t101 * t227 - t265 + t155 * t51 + (-t155 * t99 - t11) * qJD(5)) * t174 + (-t99 * t227 - t155 * t52 + (t101 * t155 - t10) * qJD(5) + t219) * t171, -t155 * t262 + t130 * t51 - t254 * t101 + (t171 * t227 - t22 + t221) * t292 + t190, -t10 * t23 - t11 * t22 + t130 * t8 + t155 * t181 - t202 * t227 + t254 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247, t76, t65, 0, 0, t188 + t263, -t167 * t93 + t187, pkin(3) * t90 - qJ(4) * t180 - t239 * t122 - t235 * t94 + t200, t122 * t91 + t189 - t263, t239 * t167 + t235 * t91 + t185, -pkin(3) * t50 - qJ(4) * t47 - t239 * t86 - t80 * t91 - t82 * t94, t12, t6, t15, t16, t248, t225 + qJ(4) * t52 + t240 * t99 + (-t259 + (t69 + t230) * t171) * t292 + t198, -qJ(4) * t51 + (t220 + t271) * t292 + t240 * t101 + (-t214 - t257) * t171 + t216, t225 + t144 * t52 + t253 * t99 + (t171 * t230 + t25) * t292 + t199, -t101 * t25 + t24 * t99 + t2 + (-t265 - t281 * t51 + (t281 * t99 - t11) * qJD(5)) * t174 + (t281 * t52 + (-t101 * t281 - t10) * qJD(5) + t219) * t171, t171 * t257 + t144 * t51 + (-t24 - t220) * t292 - t253 * t101 + t190, -t10 * t25 - t11 * t24 + t144 * t8 - t181 * t281 + t253 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t247, -t167 ^ 2 - t282, t167 * t86 + t189, 0, 0, 0, 0, 0, t286, t184 - t249, t286, t295 * t174 + (-t52 + t210) * t171, t211 * t292 + t249 - t262, -t167 * t35 - t2 + t11 * t211 + (t10 * t292 + t1) * t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266, -t99 ^ 2 + t283, -t295, t89 + (-qJD(5) + t292) * t101, -t90, -t101 * t64 + t195, t20 * t292 + t64 * t99 - t192, -t63 * t99 + t195 - t267 - 0.2e1 * t279, pkin(5) * t51 - qJ(6) * t52 + (t11 - t21) * t101 + (t10 - t238) * t99, -0.2e1 * t269 + t101 * t63 - t35 * t99 + (0.2e1 * qJD(6) - t20) * t292 + t192, -pkin(5) * t3 + qJ(6) * t1 - t10 * t21 + t11 * t238 - t35 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90 + t266, -t295, -t283 - t212, -t11 * t292 + t267 + t3;];
tauc_reg  = t7;
