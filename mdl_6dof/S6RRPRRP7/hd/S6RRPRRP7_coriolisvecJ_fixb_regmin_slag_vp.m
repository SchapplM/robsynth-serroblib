% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:06
% EndTime: 2019-03-09 12:20:16
% DurationCPUTime: 3.51s
% Computational Cost: add. (5577->409), mult. (12394->524), div. (0->0), fcn. (8041->6), ass. (0->210)
t204 = qJD(2) - qJD(4);
t148 = cos(qJ(4));
t146 = sin(qJ(2));
t205 = qJD(1) * qJD(2);
t197 = t146 * t205;
t145 = sin(qJ(4));
t149 = cos(qJ(2));
t210 = qJD(4) * t148;
t211 = qJD(4) * t145;
t212 = qJD(2) * t149;
t292 = t145 * t212 + t146 * t210 - t149 * t211;
t53 = qJD(1) * t292 - t148 * t197;
t215 = qJD(1) * t149;
t216 = qJD(1) * t146;
t95 = -t145 * t215 + t148 * t216;
t293 = t95 * t204 + t53;
t130 = pkin(7) * t216;
t106 = pkin(8) * t216 - t130;
t291 = qJD(3) - t106;
t227 = t148 * t149;
t100 = t145 * t146 + t227;
t92 = t100 * qJD(1);
t276 = qJD(5) + t92;
t144 = sin(qJ(5));
t147 = cos(qJ(5));
t97 = -qJD(1) * pkin(1) - pkin(2) * t215 - qJ(3) * t216;
t77 = pkin(3) * t215 - t97;
t40 = pkin(4) * t92 - pkin(9) * t95 + t77;
t150 = -pkin(2) - pkin(3);
t199 = t150 * qJD(2);
t80 = t199 + t291;
t131 = pkin(7) * t215;
t108 = -pkin(8) * t215 + t131;
t140 = qJD(2) * qJ(3);
t96 = t108 + t140;
t50 = t145 * t80 + t148 * t96;
t45 = -pkin(9) * t204 + t50;
t11 = t144 * t40 + t147 * t45;
t8 = qJ(6) * t276 + t11;
t290 = t276 * t8;
t185 = t147 * t204;
t64 = t144 * t95 + t185;
t289 = t276 * t64;
t66 = -t144 * t204 + t147 * t95;
t288 = t276 * t66;
t172 = qJ(3) * t148 + t145 * t150;
t105 = -pkin(9) + t172;
t242 = t105 * t53;
t275 = -t145 * qJ(3) + t148 * t150;
t81 = t148 * qJD(3) + t275 * qJD(4);
t255 = t81 * t276;
t49 = -t145 * t96 + t148 * t80;
t44 = pkin(4) * t204 - t49;
t16 = t64 * pkin(5) - t66 * qJ(6) + t44;
t279 = t276 * t16;
t287 = t279 + t242 + t255;
t262 = pkin(9) * t53;
t286 = t279 - t262;
t162 = t100 * qJD(4);
t61 = qJD(2) * t100 - t162;
t155 = t61 * qJD(1);
t208 = qJD(5) * t144;
t30 = qJD(5) * t185 - t147 * t155 + t208 * t95;
t232 = t30 * t144;
t277 = t147 * t276;
t285 = t277 * t66 - t232;
t236 = t147 * t53;
t168 = -t208 * t276 + t236;
t230 = qJD(4) * t276;
t213 = qJD(2) * t148;
t94 = t144 * t216 + t147 * t213;
t284 = t145 * (t204 * t66 + t168) + (t147 * t230 - t30) * t148 - t94 * t276;
t241 = t144 * t53;
t283 = t276 * t277 - t66 * t95 + t241;
t229 = qJD(5) * t66;
t31 = t144 * t155 + t229;
t282 = t144 * (t31 + t288) + t147 * (t30 + t289);
t281 = -0.2e1 * t205;
t58 = t106 * t148 + t108 * t145;
t280 = t58 - t81;
t246 = qJD(4) * t172 + t148 * t108 + t145 * t291;
t278 = t276 * t44;
t265 = pkin(7) - pkin(8);
t115 = t265 * t146;
t116 = t265 * t149;
t274 = t148 * t115 - t145 * t116;
t238 = t144 * t276;
t269 = t238 * t276 - t64 * t95 - t236;
t268 = t204 ^ 2;
t174 = pkin(5) * t144 - qJ(6) * t147;
t267 = t144 * qJD(6) - t276 * t174;
t266 = t66 ^ 2;
t264 = pkin(5) * t53;
t263 = pkin(5) * t95;
t10 = -t144 * t45 + t147 * t40;
t221 = qJD(6) - t10;
t7 = -pkin(5) * t276 + t221;
t261 = t7 * t95;
t260 = t8 * t95;
t259 = t10 * t95;
t258 = t11 * t276;
t257 = t11 * t95;
t256 = t66 * t64;
t254 = t276 * t95;
t252 = t95 * t92;
t251 = t267 + t50;
t250 = -t267 - t246;
t55 = t95 * pkin(4) + t92 * pkin(9);
t249 = t144 * t55 + t147 * t49;
t128 = qJ(3) * t215;
t87 = t150 * t216 + t128;
t41 = -t55 + t87;
t248 = t144 * t41 + t147 * t58;
t101 = -t145 * t149 + t146 * t148;
t113 = -t149 * pkin(2) - t146 * qJ(3) - pkin(1);
t99 = t149 * pkin(3) - t113;
t48 = pkin(4) * t100 - pkin(9) * t101 + t99;
t69 = t115 * t145 + t116 * t148;
t247 = t144 * t48 + t147 * t69;
t245 = pkin(9) * qJD(5);
t244 = qJ(6) * t95;
t243 = qJD(2) * pkin(2);
t240 = t144 * t64;
t239 = t144 * t66;
t237 = t147 * t31;
t235 = t147 * t61;
t234 = t147 * t64;
t233 = t147 * t66;
t231 = t53 * qJ(6);
t152 = qJD(1) ^ 2;
t225 = t149 * t152;
t151 = qJD(2) ^ 2;
t224 = t151 * t146;
t223 = t151 * t149;
t134 = t146 * qJD(3);
t196 = t149 * t205;
t219 = qJ(3) * t196 + qJD(1) * t134;
t218 = qJ(3) * t212 + t134;
t141 = t146 ^ 2;
t217 = -t149 ^ 2 + t141;
t214 = qJD(2) * t146;
t209 = qJD(5) * t105;
t207 = qJD(5) * t147;
t203 = t276 * t245;
t202 = t92 ^ 2 - t95 ^ 2;
t201 = t276 * t209;
t200 = t146 * t225;
t107 = t265 * t214;
t139 = qJD(2) * qJD(3);
t85 = -qJD(1) * t107 + t139;
t124 = pkin(7) * t196;
t98 = -pkin(8) * t196 + t124;
t166 = -t145 * t98 - t148 * t85 - t80 * t210 + t211 * t96;
t17 = t53 * pkin(4) + (pkin(9) * t162 + (-pkin(9) * t227 + (-pkin(9) * t145 + t150) * t146) * qJD(2)) * qJD(1) + t219;
t193 = -t144 * t166 - t147 * t17 + t45 * t207 + t40 * t208;
t190 = pkin(1) * t281;
t189 = qJD(3) - t243;
t188 = qJD(1) * t113 + t97;
t173 = -t145 * t85 + t148 * t98 - t96 * t210 - t80 * t211;
t5 = pkin(5) * t31 + qJ(6) * t30 - qJD(6) * t66 - t173;
t184 = -t5 - t203;
t183 = t146 * t199;
t182 = t5 - t201;
t180 = -t144 * t8 + t147 * t7;
t179 = t144 * t7 + t147 * t8;
t175 = pkin(5) * t147 + qJ(6) * t144;
t112 = -pkin(4) - t175;
t76 = pkin(2) * t197 - t219;
t89 = pkin(2) * t214 - t218;
t171 = -pkin(7) * t151 - qJD(1) * t89 - t76;
t170 = t16 * t66 + t193;
t169 = -t207 * t276 - t241;
t167 = -t262 + t278;
t165 = t144 * t17 - t147 * t166 + t40 * t207 - t208 * t45;
t60 = -t146 * t213 + t292;
t72 = t183 + t218;
t25 = t60 * pkin(4) - t61 * pkin(9) + t72;
t109 = qJD(2) * t116;
t36 = t274 * qJD(4) - t148 * t107 + t145 * t109;
t164 = t144 * t25 + t147 * t36 + t48 * t207 - t208 * t69;
t163 = -t242 - t278;
t161 = -t77 * t95 + t173;
t159 = t77 * t92 + t166;
t1 = qJD(6) * t276 + t165 + t231;
t2 = t193 - t264;
t158 = qJD(5) * t180 + t1 * t147 + t2 * t144;
t37 = qJD(4) * t69 - t145 * t107 - t148 * t109;
t110 = -pkin(7) * t197 + t139;
t111 = t130 + t189;
t114 = t131 + t140;
t156 = t110 * t149 + (t111 * t149 + (-t114 + t131) * t146) * qJD(2);
t91 = t144 * t213 - t147 * t216;
t153 = t91 * t276 + t64 * t211 + (-t144 * t230 - t31) * t148 + (-qJD(2) * t64 + t169) * t145;
t104 = pkin(4) - t275;
t103 = pkin(2) * t216 - t128;
t78 = -t112 - t275;
t63 = qJD(1) * t183 + t219;
t38 = pkin(5) * t66 + qJ(6) * t64;
t33 = t101 * t174 - t274;
t21 = -pkin(5) * t100 + t144 * t69 - t147 * t48;
t20 = qJ(6) * t100 + t247;
t19 = t144 * t49 - t147 * t55 - t263;
t18 = t244 + t249;
t13 = t144 * t58 - t147 * t41 + t263;
t12 = -t244 + t248;
t9 = -t30 + t289;
t6 = t174 * t61 + (qJD(5) * t175 - qJD(6) * t147) * t101 + t37;
t4 = -t60 * pkin(5) + t247 * qJD(5) + t144 * t36 - t147 * t25;
t3 = qJ(6) * t60 + qJD(6) * t100 + t164;
t14 = [0, 0, 0, 0.2e1 * t146 * t196, t217 * t281, t223, -t224, 0, -pkin(7) * t223 + t146 * t190, pkin(7) * t224 + t149 * t190, t149 * t171 + t188 * t214, t156, t146 * t171 - t188 * t212, pkin(7) * t156 + t76 * t113 + t97 * t89, t101 * t155 + t95 * t61, -t100 * t155 - t101 * t53 - t95 * t60 - t61 * t92, -t61 * t204, t60 * t204, 0, t63 * t100 + t204 * t37 + t99 * t53 + t77 * t60 + t72 * t92, t63 * t101 + t155 * t99 + t204 * t36 + t77 * t61 + t72 * t95, t61 * t233 + (-t30 * t147 - t208 * t66) * t101 (-t234 - t239) * t61 + (t232 - t237 + (-t233 + t240) * qJD(5)) * t101, -t30 * t100 + t101 * t168 + t235 * t276 + t66 * t60, -t31 * t100 + t101 * t169 - t238 * t61 - t64 * t60, t100 * t53 + t276 * t60, -t193 * t100 + t10 * t60 + t37 * t64 - t274 * t31 + ((-qJD(5) * t69 + t25) * t276 + t48 * t53 + t44 * qJD(5) * t101) * t147 + ((-qJD(5) * t48 - t36) * t276 - t69 * t53 - t173 * t101 + t44 * t61) * t144, -t164 * t276 - t247 * t53 - t165 * t100 - t11 * t60 + t37 * t66 + t274 * t30 + t44 * t235 + (-t147 * t173 - t208 * t44) * t101, t16 * t144 * t61 - t2 * t100 - t21 * t53 + t33 * t31 - t4 * t276 + t6 * t64 - t7 * t60 + (t5 * t144 + t16 * t207) * t101, -t20 * t31 - t21 * t30 - t3 * t64 + t4 * t66 + t180 * t61 + (-qJD(5) * t179 - t1 * t144 + t147 * t2) * t101, -t16 * t235 + t1 * t100 + t20 * t53 + t3 * t276 + t33 * t30 - t6 * t66 + t8 * t60 + (-t5 * t147 + t16 * t208) * t101, t1 * t20 + t16 * t6 + t2 * t21 + t3 * t8 + t33 * t5 + t4 * t7; 0, 0, 0, -t200, t217 * t152, 0, 0, 0, t152 * pkin(1) * t146, pkin(1) * t225 (t103 * t149 - t146 * t97) * qJD(1) ((t114 - t140) * t146 + (-t111 + t189) * t149) * qJD(1), 0.2e1 * t139 + (t103 * t146 + t149 * t97) * qJD(1), t110 * qJ(3) + t114 * qJD(3) - t97 * t103 + (t114 * t146 + (-t111 - t243) * t149) * qJD(1) * pkin(7), -t252, t202, 0, t293, 0, t204 * t246 - t87 * t92 - t161, -t204 * t280 - t87 * t95 - t159, -t285, t282, -t283, t269, t254, t259 + t104 * t31 + t246 * t64 + (-t173 + (-t41 - t209) * t276) * t147 + (t276 * t280 + t163) * t144, -t104 * t30 + t248 * t276 - t257 + t246 * t66 + (t173 + t201) * t144 + (t163 - t255) * t147, t13 * t276 - t287 * t144 + t182 * t147 - t250 * t64 + t31 * t78 - t261, t12 * t64 - t13 * t66 + (-t105 * t31 - t64 * t81 - t7 * t92 - t1 + (t105 * t66 - t7) * qJD(5)) * t147 + (-t105 * t30 + t66 * t81 + t8 * t92 - t2 + (t105 * t64 + t8) * qJD(5)) * t144, -t12 * t276 + t182 * t144 + t287 * t147 + t250 * t66 + t30 * t78 + t260, t5 * t78 + (t147 * t81 - t12) * t8 + (t144 * t81 - t13) * t7 - t250 * t16 + t158 * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200, 0, -t141 * t152 - t151, -qJD(2) * t114 + t216 * t97 + t124, 0, 0, 0, 0, 0, -t145 * t268 - t92 * t216, -t148 * t268 - t95 * t216, 0, 0, 0, 0, 0, t153, -t284, t153, t94 * t64 - t91 * t66 + (-t234 + t239) * t210 + (-t232 - t237 + (t233 + t240) * qJD(5)) * t145, t284, -t7 * t91 - t8 * t94 + (qJD(4) * t179 - t5) * t148 + (-t16 * t204 + t158) * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, -t202, -t92 * t204 + t155, -t293, 0, -t204 * t50 + t161, -t204 * t49 + t159, t285, -t282, t283, -t269, -t254, -pkin(4) * t31 - t259 - t50 * t64 + (t173 + (-t55 - t245) * t276) * t147 + (t276 * t49 + t167) * t144, pkin(4) * t30 + t249 * t276 + t257 - t50 * t66 + (-t173 + t203) * t144 + t167 * t147, t112 * t31 + t286 * t144 + t184 * t147 + t19 * t276 - t251 * t64 + t261, t18 * t64 - t19 * t66 + (t1 + t276 * t7 + (-t31 + t229) * pkin(9)) * t147 + (t2 - t290 + (qJD(5) * t64 - t30) * pkin(9)) * t144, t112 * t30 + t184 * t144 - t286 * t147 - t18 * t276 + t251 * t66 - t260, pkin(9) * t158 + t5 * t112 - t16 * t251 - t8 * t18 - t7 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, -t64 ^ 2 + t266, t9, -t31 + t288, t53, -t44 * t66 - t193 + t258, t10 * t276 + t44 * t64 - t165, -t38 * t64 - t170 + t258 + 0.2e1 * t264, pkin(5) * t30 - t31 * qJ(6) + (-t11 + t8) * t66 + (t7 - t221) * t64, 0.2e1 * t231 - t16 * t64 + t38 * t66 + (0.2e1 * qJD(6) - t10) * t276 + t165, -t2 * pkin(5) + t1 * qJ(6) - t7 * t11 - t16 * t38 + t221 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 + t256, t9, -t276 ^ 2 - t266, t170 - t264 - t290;];
tauc_reg  = t14;
