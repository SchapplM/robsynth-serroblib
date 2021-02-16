% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:18
% EndTime: 2021-01-16 02:21:34
% DurationCPUTime: 4.63s
% Computational Cost: add. (3253->435), mult. (7753->564), div. (0->0), fcn. (6185->14), ass. (0->231)
t145 = sin(pkin(11));
t148 = cos(pkin(11));
t152 = sin(qJ(3));
t155 = cos(qJ(3));
t105 = t145 * t155 + t148 * t152;
t303 = t105 * qJD(2);
t316 = qJD(6) + t303;
t151 = sin(qJ(6));
t318 = t151 * t316;
t154 = cos(qJ(6));
t237 = qJD(2) * qJD(3);
t220 = t152 * t237;
t173 = qJDD(2) * t105 - t145 * t220;
t219 = t155 * t237;
t55 = t148 * t219 + t173;
t50 = qJDD(6) + t55;
t47 = t154 * t50;
t192 = -t316 * t318 + t47;
t241 = t151 * qJD(3);
t254 = t148 * t155;
t223 = qJD(2) * t254;
t246 = qJD(2) * t152;
t96 = t145 * t246 - t223;
t72 = -t154 * t96 + t241;
t319 = t316 * t72;
t150 = qJ(4) + pkin(8);
t147 = sin(pkin(6));
t153 = sin(qJ(2));
t247 = qJD(1) * t153;
t225 = t147 * t247;
t210 = t150 * qJD(2) + t225;
t149 = cos(pkin(6));
t248 = qJD(1) * t149;
t67 = t152 * t248 + t155 * t210;
t60 = t145 * t67;
t66 = -t152 * t210 + t155 * t248;
t33 = t148 * t66 - t60;
t311 = -qJD(5) + t33;
t141 = qJ(3) + pkin(11);
t138 = sin(t141);
t156 = cos(qJ(2));
t264 = cos(pkin(10));
t213 = t264 * t156;
t146 = sin(pkin(10));
t259 = t146 * t153;
t93 = -t149 * t213 + t259;
t214 = t264 * t153;
t258 = t146 * t156;
t95 = t149 * t258 + t214;
t204 = g(1) * t95 + g(2) * t93;
t255 = t147 * t156;
t306 = -g(3) * t255 + t204;
t170 = t306 * t138;
t114 = t150 * t152;
t115 = t150 * t155;
t71 = -t145 * t114 + t148 * t115;
t317 = -t71 * qJDD(3) - t170;
t257 = t147 * t153;
t102 = t149 * t155 - t152 * t257;
t139 = cos(t141);
t260 = t146 * t147;
t92 = t149 * t259 - t213;
t187 = t138 * t260 - t92 * t139;
t215 = t147 * t264;
t94 = t149 * t214 + t258;
t58 = -t138 * t215 + t94 * t139;
t82 = t149 * t138 + t139 * t257;
t315 = -g(1) * t187 - g(2) * t58 - g(3) * t82;
t234 = t155 * qJDD(2);
t235 = t152 * qJDD(2);
t199 = t145 * t235 - t148 * t234;
t314 = 0.2e1 * t303 * qJD(3) + t199;
t224 = qJD(1) * t255;
t216 = qJD(3) * t150;
t89 = t155 * qJD(4) - t152 * t216;
t90 = -t152 * qJD(4) - t155 * t216;
t278 = -t105 * t224 + t145 * t89 - t148 * t90;
t304 = -t145 * t152 + t254;
t277 = t145 * t90 + t148 * t89 - t304 * t224;
t268 = t154 * t316;
t243 = qJD(6) * t151;
t313 = -t243 * t316 + t47;
t91 = t303 ^ 2;
t312 = -t96 ^ 2 - t91;
t290 = t303 * pkin(5);
t310 = t290 - t311;
t188 = t92 * t138 + t139 * t260;
t252 = t151 * t156;
t120 = t147 * t252;
t256 = t147 * t155;
t103 = t149 * t152 + t153 * t256;
t45 = -t148 * t102 + t145 * t103;
t309 = t154 * t45 + t120;
t178 = -t94 * t138 - t139 * t215;
t236 = t149 * qJDD(1);
t123 = t155 * t236;
t238 = qJD(1) * qJD(2);
t78 = qJDD(2) * pkin(8) + (qJDD(1) * t153 + t156 * t238) * t147;
t167 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t248 + t78;
t195 = t210 * qJD(3);
t20 = qJDD(3) * pkin(3) - t152 * t167 - t155 * t195 + t123;
t21 = (-t195 + t236) * t152 + t167 * t155;
t281 = t145 * t21 - t148 * t20;
t81 = t138 * t257 - t149 * t139;
t307 = -g(1) * t188 - g(2) * t178 + g(3) * t81 - t281;
t242 = qJD(6) * t154;
t272 = t151 * t50;
t305 = t242 * t316 + t272;
t135 = t155 * pkin(3) + pkin(2);
t86 = -qJD(2) * t135 + qJD(4) - t224;
t163 = -qJ(5) * t303 + t86;
t38 = t96 * pkin(4) + t163;
t302 = t303 * t38 + qJDD(5);
t70 = t148 * t114 + t145 * t115;
t301 = -t70 * qJDD(3) + t139 * t306;
t157 = qJD(3) ^ 2;
t218 = qJDD(1) * t255;
t221 = t153 * t238;
t200 = t147 * t221 - t218;
t288 = g(3) * t156;
t300 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t157 + t147 * (t221 - t288) - t200 + t204;
t245 = qJD(2) * t153;
t176 = t103 * qJD(3);
t226 = qJD(2) * t255;
t208 = t152 * t226;
t162 = -t176 - t208;
t65 = qJD(3) * t102 + t155 * t226;
t32 = t145 * t162 + t148 * t65;
t46 = t145 * t102 + t148 * t103;
t299 = t147 * (t156 * t55 - t245 * t303) + t32 * qJD(3) + t46 * qJDD(3);
t30 = t145 * t65 - t148 * t162;
t98 = t105 * qJD(3);
t54 = qJD(2) * t98 + t199;
t298 = t147 * (-t156 * t54 + t245 * t96) - t30 * qJD(3) - t45 * qJDD(3);
t136 = pkin(3) * t246;
t212 = t96 * qJ(5) + t136;
t295 = pkin(4) + pkin(9);
t142 = qJDD(3) * qJ(5);
t9 = t145 * t20 + t148 * t21;
t4 = -qJD(3) * qJD(5) - t142 - t9;
t3 = -pkin(5) * t54 - t4;
t297 = (t295 * t303 + t212) * t316 + t3 + t315;
t292 = t54 * pkin(4);
t291 = t96 * pkin(5);
t289 = pkin(3) * t152;
t287 = t145 * pkin(3);
t286 = t148 * pkin(3);
t193 = -t105 * qJ(5) - t135;
t39 = -t295 * t304 + t193;
t284 = t39 * t50;
t74 = t154 * qJD(3) + t151 * t96;
t283 = t74 * t96;
t282 = t96 * t72;
t244 = qJD(3) * t152;
t101 = qJD(3) * t254 - t145 * t244;
t280 = t101 * pkin(5) + t278;
t279 = -t98 * pkin(5) + t277;
t273 = t148 * t67;
t64 = qJD(3) * pkin(3) + t66;
t29 = t145 * t64 + t273;
t276 = -t95 * t135 - t92 * t150;
t275 = qJD(2) * pkin(2);
t24 = -qJD(3) * qJ(5) - t29;
t13 = -t24 - t291;
t274 = t13 * t151;
t270 = t151 * t98;
t26 = -qJD(6) * t241 + t154 * qJDD(3) + t151 * t54 + t96 * t242;
t267 = t26 * t154;
t266 = t55 * qJ(5);
t262 = qJDD(3) * pkin(4);
t250 = qJDD(1) - g(3);
t143 = t152 ^ 2;
t249 = -t155 ^ 2 + t143;
t233 = g(3) * t257;
t137 = pkin(3) * t244;
t229 = t154 * t255;
t228 = qJDD(5) + t281;
t134 = -pkin(4) - t286;
t227 = t147 * t245;
t222 = g(3) * (t135 * t255 + t150 * t257);
t31 = t145 * t66 + t273;
t28 = t148 * t64 - t60;
t217 = -t93 * t135 + t94 * t150;
t209 = t151 * qJDD(3) - t154 * t54;
t207 = t146 * pkin(3) * t256 + t92 * t289;
t205 = g(1) * t92 - g(2) * t94;
t203 = qJD(5) - t28;
t183 = -t101 * qJ(5) - t105 * qJD(5) + t137;
t37 = t98 * pkin(4) + t183;
t202 = -t37 + t225;
t201 = pkin(4) * t139 + qJ(5) * t138;
t12 = -t295 * qJD(3) + t203 + t290;
t25 = t295 * t96 + t163;
t6 = t151 * t12 + t154 * t25;
t196 = t102 * pkin(3);
t158 = qJD(2) ^ 2;
t194 = qJDD(2) * t156 - t153 * t158;
t191 = -g(1) * t146 + g(2) * t264;
t186 = -t151 * t45 + t229;
t181 = t194 * t147;
t179 = -t268 * t316 - t272;
t112 = -t224 - t275;
t175 = -t112 * qJD(2) - t205 - t78;
t174 = -t315 - t9;
t171 = (-t94 * t152 - t155 * t215) * pkin(3);
t169 = t30 * t303 - t32 * t96 + t45 * t55 - t46 * t54;
t166 = t31 * qJD(3) + t307;
t43 = t105 * pkin(5) + t70;
t165 = -t13 * t98 + t3 * t304 + t43 * t50 + t205;
t164 = -pkin(8) * qJDD(3) + (t112 + t224 - t275) * qJD(3);
t53 = pkin(3) * t220 - qJDD(2) * t135 + qJDD(4) + t200;
t161 = -qJD(5) * t303 - t266 + t53;
t160 = t53 - t306;
t159 = -t277 * t96 + t278 * t303 - t71 * t54 + t70 * t55 + t205 - t233;
t131 = qJ(5) + t287;
t129 = -pkin(9) + t134;
t87 = qJD(3) * t96;
t51 = -pkin(4) * t304 + t193;
t44 = pkin(5) * t304 + t71;
t40 = pkin(4) * t303 + t212;
t27 = qJD(6) * t74 + t209;
t23 = -qJD(3) * pkin(4) + t203;
t22 = t295 * t98 + t183;
t16 = t31 - t291;
t11 = t161 + t292;
t10 = t295 * t54 + t161;
t7 = t228 - t262;
t5 = t154 * t12 - t151 * t25;
t2 = pkin(5) * t55 - t295 * qJDD(3) + t228;
t1 = t154 * t2;
t8 = [t250, 0, t181, (-qJDD(2) * t153 - t156 * t158) * t147, 0, 0, 0, 0, 0, t102 * qJDD(3) + t155 * t181 + (-t176 - 0.2e1 * t208) * qJD(3), -t65 * qJD(3) - t103 * qJDD(3) + (-t152 * t194 - t156 * t219) * t147, t298, -t299, t169, -t28 * t30 + t29 * t32 + t281 * t45 + t9 * t46 - g(3) + (-t156 * t53 + t245 * t86) * t147, t169, -t298, t299, t23 * t30 - t24 * t32 - t4 * t46 + t7 * t45 - g(3) + (-t11 * t156 + t245 * t38) * t147, 0, 0, 0, 0, 0, (qJD(6) * t186 - t151 * t227 + t154 * t30) * t316 + t309 * t50 + t32 * t72 + t46 * t27, -(t309 * qJD(6) + t151 * t30 + t154 * t227) * t316 + t186 * t50 + t32 * t74 + t46 * t26; 0, qJDD(2), t218 + t306, -t250 * t257 - t205, t143 * qJDD(2) + 0.2e1 * t152 * t219, 0.2e1 * t152 * t234 - 0.2e1 * t237 * t249, qJDD(3) * t152 + t157 * t155, qJDD(3) * t155 - t157 * t152, 0, t164 * t152 + t300 * t155, -t300 * t152 + t164 * t155, -t96 * t225 - t53 * t304 - t135 * t54 + t86 * t98 + (t96 * t289 - t278) * qJD(3) + t301, -t303 * t225 + t86 * t101 + t53 * t105 - t135 * t55 + (t289 * t303 - t277) * qJD(3) + t317, -t28 * t101 + t105 * t281 - t29 * t98 + t304 * t9 + t159, t9 * t71 + t281 * t70 - t53 * t135 - g(1) * t276 - g(2) * t217 - t222 + (-t225 + t137) * t86 + t277 * t29 - t278 * t28, t23 * t101 + t7 * t105 + t24 * t98 - t304 * t4 + t159, t278 * qJD(3) + t11 * t304 + t202 * t96 - t38 * t98 - t51 * t54 - t301, t277 * qJD(3) - t38 * t101 - t11 * t105 + t202 * t303 - t51 * t55 - t317, t11 * t51 + t38 * t37 - t4 * t71 + t7 * t70 - g(1) * (-t201 * t95 + t276) - g(2) * (-t201 * t93 + t217) - t222 - t277 * t24 + t278 * t23 + (-t201 * t288 - t247 * t38) * t147, t74 * t270 - (t26 * t151 + t242 * t74) * t304, (-t151 * t72 + t154 * t74) * t98 - (-t151 * t27 + t267 + (-t151 * t74 - t154 * t72) * qJD(6)) * t304, t74 * t101 + t26 * t105 + t270 * t316 - t304 * t305, -t72 * t101 - t27 * t105 + t98 * t268 - t304 * t313, t101 * t316 + t105 * t50, t1 * t105 + t5 * t101 + t44 * t27 + t279 * t72 + (-t10 * t105 + t138 * t204 - t22 * t316 - t284) * t151 + (t280 * t316 + t165) * t154 + (t247 * t318 - g(3) * (t138 * t252 + t153 * t154)) * t147 + ((-t151 * t43 - t154 * t39) * t316 - t6 * t105 - t304 * t274) * qJD(6), -t6 * t101 + t44 * t26 + t279 * t74 + (-t284 - (qJD(6) * t12 + t10) * t105 - t13 * qJD(6) * t304 + (-qJD(6) * t43 - t22 + t225) * t316 + t170) * t154 + (-(-qJD(6) * t25 + t2) * t105 + t233 + (qJD(6) * t39 - t280) * t316 - t165) * t151; 0, 0, 0, 0, -t152 * t158 * t155, t249 * t158, t235, t234, qJDD(3), -g(3) * t102 + t152 * t175 + t191 * t256 + t123, g(3) * t103 + (-t147 * t191 - t236) * t152 + t175 * t155, -t86 * t303 + (qJDD(3) * t148 - t246 * t96) * pkin(3) + t166, t33 * qJD(3) + t86 * t96 + (-qJDD(3) * t145 - t246 * t303) * pkin(3) + t174, (t29 - t31) * t303 + (-t28 + t33) * t96 + (-t145 * t54 - t148 * t55) * pkin(3), -g(1) * t207 - g(2) * t171 - g(3) * t196 - t136 * t86 + t28 * t31 - t281 * t286 + t9 * t287 - t29 * t33, -t131 * t54 + t134 * t55 + (-t24 - t31) * t303 + (t23 + t311) * t96, t40 * t96 + (-pkin(4) + t134) * qJDD(3) - t166 + t302, t131 * qJDD(3) - t38 * t96 + t40 * t303 + t142 + (0.2e1 * qJD(5) - t33) * qJD(3) - t174, -t4 * t131 + t7 * t134 - t38 * t40 - t23 * t31 - g(1) * (pkin(4) * t188 + qJ(5) * t187 + t207) - g(2) * (pkin(4) * t178 + t58 * qJ(5) + t171) - g(3) * (-t81 * pkin(4) + t82 * qJ(5) + t196) + t311 * t24, -t318 * t74 + t267, (-t316 * t74 - t27) * t154 + (-t26 + t319) * t151, t192 + t283, t179 - t282, t316 * t96, -t16 * t268 + t131 * t27 + t5 * t96 + t310 * t72 + (t154 * t303 + t242) * t13 + t313 * t129 + t297 * t151, -t305 * t129 - t13 * t243 + t131 * t26 + t297 * t154 + t16 * t318 - t274 * t303 + t310 * t74 - t6 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, (-t96 + t223) * qJD(3) + t173, t312, t28 * t303 + t29 * t96 + t160, t312, -t314, t87 - t55, t292 - t266 - t24 * t96 + (-qJD(5) - t23) * t303 + t160, 0, 0, 0, 0, 0, t179 + t282, -t192 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87 + t55, -t303 * t96 + qJDD(3), -t91 - t157, t24 * qJD(3) - t262 + t302 - t307, 0, 0, 0, 0, 0, -qJD(3) * t72 + t192, -qJD(3) * t74 + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t72, -t72 ^ 2 + t74 ^ 2, t26 + t319, -t209 + (-qJD(6) + t316) * t74, t50, -t25 * t242 - t151 * t10 - t12 * t243 + t1 + t6 * t316 - t13 * t74 - g(1) * (-t95 * t151 - t154 * t188) - g(2) * (-t93 * t151 - t154 * t178) - g(3) * (t81 * t154 + t120), t25 * t243 - t154 * t10 - t12 * t242 - t151 * t2 + t5 * t316 + t13 * t72 - g(1) * (t151 * t188 - t95 * t154) - g(2) * (t151 * t178 - t93 * t154) - g(3) * (-t81 * t151 + t229);];
tau_reg = t8;
