% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:08
% EndTime: 2019-03-09 03:26:14
% DurationCPUTime: 2.72s
% Computational Cost: add. (4987->419), mult. (9777->518), div. (0->0), fcn. (6630->10), ass. (0->200)
t155 = sin(qJ(3));
t158 = cos(qJ(3));
t254 = sin(pkin(9));
t255 = cos(pkin(9));
t171 = t255 * t155 + t254 * t158;
t91 = t171 * qJD(1);
t284 = qJD(5) + t91;
t154 = sin(qJ(5));
t157 = cos(qJ(5));
t160 = -pkin(1) - pkin(7);
t112 = t160 * qJD(1) + qJD(2);
t237 = qJD(1) * t155;
t83 = -qJ(4) * t237 + t155 * t112;
t216 = t255 * t83;
t236 = qJD(1) * t158;
t84 = -qJ(4) * t236 + t158 * t112;
t80 = qJD(3) * pkin(3) + t84;
t45 = t254 * t80 + t216;
t38 = qJD(3) * pkin(8) + t45;
t103 = pkin(3) * t237 + qJD(1) * qJ(2) + qJD(4);
t98 = -t254 * t155 + t255 * t158;
t94 = t98 * qJD(1);
t46 = t91 * pkin(4) - t94 * pkin(8) + t103;
t17 = t154 * t46 + t157 * t38;
t8 = qJ(6) * t284 + t17;
t288 = t284 * t8;
t156 = sin(qJ(1));
t146 = g(1) * t156;
t159 = cos(qJ(1));
t283 = -g(2) * t159 + t146;
t208 = qJD(3) * t254;
t287 = qJD(1) * t208 - qJDD(1) * t255;
t209 = qJD(3) * t255;
t286 = qJD(1) * t209 + qJDD(1) * t254;
t162 = qJD(1) ^ 2;
t173 = -t162 * qJ(2) - t283;
t170 = t155 * t286 + t158 * t287;
t229 = t157 * qJD(3);
t232 = qJD(5) * t154;
t29 = -qJD(5) * t229 - t154 * qJDD(3) + t157 * t170 + t94 * t232;
t72 = t154 * qJD(3) + t157 * t94;
t93 = -t155 * t209 - t158 * t208;
t190 = -t98 * t29 + t93 * t72;
t204 = -qJD(5) * t171 - qJD(1);
t65 = t155 * t287 - t158 * t286;
t63 = -qJDD(5) + t65;
t59 = t157 * t63;
t92 = t155 * t208 - t158 * t209;
t285 = -t171 * t59 + (t154 * t204 - t157 * t92) * t284 + t190;
t148 = qJDD(1) * qJ(2);
t195 = g(1) * t159 + g(2) * t156;
t149 = qJD(1) * qJD(2);
t221 = 0.2e1 * t149;
t282 = 0.2e1 * t148 + t221 - t195;
t77 = t254 * t83;
t44 = t255 * t80 - t77;
t37 = -qJD(3) * pkin(4) - t44;
t70 = t154 * t94 - t229;
t18 = t70 * pkin(5) - t72 * qJ(6) + t37;
t126 = t254 * pkin(3) + pkin(8);
t260 = t126 * t63;
t281 = t18 * t284 + t260;
t147 = qJ(3) + pkin(9);
t136 = cos(t147);
t233 = qJD(5) * t126;
t220 = t284 * t233;
t135 = sin(t147);
t274 = g(3) * t135;
t280 = t136 * t283 + t220 - t274;
t277 = t72 ^ 2;
t275 = pkin(5) * t63;
t273 = g(3) * t136;
t272 = g(3) * t155;
t271 = g(3) * t157;
t270 = t136 * pkin(8);
t143 = t155 * pkin(3);
t269 = t17 * t284;
t268 = t70 * t91;
t267 = t72 * t70;
t266 = t72 * t94;
t265 = t94 * t70;
t231 = qJD(5) * t157;
t167 = -t157 * qJDD(3) - t154 * t170;
t30 = t72 * qJD(5) + t167;
t264 = -t154 * t30 - t70 * t231;
t111 = t160 * qJDD(1) + qJDD(2);
t100 = t158 * t111;
t223 = t158 * qJDD(1);
t226 = qJD(1) * qJD(4);
t227 = qJD(1) * qJD(3);
t235 = qJD(3) * t155;
t42 = -t158 * t226 - t112 * t235 + qJDD(3) * pkin(3) + t100 + (t155 * t227 - t223) * qJ(4);
t234 = qJD(3) * t158;
t50 = (-qJ(4) * qJD(1) + t112) * t234 + (-qJ(4) * qJDD(1) + t111 - t226) * t155;
t20 = t254 * t42 + t255 * t50;
t53 = t255 * t84 - t77;
t55 = pkin(3) * t236 + t94 * pkin(4) + t91 * pkin(8);
t263 = t154 * t55 + t157 * t53;
t246 = qJ(2) + t143;
t64 = pkin(4) * t171 - t98 * pkin(8) + t246;
t245 = qJ(4) - t160;
t104 = t245 * t155;
t105 = t245 * t158;
t68 = -t255 * t104 - t254 * t105;
t262 = t154 * t64 + t157 * t68;
t185 = pkin(5) * t154 - qJ(6) * t157;
t52 = t254 * t84 + t216;
t261 = -t154 * qJD(6) + t185 * t284 - t52;
t58 = t154 * t63;
t259 = t157 * t284;
t257 = t29 * t154;
t256 = t63 * qJ(6);
t253 = pkin(1) * qJDD(1);
t252 = qJD(5) * t98;
t251 = t156 * t154;
t250 = t156 * t157;
t249 = t159 * t154;
t248 = t159 * t157;
t16 = -t154 * t38 + t157 * t46;
t243 = qJD(6) - t16;
t242 = g(2) * t136 * t248 + t135 * t271;
t241 = t159 * pkin(1) + t156 * qJ(2);
t152 = t158 ^ 2;
t239 = t155 ^ 2 - t152;
t161 = qJD(3) ^ 2;
t238 = -t161 - t162;
t230 = t103 * qJD(1);
t228 = pkin(3) * t234 + qJD(2);
t225 = qJDD(3) * t155;
t224 = t155 * qJDD(1);
t222 = t136 * t146;
t219 = t98 * t232;
t218 = t98 * t231;
t217 = t255 * pkin(3);
t215 = t158 * t227;
t214 = -t156 * pkin(1) + t159 * qJ(2);
t213 = t284 * t72;
t15 = qJDD(3) * pkin(8) + t20;
t183 = qJDD(4) + t148 + t149 + (t215 + t224) * pkin(3);
t23 = -t65 * pkin(4) + pkin(8) * t170 + t183;
t210 = -t154 * t15 + t157 * t23 - t38 * t231 - t46 * t232;
t207 = t154 * t284;
t203 = qJDD(2) - t253;
t86 = t135 * t251 - t248;
t88 = t135 * t249 + t250;
t202 = g(1) * t88 + g(2) * t86;
t87 = t135 * t250 + t249;
t89 = t135 * t248 - t251;
t201 = -g(1) * t89 - g(2) * t87;
t177 = t157 * t15 + t154 * t23 + t46 * t231 - t38 * t232;
t1 = qJD(6) * t284 + t177 - t256;
t200 = t1 * t171 - t8 * t92;
t2 = qJDD(6) - t210 + t275;
t7 = -pkin(5) * t284 + t243;
t199 = t171 * t2 - t7 * t92;
t19 = -t254 * t50 + t255 * t42;
t14 = -qJDD(3) * pkin(4) - t19;
t3 = t30 * pkin(5) + t29 * qJ(6) - t72 * qJD(6) + t14;
t198 = t18 * t93 + t3 * t98;
t196 = t135 * pkin(4) - t270;
t193 = t14 * t98 + t37 * t93;
t192 = -t154 * t8 + t157 * t7;
t191 = -t171 * t29 - t72 * t92;
t189 = -t171 * t30 + t70 * t92;
t188 = -t284 * t93 + t63 * t98;
t81 = -t158 * qJD(4) + t245 * t235;
t82 = -qJD(3) * t105 - t155 * qJD(4);
t48 = t254 * t82 - t255 * t81;
t67 = -t254 * t104 + t255 * t105;
t184 = t231 * t284 + t91 * t259 - t58;
t182 = -t59 + (-t154 * t91 - t232) * t284;
t153 = -qJ(4) - pkin(7);
t181 = t159 * t143 + t156 * t153 + t214;
t180 = t157 * pkin(5) + t154 * qJ(6) + pkin(4);
t179 = t156 * t143 - t159 * t153 + t241;
t49 = t254 * t81 + t255 * t82;
t54 = -t92 * pkin(4) - t93 * pkin(8) + t228;
t176 = t154 * t54 + t157 * t49 + t64 * t231 - t68 * t232;
t175 = t284 * t37 + t260;
t174 = 0.2e1 * qJ(2) * t227 + qJDD(3) * t160;
t169 = g(1) * t86 - g(2) * t88 + t154 * t273 + t210;
t168 = -t171 * t20 - t19 * t98 - t44 * t93 + t45 * t92 + t283;
t166 = t18 * t72 + qJDD(6) - t169;
t165 = -t160 * t161 + t282;
t164 = -g(1) * t87 + g(2) * t89 - t136 * t271 + t177;
t163 = t171 * t58 - t98 * t30 - t93 * t70 + (t154 * t92 + t204 * t157) * t284;
t139 = qJDD(3) * t158;
t127 = -t217 - pkin(4);
t96 = -t217 - t180;
t33 = t72 * pkin(5) + t70 * qJ(6);
t31 = t185 * t98 + t67;
t25 = -pkin(5) * t171 + t154 * t68 - t157 * t64;
t24 = qJ(6) * t171 + t262;
t11 = -t94 * pkin(5) + t154 * t53 - t157 * t55;
t10 = t94 * qJ(6) + t263;
t9 = t284 * t70 - t29;
t6 = (pkin(5) * t93 + qJ(6) * t252) * t154 + (-qJ(6) * t93 + (pkin(5) * qJD(5) - qJD(6)) * t98) * t157 + t48;
t5 = t92 * pkin(5) + t262 * qJD(5) + t154 * t49 - t157 * t54;
t4 = -t92 * qJ(6) + qJD(6) * t171 + t176;
t12 = [qJDD(1), t283, t195, qJDD(2) - t283 - 0.2e1 * t253, t282, -t203 * pkin(1) - g(1) * t214 - g(2) * t241 + (t221 + t148) * qJ(2), t152 * qJDD(1) - 0.2e1 * t155 * t215, -0.2e1 * t155 * t223 + 0.2e1 * t239 * t227, -t161 * t155 + t139, -t161 * t158 - t225, 0, t155 * t165 + t158 * t174, -t155 * t174 + t158 * t165, -t170 * t67 + t48 * t94 - t49 * t91 + t68 * t65 + t168, -g(1) * t181 - g(2) * t179 + t103 * t228 + t183 * t246 - t19 * t67 + t20 * t68 - t44 * t48 + t45 * t49, t157 * t190 - t219 * t72 (-t154 * t72 - t157 * t70) * t93 + (t257 - t157 * t30 + (t154 * t70 - t157 * t72) * qJD(5)) * t98, -t157 * t188 - t219 * t284 + t191, t154 * t188 - t218 * t284 + t189, -t171 * t63 - t284 * t92, t210 * t171 - t16 * t92 + t48 * t70 + t67 * t30 + ((-qJD(5) * t68 + t54) * t284 - t64 * t63 + t37 * t252) * t157 + ((-qJD(5) * t64 - t49) * t284 + t68 * t63 + t193) * t154 + t201, t157 * t193 + t17 * t92 - t171 * t177 - t176 * t284 - t219 * t37 + t262 * t63 - t67 * t29 + t48 * t72 + t202, t154 * t198 + t18 * t218 + t25 * t63 - t284 * t5 + t31 * t30 + t6 * t70 - t199 + t201, -t24 * t30 - t25 * t29 - t4 * t70 + t5 * t72 + t192 * t93 + t195 * t136 + (-t1 * t154 + t2 * t157 + (-t154 * t7 - t157 * t8) * qJD(5)) * t98, -t157 * t198 + t18 * t219 - t24 * t63 + t284 * t4 + t31 * t29 - t6 * t72 + t200 - t202, t1 * t24 + t8 * t4 + t3 * t31 + t18 * t6 + t2 * t25 + t7 * t5 - g(1) * (t89 * pkin(5) + t88 * qJ(6) + t159 * t196 + t181) - g(2) * (t87 * pkin(5) + t86 * qJ(6) + t156 * t196 + t179); 0, 0, 0, qJDD(1), -t162, t203 + t173, 0, 0, 0, 0, 0, t238 * t155 + t139, t238 * t158 - t225, t170 * t98 + t171 * t65 + t92 * t91 - t93 * t94, -t168 - t230, 0, 0, 0, 0, 0, t163, -t285, t163 (-t204 * t72 + t189) * t157 + (-t204 * t70 + t191) * t154, t285 (-t204 * t7 + t200) * t157 + (t204 * t8 + t199) * t154 - t198 - t283; 0, 0, 0, 0, 0, 0, t158 * t162 * t155, -t239 * t162, t223, -t224, qJDD(3), t158 * t173 + t100 + t272, g(3) * t158 + (-t111 - t173) * t155 (t45 - t52) * t94 - (-t53 + t44) * t91 + (t170 * t255 + t254 * t65) * pkin(3), t44 * t52 - t45 * t53 + (t254 * t20 + t255 * t19 + t272 + (-t283 - t230) * t158) * pkin(3), t157 * t213 - t257 (-t29 - t268) * t157 - t72 * t207 + t264, t184 - t266, t182 + t265, -t284 * t94, t127 * t30 - t16 * t94 - t52 * t70 + (-t222 - t14 + (-t55 - t233) * t284) * t157 + (t284 * t53 + t175) * t154 + t242, -t127 * t29 + t263 * t284 + t17 * t94 - t52 * t72 + t175 * t157 + (t14 + t280) * t154, t11 * t284 + t96 * t30 + t7 * t94 + t261 * t70 + (-t220 - t3 - t222) * t157 + t281 * t154 + t242, -t273 + t10 * t70 - t11 * t72 - t283 * t135 + (-t126 * t30 + t7 * t91 + t1 + (t126 * t72 + t7) * qJD(5)) * t157 + (-t126 * t29 - t8 * t91 + t2 + (t126 * t70 - t8) * qJD(5)) * t154, -t10 * t284 + t96 * t29 - t8 * t94 - t261 * t72 - t281 * t157 + (-t280 - t3) * t154, t3 * t96 - t8 * t10 - t7 * t11 - g(3) * (-t143 + t270) + t261 * t18 + t180 * t274 + (qJD(5) * t192 + t1 * t157 + t2 * t154) * t126 - t283 * (pkin(3) * t158 + pkin(8) * t135 + t136 * t180); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91 ^ 2 - t94 ^ 2, t44 * t94 + t45 * t91 + t183 - t195, 0, 0, 0, 0, 0, t182 - t265, -t259 * t284 - t266 + t58, -t207 * t284 - t265 - t59 (t29 - t268) * t157 + t154 * t213 + t264, t184 + t266, -t18 * t94 + (-t2 + t288) * t157 + (t284 * t7 + t1) * t154 - t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267, -t70 ^ 2 + t277, t9, -t167 + (-qJD(5) + t284) * t72, -t63, -t37 * t72 + t169 + t269, t16 * t284 + t37 * t70 - t164, -t33 * t70 - t166 + t269 - 0.2e1 * t275, pkin(5) * t29 - t30 * qJ(6) + (-t17 + t8) * t72 + (t7 - t243) * t70, -0.2e1 * t256 - t18 * t70 + t33 * t72 + (0.2e1 * qJD(6) - t16) * t284 + t164, t1 * qJ(6) - t2 * pkin(5) - t18 * t33 - t7 * t17 - g(1) * (-t86 * pkin(5) + t87 * qJ(6)) - g(2) * (t88 * pkin(5) - t89 * qJ(6)) + t243 * t8 + t185 * t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 + t267, t9, -t284 ^ 2 - t277, t166 + t275 - t288;];
tau_reg  = t12;
