% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:31
% EndTime: 2020-01-03 12:09:39
% DurationCPUTime: 2.87s
% Computational Cost: add. (21212->302), mult. (28643->430), div. (0->0), fcn. (19226->10), ass. (0->209)
t291 = -2 * qJD(4);
t214 = sin(pkin(9));
t212 = qJD(1) + qJD(2);
t215 = cos(pkin(9));
t220 = cos(qJ(3));
t217 = sin(qJ(3));
t261 = t212 * t217;
t181 = -t212 * t215 * t220 + t214 * t261;
t183 = (t214 * t220 + t215 * t217) * t212;
t157 = t183 * t181;
t282 = qJDD(3) - t157;
t290 = t214 * t282;
t289 = t215 * t282;
t216 = sin(qJ(5));
t219 = cos(qJ(5));
t143 = t181 * t219 + t183 * t216;
t145 = -t181 * t216 + t183 * t219;
t112 = t145 * t143;
t209 = qJDD(3) + qJDD(5);
t283 = -t112 + t209;
t288 = t216 * t283;
t287 = t219 * t283;
t251 = qJD(3) * t212;
t241 = t220 * t251;
t210 = qJDD(1) + qJDD(2);
t259 = t217 * t210;
t190 = t241 + t259;
t278 = sin(qJ(1));
t279 = cos(qJ(1));
t226 = g(2) * t278 - g(3) * t279;
t196 = -qJD(1) ^ 2 * pkin(1) - t226;
t218 = sin(qJ(2));
t221 = cos(qJ(2));
t227 = -g(2) * t279 - g(3) * t278;
t224 = qJDD(1) * pkin(1) + t227;
t162 = t221 * t196 + t218 * t224;
t281 = t212 ^ 2;
t155 = -pkin(2) * t281 + pkin(7) * t210 + t162;
t264 = t155 * t217;
t118 = qJDD(3) * pkin(3) - qJ(4) * t190 - t264 + (-g(1) + (pkin(3) * t261 + qJ(4) * qJD(3)) * t212) * t220;
t246 = qJ(4) * t261;
t197 = qJD(3) * pkin(3) - t246;
t256 = t220 * t281;
t277 = t217 * g(1);
t119 = -t277 + (-t197 - t246) * qJD(3) + (-pkin(3) * t256 + qJ(4) * t210 + t155) * t220;
t230 = t118 * t215 - t119 * t214 + t183 * t291;
t242 = t217 * t251;
t257 = t220 * t210;
t228 = -t242 + t257;
t159 = t215 * t190 + t214 * t228;
t253 = qJD(3) * t181;
t232 = -t159 - t253;
t286 = t232 * pkin(8) + t230;
t211 = qJD(3) + qJD(5);
t138 = t211 * t143;
t158 = -t190 * t214 + t215 * t228;
t99 = -qJD(5) * t143 + t158 * t216 + t159 * t219;
t284 = -t138 + t99;
t139 = g(1) * t220 + t264;
t140 = t220 * t155 - t277;
t106 = t217 * t139 + t220 * t140;
t63 = t118 * t214 + t119 * t215 + t181 * t291;
t141 = t143 ^ 2;
t142 = t145 ^ 2;
t179 = t181 ^ 2;
t180 = t183 ^ 2;
t208 = t211 ^ 2;
t280 = t220 ^ 2;
t223 = pkin(4) * t282 + t286;
t171 = qJD(3) * pkin(4) - pkin(8) * t183;
t55 = -pkin(4) * t179 + pkin(8) * t158 - qJD(3) * t171 + t63;
t35 = t216 * t55 - t219 * t223;
t272 = t219 * t55;
t36 = t216 * t223 + t272;
t19 = t216 * t36 - t219 * t35;
t276 = t19 * t214;
t275 = t19 * t215;
t161 = -t218 * t196 + t221 * t224;
t154 = -pkin(2) * t210 - pkin(7) * t281 - t161;
t205 = t280 * t281;
t121 = -pkin(3) * t228 - qJ(4) * t205 + t197 * t261 + qJDD(4) + t154;
t69 = -pkin(4) * t158 - pkin(8) * t179 + t171 * t183 + t121;
t274 = t216 * t69;
t42 = t214 * t63 + t215 * t230;
t273 = t217 * t42;
t271 = t219 * t69;
t109 = t112 + t209;
t270 = t109 * t216;
t269 = t109 * t219;
t268 = t121 * t214;
t267 = t121 * t215;
t151 = qJDD(3) + t157;
t266 = t151 * t214;
t265 = t151 * t215;
t263 = t211 * t216;
t262 = t211 * t219;
t201 = t217 * t256;
t260 = t217 * (qJDD(3) + t201);
t258 = t220 * (qJDD(3) - t201);
t254 = -pkin(2) * t154 + pkin(7) * t106;
t252 = qJD(3) * t183;
t250 = qJD(3) * t214;
t249 = qJD(3) * t215;
t213 = t217 ^ 2;
t204 = t213 * t281;
t222 = qJD(3) ^ 2;
t170 = -t258 - t217 * (-t204 - t222);
t189 = 0.2e1 * t241 + t259;
t245 = -pkin(2) * t189 + pkin(7) * t170 + t154 * t217;
t169 = t220 * (-t205 - t222) - t260;
t191 = -0.2e1 * t242 + t257;
t244 = pkin(2) * t191 + pkin(7) * t169 - t154 * t220;
t20 = t216 * t35 + t219 * t36;
t15 = -pkin(4) * t69 + pkin(8) * t20;
t7 = t20 * t214 + t275;
t8 = t20 * t215 - t276;
t4 = -t217 * t7 + t220 * t8;
t243 = t217 * (-pkin(8) * t275 - qJ(4) * t7 - t15 * t214) + t220 * (-pkin(3) * t69 - pkin(8) * t276 + qJ(4) * t8 + t15 * t215) - pkin(2) * t69 + pkin(7) * t4;
t233 = -t158 * t219 + t159 * t216;
t225 = (-qJD(5) + t211) * t145 - t233;
t82 = t138 + t99;
t48 = t216 * t225 - t219 * t82;
t10 = -pkin(8) * t48 - t19;
t50 = t216 * t82 + t219 * t225;
t27 = t214 * t50 + t215 * t48;
t28 = -t214 * t48 + t215 * t50;
t14 = -t217 * t27 + t220 * t28;
t97 = -t141 - t142;
t9 = -pkin(4) * t97 + pkin(8) * t50 + t20;
t240 = t217 * (-qJ(4) * t27 + t10 * t215 - t214 * t9) + t220 * (-pkin(3) * t97 + qJ(4) * t28 + t10 * t214 + t215 * t9) - pkin(2) * t97 + pkin(7) * t14;
t43 = -t214 * t230 + t215 * t63;
t107 = -t208 - t141;
t64 = t107 * t216 + t287;
t65 = t107 * t219 - t288;
t44 = t214 * t65 + t215 * t64;
t45 = -t214 * t64 + t215 * t65;
t24 = -t217 * t44 + t220 * t45;
t77 = (qJD(5) + t211) * t145 + t233;
t37 = -pkin(4) * t77 + pkin(8) * t65 - t271;
t46 = -pkin(8) * t64 + t274;
t239 = t217 * (-qJ(4) * t44 - t214 * t37 + t215 * t46) + t220 * (-pkin(3) * t77 + qJ(4) * t45 + t214 * t46 + t215 * t37) - pkin(2) * t77 + pkin(7) * t24;
t132 = -t142 - t208;
t84 = t132 * t219 - t270;
t85 = -t132 * t216 - t269;
t52 = t214 * t85 + t215 * t84;
t53 = -t214 * t84 + t215 * t85;
t30 = -t217 * t52 + t220 * t53;
t39 = -pkin(4) * t284 + pkin(8) * t85 + t274;
t51 = -pkin(8) * t84 + t271;
t238 = t217 * (-qJ(4) * t52 - t214 * t39 + t215 * t51) + t220 * (-pkin(3) * t284 + qJ(4) * t53 + t214 * t51 + t215 * t39) - pkin(2) * t284 + pkin(7) * t30;
t129 = t158 + t252;
t100 = t129 * t214 + t215 * t232;
t101 = t129 * t215 - t214 * t232;
t127 = -t179 - t180;
t58 = -t100 * t217 + t101 * t220;
t237 = t217 * (-qJ(4) * t100 - t42) + t220 * (-pkin(3) * t127 + qJ(4) * t101 + t43) - pkin(2) * t127 + pkin(7) * t58;
t148 = -t222 - t179;
t113 = t148 * t214 + t289;
t114 = t148 * t215 - t290;
t128 = -t158 + t252;
t68 = -t113 * t217 + t114 * t220;
t236 = t217 * (-qJ(4) * t113 + t268) + t220 * (-pkin(3) * t128 + qJ(4) * t114 - t267) - pkin(2) * t128 + pkin(7) * t68;
t174 = -t180 - t222;
t122 = t174 * t215 - t266;
t123 = -t174 * t214 - t265;
t130 = t159 - t253;
t93 = -t122 * t217 + t123 * t220;
t235 = t217 * (-qJ(4) * t122 + t267) + t220 * (-pkin(3) * t130 + qJ(4) * t123 + t268) - pkin(2) * t130 + pkin(7) * t93;
t194 = (t213 + t280) * t210;
t195 = t204 + t205;
t231 = pkin(2) * t195 + pkin(7) * t194 + t106;
t22 = t220 * t43 - t273;
t229 = pkin(7) * t22 - qJ(4) * t273 - pkin(2) * t121 + t220 * (-pkin(3) * t121 + qJ(4) * t43);
t173 = -t180 + t222;
t172 = t179 - t222;
t168 = t260 + t220 * (-t204 + t222);
t167 = t217 * (t205 - t222) + t258;
t164 = (t190 + t241) * t217;
t163 = t191 * t220;
t160 = t189 * t220 + t191 * t217;
t136 = -t142 + t208;
t135 = t141 - t208;
t111 = t142 - t141;
t104 = (t217 * (-t181 * t215 + t183 * t214) + t220 * (-t181 * t214 - t183 * t215)) * qJD(3);
t103 = (-t143 * t219 + t145 * t216) * t211;
t102 = (-t143 * t216 - t145 * t219) * t211;
t98 = -qJD(5) * t145 - t233;
t96 = t217 * (t159 * t215 - t183 * t250) + t220 * (t159 * t214 + t183 * t249);
t95 = t217 * (-t158 * t214 + t181 * t249) + t220 * (t158 * t215 + t181 * t250);
t92 = t217 * (-t173 * t214 + t289) + t220 * (t173 * t215 + t290);
t91 = t217 * (t172 * t215 - t266) + t220 * (t172 * t214 + t265);
t89 = t135 * t219 - t270;
t88 = -t136 * t216 + t287;
t87 = t135 * t216 + t269;
t86 = t136 * t219 + t288;
t73 = -t145 * t263 + t219 * t99;
t72 = t145 * t262 + t216 * t99;
t71 = t143 * t262 - t216 * t98;
t70 = t143 * t263 + t219 * t98;
t57 = t217 * (-t128 * t215 - t130 * t214) + t220 * (-t128 * t214 + t130 * t215);
t49 = -t216 * t284 - t219 * t77;
t47 = -t216 * t77 + t219 * t284;
t41 = t217 * (-t102 * t214 + t103 * t215) + t220 * (t102 * t215 + t103 * t214);
t32 = t217 * (-t214 * t87 + t215 * t89) + t220 * (t214 * t89 + t215 * t87);
t31 = t217 * (-t214 * t86 + t215 * t88) + t220 * (t214 * t88 + t215 * t86);
t26 = t217 * (-t214 * t72 + t215 * t73) + t220 * (t214 * t73 + t215 * t72);
t25 = t217 * (-t214 * t70 + t215 * t71) + t220 * (t214 * t71 + t215 * t70);
t13 = t217 * (-t214 * t47 + t215 * t49) + t220 * (t214 * t49 + t215 * t47);
t1 = [0, 0, 0, 0, 0, qJDD(1), t227, t226, 0, 0, 0, 0, 0, 0, 0, t210, pkin(1) * (t210 * t221 - t218 * t281) + t161, pkin(1) * (-t210 * t218 - t221 * t281) - t162, 0, pkin(1) * (t161 * t221 + t162 * t218), t164, t160, t168, t163, t167, 0, pkin(1) * (t169 * t218 + t191 * t221) + t244, pkin(1) * (t170 * t218 - t189 * t221) + t245, pkin(1) * (t194 * t218 + t195 * t221) + t231, pkin(1) * (t106 * t218 - t154 * t221) + t254, t96, t57, t92, t95, t91, t104, pkin(1) * (-t128 * t221 + t218 * t68) + t236, pkin(1) * (-t130 * t221 + t218 * t93) + t235, pkin(1) * (-t127 * t221 + t218 * t58) + t237, pkin(1) * (-t121 * t221 + t218 * t22) + t229, t26, t13, t31, t25, t32, t41, pkin(1) * (t218 * t24 - t221 * t77) + t239, pkin(1) * (t218 * t30 - t221 * t284) + t238, pkin(1) * (t14 * t218 - t221 * t97) + t240, pkin(1) * (t218 * t4 - t221 * t69) + t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, t161, -t162, 0, 0, t164, t160, t168, t163, t167, 0, t244, t245, t231, t254, t96, t57, t92, t95, t91, t104, t236, t235, t237, t229, t26, t13, t31, t25, t32, t41, t239, t238, t240, t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201, t204 - t205, t259, t201, t257, qJDD(3), -t139, -t140, 0, 0, t157, t180 - t179, -t232, -t157, t129, qJDD(3), pkin(3) * t113 + t230, pkin(3) * t122 - t63, pkin(3) * t100, pkin(3) * t42, t112, t111, t82, -t112, t225, t209, pkin(3) * t44 + pkin(4) * t64 - t35, -t272 - t216 * t286 + pkin(3) * t52 + (-t216 * t282 + t84) * pkin(4), pkin(3) * t27 + pkin(4) * t48, pkin(3) * t7 + pkin(4) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t130, t127, t121, 0, 0, 0, 0, 0, 0, t77, t284, t97, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t111, t82, -t112, t225, t209, -t35, -t36, 0, 0;];
tauJ_reg = t1;
