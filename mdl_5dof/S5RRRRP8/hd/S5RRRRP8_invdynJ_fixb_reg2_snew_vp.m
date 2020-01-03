% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:10
% EndTime: 2019-12-31 22:02:20
% DurationCPUTime: 3.88s
% Computational Cost: add. (15466->357), mult. (31353->450), div. (0->0), fcn. (21888->8), ass. (0->236)
t215 = sin(qJ(3));
t219 = cos(qJ(3));
t216 = sin(qJ(2));
t262 = qJD(1) * t216;
t187 = t219 * qJD(2) - t215 * t262;
t188 = t215 * qJD(2) + t219 * t262;
t214 = sin(qJ(4));
t218 = cos(qJ(4));
t166 = -t218 * t187 + t188 * t214;
t168 = t187 * t214 + t188 * t218;
t128 = t168 * t166;
t256 = qJD(1) * qJD(2);
t204 = t216 * t256;
t220 = cos(qJ(2));
t255 = t220 * qJDD(1);
t192 = -t204 + t255;
t186 = -qJDD(3) + t192;
t182 = -qJDD(4) + t186;
t309 = t128 + t182;
t313 = t309 * pkin(4);
t205 = t216 * qJDD(1);
t246 = t220 * t256;
t191 = t205 + t246;
t230 = -t215 * qJDD(2) - t219 * t191;
t160 = qJD(3) * t187 - t230;
t242 = -t219 * qJDD(2) + t215 * t191;
t229 = -qJD(3) * t188 - t242;
t109 = -t166 * qJD(4) + t218 * t160 + t214 * t229;
t202 = t220 * qJD(1) - qJD(3);
t197 = -qJD(4) + t202;
t148 = t166 * t197;
t92 = -t148 + t109;
t310 = qJ(5) * t92;
t222 = qJD(1) ^ 2;
t217 = sin(qJ(1));
t221 = cos(qJ(1));
t245 = t217 * g(1) - t221 * g(2);
t180 = qJDD(1) * pkin(1) + t222 * pkin(6) + t245;
t234 = -t192 + t204;
t235 = t191 + t246;
t131 = t234 * pkin(2) - t235 * pkin(7) - t180;
t237 = t221 * g(1) + t217 * g(2);
t279 = qJDD(1) * pkin(6);
t181 = -t222 * pkin(1) - t237 + t279;
t238 = -t220 * pkin(2) - t216 * pkin(7);
t241 = t222 * t238 + t181;
t293 = t216 * g(3);
t303 = qJD(2) ^ 2;
t143 = -t303 * pkin(2) + qJDD(2) * pkin(7) + t241 * t220 - t293;
t111 = -t219 * t131 + t215 * t143;
t177 = t187 * t202;
t136 = t160 + t177;
t270 = t187 * t188;
t226 = -t186 + t270;
t72 = t226 * pkin(3) - t136 * pkin(8) - t111;
t112 = t215 * t131 + t219 * t143;
t174 = -pkin(3) * t202 - pkin(8) * t188;
t184 = t187 ^ 2;
t75 = -t184 * pkin(3) + pkin(8) * t229 + t202 * t174 + t112;
t42 = t214 * t75 - t218 * t72;
t225 = 0.2e1 * qJD(5) * t168 + t310 + t313 + t42;
t224 = -t225 - t313;
t164 = t166 ^ 2;
t196 = t197 ^ 2;
t123 = -t196 - t164;
t275 = t309 * t218;
t79 = t123 * t214 - t275;
t78 = pkin(3) * t79;
t312 = t224 + t78;
t165 = t168 ^ 2;
t140 = -t165 - t196;
t244 = t214 * t160 - t218 * t229;
t108 = -qJD(4) * t168 - t244;
t144 = -pkin(4) * t197 - qJ(5) * t168;
t43 = t214 * t72 + t218 * t75;
t27 = -t164 * pkin(4) + t108 * qJ(5) - 0.2e1 * qJD(5) * t166 + t197 * t144 + t43;
t228 = pkin(4) * t140 - t27;
t118 = -t128 + t182;
t278 = t118 * t214;
t96 = t140 * t218 + t278;
t95 = pkin(3) * t96;
t311 = t95 + t228;
t292 = t220 * g(3);
t142 = -qJDD(2) * pkin(2) - t303 * pkin(7) + t241 * t216 + t292;
t98 = -t229 * pkin(3) - t184 * pkin(8) + t188 * t174 + t142;
t276 = t309 * t214;
t46 = -t108 * pkin(4) - t164 * qJ(5) + t168 * t144 + qJDD(5) + t98;
t308 = t215 * t226;
t307 = t219 * t226;
t305 = t148 + t109;
t89 = (qJD(4) + t197) * t168 + t244;
t132 = (qJD(3) + t202) * t188 + t242;
t185 = t188 ^ 2;
t199 = t202 ^ 2;
t80 = t123 * t218 + t276;
t49 = t215 * t80 + t219 * t79;
t302 = pkin(2) * t49;
t277 = t118 * t218;
t97 = -t140 * t214 + t277;
t62 = t215 * t97 + t219 * t96;
t301 = pkin(2) * t62;
t20 = t214 * t43 - t218 * t42;
t300 = pkin(3) * t20;
t56 = -t214 * t89 - t218 * t92;
t58 = t214 * t92 - t218 * t89;
t31 = t215 * t58 + t219 * t56;
t299 = pkin(7) * t31;
t298 = pkin(7) * t49;
t297 = pkin(7) * t62;
t296 = pkin(8) * t56;
t295 = pkin(8) * t79;
t294 = pkin(8) * t96;
t114 = -t164 - t165;
t32 = -t215 * t56 + t219 * t58;
t291 = pkin(6) * (t114 * t216 + t220 * t32) - pkin(1) * t31;
t50 = -t215 * t79 + t219 * t80;
t88 = (qJD(4) - t197) * t168 + t244;
t290 = pkin(6) * (t216 * t88 + t220 * t50) - pkin(1) * t49;
t63 = -t215 * t96 + t219 * t97;
t289 = pkin(6) * (t216 * t305 + t220 * t63) - pkin(1) * t62;
t288 = -pkin(2) * t88 + pkin(7) * t50;
t287 = -pkin(2) * t305 + pkin(7) * t63;
t286 = t20 * t215;
t285 = t20 * t219;
t284 = t214 * t225;
t283 = t214 * t98;
t282 = t218 * t225;
t281 = t218 * t98;
t280 = -pkin(2) * t114 + pkin(7) * t32;
t274 = t142 * t215;
t273 = t142 * t219;
t151 = t186 + t270;
t272 = t151 * t215;
t271 = t151 * t219;
t268 = t197 * t214;
t267 = t197 * t218;
t266 = t202 * t215;
t265 = t202 * t219;
t201 = t220 * t222 * t216;
t264 = t216 * (qJDD(2) + t201);
t263 = t220 * (-t201 + qJDD(2));
t259 = -qJD(3) + t202;
t254 = t95 - t43;
t25 = pkin(4) * t225;
t9 = t214 * t27 - t282;
t253 = pkin(3) * t9 - t25;
t252 = t220 * t128;
t251 = t220 * t270;
t54 = pkin(3) * t56;
t250 = -pkin(2) * t31 - t54;
t249 = -pkin(3) * t88 + pkin(8) * t80;
t248 = -pkin(3) * t305 + pkin(8) * t97;
t247 = -pkin(3) * t114 + pkin(8) * t58;
t21 = t214 * t42 + t218 * t43;
t69 = t111 * t215 + t219 * t112;
t172 = t181 * t216 + t292;
t173 = t181 * t220 - t293;
t243 = t216 * t172 + t220 * t173;
t240 = t42 - t78;
t233 = t111 * t219 - t112 * t215;
t231 = -pkin(1) + t238;
t211 = t220 ^ 2;
t210 = t216 ^ 2;
t208 = t211 * t222;
t207 = t210 * t222;
t193 = -0.2e1 * t204 + t255;
t190 = t205 + 0.2e1 * t246;
t176 = -t185 + t199;
t175 = t184 - t199;
t170 = t185 - t184;
t169 = -t185 - t199;
t162 = -t199 - t184;
t150 = t184 + t185;
t146 = -t165 + t196;
t145 = t164 - t196;
t137 = t187 * t259 + t230;
t135 = t160 - t177;
t133 = t188 * t259 - t242;
t126 = t165 - t164;
t125 = -t169 * t215 + t271;
t124 = t169 * t219 + t272;
t122 = t162 * t219 - t308;
t121 = t162 * t215 + t307;
t116 = (t166 * t218 - t168 * t214) * t197;
t115 = (t166 * t214 + t168 * t218) * t197;
t104 = -t132 * t219 + t136 * t215;
t102 = t145 * t218 + t278;
t101 = -t146 * t214 - t275;
t100 = t145 * t214 - t277;
t99 = t146 * t218 - t276;
t85 = pkin(4) * t92;
t84 = t109 * t218 + t168 * t268;
t83 = t109 * t214 - t168 * t267;
t82 = -t108 * t214 - t166 * t267;
t81 = t108 * t218 - t166 * t268;
t76 = t115 * t219 + t116 * t215;
t73 = t216 * (-t115 * t215 + t116 * t219) + t220 * t182;
t67 = -pkin(4) * t305 + qJ(5) * t118;
t66 = t100 * t219 + t102 * t215;
t65 = t101 * t215 + t219 * t99;
t64 = t281 - t294;
t59 = t283 - t295;
t57 = -t214 * t305 - t218 * t88;
t55 = -t214 * t88 + t218 * t305;
t52 = t215 * t84 + t219 * t83;
t51 = t215 * t82 + t219 * t81;
t45 = t216 * (-t215 * t83 + t219 * t84) - t252;
t44 = t216 * (-t215 * t81 + t219 * t82) + t252;
t40 = -qJ(5) * t140 + t46;
t39 = t216 * (-t100 * t215 + t102 * t219) + t220 * t89;
t38 = t216 * (t101 * t219 - t215 * t99) - t220 * t92;
t37 = t248 + t283;
t35 = t249 - t281;
t33 = -pkin(4) * t88 + qJ(5) * t123 - t46;
t30 = t215 * t57 + t219 * t55;
t24 = t216 * (-t215 * t55 + t219 * t57) - t220 * t126;
t22 = -t214 * t67 + t218 * t40 - t294;
t19 = t225 + t310;
t18 = qJ(5) * t275 - t214 * t33 - t295;
t17 = -pkin(4) * t114 - qJ(5) * t89 + t27;
t16 = -pkin(3) * t98 + pkin(8) * t21;
t15 = t214 * t40 + t218 * t67 + t248;
t14 = qJ(5) * t276 + t218 * t33 + t249;
t13 = -pkin(4) * t46 + qJ(5) * t27;
t12 = -t20 - t296;
t11 = t21 + t247;
t10 = t218 * t27 + t284;
t8 = t21 * t219 - t286;
t7 = t21 * t215 + t285;
t6 = -t17 * t214 + t19 * t218 - t296;
t5 = t17 * t218 + t19 * t214 + t247;
t4 = t10 * t219 - t215 * t9;
t3 = t10 * t215 + t219 * t9;
t2 = -pkin(8) * t9 + qJ(5) * t282 - t13 * t214;
t1 = -pkin(3) * t46 + pkin(8) * t10 + qJ(5) * t284 + t13 * t218;
t23 = [0, 0, 0, 0, 0, qJDD(1), t245, t237, 0, 0, t235 * t216, t190 * t220 + t193 * t216, t264 + t220 * (-t207 + t303), -t234 * t220, t216 * (t208 - t303) + t263, 0, t220 * t180 + pkin(1) * t193 + pkin(6) * (t220 * (-t208 - t303) - t264), -t216 * t180 - pkin(1) * t190 + pkin(6) * (-t263 - t216 * (-t207 - t303)), pkin(1) * (t207 + t208) + (t210 + t211) * t279 + t243, pkin(1) * t180 + pkin(6) * t243, t216 * (t160 * t219 + t188 * t266) + t251, t216 * (t133 * t219 - t135 * t215) - t220 * t170, t216 * (-t176 * t215 + t307) - t220 * t136, t216 * (t187 * t265 - t215 * t229) - t251, t216 * (t175 * t219 + t272) + t220 * t132, t220 * t186 + t216 * (-t187 * t219 - t188 * t215) * t202, t216 * (-pkin(7) * t121 + t274) + t220 * (-pkin(2) * t121 + t111) - pkin(1) * t121 + pkin(6) * (t122 * t220 - t133 * t216), t216 * (-pkin(7) * t124 + t273) + t220 * (-pkin(2) * t124 + t112) - pkin(1) * t124 + pkin(6) * (t125 * t220 - t137 * t216), t216 * t233 + pkin(6) * (t104 * t220 - t150 * t216) + t231 * (-t132 * t215 - t136 * t219), pkin(6) * (t142 * t216 + t220 * t69) - t231 * t233, t45, t24, t38, t44, t39, t73, t216 * (-t215 * t35 + t219 * t59 - t298) + t220 * (t240 - t302) + t290, t216 * (-t215 * t37 + t219 * t64 - t297) + t220 * (-t254 - t301) + t289, t216 * (-t11 * t215 + t12 * t219 - t299) + t220 * t250 + t291, t216 * (-pkin(7) * t7 - pkin(8) * t285 - t16 * t215) + t220 * (-pkin(2) * t7 - t300) - pkin(1) * t7 + pkin(6) * (t216 * t98 + t220 * t8), t45, t24, t38, t44, t39, t73, t216 * (-t14 * t215 + t18 * t219 - t298) + t220 * (-t302 - t312) + t290, t216 * (-t15 * t215 + t219 * t22 - t297) + t220 * (-t301 - t311) + t289, t216 * (-t215 * t5 + t219 * t6 - t299) + t220 * (t250 + t85) + t291, t216 * (-pkin(7) * t3 - t1 * t215 + t2 * t219) + t220 * (-pkin(2) * t3 - t253) - pkin(1) * t3 + pkin(6) * (t216 * t46 + t220 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201, t207 - t208, t205, t201, t255, qJDD(2), -t172, -t173, 0, 0, t160 * t215 - t188 * t265, t133 * t215 + t135 * t219, t176 * t219 + t308, t187 * t266 + t219 * t229, t175 * t215 - t271, (-t187 * t215 + t188 * t219) * t202, pkin(2) * t133 + pkin(7) * t122 - t273, pkin(2) * t137 + pkin(7) * t125 + t274, pkin(2) * t150 + pkin(7) * t104 + t69, -pkin(2) * t142 + pkin(7) * t69, t52, t30, t65, t51, t66, t76, t215 * t59 + t219 * t35 + t288, t215 * t64 + t219 * t37 + t287, t11 * t219 + t12 * t215 + t280, -pkin(2) * t98 + pkin(7) * t8 - pkin(8) * t286 + t16 * t219, t52, t30, t65, t51, t66, t76, t14 * t219 + t18 * t215 + t288, t15 * t219 + t215 * t22 + t287, t215 * t6 + t219 * t5 + t280, -pkin(2) * t46 + pkin(7) * t4 + t1 * t219 + t2 * t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t270, t170, t136, t270, -t132, -t186, -t111, -t112, 0, 0, t128, t126, t92, -t128, -t89, -t182, -t240, t254, t54, t300, t128, t126, t92, -t128, -t89, -t182, t312, t311, -t85 + t54, t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t126, t92, -t128, -t89, -t182, -t42, -t43, 0, 0, t128, t126, t92, -t128, -t89, -t182, t224, t228, -t85, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t305, t114, t46;];
tauJ_reg = t23;
