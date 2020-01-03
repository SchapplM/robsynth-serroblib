% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRRR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:40
% EndTime: 2019-12-31 16:31:43
% DurationCPUTime: 2.71s
% Computational Cost: add. (6350->267), mult. (10213->387), div. (0->0), fcn. (7304->8), ass. (0->184)
t259 = qJD(2) + qJD(3);
t257 = t259 ^ 2;
t269 = cos(qJ(3));
t258 = qJDD(2) + qJDD(3);
t266 = sin(qJ(3));
t293 = t266 * t258;
t226 = t269 * t257 + t293;
t288 = t269 * t258;
t229 = t266 * t257 - t288;
t267 = sin(qJ(2));
t270 = cos(qJ(2));
t180 = t270 * t226 - t267 * t229;
t262 = g(3) - qJDD(1);
t213 = pkin(5) * t226 - t269 * t262;
t319 = pkin(5) * t229 - t266 * t262;
t137 = pkin(4) * t180 + t270 * t213 - t267 * t319;
t185 = t267 * t226 + t270 * t229;
t263 = sin(pkin(7));
t264 = cos(pkin(7));
t147 = t263 * t180 + t264 * t185;
t325 = pkin(4) * t185 + t267 * t213 + t270 * t319;
t332 = qJ(1) * t147 + t263 * t137 + t264 * t325;
t318 = t264 * t180 - t263 * t185;
t331 = qJ(1) * t318 + t264 * t137 - t263 * t325;
t241 = t263 * g(1) - t264 * g(2);
t242 = t264 * g(1) + t263 * g(2);
t285 = -t267 * t241 + t270 * t242;
t306 = qJD(2) ^ 2;
t191 = -t306 * pkin(2) - t285;
t273 = t270 * t241 + t267 * t242;
t272 = qJDD(2) * pkin(2) + t273;
t156 = t266 * t191 - t269 * t272;
t157 = t269 * t191 + t266 * t272;
t280 = t266 * t156 + t269 * t157;
t116 = t269 * t156 - t266 * t157;
t291 = t267 * t116;
t322 = t270 * t280 + t291;
t286 = t270 * t116;
t98 = -t267 * t280 + t286;
t87 = t263 * t98 + t264 * t322;
t328 = -t263 * t322 + t264 * t98;
t279 = -t267 * t273 - t270 * t285;
t162 = t267 * t285 - t270 * t273;
t299 = t264 * t162;
t321 = -t263 * t279 + t299;
t303 = t263 * t162;
t119 = t264 * t279 + t303;
t239 = t267 * qJDD(2) + t270 * t306;
t240 = t270 * qJDD(2) - t267 * t306;
t193 = -t263 * t239 + t264 * t240;
t217 = pkin(4) * t239 - t270 * t262;
t274 = -pkin(4) * t240 - t267 * t262;
t320 = -qJ(1) * t193 + t263 * t217 + t264 * t274;
t308 = t264 * t239 + t263 * t240;
t316 = qJ(1) * t308 + t264 * t217 - t263 * t274;
t305 = pkin(1) * t262;
t265 = sin(qJ(4));
t260 = t265 ^ 2;
t304 = t260 * t257;
t300 = t263 * t262;
t298 = t264 * t262;
t152 = -t258 * pkin(3) - t257 * pkin(6) + t156;
t297 = t265 * t152;
t268 = cos(qJ(4));
t247 = t265 * t257 * t268;
t237 = qJDD(4) + t247;
t296 = t265 * t237;
t238 = qJDD(4) - t247;
t295 = t265 * t238;
t294 = t265 * t258;
t290 = t268 * t152;
t289 = t268 * t238;
t249 = t268 * t258;
t153 = -t257 * pkin(3) + t258 * pkin(6) + t157;
t144 = t268 * t153 - t265 * t262;
t261 = t268 ^ 2;
t284 = t260 + t261;
t283 = qJD(4) * t259;
t282 = t265 * t283;
t281 = t268 * t283;
t143 = t265 * t153 + t268 * t262;
t110 = t265 * t143 + t268 * t144;
t205 = -t263 * t241 - t264 * t242;
t276 = t266 * t247;
t275 = t269 * t247;
t109 = t268 * t143 - t265 * t144;
t204 = t264 * t241 - t263 * t242;
t271 = qJD(4) ^ 2;
t250 = t261 * t257;
t246 = -t250 - t271;
t245 = t250 - t271;
t244 = -t271 - t304;
t243 = t271 - t304;
t231 = t250 - t304;
t230 = t250 + t304;
t225 = t268 * t237;
t224 = t284 * t258;
t223 = t249 - 0.2e1 * t282;
t222 = t249 - t282;
t221 = t281 + t294;
t220 = 0.2e1 * t281 + t294;
t219 = t284 * t283;
t209 = t266 * qJDD(4) + t269 * t219;
t208 = -t269 * qJDD(4) + t266 * t219;
t203 = t268 * t221 - t260 * t283;
t202 = -t265 * t222 - t261 * t283;
t201 = -t265 * t244 - t289;
t200 = -t265 * t243 + t225;
t199 = t268 * t246 - t296;
t198 = t268 * t245 - t295;
t195 = t268 * t244 - t295;
t194 = t265 * t246 + t225;
t183 = t269 * t224 - t266 * t230;
t179 = t266 * t224 + t269 * t230;
t178 = -t265 * t220 + t268 * t223;
t177 = t269 * t200 + t265 * t293;
t176 = t269 * t198 + t266 * t249;
t175 = t266 * t200 - t265 * t288;
t174 = t266 * t198 - t268 * t288;
t173 = t269 * t203 - t276;
t172 = t269 * t202 + t276;
t171 = t266 * t203 + t275;
t170 = t266 * t202 - t275;
t169 = t269 * t201 + t266 * t220;
t168 = t269 * t199 - t266 * t223;
t167 = t266 * t201 - t269 * t220;
t166 = t266 * t199 + t269 * t223;
t165 = -t267 * t208 + t270 * t209;
t164 = t270 * t208 + t267 * t209;
t161 = t269 * t178 - t266 * t231;
t158 = t266 * t178 + t269 * t231;
t155 = pkin(4) * t279 + t305;
t151 = -t267 * t179 + t270 * t183;
t150 = t270 * t179 + t267 * t183;
t141 = -t267 * t175 + t270 * t177;
t140 = -t267 * t174 + t270 * t176;
t139 = t270 * t175 + t267 * t177;
t138 = t270 * t174 + t267 * t176;
t133 = -t267 * t171 + t270 * t173;
t132 = -t267 * t170 + t270 * t172;
t131 = t270 * t171 + t267 * t173;
t130 = t270 * t170 + t267 * t172;
t129 = -t267 * t167 + t270 * t169;
t128 = -t267 * t166 + t270 * t168;
t127 = t270 * t167 + t267 * t169;
t126 = t270 * t166 + t267 * t168;
t125 = -pkin(6) * t195 + t290;
t124 = -pkin(6) * t194 + t297;
t123 = -pkin(3) * t195 + t144;
t122 = -pkin(3) * t194 + t143;
t121 = -t267 * t158 + t270 * t161;
t120 = t270 * t158 + t267 * t161;
t113 = pkin(2) * t262 + pkin(5) * t280;
t112 = -t263 * t150 + t264 * t151;
t111 = t264 * t150 + t263 * t151;
t107 = -pkin(5) * t179 + t269 * t109;
t106 = pkin(5) * t183 + t266 * t109;
t105 = -t263 * t127 + t264 * t129;
t104 = -t263 * t126 + t264 * t128;
t103 = t264 * t127 + t263 * t129;
t102 = t264 * t126 + t263 * t128;
t101 = t269 * t110 + t266 * t152;
t100 = t266 * t110 - t269 * t152;
t95 = -pkin(5) * t167 - t266 * t123 + t269 * t125;
t94 = -pkin(5) * t166 - t266 * t122 + t269 * t124;
t93 = -pkin(2) * t195 + pkin(5) * t169 + t269 * t123 + t266 * t125;
t92 = -pkin(2) * t194 + pkin(5) * t168 + t269 * t122 + t266 * t124;
t91 = -pkin(4) * t150 - t267 * t106 + t270 * t107;
t90 = pkin(4) * t151 + t270 * t106 + t267 * t107;
t89 = -t267 * t100 + t270 * t101;
t88 = t270 * t100 + t267 * t101;
t85 = pkin(4) * t98 + pkin(5) * t286 - t267 * t113;
t84 = pkin(4) * t322 + pkin(5) * t291 + t270 * t113 + t305;
t83 = -pkin(5) * t100 - (pkin(3) * t266 - pkin(6) * t269) * t109;
t82 = -pkin(4) * t127 - t267 * t93 + t270 * t95;
t81 = -pkin(4) * t126 - t267 * t92 + t270 * t94;
t80 = -pkin(1) * t195 + pkin(4) * t129 + t267 * t95 + t270 * t93;
t79 = -pkin(1) * t194 + pkin(4) * t128 + t267 * t94 + t270 * t92;
t78 = pkin(5) * t101 - (-pkin(3) * t269 - pkin(6) * t266 - pkin(2)) * t109;
t77 = -t263 * t88 + t264 * t89;
t76 = t263 * t89 + t264 * t88;
t75 = -pkin(4) * t88 - t267 * t78 + t270 * t83;
t74 = pkin(1) * t109 + pkin(4) * t89 + t267 * t83 + t270 * t78;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, 0, 0, 0, 0, 0, 0, -t308, -t193, 0, t119, 0, 0, 0, 0, 0, 0, -t318, t147, 0, t87, 0, 0, 0, 0, 0, 0, t104, t105, t112, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, 0, 0, 0, 0, 0, 0, t193, -t308, 0, -t321, 0, 0, 0, 0, 0, 0, -t147, -t318, 0, -t328, 0, 0, 0, 0, 0, 0, t102, t103, t111, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, 0, 0, 0, 0, 0, 0, t194, t195, 0, -t109; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t300, -t298, -t204, -qJ(1) * t204, 0, 0, t193, 0, -t308, 0, t320, t316, t321, pkin(4) * t299 + qJ(1) * t321 - t263 * t155, 0, 0, -t147, 0, -t318, 0, t332, t331, t328, qJ(1) * t328 - t263 * t84 + t264 * t85, -t263 * t131 + t264 * t133, -t263 * t120 + t264 * t121, -t263 * t139 + t264 * t141, -t263 * t130 + t264 * t132, -t263 * t138 + t264 * t140, -t263 * t164 + t264 * t165, -qJ(1) * t102 - t263 * t79 + t264 * t81, -qJ(1) * t103 - t263 * t80 + t264 * t82, -qJ(1) * t111 - t263 * t90 + t264 * t91, -qJ(1) * t76 - t263 * t74 + t264 * t75; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t298, -t300, t205, qJ(1) * t205, 0, 0, t308, 0, t193, 0, -t316, t320, t119, pkin(4) * t303 + qJ(1) * t119 + t264 * t155, 0, 0, t318, 0, -t147, 0, -t331, t332, t87, qJ(1) * t87 + t263 * t85 + t264 * t84, t264 * t131 + t263 * t133, t264 * t120 + t263 * t121, t264 * t139 + t263 * t141, t264 * t130 + t263 * t132, t264 * t138 + t263 * t140, t264 * t164 + t263 * t165, qJ(1) * t104 + t263 * t81 + t264 * t79, qJ(1) * t105 + t263 * t82 + t264 * t80, qJ(1) * t112 + t263 * t91 + t264 * t90, qJ(1) * t77 + t263 * t75 + t264 * t74; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t241, t242, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(1) * t240 + t273, -pkin(1) * t239 + t285, 0, -pkin(1) * t162, 0, 0, 0, 0, 0, t258, -pkin(1) * t185 - pkin(2) * t229 - t156, -pkin(1) * t180 - pkin(2) * t226 - t157, 0, -pkin(1) * t98 - pkin(2) * t116, (t221 + t281) * t265, t268 * t220 + t265 * t223, t268 * t243 + t296, (t222 - t282) * t268, t265 * t245 + t289, 0, pkin(1) * t126 + pkin(2) * t166 + pkin(3) * t223 + pkin(6) * t199 - t290, pkin(1) * t127 + pkin(2) * t167 - pkin(3) * t220 + pkin(6) * t201 + t297, pkin(1) * t150 + pkin(2) * t179 + pkin(3) * t230 + pkin(6) * t224 + t110, pkin(1) * t88 + pkin(2) * t100 - pkin(3) * t152 + pkin(6) * t110;];
tauB_reg = t1;
