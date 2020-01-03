% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPPR7
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPPR7_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:43
% EndTime: 2019-12-31 16:41:45
% DurationCPUTime: 2.09s
% Computational Cost: add. (3898->247), mult. (8658->361), div. (0->0), fcn. (5390->6), ass. (0->167)
t249 = sin(qJ(4));
t247 = sin(pkin(6));
t248 = cos(pkin(6));
t251 = cos(qJ(4));
t261 = t247 * t251 + t248 * t249;
t222 = t261 * qJD(1);
t224 = (-t247 * t249 + t248 * t251) * qJD(1);
t290 = t224 * t222;
t298 = qJDD(4) - t290;
t300 = t249 * t298;
t299 = t251 * t298;
t254 = qJD(1) ^ 2;
t250 = sin(qJ(1));
t252 = cos(qJ(1));
t234 = t250 * g(1) - t252 * g(2);
t262 = qJDD(2) - t234;
t258 = -t254 * qJ(2) + t262;
t292 = -qJ(3) - pkin(1);
t265 = -0.2e1 * qJD(1) * qJD(3) + t292 * qJDD(1) + t258;
t244 = t247 ^ 2;
t245 = t248 ^ 2;
t279 = t244 + t245;
t296 = t279 * t254;
t181 = -t248 * g(3) + t265 * t247;
t259 = t261 * qJDD(1);
t218 = t222 ^ 2;
t219 = t224 ^ 2;
t294 = pkin(3) * t254;
t293 = t247 * g(3);
t291 = qJDD(1) * pkin(1);
t172 = t293 + (-pkin(5) * qJDD(1) - t247 * t294 + t265) * t248;
t274 = qJDD(1) * t247;
t174 = -pkin(5) * t274 - t244 * t294 + t181;
t146 = -t251 * t172 + t249 * t174;
t147 = t249 * t172 + t251 * t174;
t118 = -t251 * t146 + t249 * t147;
t289 = t247 * t118;
t288 = t247 * t248;
t287 = t248 * t118;
t235 = t252 * g(1) + t250 * g(2);
t246 = qJDD(1) * qJ(2);
t260 = t235 - t246;
t275 = qJD(2) * qJD(1);
t257 = -qJDD(3) + t260 - 0.2e1 * t275;
t179 = -pkin(3) * t274 + (t279 * pkin(5) - t292) * t254 + t257;
t286 = t249 * t179;
t184 = qJDD(4) + t290;
t285 = t249 * t184;
t228 = t279 * qJDD(1);
t284 = t250 * t228;
t283 = t251 * t179;
t282 = t251 * t184;
t281 = t252 * t228;
t280 = t244 - t245;
t278 = t222 * qJD(4);
t277 = t224 * qJD(4);
t273 = qJDD(1) * t248;
t272 = t250 * qJDD(1);
t271 = t252 * qJDD(1);
t270 = t250 * t290;
t269 = t252 * t290;
t119 = t249 * t146 + t251 * t147;
t242 = 0.2e1 * t275;
t211 = -t254 * pkin(1) + t242 - t260;
t213 = -t258 + t291;
t176 = t252 * t211 - t250 * t213;
t201 = -t250 * t234 - t252 * t235;
t203 = -t292 * t254 + t257;
t266 = -t203 + t246;
t232 = -t250 * t254 + t271;
t264 = pkin(4) * t232 + t250 * g(3);
t233 = t252 * t254 + t272;
t263 = -pkin(4) * t233 + t252 * g(3);
t221 = -t249 * t274 + t251 * t273;
t180 = t248 * t265 + t293;
t150 = t248 * t180 + t247 * t181;
t151 = -t247 * t180 + t248 * t181;
t175 = t250 * t211 + t252 * t213;
t200 = t252 * t234 - t250 * t235;
t226 = t247 * t296;
t198 = -t250 * t226 + t247 * t271;
t196 = t252 * t226 + t247 * t272;
t253 = qJD(4) ^ 2;
t230 = t280 * t254;
t229 = t280 * qJDD(1);
t225 = t248 * t296;
t210 = -t219 - t253;
t209 = -t219 + t253;
t208 = t218 - t253;
t206 = t233 * t288;
t205 = t232 * t288;
t199 = -t250 * t225 + t248 * t271;
t197 = t252 * t225 + t248 * t272;
t195 = -t252 * t296 - t284;
t194 = -t250 * t296 + t281;
t192 = t219 - t218;
t191 = pkin(2) * t274 - t248 * t203;
t190 = pkin(2) * t273 + t247 * t203;
t189 = t221 - t278;
t188 = t221 - 0.2e1 * t278;
t187 = -t259 - t277;
t186 = 0.2e1 * t277 + t259;
t182 = -t253 - t218;
t178 = (-t222 * t251 + t224 * t249) * qJD(4);
t177 = (-t222 * t249 - t224 * t251) * qJD(4);
t173 = -t218 - t219;
t171 = -pkin(2) * t226 + t180;
t170 = -pkin(2) * t225 - t181;
t167 = t251 * t189 - t249 * t277;
t166 = t249 * t189 + t251 * t277;
t165 = -t249 * t187 + t251 * t278;
t164 = t251 * t187 + t249 * t278;
t163 = -t249 * t210 - t282;
t162 = -t249 * t209 + t299;
t161 = t251 * t208 - t285;
t160 = t251 * t210 - t285;
t159 = t251 * t209 + t300;
t158 = t249 * t208 + t282;
t157 = -t251 * t186 - t249 * t188;
t156 = t249 * t221 - t251 * t259;
t155 = -t249 * t186 + t251 * t188;
t154 = -t251 * t221 - t249 * t259;
t153 = t251 * t182 - t300;
t152 = t249 * t182 + t299;
t149 = -pkin(2) * t296 - t151;
t148 = -t248 * t177 - t247 * t178;
t145 = -pkin(5) * t160 - t283;
t143 = t250 * t150 - t252 * t203;
t142 = -t252 * t150 - t250 * t203;
t141 = -pkin(5) * t152 - t286;
t140 = -t248 * t166 - t247 * t167;
t139 = -t248 * t164 - t247 * t165;
t138 = -t247 * t160 + t248 * t163;
t137 = t248 * t160 + t247 * t163;
t136 = -t248 * t159 - t247 * t162;
t135 = -t248 * t158 - t247 * t161;
t134 = -pkin(3) * t188 + pkin(5) * t163 - t286;
t133 = -t247 * t154 + t248 * t156;
t132 = -t248 * t155 - t247 * t157;
t131 = t248 * t154 + t247 * t156;
t130 = -pkin(3) * t186 + pkin(5) * t153 + t283;
t129 = pkin(2) * t150 - qJ(2) * t151;
t128 = -t247 * t152 + t248 * t153;
t127 = t248 * t152 + t247 * t153;
t126 = t250 * t137 + t252 * t188;
t125 = -t252 * t137 + t250 * t188;
t124 = -pkin(2) * t203 + t292 * t151;
t123 = t250 * t127 + t252 * t186;
t122 = -t252 * t127 + t250 * t186;
t121 = t250 * t131 + t252 * t173;
t120 = -t252 * t131 + t250 * t173;
t117 = pkin(3) * t179 + pkin(5) * t119;
t116 = -pkin(5) * t154 - t118;
t115 = -pkin(3) * t173 + pkin(5) * t156 + t119;
t114 = pkin(2) * t131 + pkin(3) * t154 - qJ(2) * t133;
t113 = pkin(2) * t137 + pkin(3) * t160 - qJ(2) * t138 - t147;
t112 = t248 * t119 - t289;
t111 = t247 * t119 + t287;
t110 = t250 * t111 - t252 * t179;
t109 = -t252 * t111 - t250 * t179;
t108 = pkin(2) * t127 + pkin(3) * t152 - qJ(2) * t128 - t146;
t107 = pkin(2) * t188 - t248 * t134 + t292 * t138 - t247 * t145;
t106 = pkin(2) * t186 + t292 * t128 - t248 * t130 - t247 * t141;
t105 = pkin(2) * t173 - t248 * t115 - t247 * t116 + t292 * t133;
t104 = pkin(2) * t111 + pkin(3) * t118 - qJ(2) * t112;
t103 = -pkin(2) * t179 + pkin(5) * t289 + t292 * t112 - t248 * t117;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t233, -t232, 0, t201, 0, 0, 0, 0, 0, 0, 0, t233, t232, t176, 0, 0, 0, 0, 0, 0, t198, t199, t195, t143, 0, 0, 0, 0, 0, 0, t123, t126, t121, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t232, -t233, 0, t200, 0, 0, 0, 0, 0, 0, 0, -t232, t233, t175, 0, 0, 0, 0, 0, 0, t196, t197, t194, t142, 0, 0, 0, 0, 0, 0, t122, t125, t120, t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, 0, 0, 0, 0, 0, 0, t128, t138, t133, t112; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t232, 0, -t233, 0, -t264, -t263, -t200, -pkin(4) * t200, 0, -t232, t233, 0, 0, 0, -t175, t264, t263, -pkin(4) * t175 + (-pkin(1) * t250 + qJ(2) * t252) * g(3), t206, -t250 * t229 - t252 * t230, t199, -t206, -t198, 0, -pkin(4) * t196 + t171 * t252 - t191 * t250, -pkin(4) * t197 + t170 * t252 - t190 * t250, -pkin(2) * t281 - pkin(4) * t194 - t149 * t250, -pkin(4) * t142 - t124 * t250 + t129 * t252, -t140 * t250 + t269, -t132 * t250 + t192 * t252, -t136 * t250 + t221 * t252, -t139 * t250 - t269, -t135 * t250 - t252 * t259, qJDD(4) * t252 - t148 * t250, -pkin(4) * t122 - t106 * t250 + t108 * t252, -pkin(4) * t125 - t107 * t250 + t113 * t252, -pkin(4) * t120 - t105 * t250 + t114 * t252, -pkin(4) * t109 - t103 * t250 + t104 * t252; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t233, 0, t232, 0, t263, -t264, t201, pkin(4) * t201, 0, -t233, -t232, 0, 0, 0, t176, -t263, t264, pkin(4) * t176 + (pkin(1) * t252 + qJ(2) * t250) * g(3), -t205, t229 * t252 - t230 * t250, t197, t205, -t196, 0, pkin(4) * t198 + t171 * t250 + t191 * t252, pkin(4) * t199 + t170 * t250 + t190 * t252, -pkin(2) * t284 + pkin(4) * t195 + t149 * t252, pkin(4) * t143 + t124 * t252 + t129 * t250, t140 * t252 + t270, t132 * t252 + t192 * t250, t136 * t252 + t221 * t250, t139 * t252 - t270, t135 * t252 - t250 * t259, qJDD(4) * t250 + t148 * t252, pkin(4) * t123 + t106 * t252 + t108 * t250, pkin(4) * t126 + t107 * t252 + t113 * t250, pkin(4) * t121 + t105 * t252 + t114 * t250, pkin(4) * t110 + t103 * t252 + t104 * t250; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t234, t235, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t262 - 0.2e1 * t291, -t235 + t242 + 0.2e1 * t246, pkin(1) * t213 + qJ(2) * t211, t245 * qJDD(1), -0.2e1 * t247 * t273, 0, t244 * qJDD(1), 0, 0, -t292 * t226 + t266 * t247, -t292 * t225 + t266 * t248, -qJ(2) * t296 - t292 * t228 - t150, -qJ(2) * t203 + t292 * t150, -t166 * t247 + t167 * t248, -t155 * t247 + t157 * t248, -t159 * t247 + t162 * t248, -t164 * t247 + t165 * t248, -t158 * t247 + t161 * t248, -t177 * t247 + t178 * t248, qJ(2) * t186 + t292 * t127 - t247 * t130 + t248 * t141, qJ(2) * t188 - t247 * t134 + t292 * t137 + t248 * t145, qJ(2) * t173 - t247 * t115 + t248 * t116 + t292 * t131, -pkin(5) * t287 - qJ(2) * t179 + t292 * t111 - t247 * t117;];
tauB_reg = t1;
