% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPP5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:29
% EndTime: 2019-12-31 20:58:38
% DurationCPUTime: 3.32s
% Computational Cost: add. (5078->302), mult. (10951->324), div. (0->0), fcn. (7093->6), ass. (0->176)
t173 = sin(qJ(3));
t176 = cos(qJ(3));
t177 = cos(qJ(2));
t174 = sin(qJ(2));
t223 = qJD(1) * t174;
t140 = -qJD(1) * t176 * t177 + t173 * t223;
t169 = qJD(2) + qJD(3);
t240 = t169 * t140;
t161 = t174 * qJDD(1);
t219 = qJD(1) * qJD(2);
t215 = t177 * t219;
t148 = t161 + t215;
t162 = t177 * qJDD(1);
t216 = t174 * t219;
t149 = t162 - t216;
t202 = t176 * t148 + t173 * t149;
t83 = -qJD(3) * t140 + t202;
t268 = t83 - t240;
t228 = t174 * t176;
t142 = (t173 * t177 + t228) * qJD(1);
t109 = t142 * t140;
t168 = qJDD(2) + qJDD(3);
t266 = t109 + t168;
t277 = t176 * t266;
t139 = t142 ^ 2;
t259 = t169 ^ 2;
t99 = -t259 - t139;
t306 = t173 * t99 + t277;
t281 = t173 * t266;
t66 = t176 * t99 - t281;
t322 = pkin(6) * (t174 * t66 + t177 * t306);
t323 = pkin(1) * t268 + t322;
t319 = pkin(2) * t66;
t317 = pkin(7) * t66;
t260 = t140 ^ 2;
t264 = -t260 - t139;
t287 = pkin(1) * t264;
t269 = t240 + t83;
t45 = t176 * t269;
t208 = t173 * t148 - t149 * t176;
t55 = (qJD(3) - t169) * t142 + t208;
t29 = t173 * t55 + t45;
t249 = t173 * t269;
t32 = t176 * t55 - t249;
t321 = pkin(6) * (t174 * t29 - t177 * t32) - t287;
t320 = pkin(2) * t29;
t318 = pkin(7) * t29;
t316 = pkin(7) * t306;
t250 = t173 * t268;
t239 = t169 * t142;
t82 = qJD(3) * t142 + t208;
t270 = t239 + t82;
t303 = -t177 * (-t173 * t270 + t176 * t268) + t174 * (t176 * t270 + t250);
t265 = -t259 - t260;
t267 = -t109 + t168;
t280 = t173 * t267;
t293 = t176 * t265 - t280;
t301 = pkin(7) * t293;
t314 = -pkin(2) * t270 + t301;
t286 = pkin(2) * t264;
t313 = -pkin(7) * t32 - t286;
t263 = t260 - t259;
t310 = t174 * (-t176 * t263 + t281) - t177 * (t173 * t263 + t277);
t123 = t139 - t259;
t276 = t176 * t267;
t290 = t174 * (t123 * t173 + t276) - t177 * (t123 * t176 - t280);
t294 = t173 * t265 + t276;
t302 = pkin(2) * t294;
t300 = pkin(7) * t294;
t105 = pkin(3) * t140 - qJ(4) * t142;
t258 = 2 * qJD(4);
t179 = qJD(1) ^ 2;
t227 = t174 * t179;
t175 = sin(qJ(1));
t255 = cos(qJ(1));
t200 = g(1) * t255 + g(2) * t175;
t242 = qJDD(1) * pkin(6);
t144 = -pkin(1) * t179 - t200 + t242;
t230 = t174 * t144;
t185 = qJDD(2) * pkin(2) - t148 * pkin(7) - t230 + (pkin(2) * t227 + pkin(7) * t219 - g(3)) * t177;
t118 = -t174 * g(3) + t177 * t144;
t172 = t177 ^ 2;
t164 = t172 * t179;
t199 = qJD(2) * pkin(2) - pkin(7) * t223;
t74 = -pkin(2) * t164 + t149 * pkin(7) - qJD(2) * t199 + t118;
t37 = t173 * t185 + t176 * t74;
t204 = qJ(4) * t168 - t105 * t140 + t169 * t258 + t37;
t20 = -(t259 + t264) * pkin(3) + t204;
t291 = pkin(6) * (-t174 * t294 + t177 * t293) - pkin(1) * t270;
t299 = qJ(4) * t264;
t298 = qJ(5) * t269;
t292 = -pkin(3) * t99 + qJ(4) * t266;
t288 = 2 * qJD(5);
t257 = pkin(3) + pkin(4);
t275 = t268 * qJ(4);
t271 = pkin(3) * t267 + qJ(4) * t265;
t56 = t239 - t82;
t106 = t139 - t260;
t119 = -pkin(4) * t169 - qJ(5) * t142;
t262 = qJ(5) * t82 + t119 * t169 + t140 * t288;
t238 = t169 * t173;
t120 = t142 * t238;
t237 = t169 * t176;
t28 = t174 * (t176 * t83 - t120) + t177 * (t142 * t237 + t173 * t83);
t261 = -t82 * pkin(4) - t260 * qJ(5) + qJDD(5) + (t258 + t119) * t142;
t256 = t82 * pkin(3);
t254 = pkin(3) * t176;
t24 = -pkin(3) * t259 + t204;
t36 = t173 * t74 - t176 * t185;
t196 = -pkin(3) * t168 - qJ(4) * t259 + qJDD(4) + t36;
t241 = t142 * t105;
t26 = t196 + t241;
t253 = -pkin(3) * t26 + qJ(4) * t24;
t252 = -pkin(3) * t269 + qJ(4) * t56;
t212 = t175 * g(1) - g(2) * t255;
t198 = qJDD(1) * pkin(1) + t212;
t87 = t149 * pkin(2) - t199 * t223 + (pkin(7) * t172 + pkin(6)) * t179 + t198;
t248 = t173 * t87;
t18 = t173 * t37 - t176 * t36;
t247 = t174 * t18;
t245 = t176 * t87;
t243 = qJ(4) * t176;
t155 = t177 * t227;
t229 = t174 * (qJDD(2) + t155);
t224 = t177 * (qJDD(2) - t155);
t218 = t142 * t258;
t217 = (-t270 - t82) * pkin(3);
t214 = -qJ(4) * t173 - pkin(2);
t213 = -pkin(4) * t140 - t105;
t19 = t173 * t36 + t176 * t37;
t117 = t177 * g(3) + t230;
t209 = t174 * t117 + t118 * t177;
t193 = -t168 * pkin(4) + t196 - t298;
t13 = (-(2 * qJD(5)) - t213) * t142 + t193;
t192 = t24 + t262;
t15 = -pkin(4) * t260 + t192;
t206 = qJ(4) * t15 - t13 * t257;
t205 = qJ(4) * t55 + t257 * t269;
t195 = t177 * (-t140 * t173 - t142 * t176);
t194 = t24 + t292;
t191 = t142 * t288 - t193;
t190 = -t26 + t271;
t189 = t174 * (t140 * t237 + t173 * t82) + t177 * (t140 * t238 - t176 * t82);
t188 = -pkin(3) * t239 + t87;
t187 = (-t99 - t260) * pkin(4) + t192 + t292;
t186 = t174 * t120 + (-t140 * t228 + t195) * t169;
t184 = -t241 + (t267 - t109) * pkin(4) + t191 + t271;
t183 = t188 + t275;
t182 = t183 + t218;
t181 = t188 - t256 + 0.2e1 * t275;
t180 = t183 + t261;
t178 = qJD(2) ^ 2;
t171 = t174 ^ 2;
t163 = t171 * t179;
t150 = t162 - 0.2e1 * t216;
t147 = t161 + 0.2e1 * t215;
t143 = t179 * pkin(6) + t198;
t58 = (-qJD(3) - t169) * t140 + t202;
t38 = -qJ(4) * t270 + qJ(5) * t267;
t33 = t176 * t56 + t249;
t30 = t173 * t56 - t45;
t27 = -qJ(5) * t266 + t257 * t268;
t21 = t26 - t299;
t17 = t217 + t182;
t16 = t181 + t218;
t11 = t180 - t256;
t9 = t173 * t24 - t176 * t26;
t8 = -qJ(5) * t99 + t181 + t261;
t7 = t142 * t213 + t191 + t298 + t299;
t6 = -qJ(5) * t55 + (t260 + t264) * pkin(4) - t262 - t20;
t5 = -pkin(4) * t270 - qJ(5) * t265 + t180 + t217;
t4 = t13 * t173 + t15 * t176;
t3 = -t13 * t176 + t15 * t173;
t2 = qJ(4) * t11 - qJ(5) * t13;
t1 = -qJ(5) * t15 + t11 * t257;
t10 = [0, 0, 0, 0, 0, qJDD(1), t212, t200, 0, 0, (t148 + t215) * t174, t147 * t177 + t150 * t174, t229 + t177 * (-t163 + t178), (t149 - t216) * t177, t174 * (t164 - t178) + t224, 0, t177 * t143 + pkin(1) * t150 + pkin(6) * (t177 * (-t164 - t178) - t229), -t174 * t143 - pkin(1) * t147 + pkin(6) * (-t224 - t174 * (-t163 - t178)), pkin(1) * (t163 + t164) + (t171 + t172) * t242 + t209, pkin(1) * t143 + pkin(6) * t209, t28, -t303, t290, t189, -t310, t186, t174 * (-t248 - t300) + t177 * (t245 + t314) + t291, t174 * (-t245 - t317) + t177 * (-pkin(2) * t58 - t248 - t316) - pkin(1) * t58 - t322, t174 * (-t18 + t318) + t177 * (t19 + t313) + t321, -pkin(7) * t247 + t177 * (pkin(2) * t87 + pkin(7) * t19) + pkin(1) * t87 + pkin(6) * (t177 * t19 - t247), t28, t290, t303, t186, t310, t189, t174 * (-t17 * t173 - t243 * t270 - t300) + t177 * (t176 * t17 + t214 * t270 + t301) + t291, t174 * (-pkin(7) * t30 - t173 * t20 + t176 * t21) + t177 * (pkin(7) * t33 + t173 * t21 + t176 * t20 - t286) - t287 + pkin(6) * (-t174 * t30 + t177 * t33), t174 * (-pkin(3) * t250 + t16 * t176 + t317) + t177 * (t316 + t173 * t16 + (pkin(2) + t254) * t268) + t323, (t174 * (-pkin(3) * t173 + t243) + t177 * (-t214 + t254) + pkin(1)) * (t182 - t256) + (pkin(6) + pkin(7)) * (t177 * (t173 * t26 + t176 * t24) - t174 * t9), t28, t303, -t290, t189, -t310, (t174 * (-t140 * t176 + t142 * t173) + t195) * t169, t174 * (-t173 * t5 + t176 * t38 - t300) + t177 * (t173 * t38 + t176 * t5 + t314) + t291, t174 * (-t173 * t27 + t176 * t8 + t317) + t177 * (pkin(2) * t268 + t173 * t8 + t176 * t27 + t316) + t323, t174 * (-t173 * t6 + t176 * t7 - t318) + t177 * (t173 * t7 + t176 * t6 - t313) - t321, t174 * (-pkin(7) * t3 - t1 * t173 + t176 * t2) + t177 * (pkin(2) * t11 + pkin(7) * t4 + t1 * t176 + t173 * t2) + pkin(1) * t11 + pkin(6) * (-t174 * t3 + t177 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t164 + t163, t161, t155, t162, qJDD(2), -t117, -t118, 0, 0, t109, t106, t269, -t109, t56, t168, -t36 + t302, -t37 + t319, -t320, pkin(2) * t18, t109, t269, -t106, t168, -t56, -t109, t190 + t302, pkin(2) * t30 + t252, t194 - t319, pkin(2) * t9 + t253, t109, -t106, -t269, -t109, t56, t168, t184 + t302, t187 - t319, t205 + t320, pkin(2) * t3 + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, t106, t269, -t109, t56, t168, -t36, -t37, 0, 0, t109, t269, -t106, t168, -t56, -t109, t190, t252, t194, t253, t109, -t106, -t269, -t109, t56, t168, t184, t187, t205, t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t267, t269, t99, t26, 0, 0, 0, 0, 0, 0, -t267, t99, -t269, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t270, t268, t264, t11;];
tauJ_reg = t10;
