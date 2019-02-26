% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR13_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:50
% EndTime: 2019-02-26 21:44:52
% DurationCPUTime: 1.55s
% Computational Cost: add. (4128->141), mult. (12381->291), div. (702->12), fcn. (15714->13), ass. (0->124)
t222 = sin(qJ(2));
t223 = sin(qJ(1));
t225 = cos(qJ(2));
t226 = cos(qJ(1));
t299 = cos(pkin(6));
t255 = t226 * t299;
t209 = t222 * t223 - t225 * t255;
t221 = sin(qJ(4));
t224 = cos(qJ(4));
t219 = sin(pkin(6));
t276 = t219 * t226;
t245 = t209 * t224 + t221 * t276;
t193 = t245 ^ 2;
t277 = t219 * t225;
t241 = -t299 * t221 - t224 * t277;
t205 = 0.1e1 / t241 ^ 2;
t184 = t193 * t205 + 0.1e1;
t182 = 0.1e1 / t184;
t273 = qJD(1) * t223;
t260 = t219 * t273;
t240 = -t222 * t255 - t223 * t225;
t256 = t223 * t299;
t242 = t226 * t222 + t225 * t256;
t234 = t242 * qJD(1) - t240 * qJD(2);
t261 = t224 * t276;
t305 = -qJD(4) * t261 - t234 * t224;
t167 = (qJD(4) * t209 + t260) * t221 + t305;
t208 = -t221 * t277 + t299 * t224;
t272 = qJD(2) * t222;
t258 = t219 * t272;
t195 = t208 * qJD(4) - t224 * t258;
t204 = 0.1e1 / t241;
t283 = t245 * t205;
t248 = t167 * t204 - t195 * t283;
t151 = t248 * t182;
t185 = atan2(t245, -t241);
t180 = sin(t185);
t181 = cos(t185);
t250 = t180 * t241 + t181 * t245;
t146 = t250 * t151 - t167 * t180 + t181 * t195;
t161 = t180 * t245 - t181 * t241;
t159 = 0.1e1 / t161 ^ 2;
t306 = t146 * t159;
t270 = qJD(4) * t221;
t257 = t219 * t270;
t269 = qJD(4) * t224;
t168 = t209 * t269 + t234 * t221 + t224 * t260 + t226 * t257;
t158 = 0.1e1 / t161;
t304 = t158 * t306;
t278 = t219 * t224;
t197 = t242 * t221 + t223 * t278;
t251 = t222 * t256;
t275 = t226 * t225;
t237 = (t299 * qJD(1) + qJD(2)) * t275 - qJD(2) * t251 - t222 * t273;
t259 = qJD(1) * t276;
t170 = t197 * qJD(4) + t221 * t259 - t237 * t224;
t238 = t242 * t224;
t196 = t223 * t219 * t221 - t238;
t254 = 0.2e1 * t196 * t304;
t303 = -t159 * t170 + t254;
t171 = qJD(4) * t238 + t237 * t221 - t223 * t257 + t224 * t259;
t190 = t240 * qJD(1) - t242 * qJD(2);
t218 = sin(pkin(11));
t220 = cos(pkin(11));
t163 = t171 * t220 + t190 * t218;
t211 = -t251 + t275;
t177 = t197 * t220 + t211 * t218;
t174 = 0.1e1 / t177 ^ 2;
t302 = t163 * t174;
t301 = t195 * t205;
t279 = t219 * t222;
t262 = t245 * t279;
t243 = t204 * t240 + t205 * t262;
t300 = t224 * t243;
t173 = 0.1e1 / t177;
t192 = t196 ^ 2;
t157 = t159 * t192 + 0.1e1;
t293 = t159 * t196;
t298 = (t170 * t293 - t192 * t304) / t157 ^ 2;
t162 = t171 * t218 - t190 * t220;
t176 = t197 * t218 - t211 * t220;
t172 = t176 ^ 2;
t166 = t172 * t174 + 0.1e1;
t290 = t174 * t176;
t292 = t173 * t302;
t297 = (t162 * t290 - t172 * t292) / t166 ^ 2;
t285 = t204 * t301;
t295 = (-t167 * t283 + t193 * t285) / t184 ^ 2;
t291 = t173 * t218;
t281 = t211 * t221;
t187 = -t242 * t218 + t220 * t281;
t289 = t174 * t187;
t288 = t176 * t220;
t287 = t180 * t196;
t286 = t181 * t196;
t284 = t245 * t204;
t282 = t245 * t208;
t280 = t211 * t224;
t271 = qJD(2) * t225;
t268 = 0.2e1 * t298;
t267 = 0.2e1 * t297;
t266 = -0.2e1 * t295;
t265 = 0.2e1 * t295;
t263 = t176 * t292;
t253 = t204 * t265;
t252 = 0.2e1 * t263;
t246 = -t209 * t221 + t261;
t247 = t204 * t246 + t205 * t282;
t244 = t190 * t221 + t211 * t269;
t239 = -t180 + (t181 * t284 + t180) * t182;
t194 = t241 * qJD(4) + t221 * t258;
t191 = -qJD(1) * t251 - t223 * t272 + (qJD(2) * t299 + qJD(1)) * t275;
t186 = t218 * t281 + t242 * t220;
t179 = t218 * t240 + t220 * t246;
t178 = t218 * t246 - t220 * t240;
t164 = 0.1e1 / t166;
t155 = 0.1e1 / t157;
t154 = t182 * t300;
t153 = t247 * t182;
t148 = (-t180 * t240 - t181 * t279) * t224 + t250 * t154;
t147 = -t250 * t153 + t180 * t246 + t181 * t208;
t145 = t247 * t265 + (-0.2e1 * t282 * t285 + t168 * t204 + (t167 * t208 - t194 * t245 - t195 * t246) * t205) * t182;
t143 = t266 * t300 + (-t243 * t270 + (0.2e1 * t262 * t285 - t191 * t204 + (t195 * t240 + (-t167 * t222 + t245 * t271) * t219) * t205) * t224) * t182;
t1 = [-t196 * t253 + (t170 * t204 + t196 * t301) * t182, t143, 0, t145, 0, 0; -0.2e1 * t245 * t158 * t298 + ((-t209 * t270 - t221 * t260 - t305) * t158 - t245 * t306 - (t239 * t170 + ((-t151 * t182 * t284 + t266) * t180 + (-t245 * t253 - t151 + (t151 - t248) * t182) * t181) * t196) * t293) * t155 + (t303 * t155 + t293 * t268) * t239 * t196 (t148 * t293 + t158 * t280) * t268 + ((-t190 * t224 + t211 * t270) * t158 + t303 * t148 + (t280 * t146 - (t143 * t245 - t154 * t167 + (t222 * t270 - t224 * t271) * t219 + (t154 * t241 - t224 * t240) * t151) * t286 - (t240 * t270 + t143 * t241 - t154 * t195 + t191 * t224 + (-t154 * t245 + t222 * t278) * t151) * t287) * t159) * t155, 0 (t147 * t293 - t158 * t197) * t268 + (t147 * t254 + t171 * t158 + (-t197 * t146 - t147 * t170 - (t145 * t245 + t153 * t167 + t194 + (-t153 * t241 + t246) * t151) * t286 - (t145 * t241 + t153 * t195 - t168 + (t153 * t245 - t208) * t151) * t287) * t159) * t155, 0, 0; (-t173 * t178 + t179 * t290) * t267 + ((-t168 * t218 + t191 * t220) * t173 + t179 * t252 + (-t178 * t163 - (-t168 * t220 - t191 * t218) * t176 - t179 * t162) * t174) * t164, -0.2e1 * t186 * t173 * t297 + t176 * t267 * t289 + ((t244 * t218 + t220 * t237) * t173 - t186 * t302 - (-t218 * t237 + t244 * t220) * t290 - t162 * t289 + t187 * t252) * t164, 0 (-t174 * t288 + t291) * t196 * t267 + (-0.2e1 * t196 * t220 * t263 - t170 * t291 + (t170 * t288 + (t162 * t220 + t163 * t218) * t196) * t174) * t164, 0, 0;];
JaD_rot  = t1;
