% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:33:03
% EndTime: 2019-02-26 22:33:04
% DurationCPUTime: 1.02s
% Computational Cost: add. (4319->120), mult. (5914->254), div. (754->12), fcn. (6974->11), ass. (0->113)
t239 = sin(qJ(1));
t234 = t239 ^ 2;
t236 = qJ(2) + qJ(3);
t231 = sin(t236);
t227 = t231 ^ 2;
t232 = cos(t236);
t229 = 0.1e1 / t232 ^ 2;
t294 = t227 * t229;
t222 = t234 * t294 + 0.1e1;
t226 = t231 * t227;
t228 = 0.1e1 / t232;
t233 = qJD(2) + qJD(3);
t291 = t228 * t231;
t250 = t233 * (t226 * t228 * t229 + t291);
t242 = cos(qJ(1));
t280 = qJD(1) * t242;
t292 = t227 * t239;
t259 = t280 * t292;
t301 = (t229 * t259 + t234 * t250) / t222 ^ 2;
t312 = -0.2e1 * t301;
t266 = 0.1e1 + t294;
t311 = t239 * t266;
t219 = 0.1e1 / t222;
t267 = t231 * t280;
t288 = t233 * t239;
t270 = t229 * t288;
t186 = ((t232 * t288 + t267) * t228 + t227 * t270) * t219;
t310 = t186 - t288;
t241 = cos(qJ(4));
t309 = (-qJD(4) + qJD(6)) * t241;
t238 = sin(qJ(4));
t279 = qJD(4) * t238;
t308 = -qJD(6) * t238 + t279;
t282 = t241 * t242;
t284 = t239 * t238;
t251 = t232 * t284 + t282;
t287 = t233 * t242;
t268 = t231 * t287;
t269 = t232 * t282;
t190 = t251 * qJD(1) - qJD(4) * t269 + t238 * t268 - t239 * t279;
t261 = -qJD(1) * t232 + qJD(4);
t262 = qJD(4) * t232 - qJD(1);
t286 = t238 * t242;
t191 = -t262 * t286 + (t261 * t239 - t268) * t241;
t237 = sin(qJ(6));
t240 = cos(qJ(6));
t283 = t239 * t241;
t215 = t232 * t286 - t283;
t216 = t269 + t284;
t255 = t215 * t240 - t216 * t237;
t181 = t255 * qJD(6) - t190 * t237 + t191 * t240;
t209 = t215 * t237 + t216 * t240;
t201 = 0.1e1 / t209;
t253 = t237 * t238 + t240 * t241;
t254 = t237 * t241 - t238 * t240;
t202 = 0.1e1 / t209 ^ 2;
t298 = t202 * t255;
t307 = t254 * t201 + t253 * t298;
t285 = t239 * t231;
t221 = atan2(t285, t232);
t218 = cos(t221);
t217 = sin(t221);
t271 = t217 * t285;
t198 = t218 * t232 + t271;
t195 = 0.1e1 / t198;
t196 = 0.1e1 / t198 ^ 2;
t306 = t219 - 0.1e1;
t180 = t209 * qJD(6) + t190 * t240 + t191 * t237;
t200 = t255 ^ 2;
t185 = t200 * t202 + 0.1e1;
t203 = t201 * t202;
t302 = t181 * t203;
t305 = (-t180 * t298 - t200 * t302) / t185 ^ 2;
t235 = t242 ^ 2;
t293 = t227 * t235;
t189 = t196 * t293 + 0.1e1;
t289 = t232 * t233;
t295 = t218 * t231;
t177 = (t186 * t239 - t233) * t295 + (-t310 * t232 + t267) * t217;
t303 = t177 * t195 * t196;
t304 = (-t293 * t303 + (t231 * t235 * t289 - t259) * t196) / t189 ^ 2;
t300 = t196 * t231;
t299 = t196 * t242;
t290 = t231 * t242;
t211 = t253 * t290;
t297 = t202 * t211;
t296 = t217 * t239;
t281 = qJD(1) * t239;
t275 = 0.2e1 * t305;
t274 = -0.2e1 * t303;
t273 = -0.2e1 * t203 * t255;
t272 = t196 * t290;
t265 = -0.2e1 * t231 * t304;
t264 = t228 * t312;
t263 = t181 * t273;
t260 = t218 * t219 * t227 * t228;
t258 = t266 * t242;
t214 = -t232 * t283 + t286;
t256 = -t214 * t237 - t240 * t251;
t205 = t214 * t240 - t237 * t251;
t252 = t261 * t242;
t210 = t254 * t290;
t199 = t219 * t311;
t193 = t241 * t252 + (t231 * t233 * t241 + t262 * t238) * t239;
t192 = -t262 * t283 + (t233 * t285 + t252) * t238;
t187 = 0.1e1 / t189;
t183 = 0.1e1 / t185;
t182 = (-t306 * t231 * t217 + t239 * t260) * t242;
t179 = t232 * t296 - t295 + (-t217 * t232 + t218 * t285) * t199;
t178 = t311 * t312 + (qJD(1) * t258 + 0.2e1 * t239 * t250) * t219;
t174 = (t201 * t210 + t255 * t297) * t275 + (t180 * t297 + (t210 * t202 - t211 * t273) * t181 - t307 * t232 * t287 + (t307 * t281 + ((-t309 * t201 + t308 * t298) * t240 + (t308 * t201 + t309 * t298) * t237) * t242) * t231) * t183;
t173 = 0.2e1 * (-t179 * t300 + t195 * t232) * t242 * t304 + ((t195 * t281 + (t179 * t233 + t177) * t299) * t232 + (t195 * t287 + (t178 * t218 * t239 + t310 * t217 + (-t186 * t296 + t217 * t233 + t218 * t280) * t199) * t272 + (-t196 * t281 + t242 * t274) * t179 + ((-t178 + t280) * t217 + ((t199 * t239 - 0.1e1) * t233 + (-t199 + t239) * t186) * t218) * t232 * t299) * t231) * t187;
t1 = [t264 * t290 + (t233 * t258 - t281 * t291) * t219, t178, t178, 0, 0, 0; (t195 * t265 + (t195 * t289 + (-qJD(1) * t182 - t177) * t300) * t187) * t239 + (t196 * t265 * t182 + (((0.2e1 * t231 * t301 + t289 + (-t186 * t228 * t292 - t289) * t219) * t217 + (t264 * t292 + t186 * t231 + (t226 * t270 + (-t186 + 0.2e1 * t288) * t231) * t219) * t218) * t272 + (t196 * t289 + t231 * t274) * t182 + (t195 + ((-t234 + t235) * t260 + t306 * t271) * t196) * t231 * qJD(1)) * t187) * t242, t173, t173, 0, 0, 0; (t201 * t256 - t205 * t298) * t275 + ((t205 * qJD(6) - t192 * t240 + t193 * t237) * t201 + t205 * t263 + (t256 * t181 + (t256 * qJD(6) + t192 * t237 + t193 * t240) * t255 - t205 * t180) * t202) * t183, t174, t174 (t201 * t209 + t255 * t298) * t275 + (-t181 * t201 - t255 * t263 + (0.2e1 * t255 * t180 + t209 * t181) * t202) * t183, 0, -0.2e1 * t305 - 0.2e1 * (t180 * t202 * t183 - (-t183 * t302 - t202 * t305) * t255) * t255;];
JaD_rot  = t1;
