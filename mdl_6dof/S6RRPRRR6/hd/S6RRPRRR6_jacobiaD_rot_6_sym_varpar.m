% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:58
% EndTime: 2019-02-26 21:57:00
% DurationCPUTime: 1.42s
% Computational Cost: add. (10722->102), mult. (13073->197), div. (974->12), fcn. (16295->11), ass. (0->101)
t296 = qJ(4) + qJ(5);
t286 = sin(t296);
t287 = cos(t296);
t315 = sin(qJ(2));
t317 = cos(qJ(2));
t228 = -t317 * t286 + t315 * t287;
t316 = sin(qJ(1));
t219 = t228 * t316;
t227 = t315 * t286 + t317 * t287;
t208 = atan2(t219, t227);
t203 = sin(t208);
t204 = cos(t208);
t221 = t227 * t316;
t279 = -t203 * t227 + t204 * t219;
t217 = t219 ^ 2;
t225 = 0.1e1 / t227 ^ 2;
t207 = t217 * t225 + 0.1e1;
t205 = 0.1e1 / t207;
t224 = 0.1e1 / t227;
t298 = t219 * t228;
t271 = -t221 * t224 - t225 * t298;
t327 = t271 * t205;
t178 = -t203 * t221 + t204 * t228 + t279 * t327;
t318 = cos(qJ(1));
t223 = t228 * t318;
t345 = t178 * t223;
t191 = t203 * t219 + t204 * t227;
t188 = 0.1e1 / t191;
t189 = 0.1e1 / t191 ^ 2;
t222 = t227 * t318;
t218 = t223 ^ 2;
t187 = t218 * t189 + 0.1e1;
t295 = qJD(4) + qJD(5);
t266 = t295 * t286;
t267 = t295 * t287;
t325 = t315 * t266 + t317 * t267;
t195 = -t219 * qJD(1) + t222 * qJD(2) - t325 * t318;
t306 = t195 * t189;
t196 = -t223 * qJD(1) - t221 * qJD(2) + t325 * t316;
t324 = -t317 * t266 + t315 * t267;
t202 = -qJD(2) * t228 + t324;
t299 = t219 * t225;
t273 = -t196 * t224 - t202 * t299;
t181 = t273 * t205;
t176 = t181 * t279 - t203 * t196 + t204 * t202;
t314 = t176 * t188 * t189;
t294 = 0.2e1 * (-t218 * t314 + t223 * t306) / t187 ^ 2;
t344 = (t188 * t222 + t189 * t345) * t294;
t194 = t221 * qJD(1) + t223 * qJD(2) - t324 * t318;
t246 = sin(qJ(6));
t247 = cos(qJ(6));
t216 = t222 * t247 - t246 * t316;
t285 = qJD(1) * t318;
t192 = qJD(6) * t216 - t194 * t246 + t247 * t285;
t269 = -t222 * t246 - t247 * t316;
t335 = qJD(6) * t269;
t193 = -t194 * t247 - t246 * t285 + t335;
t209 = t269 ^ 2;
t211 = 0.1e1 / t216 ^ 2;
t200 = t209 * t211 + 0.1e1;
t198 = 0.1e1 / t200;
t210 = 0.1e1 / t216;
t301 = t211 * t269;
t272 = -t246 * t210 - t247 * t301;
t308 = t193 * t210 * t211;
t320 = -0.2e1 * t269;
t283 = t308 * t320;
t343 = (t272 * t195 - ((t246 * (-t193 - t335) - t192 * t247) * t211 + (qJD(6) * t210 + t283) * t247) * t223) * t198;
t342 = -t222 * t176 + t178 * t195;
t293 = 0.2e1 * t314;
t341 = -t194 * t188 - t293 * t345;
t201 = -t227 * qJD(2) + t325;
t339 = (t227 * t327 + t221) * t181 + t327 * t196 + t201;
t197 = qJD(1) * t222 - t219 * qJD(2) + t324 * t316;
t338 = -(-t219 * t327 - t228) * t181 + t327 * t202 + t197;
t329 = t202 * t225;
t304 = t224 * t329;
t336 = ((t196 * t228 + t201 * t219 + t202 * t221) * t225 - t197 * t224 + 0.2e1 * t298 * t304) * t205;
t300 = t219 * t224;
t323 = t205 * (-t204 * t300 + t203) - t203;
t319 = -0.2e1 * t223;
t313 = (-t192 * t301 - t209 * t308) / t200 ^ 2;
t312 = (-t196 * t299 - t217 * t304) / t207 ^ 2;
t305 = t198 * t211;
t303 = t203 * t223;
t302 = t204 * t223;
t292 = -0.2e1 * t313;
t291 = -0.2e1 * t312;
t290 = t211 * t313;
t289 = t224 * t312;
t288 = t192 * t305;
t284 = qJD(1) * t316;
t214 = -t221 * t247 - t246 * t318;
t270 = t221 * t246 - t247 * t318;
t185 = 0.1e1 / t187;
t180 = t323 * t223;
t174 = t271 * t291 + t336;
t173 = 0.2e1 * t271 * t312 - t336;
t172 = -t272 * t223 * t292 - t343;
t171 = -t344 + (((t174 * t219 - t339) * t302 + (-t174 * t227 - t338) * t303 + t342) * t189 + t341) * t185;
t1 = [t289 * t319 + (t195 * t224 - t223 * t329) * t205, t173, 0, t174, t174, 0; -t219 * t188 * t294 + (-t196 * t188 + (-t176 * t219 - t180 * t195) * t189) * t185 - ((-t180 * t293 + t323 * t306) * t185 + (-t180 * t294 + ((t181 * t205 * t300 + t291) * t303 + (0.2e1 * t219 * t289 - t181 + (t181 - t273) * t205) * t302) * t185) * t189) * t223, t344 + (((t173 * t219 + t339) * t302 + (-t173 * t227 + t338) * t303 - t342) * t189 - t341) * t185, 0, t171, t171, 0; (t290 * t320 - t288) * t214 - (-t193 * t305 + t210 * t292) * t270 + ((qJD(6) * t214 - t197 * t246 - t247 * t284) * t210 + (qJD(6) * t270 - t197 * t247 + t246 * t284) * t301 + t214 * t283) * t198, t272 * t313 * t319 + t343, 0, t172, t172, t292 + (t288 - (-t198 * t308 - t290) * t269) * t320;];
JaD_rot  = t1;
