% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR7
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

function JaD_rot = S6RRPRRR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:57:35
% EndTime: 2019-02-26 21:57:37
% DurationCPUTime: 1.30s
% Computational Cost: add. (3949->95), mult. (10113->188), div. (729->12), fcn. (12551->11), ass. (0->97)
t309 = sin(qJ(4));
t310 = sin(qJ(2));
t312 = cos(qJ(4));
t313 = cos(qJ(2));
t228 = -t313 * t309 + t310 * t312;
t311 = sin(qJ(1));
t219 = t228 * t311;
t227 = t310 * t309 + t313 * t312;
t210 = atan2(t219, t227);
t201 = sin(t210);
t202 = cos(t210);
t221 = t227 * t311;
t273 = -t201 * t227 + t202 * t219;
t217 = t219 ^ 2;
t225 = 0.1e1 / t227 ^ 2;
t205 = t217 * t225 + 0.1e1;
t203 = 0.1e1 / t205;
t224 = 0.1e1 / t227;
t292 = t219 * t228;
t268 = -t221 * t224 - t225 * t292;
t322 = t268 * t203;
t178 = -t201 * t221 + t202 * t228 + t273 * t322;
t314 = cos(qJ(1));
t223 = t228 * t314;
t339 = t178 * t223;
t193 = t201 * t219 + t202 * t227;
t190 = 0.1e1 / t193;
t191 = 0.1e1 / t193 ^ 2;
t222 = t227 * t314;
t218 = t223 ^ 2;
t187 = t218 * t191 + 0.1e1;
t329 = qJD(2) - qJD(4);
t198 = -qJD(1) * t219 + t329 * t222;
t299 = t198 * t191;
t199 = -qJD(1) * t223 - t329 * t221;
t216 = t329 * t228;
t293 = t219 * t225;
t270 = -t199 * t224 + t216 * t293;
t181 = t270 * t203;
t176 = t273 * t181 - t201 * t199 - t202 * t216;
t308 = t176 * t190 * t191;
t290 = 0.2e1 * (-t218 * t308 + t223 * t299) / t187 ^ 2;
t338 = (t190 * t222 + t191 * t339) * t290;
t249 = qJ(5) + qJ(6);
t246 = sin(t249);
t247 = cos(t249);
t248 = qJD(5) + qJD(6);
t259 = qJD(1) * t314 + t222 * t248;
t197 = t221 * qJD(1) + t329 * t223;
t274 = -t311 * t248 - t197;
t188 = t274 * t246 + t259 * t247;
t189 = -t259 * t246 + t274 * t247;
t213 = t222 * t246 + t311 * t247;
t206 = t213 ^ 2;
t214 = t222 * t247 - t311 * t246;
t208 = 0.1e1 / t214 ^ 2;
t196 = t206 * t208 + 0.1e1;
t194 = 0.1e1 / t196;
t207 = 0.1e1 / t214;
t296 = t208 * t213;
t269 = -t246 * t207 + t247 * t296;
t304 = t189 * t207 * t208;
t316 = 0.2e1 * t213;
t283 = t304 * t316;
t337 = (t269 * t198 - (t208 * (-t188 * t247 + (t213 * t248 - t189) * t246) + (t207 * t248 + t283) * t247) * t223) * t194;
t336 = -t222 * t176 + t178 * t198;
t289 = 0.2e1 * t308;
t335 = -t197 * t190 - t289 * t339;
t215 = t329 * t227;
t333 = (t227 * t322 + t221) * t181 + t322 * t199 - t215;
t200 = t222 * qJD(1) - t329 * t219;
t332 = -(-t219 * t322 - t228) * t181 - t322 * t216 + t200;
t323 = t216 * t225;
t295 = t224 * t323;
t330 = t203 * ((t199 * t228 - t215 * t219 - t216 * t221) * t225 - t200 * t224 - 0.2e1 * t292 * t295);
t294 = t219 * t224;
t319 = (-t202 * t294 + t201) * t203 - t201;
t315 = -0.2e1 * t223;
t307 = (t188 * t296 - t206 * t304) / t196 ^ 2;
t306 = (-t199 * t293 + t217 * t295) / t205 ^ 2;
t301 = t194 * t208;
t298 = t201 * t223;
t297 = t202 * t223;
t288 = -0.2e1 * t307;
t287 = -0.2e1 * t306;
t286 = t208 * t307;
t285 = t224 * t306;
t284 = t188 * t301;
t275 = -t314 * t248 - t200;
t258 = qJD(1) * t311 + t221 * t248;
t212 = -t221 * t247 - t314 * t246;
t185 = 0.1e1 / t187;
t180 = t319 * t223;
t174 = t268 * t287 + t330;
t173 = 0.2e1 * t268 * t306 - t330;
t172 = t288 + (t284 + (-t194 * t304 - t286) * t213) * t316;
t1 = [t285 * t315 + (t198 * t224 + t223 * t323) * t203, t173, 0, t174, 0, 0; -t219 * t190 * t290 + (-t199 * t190 + (-t176 * t219 - t180 * t198) * t191) * t185 - ((-t180 * t289 + t319 * t299) * t185 + (-t180 * t290 + ((t181 * t203 * t294 + t287) * t298 + (0.2e1 * t219 * t285 - t181 + (t181 - t270) * t203) * t297) * t185) * t191) * t223, t338 + (((t173 * t219 + t333) * t297 + (-t173 * t227 + t332) * t298 - t336) * t191 - t335) * t185, 0, -t338 + (((t174 * t219 - t333) * t297 + (-t174 * t227 - t332) * t298 + t336) * t191 + t335) * t185, 0, 0; (t286 * t316 - t284) * t212 + (-t189 * t301 + t207 * t288) * (-t221 * t246 + t314 * t247) + ((t275 * t246 - t258 * t247) * t207 - (t258 * t246 + t275 * t247) * t296 + t212 * t283) * t194, t269 * t307 * t315 + t337, 0, -t269 * t223 * t288 - t337, t172, t172;];
JaD_rot  = t1;
