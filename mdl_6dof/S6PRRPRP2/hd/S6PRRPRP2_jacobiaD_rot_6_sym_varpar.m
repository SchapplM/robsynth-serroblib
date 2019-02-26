% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:48
% EndTime: 2019-02-26 20:01:50
% DurationCPUTime: 1.72s
% Computational Cost: add. (11602->157), mult. (21125->305), div. (784->12), fcn. (27225->13), ass. (0->124)
t264 = qJ(3) + pkin(11);
t262 = sin(t264);
t263 = cos(t264);
t267 = sin(qJ(2));
t269 = cos(qJ(2));
t326 = cos(pkin(10));
t327 = cos(pkin(6));
t290 = t327 * t326;
t325 = sin(pkin(10));
t278 = -t267 * t290 - t269 * t325;
t265 = sin(pkin(6));
t295 = t265 * t326;
t238 = t262 * t278 - t263 * t295;
t288 = t267 * t325 - t269 * t290;
t254 = t288 * qJD(2);
t219 = qJD(3) * t238 - t254 * t263;
t268 = cos(qJ(5));
t280 = t262 * t295 + t263 * t278;
t266 = sin(qJ(5));
t282 = t288 * t266;
t229 = -t268 * t280 + t282;
t276 = t278 * qJD(2);
t201 = qJD(5) * t229 + t219 * t266 + t268 * t276;
t281 = t288 * t268;
t227 = -t266 * t280 - t281;
t222 = t227 ^ 2;
t311 = t265 * t267;
t250 = t262 * t327 + t263 * t311;
t309 = t268 * t269;
t245 = t250 * t266 + t265 * t309;
t243 = 0.1e1 / t245 ^ 2;
t213 = t222 * t243 + 0.1e1;
t211 = 0.1e1 / t213;
t249 = -t262 * t311 + t263 * t327;
t307 = qJD(2) * t265;
t296 = t269 * t307;
t236 = qJD(3) * t249 + t263 * t296;
t310 = t266 * t269;
t246 = t250 * t268 - t265 * t310;
t297 = t267 * t307;
t215 = qJD(5) * t246 + t236 * t266 - t268 * t297;
t242 = 0.1e1 / t245;
t316 = t227 * t243;
t188 = (-t201 * t242 + t215 * t316) * t211;
t214 = atan2(-t227, t245);
t208 = sin(t214);
t209 = cos(t214);
t287 = -t208 * t245 - t209 * t227;
t184 = t188 * t287 - t208 * t201 + t209 * t215;
t200 = -t208 * t227 + t209 * t245;
t197 = 0.1e1 / t200;
t198 = 0.1e1 / t200 ^ 2;
t331 = t184 * t197 * t198;
t289 = t327 * t325;
t277 = -t267 * t326 - t269 * t289;
t259 = -t267 * t289 + t326 * t269;
t294 = t265 * t325;
t279 = -t259 * t263 - t262 * t294;
t230 = -t266 * t279 + t268 * t277;
t330 = -0.2e1 * t230;
t291 = 0.2e1 * t230 * t331;
t283 = -t238 * t242 + t249 * t316;
t329 = t266 * t283;
t318 = t215 * t242 * t243;
t328 = -0.2e1 * (t201 * t316 - t222 * t318) / t213 ^ 2;
t313 = t277 * t266;
t231 = -t268 * t279 - t313;
t224 = 0.1e1 / t231;
t225 = 0.1e1 / t231 ^ 2;
t240 = -t259 * t262 + t263 * t294;
t237 = t240 ^ 2;
t315 = t237 * t225;
t210 = 0.1e1 + t315;
t255 = t277 * qJD(2);
t220 = qJD(3) * t279 - t255 * t262;
t221 = qJD(3) * t240 + t255 * t263;
t256 = t259 * qJD(2);
t204 = -qJD(5) * t230 + t221 * t268 + t256 * t266;
t321 = t204 * t224 * t225;
t299 = t237 * t321;
t317 = t225 * t240;
t324 = (t220 * t317 - t299) / t210 ^ 2;
t323 = t198 * t230;
t203 = qJD(5) * t231 + t221 * t266 - t256 * t268;
t322 = t203 * t198;
t320 = t208 * t230;
t319 = t209 * t230;
t314 = t240 * t266;
t312 = t263 * t268;
t308 = qJD(2) * t263;
t306 = qJD(3) * t262;
t305 = qJD(5) * t266;
t304 = qJD(5) * t268;
t303 = t263 * qJD(5);
t223 = t230 ^ 2;
t196 = t223 * t198 + 0.1e1;
t302 = 0.2e1 * (-t223 * t331 + t230 * t322) / t196 ^ 2;
t301 = 0.2e1 * t324;
t298 = t240 * t321;
t292 = -0.2e1 * t227 * t318;
t285 = -t229 * t242 + t246 * t316;
t232 = -t263 * t282 + t268 * t278;
t247 = (t263 * t310 - t267 * t268) * t265;
t284 = -t232 * t242 + t247 * t316;
t235 = -qJD(3) * t250 - t262 * t296;
t234 = t259 * t266 + t277 * t312;
t233 = -t259 * t268 + t263 * t313;
t218 = qJD(3) * t280 + t254 * t262;
t217 = ((-qJD(2) + t303) * t309 + (-t269 * t306 + (qJD(5) - t308) * t267) * t266) * t265;
t216 = -qJD(5) * t245 + t236 * t268 + t266 * t297;
t206 = 0.1e1 / t210;
t205 = t282 * t306 + (-t288 * t303 + t254) * t268 + (t266 * t308 - t305) * t278;
t202 = qJD(5) * t281 + t219 * t268 - t266 * t276 + t280 * t305;
t194 = 0.1e1 / t196;
t193 = t211 * t329;
t192 = t284 * t211;
t190 = t285 * t211;
t187 = (-t208 * t238 + t209 * t249) * t266 + t287 * t193;
t186 = t192 * t287 - t208 * t232 + t209 * t247;
t185 = t190 * t287 - t208 * t229 + t209 * t246;
t183 = t284 * t328 + (t247 * t292 - t205 * t242 + (t201 * t247 + t215 * t232 + t217 * t227) * t243) * t211;
t181 = t285 * t328 + (t246 * t292 - t202 * t242 + (t201 * t246 + t215 * t229 + t216 * t227) * t243) * t211;
t180 = t328 * t329 + (t283 * t304 + (t249 * t292 - t218 * t242 + (t201 * t249 + t215 * t238 + t227 * t235) * t243) * t266) * t211;
t1 = [0, t183, t180, 0, t181, 0; 0 (t186 * t323 - t197 * t233) * t302 + (t186 * t291 + (-t233 * t184 - t186 * t203 - (-t183 * t227 - t192 * t201 + t217 + (-t192 * t245 - t232) * t188) * t319 - (-t183 * t245 - t192 * t215 - t205 + (t192 * t227 - t247) * t188) * t320) * t198 + ((t277 * t303 - t255) * t268 + (qJD(5) * t259 - t256 * t263 - t277 * t306) * t266) * t197) * t194 (t187 * t323 - t197 * t314) * t302 + ((t220 * t266 + t240 * t304) * t197 + (-t322 + t291) * t187 + (-t314 * t184 - (t249 * t304 - t180 * t227 - t193 * t201 + t235 * t266 + (-t193 * t245 - t238 * t266) * t188) * t319 - (-t238 * t304 - t180 * t245 - t193 * t215 - t218 * t266 + (t193 * t227 - t249 * t266) * t188) * t320) * t198) * t194, 0 (t185 * t323 - t197 * t231) * t302 + (t185 * t291 + t204 * t197 + (-t231 * t184 - t185 * t203 - (-t181 * t227 - t190 * t201 + t216 + (-t190 * t245 - t229) * t188) * t319 - (-t181 * t245 - t190 * t215 - t202 + (t190 * t227 - t246) * t188) * t320) * t198) * t194, 0; 0 (t224 * t262 * t277 + t234 * t317) * t301 + (0.2e1 * t234 * t298 + (-qJD(3) * t263 * t277 + t256 * t262) * t224 + (-(t255 * t266 - t256 * t312 + t259 * t304) * t240 - t234 * t220 - (-t262 * t204 - (t266 * t303 + t268 * t306) * t240) * t277) * t225) * t206 (-t224 * t279 + t268 * t315) * t301 + (0.2e1 * t268 * t299 - t221 * t224 + (-0.2e1 * t220 * t240 * t268 - t204 * t279 + t237 * t305) * t225) * t206, 0, t317 * t324 * t330 + (t298 * t330 + (t203 * t240 + t220 * t230) * t225) * t206, 0;];
JaD_rot  = t1;
