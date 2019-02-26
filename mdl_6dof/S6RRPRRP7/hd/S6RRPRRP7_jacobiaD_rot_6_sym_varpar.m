% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:25
% EndTime: 2019-02-26 21:49:28
% DurationCPUTime: 2.17s
% Computational Cost: add. (5860->130), mult. (17332->273), div. (1542->14), fcn. (21924->11), ass. (0->125)
t256 = sin(qJ(5));
t257 = sin(qJ(4));
t260 = cos(qJ(2));
t355 = sin(qJ(2));
t356 = cos(qJ(4));
t232 = t355 * t257 + t260 * t356;
t258 = sin(qJ(1));
t225 = t232 * t258;
t259 = cos(qJ(5));
t357 = cos(qJ(1));
t279 = -t225 * t256 + t357 * t259;
t233 = -t260 * t257 + t355 * t356;
t326 = t233 * t256;
t207 = atan2(t279, t326);
t204 = cos(t207);
t339 = t204 * t279;
t211 = t279 ^ 2;
t230 = 0.1e1 / t233 ^ 2;
t254 = 0.1e1 / t256 ^ 2;
t328 = t230 * t254;
t208 = t211 * t328 + 0.1e1;
t205 = 0.1e1 / t208;
t224 = t233 * t258;
t229 = 0.1e1 / t233;
t253 = 0.1e1 / t256;
t327 = t232 * t253;
t309 = t279 * t327;
t275 = t224 * t229 - t230 * t309;
t360 = t275 * t205;
t203 = sin(t207);
t298 = t233 * t360 - t224;
t373 = t298 * t203 - t204 * t232;
t177 = t256 * t373 - t360 * t339;
t342 = t203 * t279;
t191 = t204 * t326 + t342;
t189 = 0.1e1 / t191 ^ 2;
t226 = t232 * t357;
t288 = -t226 * t256 - t258 * t259;
t212 = t288 ^ 2;
t187 = t212 * t189 + 0.1e1;
t227 = t233 * t357;
t365 = qJD(2) - qJD(4);
t196 = qJD(1) * t225 + t365 * t227;
t302 = qJD(1) * t357;
t321 = qJD(5) * t256;
t192 = -t196 * t256 - t258 * t321 + (qJD(5) * t226 + t302) * t259;
t349 = t192 * t189;
t188 = 0.1e1 / t191;
t210 = t365 * t232;
t320 = qJD(5) * t259;
t277 = t210 * t256 + t233 * t320;
t310 = t279 * t328;
t209 = t365 * t233;
t199 = t226 * qJD(1) - t209 * t258;
t218 = t225 * t259 + t357 * t256;
t323 = qJD(1) * t258;
t194 = t218 * qJD(5) + t199 * t256 + t259 * t323;
t329 = t229 * t253;
t312 = t194 * t329;
t179 = (-t277 * t310 - t312) * t205;
t270 = t179 * t279 + t277;
t174 = (-t179 * t326 - t194) * t203 + t270 * t204;
t372 = t174 * t189;
t368 = t188 * t372;
t354 = (-t212 * t368 - t288 * t349) / t187 ^ 2;
t319 = 0.2e1 * t354;
t350 = t189 * t288;
t291 = t319 * t350;
t381 = t177 * t291;
t358 = -0.2e1 * t288;
t283 = t358 * t368 - t349;
t380 = -t283 * t177 - (t227 * t188 + t373 * t350) * t320;
t379 = -0.2e1 * t259;
t338 = t204 * t288;
t378 = -t227 * t174 + (t298 * t179 + t209) * t338;
t197 = -t224 * qJD(1) + t365 * t226;
t377 = t197 * t379;
t198 = -t227 * qJD(1) - t365 * t225;
t304 = t254 * t320;
t336 = t210 * t230;
t335 = t229 * t336;
t370 = (t198 * t229 + t230 * (-(t209 * t253 + t232 * t304) * t279 - t194 * t327 + t210 * t224) - 0.2e1 * t309 * t335) * t205;
t367 = (-t179 * t342 - t194 * t204) * t360;
t364 = -t179 * t232 - t210 * t360 - t198;
t193 = t288 * qJD(5) - t196 * t259 - t256 * t302;
t222 = t226 * t259 - t258 * t256;
t214 = 0.1e1 / t222 ^ 2;
t361 = t193 * t214;
t195 = t279 * qJD(5) + t199 * t259 - t256 * t323;
t213 = 0.1e1 / t222;
t255 = t253 * t254;
t352 = (-t194 * t310 + (-t230 * t255 * t320 - t254 * t335) * t211) / t208 ^ 2;
t223 = t227 ^ 2;
t333 = t223 * t214;
t202 = 0.1e1 + t333;
t348 = t213 * t361;
t308 = t223 * t348;
t334 = t214 * t227;
t351 = (t197 * t334 - t308) / t202 ^ 2;
t347 = t193 * t226;
t346 = t196 * t213;
t345 = t197 * t188;
t341 = t203 * t288;
t325 = t254 * t259;
t318 = -0.2e1 * t352;
t316 = t188 * t354;
t315 = t229 * t352;
t314 = t213 * t351;
t313 = t227 * t348;
t311 = t279 * t329;
t296 = t253 * t315;
t295 = t334 * t351;
t294 = t227 * t316;
t290 = 0.2e1 * t295;
t280 = -t218 * t253 - t279 * t325;
t272 = t229 * t304 + t253 * t336;
t271 = -t203 + (-t204 * t311 + t203) * t205;
t200 = 0.1e1 / t202;
t185 = 0.1e1 / t187;
t182 = t280 * t229 * t205;
t175 = t204 * t233 * t259 - t203 * t218 + (-t203 * t326 + t339) * t182;
t173 = 0.2e1 * t275 * t352 + t370;
t172 = t275 * t318 - t370;
t170 = -0.2e1 * t280 * t315 + (-t280 * t336 + (t194 * t325 - t195 * t253 + (t218 * t325 - (-0.2e1 * t255 * t259 ^ 2 - t253) * t279) * qJD(5)) * t229) * t205;
t1 = [t296 * t358 + (-t192 * t329 - t272 * t288) * t205, t172, 0, t173, t170, 0; -0.2e1 * t279 * t316 + (-t194 * t188 - t279 * t372 + (t271 * t192 - ((t179 * t205 * t311 + t318) * t203 + (0.2e1 * t279 * t296 - t179 + (t272 * t279 + t179 + t312) * t205) * t204) * t288) * t350) * t185 - (t185 * t283 - t291) * t271 * t288, t381 + ((t172 * t339 + t367) * t350 + t380) * t185 + (0.2e1 * t294 + (-t345 + ((-t172 * t233 + t364) * t341 - t378) * t189) * t185) * t256, 0, -t381 + ((t173 * t339 - t367) * t350 - t380) * t185 + (-0.2e1 * t294 + (t345 + ((-t173 * t233 - t364) * t341 + t378) * t189) * t185) * t256 (-t175 * t350 - t188 * t222) * t319 + (t193 * t188 + t283 * t175 + (-t222 * t174 + (-t233 * t321 + t170 * t279 - t182 * t194 + t210 * t259 + (-t182 * t326 - t218) * t179) * t338 + (-t195 + (-t170 * t256 - t179 * t259) * t233 - t270 * t182) * t341) * t189) * t185, 0; -t218 * t290 + 0.2e1 * t224 * t314 + (-(-t197 * t214 + 0.2e1 * t313) * t218 + t195 * t334 + t198 * t213 + t224 * t361) * t200, -0.2e1 * t226 * t314 + (-t214 * t347 - t346) * t200 - (t259 * t290 + (t214 * t377 - (-t214 * t321 + t348 * t379) * t227) * t200) * t227, 0, 0.2e1 * (t213 * t226 + t259 * t333) * t351 + (0.2e1 * t259 * t308 + t346 + (t223 * t321 + t227 * t377 + t347) * t214) * t200, -t295 * t358 + (-t313 * t358 + (t192 * t227 - t197 * t288) * t214) * t200, 0;];
JaD_rot  = t1;
