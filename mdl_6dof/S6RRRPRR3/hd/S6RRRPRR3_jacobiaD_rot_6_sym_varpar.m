% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:08
% EndTime: 2019-02-26 22:17:09
% DurationCPUTime: 1.42s
% Computational Cost: add. (10722->106), mult. (13073->203), div. (974->12), fcn. (16295->11), ass. (0->103)
t307 = qJ(2) + qJ(3);
t297 = sin(t307);
t298 = cos(t307);
t326 = sin(qJ(5));
t328 = cos(qJ(5));
t239 = t297 * t328 - t298 * t326;
t327 = sin(qJ(1));
t230 = t239 * t327;
t238 = t297 * t326 + t298 * t328;
t219 = atan2(t230, t238);
t214 = sin(t219);
t215 = cos(t219);
t232 = t238 * t327;
t290 = -t214 * t238 + t215 * t230;
t228 = t230 ^ 2;
t236 = 0.1e1 / t238 ^ 2;
t218 = t228 * t236 + 0.1e1;
t216 = 0.1e1 / t218;
t235 = 0.1e1 / t238;
t309 = t230 * t239;
t282 = -t232 * t235 - t236 * t309;
t336 = t282 * t216;
t189 = -t214 * t232 + t215 * t239 + t290 * t336;
t329 = cos(qJ(1));
t234 = t239 * t329;
t354 = t189 * t234;
t202 = t214 * t230 + t215 * t238;
t199 = 0.1e1 / t202;
t200 = 0.1e1 / t202 ^ 2;
t233 = t238 * t329;
t229 = t234 ^ 2;
t198 = t229 * t200 + 0.1e1;
t306 = qJD(2) + qJD(3);
t278 = t306 * t298;
t266 = t329 * t278;
t277 = t306 * t297;
t267 = t326 * t277;
t206 = -t230 * qJD(1) - t233 * qJD(5) + t328 * t266 + t329 * t267;
t317 = t206 * t200;
t265 = t327 * t278;
t207 = -t234 * qJD(1) + t232 * qJD(5) - t328 * t265 - t327 * t267;
t268 = t328 * t277;
t213 = t239 * qJD(5) + t326 * t278 - t268;
t310 = t230 * t236;
t284 = -t207 * t235 - t213 * t310;
t192 = t284 * t216;
t187 = t290 * t192 - t214 * t207 + t215 * t213;
t325 = t187 * t199 * t200;
t305 = 0.2e1 * (-t229 * t325 + t234 * t317) / t198 ^ 2;
t353 = (t199 * t233 + t200 * t354) * t305;
t205 = t232 * qJD(1) - t234 * qJD(5) - t326 * t266 + t329 * t268;
t257 = sin(qJ(6));
t258 = cos(qJ(6));
t227 = t233 * t258 - t327 * t257;
t296 = qJD(1) * t329;
t203 = t227 * qJD(6) - t205 * t257 + t258 * t296;
t280 = -t233 * t257 - t327 * t258;
t344 = t280 * qJD(6);
t204 = -t205 * t258 - t257 * t296 + t344;
t220 = t280 ^ 2;
t222 = 0.1e1 / t227 ^ 2;
t211 = t220 * t222 + 0.1e1;
t209 = 0.1e1 / t211;
t221 = 0.1e1 / t227;
t312 = t222 * t280;
t283 = -t257 * t221 - t258 * t312;
t319 = t204 * t221 * t222;
t331 = -0.2e1 * t280;
t294 = t319 * t331;
t352 = (t283 * t206 - (((-t204 - t344) * t257 - t203 * t258) * t222 + (qJD(6) * t221 + t294) * t258) * t234) * t209;
t351 = -t233 * t187 + t189 * t206;
t304 = 0.2e1 * t325;
t350 = -t205 * t199 - t304 * t354;
t212 = t238 * qJD(5) - t328 * t278 - t267;
t348 = (t238 * t336 + t232) * t192 + t336 * t207 + t212;
t208 = t233 * qJD(1) + t230 * qJD(5) + t326 * t265 - t327 * t268;
t347 = -(-t230 * t336 - t239) * t192 + t336 * t213 + t208;
t338 = t213 * t236;
t315 = t235 * t338;
t345 = ((t207 * t239 + t212 * t230 + t213 * t232) * t236 - t208 * t235 + 0.2e1 * t309 * t315) * t216;
t311 = t230 * t235;
t334 = (-t215 * t311 + t214) * t216 - t214;
t330 = -0.2e1 * t234;
t324 = (-t203 * t312 - t220 * t319) / t211 ^ 2;
t323 = (-t207 * t310 - t228 * t315) / t218 ^ 2;
t316 = t209 * t222;
t314 = t214 * t234;
t313 = t215 * t234;
t303 = -0.2e1 * t324;
t302 = -0.2e1 * t323;
t301 = t222 * t324;
t300 = t235 * t323;
t299 = t203 * t316;
t295 = qJD(1) * t327;
t225 = -t232 * t258 - t329 * t257;
t281 = t232 * t257 - t329 * t258;
t196 = 0.1e1 / t198;
t191 = t334 * t234;
t185 = t282 * t302 + t345;
t184 = 0.2e1 * t282 * t323 - t345;
t183 = t283 * t324 * t330 + t352;
t182 = t353 + (((t184 * t230 + t348) * t313 + (-t184 * t238 + t347) * t314 - t351) * t200 - t350) * t196;
t1 = [t300 * t330 + (t206 * t235 - t234 * t338) * t216, t184, t184, 0, t185, 0; -t230 * t199 * t305 + (-t207 * t199 + (-t187 * t230 - t191 * t206) * t200) * t196 - ((-t191 * t304 + t334 * t317) * t196 + (-t191 * t305 + ((t192 * t216 * t311 + t302) * t314 + (0.2e1 * t230 * t300 - t192 + (t192 - t284) * t216) * t313) * t196) * t200) * t234, t182, t182, 0, -t353 + (((t185 * t230 - t348) * t313 + (-t185 * t238 - t347) * t314 + t351) * t200 + t350) * t196, 0; (t301 * t331 - t299) * t225 - (-t204 * t316 + t221 * t303) * t281 + ((t225 * qJD(6) - t208 * t257 - t258 * t295) * t221 + (t281 * qJD(6) - t208 * t258 + t257 * t295) * t312 + t225 * t294) * t209, t183, t183, 0, -t283 * t234 * t303 - t352, t303 + (t299 - (-t209 * t319 - t301) * t280) * t331;];
JaD_rot  = t1;
