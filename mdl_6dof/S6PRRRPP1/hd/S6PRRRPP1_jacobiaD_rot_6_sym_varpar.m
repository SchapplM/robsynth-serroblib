% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:45
% EndTime: 2019-02-26 20:08:47
% DurationCPUTime: 1.69s
% Computational Cost: add. (11295->157), mult. (21125->308), div. (784->12), fcn. (27225->13), ass. (0->124)
t261 = qJ(4) + pkin(11);
t259 = sin(t261);
t260 = cos(t261);
t263 = sin(qJ(3));
t265 = cos(qJ(3));
t264 = sin(qJ(2));
t266 = cos(qJ(2));
t323 = cos(pkin(10));
t324 = cos(pkin(6));
t287 = t324 * t323;
t322 = sin(pkin(10));
t275 = -t264 * t287 - t266 * t322;
t262 = sin(pkin(6));
t292 = t262 * t323;
t277 = t263 * t292 + t265 * t275;
t285 = t264 * t322 - t266 * t287;
t222 = t259 * t285 - t260 * t277;
t238 = t263 * t275 - t265 * t292;
t249 = t285 * qJD(2);
t226 = qJD(3) * t238 - t249 * t265;
t273 = qJD(2) * t275;
t198 = qJD(4) * t222 + t226 * t259 + t260 * t273;
t278 = t285 * t260;
t220 = -t259 * t277 - t278;
t215 = t220 ^ 2;
t307 = t262 * t264;
t256 = t263 * t324 + t265 * t307;
t306 = t262 * t266;
t235 = t256 * t259 + t260 * t306;
t233 = 0.1e1 / t235 ^ 2;
t207 = t215 * t233 + 0.1e1;
t205 = 0.1e1 / t207;
t236 = t256 * t260 - t259 * t306;
t255 = -t263 * t307 + t265 * t324;
t305 = qJD(2) * t262;
t293 = t266 * t305;
t243 = qJD(3) * t255 + t265 * t293;
t294 = t264 * t305;
t212 = qJD(4) * t236 + t243 * t259 - t260 * t294;
t232 = 0.1e1 / t235;
t312 = t220 * t233;
t185 = (-t198 * t232 + t212 * t312) * t205;
t208 = atan2(-t220, t235);
t203 = sin(t208);
t204 = cos(t208);
t284 = -t203 * t235 - t204 * t220;
t181 = t185 * t284 - t198 * t203 + t204 * t212;
t197 = -t203 * t220 + t204 * t235;
t194 = 0.1e1 / t197;
t195 = 0.1e1 / t197 ^ 2;
t328 = t181 * t194 * t195;
t286 = t324 * t322;
t274 = -t264 * t323 - t266 * t286;
t254 = -t264 * t286 + t323 * t266;
t291 = t262 * t322;
t276 = -t254 * t265 - t263 * t291;
t223 = -t259 * t276 + t260 * t274;
t327 = -0.2e1 * t223;
t288 = 0.2e1 * t223 * t328;
t280 = -t232 * t238 + t255 * t312;
t326 = t259 * t280;
t315 = t212 * t232 * t233;
t325 = -0.2e1 * (t198 * t312 - t215 * t315) / t207 ^ 2;
t224 = -t259 * t274 - t260 * t276;
t217 = 0.1e1 / t224;
t218 = 0.1e1 / t224 ^ 2;
t240 = -t254 * t263 + t265 * t291;
t237 = t240 ^ 2;
t314 = t218 * t237;
t211 = 0.1e1 + t314;
t250 = t274 * qJD(2);
t228 = qJD(3) * t240 + t250 * t265;
t251 = t254 * qJD(2);
t201 = -qJD(4) * t223 + t228 * t260 + t251 * t259;
t318 = t201 * t217 * t218;
t296 = t237 * t318;
t227 = qJD(3) * t276 - t250 * t263;
t311 = t227 * t240;
t321 = (t218 * t311 - t296) / t211 ^ 2;
t320 = t195 * t223;
t200 = qJD(4) * t224 + t228 * t259 - t251 * t260;
t319 = t200 * t195;
t317 = t203 * t223;
t316 = t204 * t223;
t313 = t218 * t240;
t310 = t240 * t259;
t309 = t259 * t265;
t308 = t260 * t265;
t304 = qJD(2) * t265;
t303 = qJD(3) * t263;
t302 = qJD(4) * t259;
t301 = qJD(4) * t260;
t300 = qJD(4) * t265;
t216 = t223 ^ 2;
t193 = t195 * t216 + 0.1e1;
t299 = 0.2e1 * (-t216 * t328 + t223 * t319) / t193 ^ 2;
t298 = 0.2e1 * t321;
t295 = t240 * t318;
t289 = -0.2e1 * t220 * t315;
t282 = -t222 * t232 + t236 * t312;
t279 = t265 * t285;
t229 = -t259 * t279 + t260 * t275;
t244 = (-t260 * t264 + t266 * t309) * t262;
t281 = -t229 * t232 + t244 * t312;
t242 = -qJD(3) * t256 - t263 * t293;
t231 = t254 * t259 + t274 * t308;
t230 = -t254 * t260 + t274 * t309;
t225 = qJD(3) * t277 + t249 * t263;
t214 = ((-qJD(2) + t300) * t266 * t260 + (-t266 * t303 + (qJD(4) - t304) * t264) * t259) * t262;
t213 = -qJD(4) * t235 + t243 * t260 + t259 * t294;
t209 = 0.1e1 / t211;
t202 = t249 * t260 - t275 * t302 - t279 * t301 + (t275 * t304 + t285 * t303) * t259;
t199 = qJD(4) * t278 + t226 * t260 - t259 * t273 + t277 * t302;
t191 = 0.1e1 / t193;
t190 = t205 * t326;
t188 = t281 * t205;
t187 = t282 * t205;
t184 = (-t203 * t238 + t204 * t255) * t259 + t284 * t190;
t183 = t188 * t284 - t203 * t229 + t204 * t244;
t182 = t187 * t284 - t203 * t222 + t204 * t236;
t180 = t281 * t325 + (t244 * t289 - t202 * t232 + (t198 * t244 + t212 * t229 + t214 * t220) * t233) * t205;
t178 = t282 * t325 + (t236 * t289 - t199 * t232 + (t198 * t236 + t212 * t222 + t213 * t220) * t233) * t205;
t177 = t325 * t326 + (t280 * t301 + (t255 * t289 - t225 * t232 + (t198 * t255 + t212 * t238 + t220 * t242) * t233) * t259) * t205;
t1 = [0, t180, t177, t178, 0, 0; 0 (t183 * t320 - t194 * t230) * t299 + (t183 * t288 + (-t230 * t181 - t183 * t200 - (-t180 * t220 - t188 * t198 + t214 + (-t188 * t235 - t229) * t185) * t316 - (-t180 * t235 - t188 * t212 - t202 + (t188 * t220 - t244) * t185) * t317) * t195 + ((t274 * t300 - t250) * t260 + (qJD(4) * t254 - t251 * t265 - t274 * t303) * t259) * t194) * t191 (t184 * t320 - t194 * t310) * t299 + ((t227 * t259 + t240 * t301) * t194 + (-t319 + t288) * t184 + (-t310 * t181 - (t255 * t301 - t177 * t220 - t190 * t198 + t242 * t259 + (-t190 * t235 - t238 * t259) * t185) * t316 - (-t238 * t301 - t177 * t235 - t190 * t212 - t225 * t259 + (t190 * t220 - t255 * t259) * t185) * t317) * t195) * t191 (t182 * t320 - t194 * t224) * t299 + (t182 * t288 + t201 * t194 + (-t224 * t181 - t182 * t200 - (-t178 * t220 - t187 * t198 + t213 + (-t187 * t235 - t222) * t185) * t316 - (-t178 * t235 - t187 * t212 - t199 + (t187 * t220 - t236) * t185) * t317) * t195) * t191, 0, 0; 0 (t217 * t263 * t274 + t231 * t313) * t298 + (0.2e1 * t231 * t295 + (-qJD(3) * t265 * t274 + t251 * t263) * t217 + (-(t250 * t259 - t251 * t308 + t254 * t301) * t240 - t231 * t227 - (-t263 * t201 - (t259 * t300 + t260 * t303) * t240) * t274) * t218) * t209 (-t217 * t276 + t260 * t314) * t298 + (0.2e1 * t260 * t296 - t217 * t228 + (-t201 * t276 + t237 * t302 - 0.2e1 * t260 * t311) * t218) * t209, t313 * t321 * t327 + (t295 * t327 + (t200 * t240 + t223 * t227) * t218) * t209, 0, 0;];
JaD_rot  = t1;
