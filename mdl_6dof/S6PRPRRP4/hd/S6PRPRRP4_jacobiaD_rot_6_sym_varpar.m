% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:01
% EndTime: 2019-02-26 19:52:03
% DurationCPUTime: 1.69s
% Computational Cost: add. (11602->157), mult. (21125->306), div. (784->12), fcn. (27225->13), ass. (0->123)
t261 = pkin(11) + qJ(4);
t259 = sin(t261);
t260 = cos(t261);
t264 = sin(qJ(2));
t266 = cos(qJ(2));
t322 = cos(pkin(10));
t323 = cos(pkin(6));
t287 = t323 * t322;
t321 = sin(pkin(10));
t275 = -t264 * t287 - t321 * t266;
t262 = sin(pkin(6));
t292 = t262 * t322;
t235 = t259 * t275 - t260 * t292;
t285 = t321 * t264 - t266 * t287;
t251 = t285 * qJD(2);
t216 = t235 * qJD(4) - t251 * t260;
t265 = cos(qJ(5));
t277 = t259 * t292 + t260 * t275;
t263 = sin(qJ(5));
t279 = t285 * t263;
t226 = -t265 * t277 + t279;
t273 = qJD(2) * t275;
t198 = t226 * qJD(5) + t216 * t263 + t265 * t273;
t278 = t285 * t265;
t224 = -t263 * t277 - t278;
t219 = t224 ^ 2;
t308 = t262 * t264;
t247 = t323 * t259 + t260 * t308;
t306 = t265 * t266;
t242 = t247 * t263 + t262 * t306;
t240 = 0.1e1 / t242 ^ 2;
t210 = t219 * t240 + 0.1e1;
t208 = 0.1e1 / t210;
t246 = -t259 * t308 + t323 * t260;
t304 = qJD(2) * t262;
t293 = t266 * t304;
t233 = t246 * qJD(4) + t260 * t293;
t307 = t263 * t266;
t243 = t247 * t265 - t262 * t307;
t294 = t264 * t304;
t212 = t243 * qJD(5) + t233 * t263 - t265 * t294;
t239 = 0.1e1 / t242;
t312 = t224 * t240;
t185 = (-t198 * t239 + t212 * t312) * t208;
t211 = atan2(-t224, t242);
t205 = sin(t211);
t206 = cos(t211);
t284 = -t205 * t242 - t206 * t224;
t181 = t284 * t185 - t205 * t198 + t206 * t212;
t197 = -t205 * t224 + t206 * t242;
t194 = 0.1e1 / t197;
t195 = 0.1e1 / t197 ^ 2;
t327 = t181 * t194 * t195;
t286 = t323 * t321;
t274 = -t322 * t264 - t266 * t286;
t256 = -t264 * t286 + t322 * t266;
t291 = t262 * t321;
t276 = -t256 * t260 - t259 * t291;
t227 = -t263 * t276 + t265 * t274;
t326 = -0.2e1 * t227;
t289 = 0.2e1 * t227 * t327;
t280 = -t235 * t239 + t246 * t312;
t325 = t263 * t280;
t315 = t212 * t239 * t240;
t324 = -0.2e1 * (t198 * t312 - t219 * t315) / t210 ^ 2;
t310 = t274 * t263;
t228 = -t265 * t276 - t310;
t221 = 0.1e1 / t228;
t222 = 0.1e1 / t228 ^ 2;
t237 = -t256 * t259 + t260 * t291;
t234 = t237 ^ 2;
t314 = t222 * t234;
t207 = 0.1e1 + t314;
t252 = t274 * qJD(2);
t217 = t276 * qJD(4) - t252 * t259;
t218 = t237 * qJD(4) + t252 * t260;
t253 = t256 * qJD(2);
t201 = -t227 * qJD(5) + t218 * t265 + t253 * t263;
t318 = t201 * t221 * t222;
t296 = t234 * t318;
t313 = t222 * t237;
t320 = (t217 * t313 - t296) / t207 ^ 2;
t319 = t195 * t227;
t317 = t205 * t227;
t316 = t206 * t227;
t311 = t237 * t263;
t309 = t260 * t265;
t305 = qJD(2) * t260;
t303 = qJD(4) * t259;
t302 = qJD(5) * t260;
t301 = qJD(5) * t263;
t300 = qJD(5) * t265;
t220 = t227 ^ 2;
t193 = t195 * t220 + 0.1e1;
t200 = t228 * qJD(5) + t218 * t263 - t253 * t265;
t299 = 0.2e1 * (t200 * t319 - t220 * t327) / t193 ^ 2;
t298 = 0.2e1 * t320;
t295 = t237 * t318;
t288 = -0.2e1 * t224 * t315;
t282 = -t226 * t239 + t243 * t312;
t229 = -t260 * t279 + t265 * t275;
t244 = (t260 * t307 - t264 * t265) * t262;
t281 = -t229 * t239 + t244 * t312;
t232 = -t247 * qJD(4) - t259 * t293;
t231 = t256 * t263 + t274 * t309;
t230 = -t256 * t265 + t260 * t310;
t215 = t277 * qJD(4) + t251 * t259;
t214 = ((-qJD(2) + t302) * t306 + (-t266 * t303 + (qJD(5) - t305) * t264) * t263) * t262;
t213 = -t242 * qJD(5) + t233 * t265 + t263 * t294;
t203 = 0.1e1 / t207;
t202 = -t285 * t260 * t300 + t251 * t265 + t279 * t303 + (t263 * t305 - t301) * t275;
t199 = qJD(5) * t278 + t216 * t265 - t263 * t273 + t277 * t301;
t191 = 0.1e1 / t193;
t190 = t208 * t325;
t189 = t281 * t208;
t187 = t282 * t208;
t184 = (-t205 * t235 + t206 * t246) * t263 + t284 * t190;
t183 = t284 * t189 - t205 * t229 + t206 * t244;
t182 = t284 * t187 - t205 * t226 + t206 * t243;
t180 = t281 * t324 + (t244 * t288 - t202 * t239 + (t198 * t244 + t212 * t229 + t214 * t224) * t240) * t208;
t178 = t282 * t324 + (t243 * t288 - t199 * t239 + (t198 * t243 + t212 * t226 + t213 * t224) * t240) * t208;
t177 = t324 * t325 + (t280 * t300 + (t246 * t288 - t215 * t239 + (t198 * t246 + t212 * t235 + t224 * t232) * t240) * t263) * t208;
t1 = [0, t180, 0, t177, t178, 0; 0 (t183 * t319 - t194 * t230) * t299 + (t183 * t289 + (-t230 * t181 - t183 * t200 - (-t180 * t224 - t189 * t198 + t214 + (-t189 * t242 - t229) * t185) * t316 - (-t180 * t242 - t189 * t212 - t202 + (t189 * t224 - t244) * t185) * t317) * t195 + ((t274 * t302 - t252) * t265 + (qJD(5) * t256 - t253 * t260 - t274 * t303) * t263) * t194) * t191, 0 (t184 * t319 - t194 * t311) * t299 + ((t217 * t263 + t237 * t300) * t194 + t184 * t289 + (-t184 * t200 - t311 * t181 - (t246 * t300 - t177 * t224 - t190 * t198 + t232 * t263 + (-t190 * t242 - t235 * t263) * t185) * t316 - (-t235 * t300 - t177 * t242 - t190 * t212 - t215 * t263 + (t190 * t224 - t246 * t263) * t185) * t317) * t195) * t191 (t182 * t319 - t194 * t228) * t299 + (t182 * t289 + t201 * t194 + (-t228 * t181 - t182 * t200 - (-t178 * t224 - t187 * t198 + t213 + (-t187 * t242 - t226) * t185) * t316 - (-t178 * t242 - t187 * t212 - t199 + (t187 * t224 - t243) * t185) * t317) * t195) * t191, 0; 0 (t221 * t259 * t274 + t231 * t313) * t298 + (0.2e1 * t231 * t295 + (-qJD(4) * t260 * t274 + t253 * t259) * t221 + (-(t252 * t263 - t253 * t309 + t256 * t300) * t237 - t231 * t217 - (-t259 * t201 - (t260 * t301 + t265 * t303) * t237) * t274) * t222) * t203, 0 (-t221 * t276 + t265 * t314) * t298 + (0.2e1 * t265 * t296 - t218 * t221 + (-0.2e1 * t217 * t237 * t265 - t201 * t276 + t234 * t301) * t222) * t203, t313 * t320 * t326 + (t295 * t326 + (t200 * t237 + t217 * t227) * t222) * t203, 0;];
JaD_rot  = t1;
