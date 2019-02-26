% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:57
% EndTime: 2019-02-26 19:40:58
% DurationCPUTime: 1.37s
% Computational Cost: add. (7669->100), mult. (22804->213), div. (524->12), fcn. (30024->15), ass. (0->101)
t240 = cos(pkin(6));
t237 = sin(pkin(12));
t286 = sin(pkin(11));
t265 = t286 * t237;
t239 = cos(pkin(12));
t288 = cos(pkin(11));
t267 = t288 * t239;
t235 = -t240 * t265 + t267;
t242 = sin(qJ(3));
t244 = cos(qJ(3));
t264 = t286 * t239;
t268 = t288 * t237;
t256 = t240 * t264 + t268;
t238 = sin(pkin(6));
t270 = t238 * t286;
t287 = sin(pkin(7));
t289 = cos(pkin(7));
t291 = t256 * t289 - t287 * t270;
t218 = t235 * t242 + t291 * t244;
t234 = t240 * t268 + t264;
t255 = -t240 * t267 + t265;
t271 = t238 * t287;
t252 = -t255 * t289 - t271 * t288;
t217 = t234 * t244 + t242 * t252;
t241 = sin(qJ(4));
t243 = cos(qJ(4));
t251 = -t238 * t288 * t289 + t255 * t287;
t206 = t217 * t243 + t241 * t251;
t216 = -t234 * t242 + t244 * t252;
t209 = t216 * qJD(3);
t186 = qJD(4) * t206 + t209 * t241;
t204 = t217 * t241 - t243 * t251;
t201 = t204 ^ 2;
t266 = t287 * t240;
t269 = t239 * t289;
t230 = t238 * (t237 * t244 + t242 * t269) + t242 * t266;
t233 = -t239 * t271 + t240 * t289;
t223 = t230 * t241 - t233 * t243;
t221 = 0.1e1 / t223 ^ 2;
t197 = t201 * t221 + 0.1e1;
t195 = 0.1e1 / t197;
t224 = t230 * t243 + t233 * t241;
t229 = t244 * t266 + (-t237 * t242 + t244 * t269) * t238;
t225 = t229 * qJD(3);
t199 = qJD(4) * t224 + t225 * t241;
t220 = 0.1e1 / t223;
t280 = t204 * t221;
t174 = (-t186 * t220 + t199 * t280) * t195;
t198 = atan2(-t204, t223);
t192 = sin(t198);
t193 = cos(t198);
t261 = -t192 * t223 - t193 * t204;
t171 = t174 * t261 - t192 * t186 + t193 * t199;
t185 = -t192 * t204 + t193 * t223;
t182 = 0.1e1 / t185;
t183 = 0.1e1 / t185 ^ 2;
t294 = t171 * t182 * t183;
t219 = t235 * t244 - t242 * t291;
t231 = t256 * t287 + t270 * t289;
t207 = t219 * t241 - t231 * t243;
t293 = 0.2e1 * t207 * t294;
t258 = -t216 * t220 + t229 * t280;
t292 = t241 * t258;
t281 = t199 * t220 * t221;
t290 = -0.2e1 * (t186 * t280 - t201 * t281) / t197 ^ 2;
t213 = 0.1e1 / t218;
t214 = 0.1e1 / t218 ^ 2;
t285 = t183 * t207;
t208 = t219 * t243 + t231 * t241;
t211 = t218 * qJD(3);
t188 = qJD(4) * t208 - t211 * t241;
t284 = t188 * t183;
t283 = t192 * t207;
t282 = t193 * t207;
t279 = t208 * t219;
t278 = t218 * t241;
t275 = qJD(4) * t243;
t202 = t207 ^ 2;
t180 = t202 * t183 + 0.1e1;
t274 = 0.2e1 * (-t202 * t294 + t207 * t284) / t180 ^ 2;
t189 = -qJD(4) * t207 - t211 * t243;
t203 = t208 ^ 2;
t194 = t203 * t214 + 0.1e1;
t212 = t219 * qJD(3);
t215 = t213 * t214;
t273 = 0.2e1 * (t208 * t214 * t189 - t203 * t215 * t212) / t194 ^ 2;
t263 = -0.2e1 * t204 * t281;
t259 = -t206 * t220 + t224 * t280;
t226 = t230 * qJD(3);
t210 = t217 * qJD(3);
t200 = -qJD(4) * t223 + t225 * t243;
t190 = 0.1e1 / t194;
t187 = -qJD(4) * t204 + t209 * t243;
t178 = 0.1e1 / t180;
t176 = t195 * t292;
t175 = t259 * t195;
t173 = (-t192 * t216 + t193 * t229) * t241 + t261 * t176;
t172 = t175 * t261 - t192 * t206 + t193 * t224;
t169 = t259 * t290 + (t224 * t263 - t187 * t220 + (t186 * t224 + t199 * t206 + t200 * t204) * t221) * t195;
t168 = t290 * t292 + (t258 * t275 + (t229 * t263 + t210 * t220 + (t186 * t229 + t199 * t216 - t204 * t226) * t221) * t241) * t195;
t1 = [0, 0, t168, t169, 0, 0; 0, 0 (t173 * t285 + t182 * t278) * t274 + ((-t212 * t241 - t218 * t275) * t182 + (-t284 + t293) * t173 + (t278 * t171 - (t229 * t275 - t168 * t204 - t176 * t186 - t226 * t241 + (-t176 * t223 - t216 * t241) * t174) * t282 - (-t216 * t275 - t168 * t223 - t176 * t199 + t210 * t241 + (t176 * t204 - t229 * t241) * t174) * t283) * t183) * t178 (t172 * t285 - t182 * t208) * t274 + (t172 * t293 + t189 * t182 + (-t208 * t171 - t172 * t188 - (-t169 * t204 - t175 * t186 + t200 + (-t175 * t223 - t206) * t174) * t282 - (-t169 * t223 - t175 * t199 - t187 + (t175 * t204 - t224) * t174) * t283) * t183) * t178, 0, 0; 0, 0 (t213 * t218 * t243 + t214 * t279) * t273 + (qJD(4) * t213 * t278 + (-t189 * t219 + t208 * t211) * t214 + (0.2e1 * t215 * t279 + (t214 * t218 - t213) * t243) * t212) * t190, t207 * t213 * t273 + (t207 * t212 * t214 - t188 * t213) * t190, 0, 0;];
JaD_rot  = t1;
