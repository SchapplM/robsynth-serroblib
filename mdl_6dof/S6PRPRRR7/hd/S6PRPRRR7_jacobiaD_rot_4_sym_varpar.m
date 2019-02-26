% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR7_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiaD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:30
% EndTime: 2019-02-26 19:57:32
% DurationCPUTime: 1.02s
% Computational Cost: add. (4236->99), mult. (12947->219), div. (275->12), fcn. (16774->17), ass. (0->102)
t277 = sin(pkin(13));
t282 = cos(pkin(13));
t289 = cos(qJ(2));
t285 = cos(pkin(6));
t287 = sin(qJ(2));
t308 = t285 * t287;
t272 = t277 * t289 + t282 * t308;
t276 = sin(pkin(14));
t278 = sin(pkin(8));
t279 = sin(pkin(7));
t281 = cos(pkin(14));
t283 = cos(pkin(8));
t284 = cos(pkin(7));
t307 = t285 * t289;
t298 = -t277 * t287 + t282 * t307;
t280 = sin(pkin(6));
t311 = t280 * t284;
t314 = t279 * t280;
t242 = (-t272 * t276 + (-t282 * t314 + t298 * t284) * t281) * t278 - (-t298 * t279 - t282 * t311) * t283;
t309 = t284 * t281;
t299 = -t276 * t287 + t289 * t309;
t312 = t279 * t289;
t253 = -(t285 * t279 * t281 + t299 * t280) * t278 + (-t280 * t312 + t285 * t284) * t283;
t237 = atan2(t242, t253);
t232 = sin(t237);
t233 = cos(t237);
t219 = t232 * t242 + t233 * t253;
t216 = 0.1e1 / t219;
t296 = t277 * t308 - t282 * t289;
t297 = t277 * t307 + t282 * t287;
t300 = t277 * t314 - t284 * t297;
t255 = t300 * t276 - t281 * t296;
t286 = sin(qJ(4));
t288 = cos(qJ(4));
t254 = t276 * t296 + t300 * t281;
t267 = t277 * t311 + t279 * t297;
t304 = t254 * t283 + t267 * t278;
t231 = t255 * t288 + t304 * t286;
t227 = 0.1e1 / t231;
t250 = 0.1e1 / t253;
t217 = 0.1e1 / t219 ^ 2;
t228 = 0.1e1 / t231 ^ 2;
t251 = 0.1e1 / t253 ^ 2;
t240 = t242 ^ 2;
t236 = t240 * t251 + 0.1e1;
t234 = 0.1e1 / t236;
t268 = t298 * qJD(2);
t269 = t272 * qJD(2);
t313 = t279 * t283;
t244 = (-t268 * t276 - t269 * t309) * t278 - t269 * t313;
t264 = (-(-t276 * t289 - t287 * t309) * t278 + t287 * t313) * t280;
t263 = qJD(2) * t264;
t322 = t242 * t251;
t211 = (t244 * t250 - t263 * t322) * t234;
t305 = -t232 * t253 + t233 * t242;
t208 = t305 * t211 + t232 * t244 + t233 * t263;
t328 = t208 * t216 * t217;
t310 = t283 * t288;
t315 = t278 * t288;
t319 = t255 * t286;
t230 = -t254 * t310 - t267 * t315 + t319;
t226 = t230 ^ 2;
t223 = t226 * t228 + 0.1e1;
t270 = t297 * qJD(2);
t271 = t296 * qJD(2);
t317 = t276 * t284;
t257 = -t270 * t281 + t271 * t317;
t256 = t270 * t276 + t271 * t309;
t316 = t278 * t279;
t302 = -t256 * t283 + t271 * t316;
t224 = t231 * qJD(4) + t257 * t286 + t302 * t288;
t324 = t228 * t230;
t225 = t257 * t288 - t302 * t286 + (t304 * t288 - t319) * qJD(4);
t325 = t225 * t227 * t228;
t327 = (t224 * t324 - t226 * t325) / t223 ^ 2;
t243 = -t254 * t278 + t267 * t283;
t326 = t217 * t243;
t261 = -t281 * t297 + t296 * t317;
t260 = t276 * t297 + t296 * t309;
t301 = t260 * t283 - t296 * t316;
t239 = t261 * t288 + t301 * t286;
t323 = t230 * t239;
t321 = t242 * t264;
t320 = t250 * t251 * t263;
t318 = t261 * t286;
t306 = t279 * t315;
t247 = (-t272 * t309 - t298 * t276) * t278 - t272 * t313;
t303 = -t247 * t250 + t251 * t321;
t262 = (t299 * t278 + t283 * t312) * t280 * qJD(2);
t259 = t270 * t317 + t271 * t281;
t258 = t270 * t309 - t271 * t276;
t248 = -t260 * t278 - t296 * t313;
t246 = -t256 * t278 - t271 * t313;
t245 = (-t268 * t309 + t269 * t276) * t278 - t268 * t313;
t241 = t243 ^ 2;
t238 = -t301 * t288 + t318;
t221 = 0.1e1 / t223;
t215 = t241 * t217 + 0.1e1;
t212 = t303 * t234;
t209 = -t305 * t212 + t232 * t247 + t233 * t264;
t207 = 0.2e1 * t303 / t236 ^ 2 * (-t240 * t320 + t244 * t322) + (0.2e1 * t320 * t321 + t245 * t250 + (-t242 * t262 - t244 * t264 - t247 * t263) * t251) * t234;
t1 = [0, t207, 0, 0, 0, 0; 0, 0.2e1 * (t209 * t326 - t216 * t248) / t215 ^ 2 * (-t241 * t328 + t246 * t326) + ((-t258 * t278 - t270 * t313) * t216 + (-t248 * t208 - t209 * t246) * t217 + (0.2e1 * t209 * t328 + (-(t207 * t242 - t212 * t244 + t262 + (t212 * t253 + t247) * t211) * t233 - (-t207 * t253 + t212 * t263 + t245 + (t212 * t242 - t264) * t211) * t232) * t217) * t243) / t215, 0, 0, 0, 0; 0, 0.2e1 * (-t227 * t238 + t228 * t323) * t327 + ((-t258 * t310 + t259 * t286 + t270 * t306) * t227 + 0.2e1 * t323 * t325 + (-t238 * t225 - (t259 * t288 + (t258 * t283 - t270 * t316) * t286) * t230 - t239 * t224) * t228 + (t239 * t227 - (t260 * t310 - t296 * t306 - t318) * t324) * qJD(4)) * t221, 0, -0.2e1 * t327 + 0.2e1 * (t224 * t228 * t221 + (-t221 * t325 - t228 * t327) * t230) * t230, 0, 0;];
JaD_rot  = t1;
