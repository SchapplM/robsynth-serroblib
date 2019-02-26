% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:12
% EndTime: 2019-02-26 21:18:12
% DurationCPUTime: 0.34s
% Computational Cost: add. (580->75), mult. (536->106), div. (0->0), fcn. (394->10), ass. (0->67)
t322 = pkin(10) + r_i_i_C(3);
t273 = cos(qJ(6));
t317 = r_i_i_C(1) * t273;
t329 = pkin(5) + t317;
t269 = qJ(3) + qJ(4);
t264 = sin(t269);
t271 = sin(qJ(3));
t266 = qJ(5) + t269;
t262 = cos(t266);
t302 = t322 * t262;
t261 = sin(t266);
t319 = pkin(5) * t261;
t327 = -t271 * pkin(3) - pkin(4) * t264 - qJ(2) + t302 - t319;
t275 = cos(qJ(1));
t288 = qJD(1) * t261 + qJD(6);
t326 = t288 * t275;
t325 = t322 * t261;
t265 = cos(t269);
t270 = sin(qJ(6));
t316 = r_i_i_C(2) * t270;
t304 = t262 * t316;
t324 = pkin(4) * t265 - t304;
t272 = sin(qJ(1));
t267 = qJD(3) + qJD(4);
t263 = qJD(5) + t267;
t311 = t263 * t275;
t323 = -t262 * t311 + t288 * t272;
t320 = pkin(4) * t267;
t318 = pkin(5) * t262;
t315 = -pkin(1) - pkin(9) - pkin(8) - pkin(7);
t314 = pkin(3) * qJD(3);
t313 = t261 * t270;
t312 = t263 * t273;
t310 = qJD(1) * t272;
t309 = qJD(1) * t275;
t308 = qJD(6) * t270;
t307 = qJD(6) * t273;
t306 = t264 * t320;
t305 = t265 * t320;
t274 = cos(qJ(3));
t303 = t274 * t314;
t301 = t263 * t313;
t299 = t262 * t263 * t272;
t295 = t262 * t309;
t294 = t262 * t308;
t293 = t262 * t307;
t292 = r_i_i_C(2) * t301;
t291 = qJD(1) * t262 * t317;
t290 = r_i_i_C(2) * t293;
t289 = -qJD(6) * t261 - qJD(1);
t286 = t289 * t275;
t285 = qJD(1) * (t274 * pkin(3) + t324);
t284 = pkin(5) * t295 + t272 * t292 + t275 * t291 + t322 * (t261 * t309 + t299);
t283 = qJD(1) * t324;
t282 = t261 * t312 + t294;
t281 = t272 * t291 + t310 * t318 + (r_i_i_C(1) * t294 + t290) * t275 + (t322 * t310 + t329 * t311) * t261;
t280 = (-r_i_i_C(2) * t313 - t302) * t263;
t279 = qJD(2) + t303 + t305 + (t318 + t325) * t263;
t278 = (-t329 * t262 + t304 - t325) * t263 + (r_i_i_C(1) * t308 + r_i_i_C(2) * t307) * t261;
t277 = -t282 * r_i_i_C(1) - t263 * t319 - t290;
t276 = t278 - t305;
t241 = -t271 * t314 - t306;
t238 = t273 * t326 + (t262 * t312 + t289 * t270) * t272;
t237 = t289 * t273 * t272 + (-t299 - t326) * t270;
t236 = t270 * t286 - t323 * t273;
t235 = t323 * t270 + t273 * t286;
t1 = [t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t279 * t275 + (t327 * t272 + t315 * t275) * qJD(1), t309, t275 * t285 + (t241 + t277) * t272 + t284, t275 * t283 + (t277 - t306) * t272 + t284, t277 * t272 - t295 * t316 + t284, t237 * r_i_i_C(1) - t238 * r_i_i_C(2); t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t279 * t272 + (t315 * t272 - t327 * t275) * qJD(1), t310, t272 * t285 + (-t241 + t280) * t275 + t281, t272 * t283 + (t280 + t306) * t275 + t281, -t275 * t292 + (-t310 * t316 - t322 * t311) * t262 + t281, -t235 * r_i_i_C(1) + t236 * r_i_i_C(2); 0, 0, t276 - t303, t276, t278, t282 * r_i_i_C(2) + (-t293 + t301) * r_i_i_C(1);];
JaD_transl  = t1;
