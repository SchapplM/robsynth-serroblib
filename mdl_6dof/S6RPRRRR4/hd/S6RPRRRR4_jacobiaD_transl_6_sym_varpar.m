% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:27
% EndTime: 2019-02-26 21:16:27
% DurationCPUTime: 0.30s
% Computational Cost: add. (738->75), mult. (530->107), div. (0->0), fcn. (392->11), ass. (0->65)
t275 = pkin(11) + qJ(3);
t273 = qJ(4) + t275;
t269 = qJ(5) + t273;
t265 = sin(t269);
t279 = cos(qJ(6));
t333 = r_i_i_C(1) * t279 + pkin(5);
t294 = t333 * t265;
t277 = sin(qJ(6));
t315 = qJD(6) * t279;
t266 = cos(t269);
t276 = qJD(3) + qJD(4);
t272 = qJD(5) + t276;
t323 = t266 * t272;
t335 = t265 * t315 + t277 * t323;
t330 = pkin(10) + r_i_i_C(3);
t308 = t330 * t266;
t270 = sin(t275);
t324 = pkin(3) * qJD(3);
t311 = t270 * t324;
t267 = sin(t273);
t328 = pkin(4) * t276;
t314 = t267 * t328;
t327 = pkin(5) * t265;
t334 = (t308 - t327) * t272 - t311 - t314;
t316 = qJD(6) * t277;
t301 = t265 * t316;
t331 = r_i_i_C(1) * t301 + t335 * r_i_i_C(2);
t268 = cos(t273);
t293 = t333 * t266;
t309 = t330 * t265;
t283 = (-t293 - t309) * t272 - t268 * t328;
t329 = pkin(4) * t267;
t325 = r_i_i_C(2) * t277;
t278 = sin(qJ(1));
t322 = t272 * t278;
t321 = t272 * t279;
t280 = cos(qJ(1));
t320 = t272 * t280;
t319 = t279 * t280;
t318 = qJD(1) * t278;
t317 = qJD(1) * t280;
t312 = t265 * t325;
t310 = qJD(1) * t325;
t307 = t330 * t278;
t306 = t265 * t321;
t296 = qJD(6) * t266 - qJD(1);
t295 = qJD(1) * t266 - qJD(6);
t292 = t333 * t280;
t291 = t331 * t280 + t318 * t294;
t290 = t296 * t277;
t289 = t280 * t265 * t310 + t331 * t278 + t317 * t308;
t288 = -t308 - t312;
t271 = cos(t275);
t287 = -pkin(5) * t266 - pkin(4) * t268 - pkin(3) * t271 - cos(pkin(11)) * pkin(2) - pkin(1) - t309;
t286 = t265 * t320 + t295 * t278;
t284 = -t271 * t324 + t283;
t282 = -t266 * r_i_i_C(2) * t315 + (-t266 * t316 - t306) * r_i_i_C(1) + t330 * t323 + (-t327 + t312) * t272;
t281 = t282 - t314;
t274 = -pkin(9) - pkin(8) - pkin(7) - qJ(2);
t259 = -pkin(3) * t270 - t329;
t246 = -t295 * t319 + (t290 + t306) * t278;
t245 = t296 * t279 * t278 + (-t265 * t322 + t295 * t280) * t277;
t244 = t286 * t279 + t280 * t290;
t243 = t286 * t277 - t296 * t319;
t1 = [t246 * r_i_i_C(1) + t245 * r_i_i_C(2) + t280 * qJD(2) - t334 * t278 + (t274 * t278 + t287 * t280) * qJD(1), t317 (-t259 + t288) * t318 + t284 * t280 + t291 (t288 + t329) * t318 + t283 * t280 + t291 (-t278 * t310 - t330 * t320) * t265 + (-qJD(1) * t307 - t272 * t292) * t266 + t291, t243 * r_i_i_C(1) + t244 * r_i_i_C(2); -t244 * r_i_i_C(1) + t243 * r_i_i_C(2) + t278 * qJD(2) + t334 * t280 + (-t274 * t280 + t287 * t278) * qJD(1), t318 (t259 - t294) * t317 + t284 * t278 + t289 (-t294 - t329) * t317 + t283 * t278 + t289, -t293 * t322 + (-qJD(1) * t292 - t272 * t307) * t265 + t289, -t245 * r_i_i_C(1) + t246 * r_i_i_C(2); 0, 0, t281 - t311, t281, t282 (-t266 * t321 + t301) * r_i_i_C(2) - t335 * r_i_i_C(1);];
JaD_transl  = t1;
