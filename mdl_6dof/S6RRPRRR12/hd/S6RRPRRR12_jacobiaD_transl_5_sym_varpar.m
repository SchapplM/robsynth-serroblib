% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR12_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR12_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:00:25
% EndTime: 2019-02-26 22:00:25
% DurationCPUTime: 0.25s
% Computational Cost: add. (378->75), mult. (797->120), div. (0->0), fcn. (746->10), ass. (0->56)
t285 = sin(qJ(4));
t326 = t285 * pkin(4);
t288 = cos(qJ(4));
t325 = t288 * pkin(4) + pkin(3) + pkin(8);
t284 = cos(pkin(6));
t289 = cos(qJ(2));
t290 = cos(qJ(1));
t317 = t290 * t289;
t306 = t284 * t317;
t286 = sin(qJ(2));
t287 = sin(qJ(1));
t320 = t287 * t286;
t269 = -t306 + t320;
t281 = qJD(4) + qJD(5);
t324 = t269 * t281;
t283 = sin(pkin(6));
t323 = t283 * t287;
t322 = t283 * t289;
t321 = t283 * t290;
t319 = t287 * t289;
t318 = t290 * t286;
t282 = qJ(4) + qJ(5);
t279 = sin(t282);
t280 = cos(t282);
t313 = qJD(1) * t287;
t305 = t283 * t313;
t296 = -t305 - t324;
t271 = t284 * t319 + t318;
t297 = t284 * t318 + t319;
t267 = t271 * qJD(1) + t297 * qJD(2);
t298 = t281 * t321 + t267;
t316 = (t296 * t279 + t298 * t280) * r_i_i_C(1) + (-t298 * t279 + t296 * t280) * r_i_i_C(2);
t312 = qJD(1) * t290;
t304 = t283 * t312;
t295 = t271 * t281 + t304;
t302 = qJD(2) * t284 + qJD(1);
t265 = -qJD(1) * t306 - qJD(2) * t317 + t302 * t320;
t299 = -t281 * t323 - t265;
t261 = -t295 * t279 + t299 * t280;
t262 = t299 * t279 + t295 * t280;
t315 = t261 * r_i_i_C(1) - t262 * r_i_i_C(2);
t311 = qJD(2) * t286;
t303 = t283 * t311;
t294 = -t281 * t284 + t303;
t308 = t281 * t322;
t314 = (t279 * t308 + t294 * t280) * r_i_i_C(1) + (-t294 * t279 + t280 * t308) * r_i_i_C(2);
t310 = qJD(4) * t288;
t309 = -r_i_i_C(3) - pkin(10) - pkin(9) - pkin(2);
t307 = t284 * t320;
t301 = r_i_i_C(1) * t280 - r_i_i_C(2) * t279;
t300 = -t279 * r_i_i_C(1) - t280 * r_i_i_C(2);
t293 = qJ(3) - t300 + t326;
t292 = pkin(4) * t310 + t301 * t281 + qJD(3);
t268 = -qJD(1) * t307 - t287 * t311 + t302 * t317;
t266 = t297 * qJD(1) + t271 * qJD(2);
t1 = [(-t267 * t279 - t280 * t324) * r_i_i_C(1) + (-t267 * t280 + t279 * t324) * r_i_i_C(2) - t267 * qJ(3) - t269 * qJD(3) - pkin(1) * t312 + (-t267 * t285 - t269 * t310) * pkin(4) + t309 * t268 + ((-qJD(4) * t326 + t300 * t281) * t290 + (-t301 - t325) * t313) * t283, -t309 * t265 + t292 * (-t307 + t317) - t293 * t266, -t265 (-t285 * t304 - t265 * t288 + (-t271 * t285 - t288 * t323) * qJD(4)) * pkin(4) + t315, t315, 0; t262 * r_i_i_C(1) + t261 * r_i_i_C(2) - t265 * qJ(3) + t271 * qJD(3) + t309 * t266 + (-pkin(1) * t287 + t325 * t321) * qJD(1) + (-t265 * t285 + (t271 * t288 - t285 * t323) * qJD(4)) * pkin(4), t309 * t267 + t293 * t268 + t292 * t297, t267 (-t285 * t305 + t267 * t288 + (-t269 * t285 + t288 * t321) * qJD(4)) * pkin(4) + t316, t316, 0; 0 (t292 * t286 + (t309 * t286 + t293 * t289) * qJD(2)) * t283, t303 (t288 * t303 + (-t284 * t288 + t285 * t322) * qJD(4)) * pkin(4) + t314, t314, 0;];
JaD_transl  = t1;
