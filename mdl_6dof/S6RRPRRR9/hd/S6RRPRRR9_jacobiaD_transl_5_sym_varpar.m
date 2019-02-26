% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:51
% EndTime: 2019-02-26 21:58:51
% DurationCPUTime: 0.30s
% Computational Cost: add. (456->78), mult. (710->123), div. (0->0), fcn. (669->12), ass. (0->61)
t292 = sin(qJ(1));
t290 = cos(pkin(6));
t303 = qJD(2) * t290 + qJD(1);
t291 = sin(qJ(2));
t323 = t292 * t291;
t309 = t290 * t323;
t314 = qJD(2) * t291;
t293 = cos(qJ(2));
t294 = cos(qJ(1));
t320 = t294 * t293;
t268 = -qJD(1) * t309 - t292 * t314 + t303 * t320;
t288 = qJD(4) + qJD(5);
t289 = sin(pkin(6));
t324 = t289 * t294;
t333 = t288 * t324 - t268;
t287 = pkin(12) + qJ(4);
t285 = qJ(5) + t287;
t281 = sin(t285);
t282 = cos(t285);
t302 = r_i_i_C(1) * t281 + r_i_i_C(2) * t282;
t283 = sin(t287);
t329 = pkin(4) * qJD(4);
t312 = t283 * t329;
t332 = t302 * t288 + t312;
t331 = pkin(8) + pkin(4) * t283 + sin(pkin(12)) * pkin(3);
t330 = r_i_i_C(3) + pkin(10) + pkin(9) + qJ(3);
t321 = t294 * t291;
t322 = t292 * t293;
t271 = t290 * t321 + t322;
t328 = t271 * t288;
t327 = t281 * t288;
t326 = t289 * t291;
t325 = t289 * t292;
t273 = -t309 + t320;
t315 = qJD(1) * t294;
t306 = t289 * t315;
t297 = -t273 * t288 + t306;
t299 = t290 * t322 + t321;
t266 = t271 * qJD(1) + t299 * qJD(2);
t301 = t288 * t325 - t266;
t261 = -t301 * t281 + t297 * t282;
t262 = t297 * t281 + t301 * t282;
t319 = t261 * r_i_i_C(1) - t262 * r_i_i_C(2);
t316 = qJD(1) * t292;
t307 = t289 * t316;
t298 = t307 - t328;
t304 = t333 * t282;
t318 = (t333 * t281 + t298 * t282) * r_i_i_C(1) + (-t298 * t281 + t304) * r_i_i_C(2);
t313 = qJD(2) * t293;
t305 = t289 * t313;
t296 = -t288 * t290 - t305;
t311 = t288 * t326;
t317 = (t296 * t281 - t282 * t311) * r_i_i_C(1) + (t281 * t311 + t296 * t282) * r_i_i_C(2);
t308 = t290 * t320;
t284 = cos(t287);
t274 = pkin(4) * t284 + cos(pkin(12)) * pkin(3) + pkin(2);
t300 = t282 * r_i_i_C(1) - t281 * r_i_i_C(2) + t274;
t270 = -t308 + t323;
t267 = t299 * qJD(1) + t271 * qJD(2);
t265 = -qJD(1) * t308 - t294 * t313 + t303 * t323;
t1 = [(t271 * t327 + t304) * r_i_i_C(1) + (t268 * t281 + t282 * t328) * r_i_i_C(2) - t268 * t274 + t271 * t312 - t270 * qJD(3) - pkin(1) * t315 - t330 * t267 + ((-r_i_i_C(2) * t327 + t284 * t329) * t294 + (-t302 - t331) * t316) * t289, t273 * qJD(3) + t300 * t265 - t330 * t266 + t299 * t332, -t265 (t284 * t306 + t266 * t283 + (-t273 * t284 - t283 * t325) * qJD(4)) * pkin(4) + t319, t319, 0; t262 * r_i_i_C(1) + t261 * r_i_i_C(2) + t299 * qJD(3) - t266 * t274 - t330 * t265 + (-t273 * t283 + t284 * t325) * t329 + (-pkin(1) * t292 + t331 * t324) * qJD(1), t271 * qJD(3) - t300 * t267 + t330 * t268 + t332 * t270, t267 (t284 * t307 - t268 * t283 + (-t271 * t284 + t283 * t324) * qJD(4)) * pkin(4) + t318, t318, 0; 0 (qJD(3) * t291 - t332 * t293 + (-t300 * t291 + t330 * t293) * qJD(2)) * t289, t289 * t314 (-t283 * t305 + (-t283 * t290 - t284 * t326) * qJD(4)) * pkin(4) + t317, t317, 0;];
JaD_transl  = t1;
