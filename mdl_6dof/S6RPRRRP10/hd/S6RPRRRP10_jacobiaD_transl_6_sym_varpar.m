% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:06
% EndTime: 2019-02-26 21:13:06
% DurationCPUTime: 0.34s
% Computational Cost: add. (527->69), mult. (746->98), div. (0->0), fcn. (626->8), ass. (0->55)
t285 = qJ(4) + qJ(5);
t282 = sin(t285);
t284 = qJD(4) + qJD(5);
t286 = sin(qJ(4));
t329 = pkin(4) * qJD(4);
t331 = r_i_i_C(2) + pkin(9) + pkin(8);
t295 = t331 * qJD(3) + qJD(6) * t282 - t286 * t329;
t283 = cos(t285);
t330 = r_i_i_C(3) + qJ(6);
t309 = t330 * t283;
t333 = r_i_i_C(1) + pkin(5);
t341 = (-t333 * t282 + t309) * t284 + t295;
t289 = cos(qJ(4));
t281 = t289 * pkin(4) + pkin(3);
t310 = t330 * t282;
t301 = -t333 * t283 - t310;
t340 = -t281 + t301;
t287 = sin(qJ(3));
t290 = cos(qJ(3));
t316 = t289 * t329;
t317 = t283 * qJD(6);
t338 = (-t281 * t287 + t331 * t290 - qJ(2)) * qJD(1) - t316 + t317;
t337 = t289 * (qJD(4) * t287 + qJD(1));
t332 = pkin(4) * t286;
t328 = t284 * t287;
t327 = t284 * t290;
t288 = sin(qJ(1));
t326 = t288 * t282;
t325 = t288 * t283;
t291 = cos(qJ(1));
t324 = t291 * t282;
t323 = qJD(1) * t287;
t322 = qJD(1) * t288;
t321 = qJD(1) * t291;
t320 = qJD(3) * t287;
t319 = qJD(3) * t290;
t318 = qJD(3) * t291;
t315 = t291 * t287 * t283;
t313 = t282 * t320;
t314 = t290 * t317 + t333 * t313;
t312 = t288 * t319;
t311 = t290 * t318;
t307 = t284 + t323;
t305 = -qJD(4) - t323;
t302 = t287 * t325 + t324;
t266 = -t282 * t311 - t283 * t321 - t284 * t315 + t307 * t326;
t267 = -t283 * t311 + (t287 * t324 + t325) * t284 + t302 * qJD(1);
t300 = -(t315 - t326) * qJD(6) + t330 * t267 - t333 * t266;
t297 = t307 * t291 + t312;
t268 = (qJD(1) + t328) * t325 + t297 * t282;
t269 = -t282 * t322 + t297 * t283 - t326 * t328;
t299 = t302 * qJD(6) - t333 * t268 + t330 * t269;
t296 = qJD(3) * t340;
t293 = t281 * t319 + qJD(2) + (-pkin(1) - pkin(7) - t332) * qJD(1) + t295 * t287;
t1 = [-t330 * t266 - t333 * t267 + t338 * t288 + t293 * t291, t321 (t331 * t287 - t290 * t340) * t321 + (t287 * t296 + t341 * t290) * t288 (-t288 * t337 + (t305 * t291 - t312) * t286) * pkin(4) + t299, t299, t268; t330 * t268 + t333 * t269 + t293 * t288 - t338 * t291, t322 (-t318 * t340 + t331 * t322) * t287 + (-t291 * t341 - t322 * t340) * t290 (t291 * t337 + (t305 * t288 + t311) * t286) * pkin(4) + t300, t300, t266; 0, 0, -t287 * t341 + t290 * t296 (-t309 + t332) * t320 + (t301 * t284 - t316) * t290 + t314, -t310 * t327 + (-t330 * t320 - t333 * t327) * t283 + t314, t283 * t327 - t313;];
JaD_transl  = t1;
