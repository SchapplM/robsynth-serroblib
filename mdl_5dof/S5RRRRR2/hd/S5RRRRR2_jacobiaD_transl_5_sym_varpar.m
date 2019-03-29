% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiaD_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_jacobiaD_transl_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiaD_transl_5_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:52
% EndTime: 2019-03-29 15:26:52
% DurationCPUTime: 0.26s
% Computational Cost: add. (400->53), mult. (394->92), div. (0->0), fcn. (313->10), ass. (0->56)
t274 = qJ(3) + qJ(4);
t268 = sin(t274);
t270 = cos(t274);
t278 = cos(qJ(5));
t304 = qJD(5) * t278;
t272 = qJD(3) + qJD(4);
t276 = sin(qJ(5));
t309 = t272 * t276;
t322 = t268 * t304 + t270 * t309;
t305 = qJD(5) * t276;
t296 = t268 * t305;
t321 = r_i_i_C(1) * t296 + r_i_i_C(2) * t322;
t273 = qJD(1) + qJD(2);
t311 = t270 * t273;
t292 = -qJD(5) + t311;
t320 = t278 * t292;
t291 = qJD(5) * t270 - t273;
t300 = t268 * t309;
t319 = t278 * t291 - t300;
t277 = sin(qJ(3));
t318 = pkin(2) * t277;
t317 = r_i_i_C(1) * t278;
t316 = r_i_i_C(3) * t270;
t315 = pkin(1) * qJD(1);
t314 = t268 * t272;
t275 = qJ(1) + qJ(2);
t269 = sin(t275);
t313 = t269 * t273;
t312 = t270 * t272;
t271 = cos(t275);
t310 = t271 * t273;
t308 = t272 * t278;
t279 = cos(qJ(3));
t307 = t273 * t279;
t306 = qJD(3) * t277;
t303 = r_i_i_C(2) * t268 * t276;
t267 = r_i_i_C(3) * t312;
t302 = pkin(2) * t306;
t301 = t268 * t313;
t299 = t268 * t308;
t297 = t270 * t308;
t294 = t273 * t303;
t288 = t271 * t321 + t301 * t317;
t287 = t269 * t321 + t271 * t294 + t310 * t316;
t286 = t292 * t276;
t285 = -t268 * t310 - t269 * t312;
t284 = t276 * t291 + t299;
t283 = -pkin(2) * qJD(3) * t279 + (-r_i_i_C(3) * t268 - t270 * t317) * t272;
t282 = t267 + (-t270 * t305 - t299) * r_i_i_C(1) + (-t270 * t304 + t300) * r_i_i_C(2);
t253 = t269 * t319 + t271 * t286;
t254 = t269 * t284 - t271 * t320;
t281 = -pkin(2) * t271 * t307 + r_i_i_C(1) * t254 + r_i_i_C(2) * t253 + r_i_i_C(3) * t285 + t269 * t302;
t251 = t269 * t286 - t271 * t319;
t252 = t269 * t320 + t271 * t284;
t280 = -r_i_i_C(3) * t301 + t251 * r_i_i_C(2) - t252 * r_i_i_C(1) + t271 * t267 + (-t269 * t307 - t271 * t306) * pkin(2);
t1 = [-cos(qJ(1)) * t315 + t281, t281, t283 * t271 + (-t303 - t316 + t318) * t313 + t288, -t271 * r_i_i_C(1) * t297 - t269 * t294 + (-t269 * t311 - t271 * t314) * r_i_i_C(3) + t288, r_i_i_C(1) * t251 + r_i_i_C(2) * t252; -sin(qJ(1)) * t315 + t280, t280 (-t268 * t317 - t318) * t310 + t283 * t269 + t287, -r_i_i_C(3) * t269 * t314 + t285 * t317 + t287, -r_i_i_C(1) * t253 + r_i_i_C(2) * t254; 0, 0, t282 - t302, t282 (t296 - t297) * r_i_i_C(2) - t322 * r_i_i_C(1);];
JaD_transl  = t1;
