% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:49:27
% EndTime: 2019-02-26 22:49:27
% DurationCPUTime: 0.34s
% Computational Cost: add. (649->85), mult. (874->128), div. (0->0), fcn. (805->12), ass. (0->66)
t295 = cos(pkin(6));
t297 = sin(qJ(2));
t298 = sin(qJ(1));
t327 = t298 * t297;
t316 = t295 * t327;
t300 = cos(qJ(2));
t301 = cos(qJ(1));
t325 = t300 * t301;
t270 = -qJD(1) * t316 - qJD(2) * t327 + (qJD(2) * t295 + qJD(1)) * t325;
t291 = qJD(3) + qJD(4);
t286 = qJD(5) + t291;
t294 = sin(pkin(6));
t329 = t294 * t301;
t339 = t286 * t329 - t270;
t293 = qJ(3) + qJ(4);
t287 = sin(t293);
t338 = pkin(4) * t287;
t296 = sin(qJ(3));
t279 = pkin(3) * t296 + t338;
t337 = pkin(8) + t279;
t336 = r_i_i_C(3) + pkin(11) + pkin(10) + pkin(9);
t335 = pkin(3) * qJD(3);
t326 = t298 * t300;
t328 = t297 * t301;
t275 = t295 * t328 + t326;
t334 = t275 * t286;
t289 = qJ(5) + t293;
t284 = sin(t289);
t333 = t284 * t286;
t288 = cos(t293);
t332 = t288 * t291;
t331 = t294 * t297;
t330 = t294 * t298;
t285 = cos(t289);
t306 = t316 - t325;
t320 = qJD(1) * t301;
t314 = t294 * t320;
t304 = t286 * t306 + t314;
t307 = t295 * t326 + t328;
t268 = t275 * qJD(1) + t307 * qJD(2);
t310 = t286 * t330 - t268;
t263 = -t310 * t284 + t304 * t285;
t264 = t304 * t284 + t310 * t285;
t324 = t263 * r_i_i_C(1) - t264 * r_i_i_C(2);
t321 = qJD(1) * t298;
t315 = t294 * t321;
t305 = t315 - t334;
t312 = t339 * t285;
t323 = (t339 * t284 + t305 * t285) * r_i_i_C(1) + (-t305 * t284 + t312) * r_i_i_C(2);
t319 = qJD(2) * t300;
t313 = t294 * t319;
t303 = -t286 * t295 - t313;
t318 = t286 * t331;
t322 = (t303 * t284 - t285 * t318) * r_i_i_C(1) + (t284 * t318 + t303 * t285) * r_i_i_C(2);
t299 = cos(qJ(3));
t280 = t299 * pkin(3) + pkin(4) * t288;
t311 = -r_i_i_C(1) * t284 - r_i_i_C(2) * t285;
t278 = pkin(2) + t280;
t309 = r_i_i_C(1) * t285 - r_i_i_C(2) * t284 + t278;
t308 = t295 * t325 - t327;
t272 = -t291 * t338 - t296 * t335;
t302 = t311 * t286 + t272;
t273 = pkin(4) * t332 + t299 * t335;
t269 = t307 * qJD(1) + t275 * qJD(2);
t267 = -t308 * qJD(1) + t306 * qJD(2);
t1 = [(t275 * t333 + t312) * r_i_i_C(1) + (t270 * t284 + t285 * t334) * r_i_i_C(2) - t270 * t278 - t275 * t272 - pkin(1) * t320 - t336 * t269 + ((-r_i_i_C(2) * t333 + t273) * t301 + (t311 - t337) * t321) * t294, t309 * t267 - t336 * t268 - t302 * t307, t268 * t279 + t306 * t273 + (t272 * t298 + t280 * t320) * t294 + t324 ((t291 * t306 + t314) * t288 + (-t291 * t330 + t268) * t287) * pkin(4) + t324, t324, 0; t273 * t330 + t264 * r_i_i_C(1) + t263 * r_i_i_C(2) - t268 * t278 - t306 * t272 - t336 * t267 + (-pkin(1) * t298 + t337 * t329) * qJD(1), -t309 * t269 + t336 * t270 + t302 * t308, -t270 * t279 - t275 * t273 + (-t272 * t301 + t280 * t321) * t294 + t323 ((-t275 * t291 + t315) * t288 + (t291 * t329 - t270) * t287) * pkin(4) + t323, t323, 0; 0 (t302 * t300 + (-t309 * t297 + t336 * t300) * qJD(2)) * t294, t295 * t272 + (-t273 * t297 - t279 * t319) * t294 + t322 (-t331 * t332 + (-t291 * t295 - t313) * t287) * pkin(4) + t322, t322, 0;];
JaD_transl  = t1;
