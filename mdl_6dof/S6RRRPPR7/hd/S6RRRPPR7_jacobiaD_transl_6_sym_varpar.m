% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:02
% EndTime: 2019-02-26 22:07:03
% DurationCPUTime: 0.54s
% Computational Cost: add. (553->96), mult. (1296->147), div. (0->0), fcn. (1206->10), ass. (0->62)
t299 = sin(qJ(2));
t298 = sin(qJ(3));
t301 = cos(qJ(3));
t295 = pkin(10) + qJ(6);
t294 = cos(t295);
t293 = sin(t295);
t322 = sin(pkin(10)) * pkin(5) + qJ(4);
t311 = -t293 * r_i_i_C(1) - t322;
t307 = t294 * r_i_i_C(2) - t311;
t342 = -pkin(3) - cos(pkin(10)) * pkin(5) - pkin(4);
t319 = t293 * r_i_i_C(2) + t342;
t310 = t294 * r_i_i_C(1) - t319;
t344 = t298 * t307 + t301 * t310 + pkin(2);
t302 = cos(qJ(2));
t328 = -r_i_i_C(3) - pkin(9) - qJ(5) + pkin(8);
t350 = t328 * t302;
t356 = t344 * t299 - t350;
t355 = (-t299 * pkin(2) + t350) * qJD(2) - t299 * qJD(5);
t354 = t299 * (qJD(3) - qJD(6));
t303 = cos(qJ(1));
t338 = t303 * t301;
t300 = sin(qJ(1));
t340 = t300 * t302;
t276 = t298 * t340 + t338;
t339 = t303 * t298;
t277 = t301 * t340 - t339;
t316 = t276 * t293 + t277 * t294;
t317 = t276 * t294 - t277 * t293;
t352 = (r_i_i_C(1) * t316 + r_i_i_C(2) * t317) * qJD(6);
t312 = t293 * t298 + t294 * t301;
t313 = t293 * t301 - t294 * t298;
t351 = (t298 * t310 - t301 * t307) * qJD(3) + (r_i_i_C(1) * t313 + r_i_i_C(2) * t312) * qJD(6) - t328 * qJD(2) - t298 * qJD(4);
t332 = qJD(3) * t301;
t336 = qJD(1) * t303;
t309 = t298 * t336 + t300 * t332;
t331 = qJD(3) * t303;
t324 = t298 * t331;
t335 = qJD(2) * t299;
t327 = t300 * t335;
t337 = qJD(1) * t300;
t274 = -t298 * t327 - t301 * t337 + t302 * t309 - t324;
t271 = t274 * t294;
t341 = t300 * t298;
t334 = qJD(2) * t302;
t333 = qJD(2) * t303;
t326 = qJD(3) * t341;
t325 = t299 * t333;
t323 = t301 * t331;
t318 = -(-t312 * t354 + t313 * t334) * r_i_i_C(1) - (t312 * t334 + t313 * t354) * r_i_i_C(2);
t278 = -t300 * t301 + t302 * t339;
t279 = t302 * t338 + t341;
t315 = t278 * t294 - t279 * t293;
t314 = t278 * t293 + t279 * t294;
t308 = -pkin(2) * t302 - t299 * t328 - pkin(1);
t305 = -qJD(2) * t344 - qJD(5);
t304 = t351 * t299 + t305 * t302;
t275 = qJD(1) * t279 - t301 * t327 - t302 * t326 - t323;
t273 = t302 * t324 + (t302 * t337 + t325) * t301 - t309;
t272 = qJD(1) * t276 + t298 * t325 - t302 * t323 - t326;
t268 = qJD(6) * t315 - t272 * t293 - t273 * t294;
t267 = -qJD(6) * t314 - t272 * t294 + t273 * t293;
t1 = [-t271 * r_i_i_C(2) - t276 * qJD(4) - t310 * t275 + t311 * t274 + (-r_i_i_C(1) * t317 + r_i_i_C(2) * t316) * qJD(6) - t355 * t300 + (-t300 * pkin(7) + t303 * t308) * qJD(1), t304 * t303 + t356 * t337, t279 * qJD(4) - t307 * t273 + t310 * t272 + (r_i_i_C(1) * t314 + r_i_i_C(2) * t315) * qJD(6), -t272, t299 * t337 - t302 * t333, t267 * r_i_i_C(1) - t268 * r_i_i_C(2); t268 * r_i_i_C(1) + t267 * r_i_i_C(2) + t278 * qJD(4) + t342 * t273 - t322 * t272 + t355 * t303 + (pkin(7) * t303 + t300 * t308) * qJD(1), t304 * t300 - t356 * t336, -t271 * r_i_i_C(1) + t277 * qJD(4) + t319 * t274 + t307 * t275 + t352, t274, -t299 * t336 - t300 * t334 (-t275 * t293 + t271) * r_i_i_C(1) + (-t274 * t293 - t275 * t294) * r_i_i_C(2) - t352; 0, t305 * t299 - t351 * t302 (t298 * t342 + t301 * t322) * t334 + (qJD(4) * t301 + (-t322 * t298 + t301 * t342) * qJD(3)) * t299 - t318, t298 * t334 + t299 * t332, -t335, t318;];
JaD_transl  = t1;
