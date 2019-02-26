% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:10
% EndTime: 2019-02-26 22:05:11
% DurationCPUTime: 0.57s
% Computational Cost: add. (772->106), mult. (1258->159), div. (0->0), fcn. (1161->10), ass. (0->69)
t304 = sin(qJ(2));
t307 = cos(qJ(3));
t351 = t307 * pkin(3);
t297 = pkin(2) + t351;
t300 = qJ(3) + pkin(10);
t298 = sin(t300);
t299 = cos(t300);
t306 = cos(qJ(6));
t302 = sin(qJ(6));
t353 = -pkin(4) - pkin(5);
t330 = t302 * r_i_i_C(2) + t353;
t315 = t306 * r_i_i_C(1) - t330;
t331 = -t302 * r_i_i_C(1) - qJ(5);
t316 = t306 * r_i_i_C(2) - t331;
t355 = t298 * t316 + t299 * t315 + t297;
t308 = cos(qJ(2));
t337 = -r_i_i_C(3) - pkin(9) + qJ(4) + pkin(8);
t361 = t337 * t308;
t366 = t355 * t304 - t361;
t329 = qJD(3) * t308 - qJD(1);
t303 = sin(qJ(3));
t352 = t303 * pkin(3);
t365 = qJD(1) * pkin(7) + (-t304 * t297 + t361) * qJD(2) + t304 * qJD(4) - t329 * t352;
t309 = cos(qJ(1));
t347 = t309 * t299;
t305 = sin(qJ(1));
t349 = t305 * t308;
t281 = t298 * t349 + t347;
t348 = t309 * t298;
t282 = t299 * t349 - t348;
t322 = t281 * t302 + t282 * t306;
t323 = t281 * t306 - t282 * t302;
t363 = (r_i_i_C(1) * t322 + r_i_i_C(2) * t323) * qJD(6);
t318 = t298 * t302 + t299 * t306;
t319 = t298 * t306 - t299 * t302;
t362 = (t298 * t315 - t299 * t316 + t352) * qJD(3) - (r_i_i_C(1) * t319 - r_i_i_C(2) * t318) * qJD(6) - t337 * qJD(2) - t298 * qJD(5);
t358 = t304 * (-qJD(3) + qJD(6));
t340 = qJD(3) * t305;
t344 = qJD(1) * t309;
t314 = t298 * t344 + t299 * t340;
t339 = qJD(3) * t309;
t333 = t298 * t339;
t343 = qJD(2) * t304;
t336 = t305 * t343;
t346 = qJD(1) * t305;
t279 = -t298 * t336 - t299 * t346 + t308 * t314 - t333;
t276 = t279 * t306;
t345 = qJD(1) * t308;
t342 = qJD(2) * t308;
t341 = qJD(2) * t309;
t335 = t304 * t341;
t334 = t298 * t340;
t332 = t299 * t339;
t328 = -qJD(3) + t345;
t325 = -(t318 * t358 - t319 * t342) * r_i_i_C(1) - (t318 * t342 + t319 * t358) * r_i_i_C(2);
t283 = -t305 * t299 + t308 * t348;
t284 = t305 * t298 + t308 * t347;
t321 = t283 * t306 - t284 * t302;
t320 = t283 * t302 + t284 * t306;
t317 = t329 * t307;
t313 = qJD(3) * t351 + (-t297 * t308 - t304 * t337 - pkin(1)) * qJD(1);
t311 = -qJD(2) * t355 + qJD(4);
t310 = t362 * t304 + t311 * t308;
t280 = qJD(1) * t284 - t299 * t336 - t308 * t334 - t332;
t278 = t308 * t333 + (t305 * t345 + t335) * t299 - t314;
t277 = qJD(1) * t281 + t298 * t335 - t308 * t332 - t334;
t273 = qJD(6) * t321 - t277 * t302 - t278 * t306;
t272 = -qJD(6) * t320 - t277 * t306 + t278 * t302;
t1 = [-t276 * r_i_i_C(2) - t281 * qJD(5) + t331 * t279 - t315 * t280 + (-r_i_i_C(1) * t323 + r_i_i_C(2) * t322) * qJD(6) + t313 * t309 - t365 * t305, t310 * t309 + t366 * t346, t284 * qJD(5) - t316 * t278 + t315 * t277 + (r_i_i_C(1) * t320 + r_i_i_C(2) * t321) * qJD(6) + (-t309 * t317 + (t305 * t328 + t335) * t303) * pkin(3), -t304 * t346 + t308 * t341, -t277, t272 * r_i_i_C(1) - t273 * r_i_i_C(2); t273 * r_i_i_C(1) + t272 * r_i_i_C(2) - t277 * qJ(5) + t283 * qJD(5) + t353 * t278 + t313 * t305 + t365 * t309, t310 * t305 - t366 * t344, -t276 * r_i_i_C(1) + t282 * qJD(5) + t316 * t280 + t330 * t279 + t363 + (-t305 * t317 + (-t309 * t328 + t336) * t303) * pkin(3), t304 * t344 + t305 * t342, t279 (-t280 * t302 + t276) * r_i_i_C(1) + (-t279 * t302 - t280 * t306) * r_i_i_C(2) - t363; 0, t311 * t304 - t362 * t308 (t299 * qJ(5) + t298 * t353 - t352) * t342 + (qJD(5) * t299 + (-t298 * qJ(5) + t299 * t353 - t351) * qJD(3)) * t304 - t325, t343, t304 * qJD(3) * t299 + t298 * t342, t325;];
JaD_transl  = t1;
