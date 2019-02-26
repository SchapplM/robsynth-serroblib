% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:42:26
% EndTime: 2019-02-26 22:42:26
% DurationCPUTime: 0.45s
% Computational Cost: add. (1099->89), mult. (968->122), div. (0->0), fcn. (812->10), ass. (0->69)
t315 = qJ(3) + qJ(4);
t311 = qJ(5) + t315;
t307 = cos(t311);
t363 = r_i_i_C(3) + qJ(6);
t372 = t363 * t307;
t313 = qJD(3) + qJD(4);
t308 = qJD(5) + t313;
t309 = sin(t315);
t316 = sin(qJ(3));
t362 = pkin(3) * qJD(3);
t366 = pkin(4) * t313;
t289 = -t309 * t366 - t316 * t362;
t306 = sin(t311);
t364 = r_i_i_C(2) + pkin(10) + pkin(9) + pkin(8);
t329 = t364 * qJD(2) + qJD(6) * t306 + t289;
t368 = pkin(5) + r_i_i_C(1);
t349 = t368 * t306;
t322 = (t349 - t372) * t308 - t329;
t310 = cos(t315);
t319 = cos(qJ(3));
t301 = t319 * pkin(3) + pkin(4) * t310;
t297 = pkin(2) + t301;
t320 = cos(qJ(2));
t317 = sin(qJ(2));
t353 = qJD(2) * t317;
t371 = -t297 * t353 + t329 * t320;
t331 = -t363 * t306 - t368 * t307;
t326 = -t297 + t331;
t367 = pkin(4) * t309;
t300 = pkin(3) * t316 + t367;
t365 = pkin(7) + t300;
t321 = cos(qJ(1));
t361 = t307 * t321;
t360 = t308 * t317;
t359 = t308 * t321;
t318 = sin(qJ(1));
t358 = t318 * t306;
t357 = t318 * t320;
t356 = qJD(1) * t318;
t355 = qJD(1) * t320;
t354 = qJD(1) * t321;
t352 = qJD(2) * t320;
t351 = qJD(6) * t307;
t350 = t310 * t366;
t348 = t307 * t360;
t347 = t306 * t357;
t346 = t306 * t359;
t345 = t307 * t359;
t344 = t317 * t351 + t352 * t372;
t341 = t318 * t353;
t340 = t321 * t353;
t339 = -t313 + t355;
t290 = t319 * t362 + t350;
t337 = t290 - t351;
t336 = t300 * t355 + t289;
t335 = t310 * (-t313 * t320 + qJD(1));
t334 = t320 * t361 + t358;
t333 = -t297 * t320 - t364 * t317 - pkin(1);
t332 = t318 * t308 * t307 + t306 * t354;
t281 = t306 * t340 - t320 * t345 - t308 * t358 + (t347 + t361) * qJD(1);
t282 = t320 * t346 + (t318 * t355 + t340) * t307 - t332;
t328 = t334 * qJD(6) + t368 * t281 - t363 * t282;
t283 = -t306 * t341 - t307 * t356 + t332 * t320 - t346;
t284 = t334 * qJD(1) - t307 * t341 - t308 * t347 - t345;
t327 = -(t306 * t321 - t307 * t357) * qJD(6) + t363 * t284 - t368 * t283;
t325 = t331 * t308;
t324 = qJD(1) * t301 - t290 * t320 + t300 * t353;
t323 = qJD(2) * t326;
t1 = [t337 * t321 - t368 * t284 - t363 * t283 - t371 * t318 + (-t365 * t318 + t333 * t321) * qJD(1) (t321 * t323 - t364 * t356) * t320 + (t322 * t321 - t326 * t356) * t317, t336 * t318 + t324 * t321 + t328 (t321 * t335 + (t339 * t318 + t340) * t309) * pkin(4) + t328, t328, -t281; t337 * t318 - t368 * t282 - t363 * t281 + t371 * t321 + (t333 * t318 + t365 * t321) * qJD(1) (t318 * t323 + t364 * t354) * t320 + (t322 * t318 + t326 * t354) * t317, t324 * t318 - t336 * t321 + t327 (t318 * t335 + (-t339 * t321 + t341) * t309) * pkin(4) + t327, t327, t283; 0, t317 * t323 - t322 * t320 (-t300 - t349) * t352 + (-t290 + t325) * t317 + t344 (-t349 - t367) * t352 + (t325 - t350) * t317 + t344, -t368 * t348 + (-t368 * t352 - t363 * t360) * t306 + t344, t306 * t352 + t348;];
JaD_transl  = t1;
