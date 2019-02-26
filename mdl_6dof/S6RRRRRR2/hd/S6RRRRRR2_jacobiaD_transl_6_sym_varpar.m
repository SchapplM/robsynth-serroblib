% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:47:50
% EndTime: 2019-02-26 22:47:51
% DurationCPUTime: 0.47s
% Computational Cost: add. (912->91), mult. (706->119), div. (0->0), fcn. (533->12), ass. (0->79)
t307 = qJ(2) + qJ(3);
t302 = qJ(4) + t307;
t294 = sin(t302);
t306 = qJ(5) + qJ(6);
t298 = sin(t306);
t300 = cos(t306);
t303 = qJD(5) + qJD(6);
t359 = t300 * t303;
t295 = cos(t302);
t304 = qJD(2) + qJD(3);
t297 = qJD(4) + t304;
t363 = t295 * t297;
t384 = t294 * t359 + t298 * t363;
t314 = -pkin(11) - pkin(10);
t308 = sin(qJ(5));
t366 = pkin(5) * qJD(5);
t350 = t308 * t366;
t383 = t297 * t314 + t350;
t311 = cos(qJ(5));
t296 = pkin(5) * t311 + pkin(4);
t381 = r_i_i_C(1) * t300 + t296;
t382 = t294 * t381 + t295 * t314;
t309 = sin(qJ(2));
t367 = pkin(2) * qJD(2);
t347 = t309 * t367;
t299 = sin(t307);
t373 = pkin(3) * t304;
t353 = t299 * t373;
t365 = t294 * t297;
t368 = r_i_i_C(3) - t314;
t372 = pkin(5) * t308;
t380 = (t297 * t368 - t350) * t295 + (pkin(9) + pkin(8) + pkin(7) + t372) * qJD(1) - t296 * t365 - t347 - t353;
t364 = t294 * t298;
t343 = t303 * t364;
t379 = r_i_i_C(1) * t343 + r_i_i_C(2) * t384 + t294 * t383;
t313 = cos(qJ(1));
t334 = t295 * t303 - qJD(1);
t378 = t313 * t334;
t301 = cos(t307);
t321 = (-r_i_i_C(3) * t294 - t295 * t381) * t297;
t318 = -t301 * t373 + t321;
t356 = qJD(1) * t295;
t333 = -t303 + t356;
t310 = sin(qJ(1));
t346 = t310 * t365;
t376 = t313 * t333 - t346;
t374 = pkin(3) * t299;
t370 = r_i_i_C(2) * t300;
t369 = r_i_i_C(3) * t295;
t361 = t297 * t313;
t345 = t294 * t361;
t320 = t310 * t333 + t345;
t268 = t298 * t320 - t300 * t378;
t269 = t298 * t378 + t300 * t320;
t358 = r_i_i_C(1) * t268 + r_i_i_C(2) * t269;
t327 = t334 * t310;
t270 = t298 * t376 + t300 * t327;
t271 = t298 * t327 - t300 * t376;
t357 = -r_i_i_C(1) * t270 + r_i_i_C(2) * t271;
t355 = qJD(1) * t310;
t354 = qJD(1) * t313;
t351 = r_i_i_C(2) * t364;
t349 = t311 * t366;
t348 = r_i_i_C(2) * qJD(1) * t298;
t331 = -qJD(5) + t356;
t330 = -r_i_i_C(1) * t298 - t370;
t328 = t381 * t297;
t326 = (-qJD(5) * t295 + qJD(1)) * t311;
t325 = -t351 - t369;
t324 = t313 * t294 * t348 + t310 * t379 + t354 * t369;
t322 = t379 * t313 + t355 * t382;
t312 = cos(qJ(2));
t319 = -t312 * t367 + t318;
t317 = t349 + (-pkin(2) * t312 - pkin(3) * t301 - t294 * t368 - t295 * t296 - pkin(1)) * qJD(1);
t316 = t297 * t351 + r_i_i_C(3) * t363 - t294 * t328 + (t303 * t330 - t383) * t295;
t315 = t316 - t353;
t293 = -pkin(2) * t309 - t374;
t286 = r_i_i_C(2) * t343;
t1 = [t271 * r_i_i_C(1) + t270 * r_i_i_C(2) - t310 * t380 + t313 * t317 (-t293 + t325) * t355 + t319 * t313 + t322 (t325 + t374) * t355 + t318 * t313 + t322 (-r_i_i_C(3) * t361 - t310 * t348) * t294 + (-r_i_i_C(3) * t355 - t313 * t328) * t295 + t322 (t313 * t326 + (t310 * t331 + t345) * t308) * pkin(5) + t358, t358; -t269 * r_i_i_C(1) + t268 * r_i_i_C(2) + t310 * t317 + t313 * t380, t319 * t310 + (t293 - t382) * t354 + t324, t318 * t310 + (-t382 - t374) * t354 + t324, t310 * t321 - t354 * t382 + t324 (t310 * t326 + (-t313 * t331 + t346) * t308) * pkin(5) + t357, t357; 0, t315 - t347, t315, t316, t286 + (-r_i_i_C(1) * t359 - t349) * t294 + (t330 - t372) * t363, -r_i_i_C(1) * t384 - t363 * t370 + t286;];
JaD_transl  = t1;
