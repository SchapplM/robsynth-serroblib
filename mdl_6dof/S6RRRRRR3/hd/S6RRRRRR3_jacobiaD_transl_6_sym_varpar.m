% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR3
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
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:19
% EndTime: 2019-02-26 22:48:20
% DurationCPUTime: 0.46s
% Computational Cost: add. (920->96), mult. (723->123), div. (0->0), fcn. (559->12), ass. (0->77)
t304 = qJ(4) + qJ(5);
t299 = qJ(6) + t304;
t291 = sin(t299);
t292 = cos(t299);
t305 = qJ(2) + qJ(3);
t298 = cos(t305);
t302 = qJD(2) + qJD(3);
t352 = t298 * t302;
t301 = qJD(4) + qJD(5);
t294 = qJD(6) + t301;
t296 = sin(t305);
t354 = t294 * t296;
t373 = t291 * t352 + t292 * t354;
t303 = -pkin(11) - pkin(10) - pkin(9);
t297 = cos(t304);
t309 = cos(qJ(4));
t287 = t309 * pkin(4) + pkin(5) * t297;
t285 = pkin(3) + t287;
t360 = r_i_i_C(1) * t292;
t371 = t285 + t360;
t372 = t296 * t371 + t298 * t303;
t342 = r_i_i_C(2) * t291 * t296;
t370 = r_i_i_C(3) * t298 + t342;
t340 = t291 * t354;
t369 = r_i_i_C(1) * t340 + t373 * r_i_i_C(2);
t295 = sin(t304);
t306 = sin(qJ(4));
t355 = pkin(4) * qJD(4);
t361 = pkin(5) * t301;
t277 = -t295 * t361 - t306 * t355;
t362 = pkin(5) * t295;
t286 = t306 * pkin(4) + t362;
t307 = sin(qJ(2));
t356 = pkin(2) * qJD(2);
t341 = t307 * t356;
t353 = t296 * t302;
t357 = r_i_i_C(3) - t303;
t368 = (t357 * t302 + t277) * t298 + (t286 + pkin(8) + pkin(7)) * qJD(1) - t285 * t353 - t341;
t311 = cos(qJ(1));
t330 = t294 * t298 - qJD(1);
t367 = t311 * t330;
t347 = qJD(1) * t298;
t329 = -t294 + t347;
t308 = sin(qJ(1));
t337 = t308 * t353;
t365 = t329 * t311 - t337;
t363 = pkin(2) * t307;
t359 = r_i_i_C(2) * t292;
t336 = t311 * t353;
t317 = t329 * t308 + t336;
t264 = t317 * t291 - t292 * t367;
t265 = t291 * t367 + t317 * t292;
t350 = t264 * r_i_i_C(1) + t265 * r_i_i_C(2);
t323 = t330 * t308;
t266 = t365 * t291 + t292 * t323;
t267 = t291 * t323 - t365 * t292;
t349 = -t266 * r_i_i_C(1) + t267 * r_i_i_C(2);
t346 = qJD(1) * t308;
t345 = qJD(1) * t311;
t344 = t297 * t361;
t343 = t294 * t360;
t328 = -t301 + t347;
t327 = -r_i_i_C(1) * t291 - t359;
t326 = t286 * t347 + t277;
t325 = t371 * t302;
t324 = t297 * (-t298 * t301 + qJD(1));
t321 = t303 * t337 + t369 * t308 + t370 * t345;
t320 = t303 * t336 + t369 * t311 + t372 * t346;
t278 = t309 * t355 + t344;
t318 = qJD(1) * t287 - t278 * t298 + t286 * t353;
t310 = cos(qJ(2));
t316 = t278 + (-t310 * pkin(2) - t285 * t298 - t357 * t296 - pkin(1)) * qJD(1);
t315 = -t277 * t296 + (-r_i_i_C(3) * t296 - t298 * t371) * t302;
t314 = -t310 * t356 + t315;
t313 = r_i_i_C(3) * t352 - t296 * t325 + t302 * t342 + (t327 * t294 - t302 * t303 + t277) * t298;
t276 = r_i_i_C(2) * t340;
t1 = [t267 * r_i_i_C(1) + t266 * r_i_i_C(2) - t368 * t308 + t316 * t311 (-t370 + t363) * t346 + t314 * t311 + t320, t315 * t311 - t346 * t370 + t320, t326 * t308 + t318 * t311 + t350 (t311 * t324 + (t328 * t308 + t336) * t295) * pkin(5) + t350, t350; -t265 * r_i_i_C(1) + t264 * r_i_i_C(2) + t316 * t308 + t368 * t311 (-t363 - t372) * t345 + t314 * t308 + t321 (-t303 * t345 - t308 * t325) * t298 + ((-r_i_i_C(3) * t302 - t277) * t308 - t371 * t345) * t296 + t321, t318 * t308 - t326 * t311 + t349 (t308 * t324 + (-t328 * t311 + t337) * t295) * pkin(5) + t349, t349; 0, t313 - t341, t313, t276 + (-t278 - t343) * t296 + (-t286 + t327) * t352, t276 + (-t343 - t344) * t296 + (t327 - t362) * t352, -t373 * r_i_i_C(1) - t352 * t359 + t276;];
JaD_transl  = t1;
