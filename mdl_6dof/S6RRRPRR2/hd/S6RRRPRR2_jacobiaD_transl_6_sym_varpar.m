% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:34
% EndTime: 2019-02-26 22:16:34
% DurationCPUTime: 0.46s
% Computational Cost: add. (698->85), mult. (592->106), div. (0->0), fcn. (453->12), ass. (0->74)
t306 = qJ(2) + qJ(3);
t297 = pkin(11) + t306;
t295 = cos(t297);
t313 = -pkin(10) - pkin(9);
t294 = sin(t297);
t310 = cos(qJ(5));
t296 = pkin(5) * t310 + pkin(4);
t305 = qJ(5) + qJ(6);
t300 = cos(t305);
t377 = r_i_i_C(1) * t300 + t296;
t326 = t294 * t377;
t380 = t295 * t313 + t326;
t298 = sin(t305);
t303 = qJD(5) + qJD(6);
t355 = t300 * t303;
t304 = qJD(2) + qJD(3);
t357 = t295 * t304;
t379 = t294 * t355 + t298 * t357;
t307 = sin(qJ(5));
t360 = pkin(5) * qJD(5);
t347 = t307 * t360;
t378 = t304 * t313 + t347;
t359 = t294 * t298;
t348 = r_i_i_C(2) * t359;
t376 = r_i_i_C(3) * t295 + t348;
t299 = sin(t306);
t308 = sin(qJ(2));
t361 = pkin(2) * qJD(2);
t345 = t308 * t361;
t358 = t294 * t304;
t362 = r_i_i_C(3) - t313;
t366 = pkin(5) * t307;
t367 = pkin(3) * t304;
t375 = (t362 * t304 - t347) * t295 + (qJ(4) + pkin(8) + pkin(7) + t366) * qJD(1) - t296 * t358 - t299 * t367 - t345;
t344 = t303 * t359;
t374 = r_i_i_C(1) * t344 + t379 * r_i_i_C(2) + t378 * t294;
t312 = cos(qJ(1));
t332 = t295 * t303 - qJD(1);
t373 = t312 * t332;
t351 = qJD(1) * t295;
t331 = -t303 + t351;
t309 = sin(qJ(1));
t342 = t309 * t358;
t371 = t331 * t312 - t342;
t369 = pkin(3) * t299;
t301 = cos(t306);
t368 = pkin(3) * t301;
t364 = r_i_i_C(2) * t300;
t341 = t312 * t358;
t318 = t331 * t309 + t341;
t268 = t318 * t298 - t300 * t373;
t269 = t298 * t373 + t318 * t300;
t353 = t268 * r_i_i_C(1) + t269 * r_i_i_C(2);
t325 = t332 * t309;
t270 = t371 * t298 + t300 * t325;
t271 = t298 * t325 - t371 * t300;
t352 = -t270 * r_i_i_C(1) + t271 * r_i_i_C(2);
t350 = qJD(1) * t309;
t349 = qJD(1) * t312;
t346 = t310 * t360;
t329 = -qJD(5) + t351;
t328 = -r_i_i_C(1) * t298 - t364;
t324 = (-qJD(5) * t295 + qJD(1)) * t310;
t322 = t374 * t309 + t376 * t349;
t321 = -r_i_i_C(3) * t294 - t295 * t377;
t319 = t374 * t312 + t380 * t350;
t311 = cos(qJ(2));
t317 = -t301 * t367 + t321 * t304 - t311 * t361;
t316 = t304 * (t321 - t368);
t315 = t346 + qJD(4) + (-pkin(2) * t311 - t362 * t294 - t295 * t296 - pkin(1) - t368) * qJD(1);
t314 = r_i_i_C(3) * t357 + (t328 * t303 - t378) * t295 + (t348 - t326 - t369) * t304;
t293 = -pkin(2) * t308 - t369;
t283 = r_i_i_C(2) * t344;
t1 = [t271 * r_i_i_C(1) + t270 * r_i_i_C(2) - t375 * t309 + t315 * t312 (-t293 - t376) * t350 + t317 * t312 + t319 (-t376 + t369) * t350 + t312 * t316 + t319, t349 (t312 * t324 + (t329 * t309 + t341) * t307) * pkin(5) + t353, t353; -t269 * r_i_i_C(1) + t268 * r_i_i_C(2) + t315 * t309 + t375 * t312, t317 * t309 + (t293 - t380) * t349 + t322, t309 * t316 + (-t380 - t369) * t349 + t322, t350 (t309 * t324 + (-t329 * t312 + t342) * t307) * pkin(5) + t352, t352; 0, t314 - t345, t314, 0, t283 + (-r_i_i_C(1) * t355 - t346) * t294 + (t328 - t366) * t357, -t379 * r_i_i_C(1) - t357 * t364 + t283;];
JaD_transl  = t1;
