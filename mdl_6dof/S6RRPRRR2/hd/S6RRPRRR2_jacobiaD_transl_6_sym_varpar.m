% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:41
% EndTime: 2019-02-26 21:54:42
% DurationCPUTime: 0.44s
% Computational Cost: add. (680->82), mult. (578->111), div. (0->0), fcn. (444->12), ass. (0->68)
t302 = qJ(2) + pkin(11);
t296 = qJ(4) + t302;
t291 = sin(t296);
t303 = qJ(5) + qJ(6);
t297 = sin(t303);
t298 = cos(t303);
t300 = qJD(5) + qJD(6);
t352 = t298 * t300;
t292 = cos(t296);
t301 = qJD(2) + qJD(4);
t354 = t292 * t301;
t372 = t291 * t352 + t297 * t354;
t310 = -pkin(10) - pkin(9);
t304 = sin(qJ(5));
t357 = pkin(5) * qJD(5);
t343 = t304 * t357;
t371 = t301 * t310 + t343;
t307 = cos(qJ(5));
t293 = t307 * pkin(5) + pkin(4);
t369 = r_i_i_C(1) * t298 + t293;
t370 = t291 * t369 + t292 * t310;
t289 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t302);
t284 = t289 * qJD(2);
t355 = t291 * t301;
t358 = r_i_i_C(3) - t310;
t362 = pkin(5) * t304;
t368 = (t358 * t301 - t343) * t292 + (pkin(8) + qJ(3) + pkin(7) + t362) * qJD(1) - t293 * t355 + t284;
t356 = t291 * t297;
t340 = t300 * t356;
t367 = r_i_i_C(1) * t340 + t372 * r_i_i_C(2) + t371 * t291;
t309 = cos(qJ(1));
t328 = t292 * t300 - qJD(1);
t366 = t309 * t328;
t347 = qJD(1) * t292;
t327 = -t300 + t347;
t306 = sin(qJ(1));
t338 = t306 * t355;
t364 = t327 * t309 - t338;
t360 = r_i_i_C(2) * t298;
t359 = r_i_i_C(3) * t292;
t351 = t301 * t309;
t337 = t291 * t351;
t314 = t327 * t306 + t337;
t265 = t314 * t297 - t298 * t366;
t266 = t297 * t366 + t314 * t298;
t349 = t265 * r_i_i_C(1) + t266 * r_i_i_C(2);
t320 = t328 * t306;
t267 = t364 * t297 + t298 * t320;
t268 = t297 * t320 - t364 * t298;
t348 = -t267 * r_i_i_C(1) + t268 * r_i_i_C(2);
t346 = qJD(1) * t306;
t345 = qJD(1) * t309;
t344 = r_i_i_C(2) * t356;
t342 = t307 * t357;
t341 = r_i_i_C(2) * qJD(1) * t297;
t325 = -qJD(5) + t347;
t324 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t302);
t323 = -r_i_i_C(1) * t297 - t360;
t321 = t369 * t301;
t319 = (-qJD(5) * t292 + qJD(1)) * t307;
t318 = t309 * t291 * t341 + t367 * t306 + t345 * t359;
t316 = t367 * t309 + t370 * t346;
t315 = (-r_i_i_C(3) * t291 - t292 * t369) * t301;
t313 = t324 * qJD(2) + t315;
t312 = t342 + qJD(3) + (-t358 * t291 - t292 * t293 - pkin(1) + t324) * qJD(1);
t311 = t301 * t344 + r_i_i_C(3) * t354 - t291 * t321 + (t323 * t300 - t371) * t292;
t278 = r_i_i_C(2) * t340;
t1 = [t268 * r_i_i_C(1) + t267 * r_i_i_C(2) - t368 * t306 + t312 * t309 (-t289 - t344 - t359) * t346 + t313 * t309 + t316, t345 (-r_i_i_C(3) * t351 - t306 * t341) * t291 + (-r_i_i_C(3) * t346 - t309 * t321) * t292 + t316 (t309 * t319 + (t325 * t306 + t337) * t304) * pkin(5) + t349, t349; -t266 * r_i_i_C(1) + t265 * r_i_i_C(2) + t312 * t306 + t368 * t309, t313 * t306 + (t289 - t370) * t345 + t318, t346, t306 * t315 - t345 * t370 + t318 (t306 * t319 + (-t325 * t309 + t338) * t304) * pkin(5) + t348, t348; 0, t284 + t311, 0, t311, t278 + (-r_i_i_C(1) * t352 - t342) * t291 + (t323 - t362) * t354, -t372 * r_i_i_C(1) - t354 * t360 + t278;];
JaD_transl  = t1;
