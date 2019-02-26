% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP8
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
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:12:05
% EndTime: 2019-02-26 21:12:05
% DurationCPUTime: 0.46s
% Computational Cost: add. (478->69), mult. (714->101), div. (0->0), fcn. (573->8), ass. (0->60)
t306 = sin(qJ(5));
t366 = -r_i_i_C(1) - pkin(5);
t381 = t306 * t366;
t380 = qJD(5) * t381;
t309 = cos(qJ(5));
t364 = r_i_i_C(3) + qJ(6);
t321 = -t306 * t364 + t309 * t366;
t319 = -pkin(4) + t321;
t305 = qJ(3) + qJ(4);
t303 = cos(t305);
t379 = t303 * t380;
t304 = qJD(3) + qJD(4);
t352 = qJD(6) * t306;
t367 = pkin(9) + r_i_i_C(2);
t378 = t304 * t367 + t352;
t335 = t364 * t309;
t376 = -qJD(5) * t335 - t378;
t374 = t303 * t304;
t302 = sin(t305);
t307 = sin(qJ(3));
t368 = -pkin(3) * t307 - pkin(4) * t302 + t303 * t367 - qJ(2);
t365 = -pkin(1) - pkin(8) - pkin(7);
t363 = pkin(3) * qJD(3);
t362 = t302 * t304;
t308 = sin(qJ(1));
t361 = t303 * t308;
t360 = t303 * t309;
t311 = cos(qJ(1));
t359 = t304 * t311;
t358 = t308 * t306;
t357 = t308 * t309;
t356 = t311 * t306;
t355 = qJD(1) * t308;
t354 = qJD(1) * t311;
t351 = t309 * qJD(6);
t310 = cos(qJ(3));
t350 = pkin(3) * qJD(1) * t310;
t349 = t307 * t363;
t348 = t310 * t363;
t347 = t304 * t361;
t346 = t311 * t302 * t309;
t344 = t303 * t359;
t340 = t303 * t354;
t339 = t309 * t354;
t336 = qJD(5) * t360;
t334 = qJD(5) * t302 + qJD(1);
t333 = qJD(1) * t302 + qJD(5);
t324 = t333 * t311;
t322 = t302 * t357 + t356;
t318 = -t366 * t303 * t339 + pkin(4) * t340 + t352 * t361 + t364 * (t306 * t340 + t308 * t336) + t367 * (t302 * t354 + t347);
t317 = t376 * t303;
t316 = t367 * t302 * t355 - t311 * t379 - t319 * (t302 * t359 + t303 * t355);
t315 = pkin(4) * t374 + t378 * t302 + qJD(2) + t348;
t314 = t319 * t362 + t379;
t313 = t319 * t374 + (t376 - t380) * t302;
t269 = t309 * t324 + (t304 * t360 - t306 * t334) * t308;
t268 = t334 * t357 + (t324 + t347) * t306;
t267 = -t309 * t344 + (t302 * t356 + t357) * qJD(5) + t322 * qJD(1);
t266 = -qJD(5) * t346 - t306 * t344 + t333 * t358 - t339;
t1 = [t308 * t351 + t366 * t267 - t364 * t266 + t315 * t311 + (t368 * t308 + t365 * t311) * qJD(1), t354, t318 + t311 * t350 + (t314 - t349) * t308, t308 * t314 + t318, qJD(6) * t322 + t268 * t366 + t269 * t364, t268; -t311 * t351 - t366 * t269 + t364 * t268 + t315 * t308 + (t365 * t308 - t368 * t311) * qJD(1), t355, t316 + t308 * t350 + (t317 + t349) * t311, t311 * t317 + t316 -(t346 - t358) * qJD(6) + t364 * t267 + t366 * t266, t266; 0, 0, t313 - t348, t313 (-t335 - t381) * t362 + (qJD(5) * t321 + t351) * t303, -t306 * t362 + t336;];
JaD_transl  = t1;
