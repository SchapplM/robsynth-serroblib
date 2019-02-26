% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:24
% EndTime: 2019-02-26 22:21:24
% DurationCPUTime: 0.62s
% Computational Cost: add. (752->119), mult. (1664->184), div. (0->0), fcn. (1558->10), ass. (0->72)
t323 = qJD(5) + qJD(6);
t327 = sin(qJ(2));
t392 = (qJD(3) - t323) * t327;
t326 = sin(qJ(3));
t330 = cos(qJ(3));
t324 = qJ(5) + qJ(6);
t322 = cos(t324);
t321 = sin(t324);
t355 = -t321 * r_i_i_C(1) - qJ(4);
t341 = t322 * r_i_i_C(2) - t355;
t325 = sin(qJ(5));
t383 = t325 * pkin(5);
t337 = t341 + t383;
t329 = cos(qJ(5));
t381 = -t329 * pkin(5) - pkin(3) - pkin(4);
t350 = t321 * r_i_i_C(2) + t381;
t340 = t322 * r_i_i_C(1) - t350;
t343 = t325 * t330 - t326 * t329;
t344 = t321 * t326 + t322 * t330;
t345 = t321 * t330 - t322 * t326;
t362 = -r_i_i_C(3) - pkin(10) - pkin(9) + pkin(8);
t380 = pkin(5) * qJD(5);
t334 = -t362 * qJD(2) + (t326 * t340 - t330 * t337) * qJD(3) - t326 * qJD(4) + (r_i_i_C(1) * t345 + r_i_i_C(2) * t344) * t323 + t343 * t380;
t331 = cos(qJ(2));
t367 = qJD(2) * t331;
t371 = -(-t344 * t392 + t345 * t367) * r_i_i_C(1) - (t344 * t367 + t345 * t392) * r_i_i_C(2);
t390 = -t327 * pkin(2) + t362 * t331;
t328 = sin(qJ(1));
t365 = qJD(3) * t330;
t332 = cos(qJ(1));
t369 = qJD(1) * t332;
t339 = t326 * t369 + t328 * t365;
t364 = qJD(3) * t332;
t358 = t326 * t364;
t366 = qJD(2) * t332;
t359 = t327 * t366;
t370 = qJD(1) * t328;
t301 = t331 * t358 + (t331 * t370 + t359) * t330 - t339;
t375 = t332 * t326;
t306 = -t328 * t330 + t331 * t375;
t374 = t332 * t330;
t377 = t328 * t326;
t307 = t331 * t374 + t377;
t389 = -t301 * t325 + (t306 * t325 + t307 * t329) * qJD(5);
t357 = t330 * t364;
t360 = qJD(3) * t377;
t368 = qJD(2) * t328;
t361 = t327 * t368;
t303 = qJD(1) * t307 - t330 * t361 - t331 * t360 - t357;
t376 = t328 * t331;
t304 = t326 * t376 + t374;
t305 = t330 * t376 - t375;
t388 = t303 * t325 + (t304 * t325 + t305 * t329) * qJD(5);
t385 = t326 * t337 + t330 * t340 + pkin(2);
t353 = t306 * t323 - t301;
t300 = qJD(1) * t304 + t326 * t359 - t331 * t357 - t360;
t354 = -t307 * t323 - t300;
t293 = -t321 * t353 + t322 * t354;
t294 = t321 * t354 + t322 * t353;
t373 = t293 * r_i_i_C(1) - t294 * r_i_i_C(2);
t351 = -t304 * t323 - t303;
t302 = -t326 * t361 - t370 * t330 + t331 * t339 - t358;
t352 = t305 * t323 - t302;
t372 = (t321 * t351 - t322 * t352) * r_i_i_C(1) + (t321 * t352 + t322 * t351) * r_i_i_C(2);
t356 = qJ(4) + t383;
t349 = t304 * t322 - t305 * t321;
t348 = t304 * t321 + t305 * t322;
t342 = t325 * t326 + t329 * t330;
t338 = -pkin(2) * t331 - t327 * t362 - pkin(1);
t335 = qJD(2) * t385;
t299 = t302 * t322;
t1 = [-t299 * r_i_i_C(2) - t304 * qJD(4) + t355 * t302 + (-r_i_i_C(1) * t349 + r_i_i_C(2) * t348) * t323 - t340 * t303 + (-t302 * t325 + (-t304 * t329 + t305 * t325) * qJD(5)) * pkin(5) - t390 * t368 + (-t328 * pkin(7) + t332 * t338) * qJD(1) (-t332 * t335 - t362 * t370) * t331 + (t334 * t332 + t385 * t370) * t327, t307 * qJD(4) - t341 * t301 + ((t306 * t321 + t307 * t322) * r_i_i_C(1) + (t306 * t322 - t307 * t321) * r_i_i_C(2)) * t323 + t340 * t300 + t389 * pkin(5), -t300 (-t300 * t329 - t389) * pkin(5) + t373, t373; t294 * r_i_i_C(1) + t293 * r_i_i_C(2) - t300 * qJ(4) + t306 * qJD(4) + t381 * t301 + (-t300 * t325 + (t306 * t329 - t307 * t325) * qJD(5)) * pkin(5) + t390 * t366 + (pkin(7) * t332 + t328 * t338) * qJD(1) (-t328 * t335 + t362 * t369) * t331 + (t328 * t334 - t369 * t385) * t327, -t299 * r_i_i_C(1) + t305 * qJD(4) + t341 * t303 + (r_i_i_C(1) * t348 + r_i_i_C(2) * t349) * t323 + t350 * t302 + t388 * pkin(5), t302 (t302 * t329 - t388) * pkin(5) + t372, t372; 0, -t327 * t335 - t334 * t331 (t326 * t381 + t330 * t356) * t367 + (qJD(4) * t330 + t342 * t380 + (-t326 * t356 + t330 * t381) * qJD(3)) * t327 - t371, t326 * t367 + t327 * t365 (-t343 * t367 + (qJD(3) - qJD(5)) * t327 * t342) * pkin(5) + t371, t371;];
JaD_transl  = t1;
