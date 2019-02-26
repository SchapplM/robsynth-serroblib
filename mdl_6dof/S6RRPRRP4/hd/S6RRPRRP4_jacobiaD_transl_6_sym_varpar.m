% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:47
% EndTime: 2019-02-26 21:47:47
% DurationCPUTime: 0.36s
% Computational Cost: add. (715->72), mult. (764->97), div. (0->0), fcn. (639->10), ass. (0->59)
t321 = cos(qJ(4));
t308 = t321 * pkin(4) + pkin(3);
t315 = qJ(2) + pkin(10);
t310 = sin(t315);
t311 = cos(t315);
t366 = r_i_i_C(2) + pkin(9) + pkin(8);
t339 = t366 * t311 - sin(qJ(2)) * pkin(2);
t342 = qJD(4) * t311 - qJD(1);
t316 = qJ(4) + qJ(5);
t312 = sin(t316);
t354 = qJD(6) * t312;
t318 = sin(qJ(4));
t368 = pkin(4) * t318;
t379 = (-t308 * t310 + t339) * qJD(2) - t342 * t368 - qJD(1) * (-qJ(3) - pkin(7)) + t311 * t354;
t313 = cos(t316);
t365 = r_i_i_C(3) + qJ(6);
t378 = t365 * t313;
t370 = pkin(5) + r_i_i_C(1);
t332 = -t365 * t312 - t370 * t313;
t329 = -t308 + t332;
t327 = t310 * t329 + t339;
t314 = qJD(4) + qJD(5);
t351 = t370 * t312;
t364 = pkin(4) * qJD(4);
t376 = (t351 - t378) * t314 + t318 * t364 - t354;
t372 = -t366 * t310 - cos(qJ(2)) * pkin(2);
t320 = sin(qJ(1));
t362 = t320 * t312;
t361 = t320 * t313;
t323 = cos(qJ(1));
t360 = t323 * t312;
t359 = t323 * t313;
t358 = qJD(1) * t320;
t357 = qJD(1) * t323;
t356 = qJD(2) * t310;
t355 = qJD(2) * t311;
t353 = t313 * qJD(6);
t352 = t321 * t364;
t350 = t314 * t362;
t349 = t314 * t359;
t348 = t310 * t353 + t355 * t378;
t347 = t312 * t355;
t345 = t320 * t356;
t344 = t323 * t356;
t341 = qJD(1) * t311 - qJD(4);
t340 = t342 * t321;
t336 = t311 * t359 + t362;
t334 = t313 * t358 + t314 * t360;
t333 = t312 * t357 + t314 * t361;
t289 = t312 * t344 - t311 * t349 - t350 + (t311 * t362 + t359) * qJD(1);
t290 = t311 * t334 + t313 * t344 - t333;
t331 = t336 * qJD(6) + t370 * t289 - t365 * t290;
t291 = t311 * t333 - t312 * t345 - t334;
t292 = qJD(1) * t336 - t311 * t350 - t313 * t345 - t349;
t330 = -(-t311 * t361 + t360) * qJD(6) + t365 * t292 - t370 * t291;
t328 = t332 * t314;
t326 = t352 - t353 + qJD(3) + (-t308 * t311 - pkin(1) + t372) * qJD(1);
t325 = t376 * t310 + (t311 * t329 + t372) * qJD(2);
t1 = [-t365 * t291 - t370 * t292 - t379 * t320 + t326 * t323, t325 * t323 - t327 * t358, t357 (-t323 * t340 + (t320 * t341 + t344) * t318) * pkin(4) + t331, t331, -t289; -t365 * t289 - t370 * t290 + t326 * t320 + t379 * t323, t320 * t325 + t327 * t357, t358 (-t320 * t340 + (-t323 * t341 + t345) * t318) * pkin(4) + t330, t330, t291; 0, t327 * qJD(2) - t376 * t311, 0 (-t351 - t368) * t355 + (t328 - t352) * t310 + t348, t310 * t328 - t370 * t347 + t348, t310 * t314 * t313 + t347;];
JaD_transl  = t1;
