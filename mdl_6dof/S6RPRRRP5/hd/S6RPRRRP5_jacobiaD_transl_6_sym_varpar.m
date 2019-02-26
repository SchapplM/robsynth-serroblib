% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:20
% EndTime: 2019-02-26 21:10:20
% DurationCPUTime: 0.36s
% Computational Cost: add. (670->74), mult. (708->100), div. (0->0), fcn. (571->9), ass. (0->58)
t320 = sin(qJ(5));
t322 = cos(qJ(5));
t323 = cos(qJ(1));
t359 = qJD(5) * t323;
t321 = sin(qJ(1));
t362 = qJD(1) * t321;
t331 = t320 * t359 + t322 * t362;
t371 = r_i_i_C(1) + pkin(5);
t369 = r_i_i_C(3) + qJ(6);
t376 = t369 * t320;
t318 = pkin(10) + qJ(3);
t316 = qJ(4) + t318;
t313 = cos(t316);
t372 = pkin(9) + r_i_i_C(2);
t353 = t372 * t313;
t319 = qJD(3) + qJD(4);
t352 = t372 * t319;
t314 = sin(t318);
t368 = pkin(3) * qJD(3);
t356 = t314 * t368;
t358 = qJD(6) * t320;
t312 = sin(t316);
t367 = t312 * t319;
t375 = -pkin(4) * t367 + (t352 + t358) * t313 - t356;
t360 = qJD(5) * t322;
t329 = -t369 * t360 - t358;
t370 = pkin(3) * t314;
t366 = t313 * t319;
t365 = t313 * t321;
t364 = t321 * t320;
t363 = t323 * t322;
t361 = qJD(1) * t323;
t357 = t322 * qJD(6);
t355 = t371 * t320;
t354 = t372 * t312;
t351 = t321 * t367;
t350 = t323 * t367;
t345 = qJD(5) * t364;
t343 = t322 * t359;
t342 = qJD(2) - t357;
t341 = t371 * t312 * t345 + t361 * t353;
t336 = t313 * t363 + t364;
t315 = cos(t318);
t334 = -pkin(4) * t313 - pkin(3) * t315 - cos(pkin(10)) * pkin(2) - pkin(1) - t354;
t333 = (t371 * t331 + (pkin(4) + t376) * t362) * t312;
t332 = -t371 * t322 - t376;
t330 = t320 * t361 + t321 * t360;
t328 = -pkin(4) + t332;
t327 = t328 * t319;
t326 = t312 * t327 + t372 * t366 + (-qJD(5) * t355 - t329) * t313;
t325 = t329 * t312 + (t328 * t313 - t354) * t319;
t324 = -t315 * t368 + t325;
t317 = -pkin(8) - pkin(7) - qJ(2);
t284 = t336 * qJD(1) - t313 * t345 - t322 * t351 - t343;
t283 = t330 * t313 - t320 * t351 - t331;
t282 = t331 * t313 + t322 * t350 - t330;
t281 = t320 * t350 - t313 * t343 - t345 + (t313 * t364 + t363) * qJD(1);
t1 = [t342 * t323 - t371 * t284 - t369 * t283 - t375 * t321 + (t321 * t317 + t334 * t323) * qJD(1), t361 (-t353 + t370) * t362 + t324 * t323 + t333, t325 * t323 - t353 * t362 + t333, t336 * qJD(6) + t371 * t281 - t369 * t282, -t281; t342 * t321 - t371 * t282 - t369 * t281 + t375 * t323 + (-t323 * t317 + t334 * t321) * qJD(1), t362 (t328 * t312 - t370) * t361 + t324 * t321 + t341, t327 * t365 + ((-t352 + t329) * t321 + t328 * t361) * t312 + t341 -(t323 * t320 - t322 * t365) * qJD(6) + t369 * t284 - t371 * t283, t283; 0, 0, t326 - t356, t326 (t369 * t322 - t355) * t366 + (t332 * qJD(5) + t357) * t312, t312 * t360 + t320 * t366;];
JaD_transl  = t1;
