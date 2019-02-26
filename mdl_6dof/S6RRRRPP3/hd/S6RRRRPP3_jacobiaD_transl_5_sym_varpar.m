% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:23
% EndTime: 2019-02-26 22:26:23
% DurationCPUTime: 0.44s
% Computational Cost: add. (474->75), mult. (702->108), div. (0->0), fcn. (565->8), ass. (0->61)
t315 = sin(qJ(4));
t366 = r_i_i_C(3) + qJ(5);
t341 = t366 * t315;
t377 = pkin(3) + t341;
t314 = qJ(2) + qJ(3);
t312 = cos(t314);
t316 = sin(qJ(2));
t365 = pkin(2) * qJD(2);
t351 = t316 * t365;
t311 = sin(t314);
t313 = qJD(2) + qJD(3);
t364 = t311 * t313;
t355 = qJD(5) * t315;
t371 = pkin(9) + r_i_i_C(1);
t372 = t371 * t313 + t355;
t376 = -pkin(3) * t364 + t372 * t312 - t351;
t318 = cos(qJ(4));
t367 = r_i_i_C(2) * t318;
t375 = t311 * t367 + t371 * t312;
t370 = -r_i_i_C(2) + pkin(4);
t352 = t313 * t367;
t358 = qJD(4) * t315;
t368 = pkin(4) * t311;
t374 = t312 * t352 + t358 * t368;
t373 = t366 * t318;
t369 = pkin(2) * t316;
t363 = t312 * t313;
t317 = sin(qJ(1));
t362 = t317 * t315;
t320 = cos(qJ(1));
t361 = t320 * t318;
t360 = qJD(1) * t317;
t359 = qJD(1) * t320;
t357 = qJD(4) * t318;
t356 = qJD(4) * t320;
t354 = t318 * qJD(5);
t350 = t317 * t364;
t349 = t320 * t364;
t347 = t318 * t360;
t343 = t317 * t358;
t342 = t318 * t356;
t337 = t374 * t317 + t375 * t359;
t336 = t377 * t311 * t360 + t374 * t320 + t347 * t368;
t335 = t312 * t361 + t362;
t319 = cos(qJ(2));
t331 = -t319 * pkin(2) - pkin(3) * t312 - t371 * t311 - pkin(1);
t330 = -pkin(4) * t318 - t377;
t329 = t315 * t356 + t347;
t328 = t315 * t359 + t317 * t357;
t327 = t330 * t313;
t326 = t312 * t327;
t325 = (-r_i_i_C(2) * t315 - t373) * qJD(4) - t372;
t324 = t371 * t363 + (t327 + t352) * t311 + (t366 * t357 - t370 * t358 + t355) * t312;
t323 = t325 * t311 + t326;
t322 = -t319 * t365 + t323;
t321 = -pkin(8) - pkin(7);
t281 = t335 * qJD(1) - t312 * t343 - t318 * t350 - t342;
t280 = t328 * t312 - t315 * t350 - t329;
t279 = t329 * t312 + t318 * t349 - t328;
t278 = t315 * t349 - t312 * t342 - t343 + (t312 * t362 + t361) * qJD(1);
t1 = [-t320 * t354 - t370 * t281 - t366 * t280 - t376 * t317 + (t317 * t321 + t331 * t320) * qJD(1) (-t375 + t369) * t360 + t322 * t320 + t336, t323 * t320 - t360 * t375 + t336, t335 * qJD(5) + t370 * t278 - t366 * t279, -t278, 0; -t317 * t354 - t370 * t279 - t366 * t278 + t376 * t320 + (t331 * t317 - t320 * t321) * qJD(1) (t330 * t311 - t369) * t359 + t322 * t317 + t337, t317 * t326 + (t325 * t317 + t330 * t359) * t311 + t337 -(-t317 * t312 * t318 + t320 * t315) * qJD(5) + t366 * t281 - t370 * t280, t280, 0; 0, t324 - t351, t324 (-t370 * t315 + t373) * t363 + (t354 + (-t370 * t318 - t341) * qJD(4)) * t311, t311 * t357 + t315 * t363, 0;];
JaD_transl  = t1;
