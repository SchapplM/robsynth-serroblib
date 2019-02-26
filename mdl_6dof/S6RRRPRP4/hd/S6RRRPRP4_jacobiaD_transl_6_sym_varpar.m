% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:06
% EndTime: 2019-02-26 22:11:07
% DurationCPUTime: 0.40s
% Computational Cost: add. (538->80), mult. (769->116), div. (0->0), fcn. (616->8), ass. (0->65)
t382 = r_i_i_C(1) + pkin(5);
t379 = r_i_i_C(3) + qJ(6);
t362 = pkin(3) + pkin(9) + r_i_i_C(2);
t315 = qJ(2) + qJ(3);
t313 = cos(t315);
t314 = qJD(2) + qJD(3);
t377 = t313 * t314;
t306 = qJ(4) * t377;
t312 = sin(t315);
t319 = cos(qJ(5));
t364 = qJD(6) * t319;
t331 = -t362 * t314 - t364;
t317 = sin(qJ(2));
t378 = pkin(2) * qJD(2);
t361 = t317 * t378;
t386 = (qJD(4) + t331) * t312 + (pkin(4) + pkin(8) + pkin(7)) * qJD(1) + t306 - t361;
t316 = sin(qJ(5));
t385 = t382 * t316;
t381 = pkin(2) * t317;
t318 = sin(qJ(1));
t376 = t313 * t318;
t375 = t314 * t312;
t321 = cos(qJ(1));
t374 = t314 * t321;
t373 = t316 * t318;
t372 = t316 * t321;
t371 = t318 * t319;
t370 = t319 * t321;
t369 = qJD(1) * t318;
t368 = qJD(1) * t321;
t367 = qJD(4) * t313;
t366 = qJD(5) * t316;
t365 = qJD(5) * t319;
t363 = t316 * qJD(6);
t359 = t314 * t376;
t358 = t313 * t374;
t357 = t314 * t371;
t356 = t314 * t370;
t355 = t312 * t369;
t354 = t312 * t368;
t353 = t313 * t368;
t350 = t313 * t366;
t349 = qJD(5) * t376;
t348 = t321 * t365;
t347 = t379 * t319;
t346 = t362 * t313;
t345 = qJD(5) * t312 + qJD(1);
t344 = qJD(1) * t312 + qJD(5);
t335 = -qJ(4) - t385;
t334 = t345 * t319;
t333 = t335 * t312;
t332 = t312 * t372 + t371;
t329 = -t362 * t312 - t313 * t347;
t328 = qJ(4) * t353 + t318 * t367 + t382 * (t316 * t353 + t319 * t349) + t379 * (t312 * t357 + t316 * t349);
t327 = t382 * t313 * t348 + t321 * t367 + t362 * t355 + t379 * (t313 * t319 * t369 + t312 * t356 + t321 * t350);
t320 = cos(qJ(2));
t326 = t363 + (-pkin(2) * t320 - qJ(4) * t312 - pkin(1) - t346) * qJD(1);
t325 = -t313 * t364 + (-t346 + t333) * t314;
t324 = -t320 * t378 + t325;
t323 = t329 * t314 + t306 + t382 * (t312 * t365 + t316 * t377) + (t379 * t366 + qJD(4) - t364) * t312;
t276 = t321 * t334 + (-t344 * t318 + t358) * t316;
t275 = -t313 * t356 + t332 * qJD(5) + (t312 * t371 + t372) * qJD(1);
t274 = t318 * t334 + (t344 * t321 + t359) * t316;
t273 = -t313 * t357 - t319 * t354 + t345 * t373 - t348;
t1 = [-t379 * t273 - t382 * t274 - t386 * t318 + t326 * t321 (t335 * t313 + t381) * t369 + t327 + t324 * t321, t327 + t333 * t374 + (t331 * t321 + t335 * t369) * t313, -t355 + t358, t332 * qJD(6) - t382 * t275 + t379 * t276, t275; t379 * t275 + t382 * t276 + t326 * t318 + t386 * t321 (t329 - t381) * t368 + t324 * t318 + t328, t325 * t318 + t329 * t368 + t328, t354 + t359 -(-t312 * t373 + t370) * qJD(6) + t379 * t274 - t382 * t273, t273; 0, t323 - t361, t323, t375 (t379 * t316 + t382 * t319) * t375 + (-t363 + (-t347 + t385) * qJD(5)) * t313, -t319 * t375 - t350;];
JaD_transl  = t1;
