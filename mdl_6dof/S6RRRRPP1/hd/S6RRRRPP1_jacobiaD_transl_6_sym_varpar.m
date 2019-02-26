% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:54
% EndTime: 2019-02-26 22:24:55
% DurationCPUTime: 0.44s
% Computational Cost: add. (750->98), mult. (823->132), div. (0->0), fcn. (662->10), ass. (0->70)
t388 = pkin(5) + r_i_i_C(1);
t323 = qJ(4) + pkin(10);
t318 = sin(t323);
t382 = r_i_i_C(3) + qJ(6);
t396 = t382 * t318;
t319 = cos(t323);
t326 = sin(qJ(4));
t386 = pkin(4) * t326;
t395 = -t388 * t318 + t382 * t319 - t386;
t324 = qJ(2) + qJ(3);
t320 = sin(t324);
t321 = cos(t324);
t380 = pkin(4) * qJD(4);
t364 = t326 * t380;
t393 = qJD(5) * t321 + t320 * t364;
t328 = sin(qJ(1));
t374 = qJD(1) * t328;
t361 = t321 * t374;
t322 = qJD(2) + qJD(3);
t331 = cos(qJ(1));
t377 = t322 * t331;
t362 = t320 * t377;
t392 = t361 + t362;
t329 = cos(qJ(4));
t384 = t329 * pkin(4);
t316 = pkin(3) + t384;
t327 = sin(qJ(2));
t381 = pkin(2) * qJD(2);
t365 = t327 * t381;
t367 = t318 * qJD(6);
t325 = -qJ(5) - pkin(9);
t383 = r_i_i_C(2) - t325;
t391 = (t383 * t322 - t364 + t367) * t321 + (pkin(8) + pkin(7) + t386) * qJD(1) - (t316 * t322 - qJD(5)) * t320 - t365;
t387 = pkin(2) * t327;
t385 = r_i_i_C(2) * t321;
t379 = t321 * t322;
t378 = t322 * t328;
t376 = t328 * t318;
t375 = t331 * t319;
t373 = qJD(1) * t331;
t372 = qJD(4) * t319;
t371 = qJD(4) * t321;
t370 = qJD(4) * t328;
t369 = qJD(4) * t331;
t366 = t319 * qJD(6);
t363 = t320 * t378;
t360 = t320 * t374;
t358 = t318 * t370;
t357 = t318 * t369;
t356 = t319 * t369;
t352 = qJD(1) * t321 - qJD(4);
t347 = (qJD(1) - t371) * t329;
t345 = t388 * t320 * t358 + t325 * t363 + t393 * t328 + t373 * t385;
t343 = t321 * t375 + t376;
t342 = -t388 * t319 - t396;
t341 = t318 * t373 + t319 * t370;
t340 = -t382 * t372 - t367;
t339 = -t316 + t342;
t338 = t392 * t325 + t393 * t331 + t388 * (t319 * t360 + t320 * t357) + (t316 + t396) * t360;
t337 = t339 * t320 - t321 * t325;
t330 = cos(qJ(2));
t336 = t329 * t380 - t366 + (-t330 * pkin(2) - t316 * t321 - t383 * t320 - pkin(1)) * qJD(1);
t335 = t340 * t320 + (-r_i_i_C(2) * t320 + t339 * t321) * t322;
t334 = -t330 * t381 + t335;
t333 = r_i_i_C(2) * t379 + t320 * qJD(5) + t321 * t367 + t337 * t322 + t395 * t371;
t283 = t343 * qJD(1) - t319 * t363 - t321 * t358 - t356;
t282 = -t318 * t363 - t319 * t374 + t341 * t321 - t357;
t281 = t392 * t319 + t321 * t357 - t341;
t280 = t318 * t362 - t321 * t356 - t358 + (t321 * t376 + t375) * qJD(1);
t1 = [-t382 * t282 - t388 * t283 - t391 * t328 + t336 * t331, t338 + t334 * t331 + (-t385 + t387) * t374, -r_i_i_C(2) * t361 + t331 * t335 + t338, t343 * qJD(6) - t382 * t281 + t388 * t280 + (t331 * t347 + (t328 * t352 + t362) * t326) * pkin(4), t321 * t377 - t360, -t280; -t382 * t280 - t388 * t281 + t336 * t328 + t391 * t331 (t337 - t387) * t373 + t334 * t328 + t345 (-t325 * t373 + t339 * t378) * t321 + ((-r_i_i_C(2) * t322 + t340) * t328 + t339 * t373) * t320 + t345 -(-t328 * t321 * t319 + t331 * t318) * qJD(6) + t382 * t283 - t388 * t282 + (t328 * t347 + (-t331 * t352 + t363) * t326) * pkin(4), t320 * t373 + t321 * t378, t282; 0, t333 - t365, t333, t395 * t379 + (t366 + (t342 - t384) * qJD(4)) * t320, t322 * t320, t318 * t379 + t320 * t372;];
JaD_transl  = t1;
