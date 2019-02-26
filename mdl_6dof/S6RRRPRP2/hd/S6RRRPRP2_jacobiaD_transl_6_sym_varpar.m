% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:10
% EndTime: 2019-02-26 22:10:11
% DurationCPUTime: 0.38s
% Computational Cost: add. (695->79), mult. (734->102), div. (0->0), fcn. (587->10), ass. (0->64)
t329 = sin(qJ(5));
t332 = cos(qJ(5));
t334 = cos(qJ(1));
t371 = qJD(5) * t334;
t331 = sin(qJ(1));
t374 = qJD(1) * t331;
t344 = t329 * t371 + t332 * t374;
t384 = r_i_i_C(1) + pkin(5);
t380 = r_i_i_C(3) + qJ(6);
t390 = t329 * t380;
t328 = qJ(2) + qJ(3);
t323 = pkin(10) + t328;
t322 = cos(t323);
t385 = pkin(9) + r_i_i_C(2);
t365 = t385 * t322;
t370 = qJD(6) * t329;
t372 = qJD(5) * t332;
t389 = -t372 * t380 - t370;
t324 = sin(t328);
t327 = qJD(2) + qJD(3);
t330 = sin(qJ(2));
t379 = pkin(2) * qJD(2);
t368 = t330 * t379;
t321 = sin(t323);
t378 = t321 * t327;
t381 = pkin(3) * t327;
t387 = -pkin(4) * t378 + (t327 * t385 + t370) * t322 - t324 * t381 - t368;
t383 = pkin(3) * t324;
t325 = cos(t328);
t382 = pkin(3) * t325;
t377 = t322 * t327;
t376 = t331 * t329;
t375 = t334 * t332;
t373 = qJD(1) * t334;
t369 = t332 * qJD(6);
t367 = t384 * t329;
t366 = t385 * t321;
t364 = t331 * t378;
t363 = t334 * t378;
t358 = qJD(5) * t376;
t356 = t332 * t371;
t355 = qJD(4) - t369;
t354 = t384 * t321 * t358 + t373 * t365;
t348 = t322 * t375 + t376;
t333 = cos(qJ(2));
t347 = -t333 * pkin(2) - pkin(4) * t322 - pkin(1) - t366 - t382;
t346 = (t384 * t344 + (pkin(4) + t390) * t374) * t321;
t345 = -t332 * t384 - t390;
t343 = t329 * t373 + t331 * t372;
t342 = -pkin(4) + t345;
t341 = t389 * t321;
t340 = t342 * t321;
t339 = t340 - t383;
t338 = t322 * t342 - t366;
t337 = t327 * t339 + t385 * t377 + (-qJD(5) * t367 - t389) * t322;
t336 = -t325 * t381 + t327 * t338 - t333 * t379 + t341;
t335 = t341 + (t338 - t382) * t327;
t326 = -qJ(4) - pkin(8) - pkin(7);
t313 = -t330 * pkin(2) - t383;
t290 = qJD(1) * t348 - t322 * t358 - t332 * t364 - t356;
t289 = t322 * t343 - t329 * t364 - t344;
t288 = t322 * t344 + t332 * t363 - t343;
t287 = t329 * t363 - t322 * t356 - t358 + (t322 * t376 + t375) * qJD(1);
t1 = [t355 * t334 - t384 * t290 - t380 * t289 - t387 * t331 + (t331 * t326 + t334 * t347) * qJD(1) (-t313 - t365) * t374 + t336 * t334 + t346 (-t365 + t383) * t374 + t335 * t334 + t346, t373, qJD(6) * t348 + t287 * t384 - t288 * t380, -t287; t355 * t331 - t384 * t288 - t380 * t287 + t387 * t334 + (-t334 * t326 + t331 * t347) * qJD(1) (t313 + t340) * t373 + t336 * t331 + t354, t331 * t335 + t339 * t373 + t354, t374 -(-t331 * t322 * t332 + t334 * t329) * qJD(6) + t380 * t290 - t384 * t289, t289; 0, t337 - t368, t337, 0 (t332 * t380 - t367) * t377 + (qJD(5) * t345 + t369) * t321, t321 * t372 + t329 * t377;];
JaD_transl  = t1;
