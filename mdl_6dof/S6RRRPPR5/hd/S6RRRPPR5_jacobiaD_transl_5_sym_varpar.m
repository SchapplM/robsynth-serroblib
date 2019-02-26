% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:47
% EndTime: 2019-02-26 22:05:48
% DurationCPUTime: 0.47s
% Computational Cost: add. (550->95), mult. (1188->149), div. (0->0), fcn. (1158->12), ass. (0->56)
t345 = sin(qJ(1));
t341 = cos(pkin(6));
t358 = qJD(2) * t341 + qJD(1);
t344 = sin(qJ(2));
t374 = t345 * t344;
t365 = t341 * t374;
t369 = qJD(2) * t344;
t347 = cos(qJ(2));
t348 = cos(qJ(1));
t371 = t348 * t347;
t318 = -qJD(1) * t365 - t345 * t369 + t358 * t371;
t339 = sin(pkin(6));
t375 = t339 * t348;
t386 = -qJD(3) * t375 + t318;
t372 = t348 * t344;
t373 = t345 * t347;
t322 = t341 * t372 + t373;
t370 = qJD(1) * t339;
t363 = t345 * t370;
t385 = -qJD(3) * t322 + t363;
t337 = qJ(3) + pkin(11);
t335 = sin(t337);
t336 = cos(t337);
t343 = sin(qJ(3));
t338 = sin(pkin(12));
t340 = cos(pkin(12));
t357 = t340 * r_i_i_C(1) - t338 * r_i_i_C(2) + pkin(4);
t379 = r_i_i_C(3) + qJ(5);
t384 = (pkin(3) * t343 + t357 * t335 - t379 * t336) * qJD(3) - t335 * qJD(5);
t383 = t385 * t335 + t386 * t336;
t346 = cos(qJ(3));
t334 = t346 * pkin(3) + pkin(2);
t382 = t379 * t335 + t357 * t336 + t334;
t378 = t339 * t344;
t377 = t339 * t345;
t376 = t339 * t346;
t368 = qJD(2) * t347;
t364 = t341 * t371;
t362 = t348 * t370;
t361 = t339 * t368;
t342 = -qJ(4) - pkin(9);
t356 = t338 * r_i_i_C(1) + t340 * r_i_i_C(2) - t342;
t324 = -t365 + t371;
t355 = -t324 * t335 + t336 * t377;
t354 = t324 * t336 + t335 * t377;
t353 = t341 * t335 + t336 * t378;
t352 = t341 * t373 + t372;
t311 = t386 * t335 - t385 * t336;
t321 = -t364 + t374;
t319 = t353 * qJD(3) + t335 * t361;
t317 = t352 * qJD(1) + t322 * qJD(2);
t316 = t322 * qJD(1) + t352 * qJD(2);
t315 = -qJD(1) * t364 - t348 * t368 + t358 * t374;
t310 = t355 * qJD(3) - t316 * t336 + t335 * t362;
t309 = t354 * qJD(3) - t316 * t335 - t336 * t362;
t1 = [(-t317 * t338 - t340 * t383) * r_i_i_C(1) + (-t317 * t340 + t338 * t383) * r_i_i_C(2) - t383 * pkin(4) - (t322 * t335 + t336 * t375) * qJD(5) - t318 * t334 + t317 * t342 - t321 * qJD(4) - t379 * t311 + (-t348 * pkin(1) - pkin(8) * t377) * qJD(1) + (-t343 * t363 + (t322 * t343 + t346 * t375) * qJD(3)) * pkin(3), t324 * qJD(4) + t382 * t315 - t356 * t316 + t352 * t384, t354 * qJD(5) + t379 * t310 - t357 * t309 + (t346 * t362 + t316 * t343 + (-t324 * t346 - t343 * t377) * qJD(3)) * pkin(3), -t315, t309, 0; (t310 * t340 - t315 * t338) * r_i_i_C(1) + (-t310 * t338 - t315 * t340) * r_i_i_C(2) + t310 * pkin(4) - t355 * qJD(5) - t316 * t334 + t315 * t342 + t352 * qJD(4) + t379 * t309 + (-t345 * pkin(1) + pkin(8) * t375) * qJD(1) + (t343 * t362 + (-t324 * t343 + t345 * t376) * qJD(3)) * pkin(3), t322 * qJD(4) - t317 * t382 + t356 * t318 + t384 * t321 -(-t322 * t336 + t335 * t375) * qJD(5) + t379 * t383 - t357 * t311 + (t346 * t363 - t318 * t343 + (-t322 * t346 + t343 * t375) * qJD(3)) * pkin(3), t317, t311, 0; 0 ((-qJD(2) * t382 + qJD(4)) * t344 + (t356 * qJD(2) - t384) * t347) * t339, t353 * qJD(5) + t379 * (t336 * t361 + (-t335 * t378 + t336 * t341) * qJD(3)) - t357 * t319 + (-t343 * t361 + (-t341 * t343 - t344 * t376) * qJD(3)) * pkin(3), t339 * t369, t319, 0;];
JaD_transl  = t1;
