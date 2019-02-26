% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:29
% EndTime: 2019-02-26 20:06:29
% DurationCPUTime: 0.40s
% Computational Cost: add. (654->81), mult. (1172->133), div. (0->0), fcn. (1186->14), ass. (0->65)
t343 = sin(qJ(3));
t345 = cos(qJ(3));
t339 = pkin(12) + qJ(5);
t336 = cos(t339);
t337 = qJ(6) + t339;
t333 = sin(t337);
t334 = cos(t337);
t361 = t334 * r_i_i_C(1) - t333 * r_i_i_C(2);
t357 = pkin(5) * t336 + cos(pkin(12)) * pkin(4) + pkin(3) + t361;
t386 = r_i_i_C(3) + pkin(10) + pkin(9) + qJ(4);
t335 = sin(t339);
t340 = qJD(5) + qJD(6);
t360 = t333 * r_i_i_C(1) + t334 * r_i_i_C(2);
t385 = pkin(5) * qJD(5);
t388 = t335 * t385 + t360 * t340;
t347 = t388 * t345 + (t357 * t343 - t386 * t345) * qJD(3) - t343 * qJD(4);
t341 = sin(pkin(11));
t344 = sin(qJ(2));
t346 = cos(qJ(2));
t383 = cos(pkin(11));
t384 = cos(pkin(6));
t358 = t384 * t383;
t322 = t341 * t346 + t344 * t358;
t342 = sin(pkin(6));
t369 = t342 * t383;
t312 = t322 * t345 - t343 * t369;
t370 = t341 * t384;
t324 = -t344 * t370 + t383 * t346;
t381 = t342 * t343;
t380 = t342 * t345;
t379 = t342 * t346;
t318 = t322 * qJD(2);
t365 = t312 * t340 - t318;
t356 = t346 * t358;
t375 = qJD(2) * t344;
t317 = -qJD(2) * t356 + t341 * t375;
t352 = -t322 * t343 - t345 * t369;
t308 = t352 * qJD(3) - t317 * t345;
t321 = t341 * t344 - t356;
t367 = -t321 * t340 - t308;
t378 = (t367 * t333 - t365 * t334) * r_i_i_C(1) + (t365 * t333 + t367 * t334) * r_i_i_C(2);
t314 = t324 * t345 + t341 * t381;
t320 = t324 * qJD(2);
t364 = t314 * t340 - t320;
t323 = t383 * t344 + t346 * t370;
t319 = t323 * qJD(2);
t355 = -t324 * t343 + t341 * t380;
t310 = t355 * qJD(3) - t319 * t345;
t366 = -t323 * t340 - t310;
t377 = (t366 * t333 - t364 * t334) * r_i_i_C(1) + (t364 * t333 + t366 * t334) * r_i_i_C(2);
t326 = t384 * t343 + t344 * t380;
t372 = t342 * t375;
t354 = -t326 * t340 + t372;
t351 = -t344 * t381 + t384 * t345;
t371 = qJD(2) * t379;
t316 = t351 * qJD(3) + t345 * t371;
t359 = t340 * t379 - t316;
t376 = (t359 * t333 + t354 * t334) * r_i_i_C(1) + (-t354 * t333 + t359 * t334) * r_i_i_C(2);
t353 = pkin(8) + pkin(5) * t335 + sin(pkin(12)) * pkin(4) + t360;
t349 = t336 * t385 + t361 * t340;
t348 = -t386 * t343 - t357 * t345 - pkin(2);
t315 = t326 * qJD(3) + t343 * t371;
t309 = t314 * qJD(3) - t319 * t343;
t307 = t312 * qJD(3) - t317 * t343;
t1 = [0, -t353 * t319 + t348 * t320 + t347 * t323 + t349 * t324, t314 * qJD(4) - t357 * t309 + t386 * t310 - t355 * t388, t309 (-t310 * t335 + t320 * t336 + (-t314 * t336 - t323 * t335) * qJD(5)) * pkin(5) + t377, t377; 0, -t353 * t317 + t348 * t318 + t347 * t321 + t349 * t322, t312 * qJD(4) - t357 * t307 + t386 * t308 - t352 * t388, t307 (-t308 * t335 + t318 * t336 + (-t312 * t336 - t321 * t335) * qJD(5)) * pkin(5) + t378, t378; 0 ((t348 * qJD(2) + t349) * t344 + (t353 * qJD(2) - t347) * t346) * t342, t326 * qJD(4) - t357 * t315 + t386 * t316 - t351 * t388, t315 (t336 * t372 - t316 * t335 + (-t326 * t336 + t335 * t379) * qJD(5)) * pkin(5) + t376, t376;];
JaD_transl  = t1;
