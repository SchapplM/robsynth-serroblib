% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR11_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR11_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:51
% EndTime: 2019-02-26 22:21:52
% DurationCPUTime: 0.50s
% Computational Cost: add. (623->93), mult. (1848->153), div. (0->0), fcn. (1896->10), ass. (0->59)
t338 = cos(pkin(6));
t342 = sin(qJ(1));
t341 = sin(qJ(2));
t375 = t342 * t341;
t367 = t338 * t375;
t370 = qJD(2) * t341;
t345 = cos(qJ(2));
t346 = cos(qJ(1));
t372 = t346 * t345;
t313 = -qJD(1) * t367 - t342 * t370 + (qJD(2) * t338 + qJD(1)) * t372;
t337 = sin(pkin(6));
t376 = t337 * t346;
t388 = -qJD(3) * t376 + t313;
t373 = t346 * t341;
t374 = t342 * t345;
t325 = t338 * t373 + t374;
t371 = qJD(1) * t337;
t387 = -qJD(3) * t325 + t342 * t371;
t340 = sin(qJ(3));
t344 = cos(qJ(3));
t316 = t325 * t340 + t344 * t376;
t317 = t325 * t344 - t340 * t376;
t339 = sin(qJ(5));
t343 = cos(qJ(5));
t383 = t316 * t339 + t317 * t343;
t384 = t316 * t343 - t317 * t339;
t386 = (t383 * r_i_i_C(1) + t384 * r_i_i_C(2)) * qJD(5);
t322 = t337 * t341 * t340 - t338 * t344;
t377 = t337 * t344;
t323 = t338 * t340 + t341 * t377;
t385 = ((t322 * t339 + t323 * t343) * r_i_i_C(1) + (t322 * t343 - t323 * t339) * r_i_i_C(2)) * qJD(5);
t307 = t387 * t340 + t388 * t344;
t379 = pkin(3) + pkin(4);
t354 = t343 * r_i_i_C(1) - t339 * r_i_i_C(2) + t379;
t355 = t339 * r_i_i_C(1) + t343 * r_i_i_C(2) + qJ(4);
t380 = t355 * t340 + t354 * t344 + pkin(2);
t378 = t337 * t342;
t368 = r_i_i_C(3) + pkin(10) - pkin(9);
t365 = t346 * t371;
t364 = qJD(2) * t337 * t345;
t350 = t367 - t372;
t321 = t340 * t378 - t344 * t350;
t353 = t340 * t350 + t342 * t377;
t359 = -t321 * t339 - t343 * t353;
t358 = t321 * t343 - t339 * t353;
t352 = t338 * t372 - t375;
t351 = t338 * t374 + t373;
t306 = t388 * t340 - t387 * t344;
t347 = t340 * qJD(4) + ((-t339 * t344 + t340 * t343) * r_i_i_C(1) + (-t339 * t340 - t343 * t344) * r_i_i_C(2)) * qJD(5) + (-t354 * t340 + t355 * t344) * qJD(3);
t315 = -t322 * qJD(3) + t344 * t364;
t314 = t323 * qJD(3) + t340 * t364;
t312 = t351 * qJD(1) + t325 * qJD(2);
t311 = t325 * qJD(1) + t351 * qJD(2);
t310 = -t352 * qJD(1) + t350 * qJD(2);
t305 = t353 * qJD(3) - t311 * t344 + t340 * t365;
t304 = t321 * qJD(3) - t311 * t340 - t344 * t365;
t303 = t359 * qJD(5) + t304 * t339 + t305 * t343;
t302 = -t358 * qJD(5) + t304 * t343 - t305 * t339;
t1 = [-t313 * pkin(2) - t316 * qJD(4) - t355 * t306 - t354 * t307 + (-t384 * r_i_i_C(1) + t383 * r_i_i_C(2)) * qJD(5) + (-t346 * pkin(1) - pkin(8) * t378) * qJD(1) + t368 * t312, t310 * t380 + t368 * t311 - t347 * t351, t321 * qJD(4) + t355 * t305 - t354 * t304 + (t358 * r_i_i_C(1) + t359 * r_i_i_C(2)) * qJD(5), t304, t302 * r_i_i_C(1) - t303 * r_i_i_C(2), 0; -t311 * pkin(2) + t303 * r_i_i_C(1) + t302 * r_i_i_C(2) + t304 * qJ(4) - t353 * qJD(4) + t379 * t305 + (-pkin(1) * t342 + pkin(8) * t376) * qJD(1) + t368 * t310, -t312 * t380 - t368 * t313 + t347 * t352, t317 * qJD(4) - t354 * t306 + t355 * t307 + t386, t306 (t306 * t343 - t307 * t339) * r_i_i_C(1) + (-t306 * t339 - t307 * t343) * r_i_i_C(2) - t386, 0; 0 (-t380 * t370 + (-t368 * qJD(2) + t347) * t345) * t337, t323 * qJD(4) - t354 * t314 + t355 * t315 + t385, t314 (t314 * t343 - t315 * t339) * r_i_i_C(1) + (-t314 * t339 - t315 * t343) * r_i_i_C(2) - t385, 0;];
JaD_transl  = t1;
