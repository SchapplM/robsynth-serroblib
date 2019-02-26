% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:42:34
% EndTime: 2019-02-26 21:42:35
% DurationCPUTime: 0.34s
% Computational Cost: add. (520->78), mult. (1089->122), div. (0->0), fcn. (1071->12), ass. (0->53)
t344 = sin(qJ(1));
t341 = cos(pkin(6));
t356 = qJD(2) * t341 + qJD(1);
t343 = sin(qJ(2));
t373 = t344 * t343;
t364 = t341 * t373;
t368 = qJD(2) * t343;
t345 = cos(qJ(2));
t346 = cos(qJ(1));
t370 = t346 * t345;
t317 = -qJD(1) * t364 - t344 * t368 + t356 * t370;
t339 = sin(pkin(6));
t374 = t339 * t346;
t383 = -qJD(4) * t374 + t317;
t371 = t346 * t343;
t372 = t344 * t345;
t321 = t341 * t371 + t372;
t369 = qJD(1) * t339;
t382 = -qJD(4) * t321 + t344 * t369;
t336 = pkin(11) + qJ(4);
t334 = sin(t336);
t335 = cos(t336);
t337 = sin(pkin(12));
t340 = cos(pkin(12));
t355 = t340 * r_i_i_C(1) - t337 * r_i_i_C(2) + pkin(4);
t377 = r_i_i_C(3) + qJ(5);
t381 = (t355 * t334 - t377 * t335) * qJD(4) - t334 * qJD(5);
t380 = t382 * t334 + t383 * t335;
t333 = cos(pkin(11)) * pkin(3) + pkin(2);
t378 = t377 * t334 + t355 * t335 + t333;
t376 = t339 * t343;
t375 = t339 * t344;
t367 = qJD(2) * t345;
t363 = t341 * t370;
t362 = pkin(3) * sin(pkin(11)) + pkin(8);
t360 = t346 * t369;
t359 = t339 * t367;
t342 = -pkin(9) - qJ(3);
t354 = t337 * r_i_i_C(1) + t340 * r_i_i_C(2) - t342;
t323 = -t364 + t370;
t353 = -t323 * t334 + t335 * t375;
t352 = t323 * t335 + t334 * t375;
t351 = t341 * t334 + t335 * t376;
t350 = t341 * t372 + t371;
t310 = t383 * t334 - t382 * t335;
t320 = -t363 + t373;
t318 = t351 * qJD(4) + t334 * t359;
t316 = t350 * qJD(1) + t321 * qJD(2);
t315 = t321 * qJD(1) + t350 * qJD(2);
t314 = -qJD(1) * t363 - t346 * t367 + t356 * t373;
t309 = t353 * qJD(4) - t315 * t335 + t334 * t360;
t308 = t352 * qJD(4) - t315 * t334 - t335 * t360;
t1 = [(-t316 * t337 - t340 * t380) * r_i_i_C(1) + (-t316 * t340 + t337 * t380) * r_i_i_C(2) - t380 * pkin(4) - (t321 * t334 + t335 * t374) * qJD(5) - t317 * t333 + t316 * t342 - t320 * qJD(3) - t377 * t310 + (-t346 * pkin(1) - t362 * t375) * qJD(1), t323 * qJD(3) + t378 * t314 - t354 * t315 + t350 * t381, -t314, t352 * qJD(5) - t355 * t308 + t377 * t309, t308, 0; (t309 * t340 - t314 * t337) * r_i_i_C(1) + (-t309 * t337 - t314 * t340) * r_i_i_C(2) + t309 * pkin(4) - t353 * qJD(5) - t315 * t333 + t314 * t342 + t350 * qJD(3) + t377 * t308 + (-t344 * pkin(1) + t362 * t374) * qJD(1), t321 * qJD(3) - t316 * t378 + t354 * t317 + t381 * t320, t316 -(-t321 * t335 + t334 * t374) * qJD(5) + t377 * t380 - t355 * t310, t310, 0; 0 (t343 * qJD(3) - t381 * t345 + (-t343 * t378 + t354 * t345) * qJD(2)) * t339, t339 * t368, t351 * qJD(5) + t377 * (t335 * t359 + (-t334 * t376 + t335 * t341) * qJD(4)) - t355 * t318, t318, 0;];
JaD_transl  = t1;
