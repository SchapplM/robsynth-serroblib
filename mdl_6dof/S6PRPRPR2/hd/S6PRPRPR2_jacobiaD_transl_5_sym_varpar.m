% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:05
% EndTime: 2019-02-26 19:47:05
% DurationCPUTime: 0.24s
% Computational Cost: add. (297->50), mult. (960->98), div. (0->0), fcn. (1004->12), ass. (0->46)
t334 = sin(pkin(11));
t338 = cos(pkin(11));
t344 = cos(qJ(2));
t358 = qJD(2) * t344;
t342 = sin(qJ(2));
t359 = qJD(2) * t342;
t367 = t334 * t359 - t338 * t358;
t340 = cos(pkin(6));
t352 = t334 * t344 + t338 * t342;
t323 = t352 * t340;
t335 = sin(pkin(10));
t339 = cos(pkin(10));
t351 = t334 * t342 - t338 * t344;
t312 = t323 * t339 - t335 * t351;
t343 = cos(qJ(4));
t336 = sin(pkin(6));
t341 = sin(qJ(4));
t362 = t336 * t341;
t366 = -t312 * t343 + t339 * t362;
t365 = r_i_i_C(3) + qJ(5);
t364 = pkin(2) * qJD(2);
t361 = t336 * t343;
t360 = t340 * t342;
t318 = t367 * t340;
t325 = -t334 * t358 - t338 * t359;
t354 = t318 * t339 - t325 * t335;
t309 = t318 * t335 + t325 * t339;
t321 = t352 * t336;
t353 = t321 * t343 + t340 * t341;
t333 = sin(pkin(12));
t337 = cos(pkin(12));
t350 = r_i_i_C(1) * t337 - r_i_i_C(2) * t333 + pkin(4);
t349 = -r_i_i_C(1) * t333 - r_i_i_C(2) * t337 - pkin(8);
t314 = -t323 * t335 - t339 * t351;
t348 = t314 * t343 + t335 * t362;
t347 = qJD(2) * t352;
t346 = t341 * t365 + t343 * t350 + pkin(3);
t345 = t341 * qJD(5) + (-t341 * t350 + t343 * t365) * qJD(4);
t324 = t351 * qJD(2);
t322 = t351 * t340;
t319 = t340 * t347;
t316 = t367 * t336;
t303 = qJD(4) * t353 - t316 * t341;
t301 = qJD(4) * t348 + t309 * t341;
t299 = -qJD(4) * t366 - t341 * t354;
t1 = [0 (t335 * t360 - t339 * t344) * t364 - t349 * t309 + t345 * (t322 * t335 - t339 * t352) + t346 * (t319 * t335 + t324 * t339) 0, t348 * qJD(5) + t365 * (t309 * t343 + (-t314 * t341 + t335 * t361) * qJD(4)) - t350 * t301, t301, 0; 0 (-t335 * t344 - t339 * t360) * t364 + t349 * t354 + t345 * (-t322 * t339 - t335 * t352) + t346 * (-t319 * t339 + t324 * t335) 0, -t366 * qJD(5) + t365 * (-t354 * t343 + (-t312 * t341 - t339 * t361) * qJD(4)) - t350 * t299, t299, 0; 0, t349 * t316 + (-pkin(2) * t359 - t345 * t351 - t346 * t347) * t336, 0, t353 * qJD(5) + t365 * (-t316 * t343 + (-t321 * t341 + t340 * t343) * qJD(4)) - t350 * t303, t303, 0;];
JaD_transl  = t1;
