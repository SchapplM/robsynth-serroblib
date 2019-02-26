% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:12
% EndTime: 2019-02-26 21:30:12
% DurationCPUTime: 0.36s
% Computational Cost: add. (366->80), mult. (1084->138), div. (0->0), fcn. (1136->10), ass. (0->54)
t318 = cos(pkin(6));
t315 = sin(pkin(11));
t317 = cos(pkin(11));
t320 = sin(qJ(2));
t323 = cos(qJ(2));
t334 = t323 * t315 + t320 * t317;
t304 = t334 * t318;
t308 = t320 * t315 - t323 * t317;
t327 = t308 * qJD(2);
t344 = qJD(2) * t323;
t339 = t318 * t344;
t345 = qJD(2) * t320;
t298 = t318 * t315 * t345 - t317 * t339;
t306 = -t315 * t344 - t317 * t345;
t321 = sin(qJ(1));
t324 = cos(qJ(1));
t346 = qJD(1) * t321;
t355 = t304 * t346 - t321 * t306 + (qJD(1) * t308 + t298) * t324;
t335 = t324 * t304 - t321 * t308;
t354 = t335 * qJD(1) - t321 * t298 - t324 * t306;
t319 = sin(qJ(5));
t322 = cos(qJ(5));
t337 = t322 * r_i_i_C(1) - t319 * r_i_i_C(2);
t353 = t337 * qJD(5) + qJD(4);
t316 = sin(pkin(6));
t352 = t316 * t321;
t351 = t316 * t324;
t350 = t320 * t321;
t349 = t320 * t324;
t348 = t321 * t323;
t347 = t323 * t324;
t343 = -r_i_i_C(3) - pkin(9) - pkin(3);
t342 = pkin(2) * t345;
t341 = t316 * t346;
t340 = qJD(1) * t351;
t333 = t319 * r_i_i_C(1) + t322 * r_i_i_C(2) + qJ(4);
t328 = t308 * t318;
t288 = -t321 * t334 - t324 * t328;
t332 = t288 * t319 + t322 * t351;
t331 = t288 * t322 - t319 * t351;
t329 = t334 * t316;
t326 = t321 * t328;
t314 = t323 * pkin(2) + pkin(1);
t307 = pkin(2) * t339 - t316 * qJD(3);
t305 = t318 * t320 * pkin(2) + (-pkin(8) - qJ(3)) * t316;
t302 = t308 * t316;
t299 = qJD(2) * t304;
t296 = qJD(2) * t329;
t290 = -t324 * t334 + t326;
t285 = qJD(1) * t326 + (-qJD(1) * t334 - t299) * t324 + t321 * t327;
t282 = t288 * qJD(1) - t321 * t299 - t324 * t327;
t280 = t322 * t340 + t282 * t319 + (-t290 * t322 - t319 * t352) * qJD(5);
t279 = -t319 * t340 + t282 * t322 + (t290 * t319 - t322 * t352) * qJD(5);
t1 = [t321 * t342 + t288 * qJD(4) - t324 * t307 + t333 * t285 + (t331 * r_i_i_C(1) - t332 * r_i_i_C(2)) * qJD(5) - t343 * t355 + (-t324 * t314 + (t305 + (-pkin(4) - t337) * t316) * t321) * qJD(1), -t353 * (t321 * t304 + t324 * t308) + t343 * t282 - t333 * t354 + ((t318 * t350 - t347) * qJD(2) + (-t318 * t347 + t350) * qJD(1)) * pkin(2), t340, t282, t279 * r_i_i_C(1) - t280 * r_i_i_C(2), 0; -t324 * t342 + t280 * r_i_i_C(1) + t279 * r_i_i_C(2) + t282 * qJ(4) - t290 * qJD(4) - t321 * t307 + t343 * t354 + (-t314 * t321 + (pkin(4) * t316 - t305) * t324) * qJD(1), t353 * t335 - t343 * t285 - t333 * t355 + ((-t318 * t349 - t348) * qJD(2) + (-t318 * t348 - t349) * qJD(1)) * pkin(2), t341, -t285 (-t285 * t322 - t319 * t341) * r_i_i_C(1) + (t285 * t319 - t322 * t341) * r_i_i_C(2) + (t332 * r_i_i_C(1) + t331 * r_i_i_C(2)) * qJD(5), 0; 0, t353 * t329 + t343 * t296 + (-t333 * t327 - t342) * t316, 0, t296, t337 * t296 + ((-t302 * t319 - t318 * t322) * r_i_i_C(1) + (-t302 * t322 + t318 * t319) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
