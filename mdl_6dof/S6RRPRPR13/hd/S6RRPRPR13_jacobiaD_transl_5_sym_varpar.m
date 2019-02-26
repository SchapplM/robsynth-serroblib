% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR13_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR13_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:52
% EndTime: 2019-02-26 21:44:52
% DurationCPUTime: 0.27s
% Computational Cost: add. (378->66), mult. (1129->103), div. (0->0), fcn. (1106->10), ass. (0->48)
t334 = cos(pkin(6));
t335 = sin(qJ(4));
t338 = cos(qJ(4));
t332 = sin(pkin(6));
t339 = cos(qJ(2));
t365 = t332 * t339;
t370 = -t334 * t338 + t335 * t365;
t336 = sin(qJ(2));
t340 = cos(qJ(1));
t360 = t340 * t336;
t337 = sin(qJ(1));
t361 = t337 * t339;
t323 = t334 * t361 + t360;
t345 = t334 * t360 + t361;
t317 = qJD(1) * t323 + qJD(2) * t345;
t359 = t340 * t339;
t362 = t337 * t336;
t321 = -t334 * t359 + t362;
t364 = t332 * t340;
t348 = t321 * t338 + t335 * t364;
t358 = qJD(1) * t332;
t353 = t337 * t358;
t369 = qJD(4) * t348 + t317 * t335 + t338 * t353;
t331 = sin(pkin(11));
t333 = cos(pkin(11));
t350 = r_i_i_C(1) * t333 - r_i_i_C(2) * t331 + pkin(4);
t367 = r_i_i_C(3) + qJ(5);
t342 = t335 * t350 - t338 * t367 + qJ(3);
t368 = pkin(3) + pkin(8);
t366 = t332 * t337;
t357 = qJD(2) * t336;
t355 = t334 * t362;
t354 = t338 * t364;
t352 = t340 * t358;
t351 = t332 * t357;
t349 = r_i_i_C(1) * t331 + r_i_i_C(2) * t333 + pkin(2) + pkin(9);
t347 = t323 * t335 + t338 * t366;
t346 = t323 * t338 - t335 * t366;
t344 = t355 - t359;
t307 = -t317 * t338 - qJD(4) * t354 + (qJD(4) * t321 + t353) * t335;
t341 = -t338 * qJD(5) + qJD(3) + (t335 * t367 + t338 * t350) * qJD(4);
t320 = -qJD(4) * t370 - t338 * t351;
t318 = -qJD(1) * t355 - t337 * t357 + (qJD(2) * t334 + qJD(1)) * t359;
t316 = qJD(1) * t345 + qJD(2) * t323;
t315 = qJD(1) * t321 + qJD(2) * t344;
t312 = qJD(4) * t346 - t315 * t335 + t338 * t352;
t311 = qJD(4) * t347 + t315 * t338 + t335 * t352;
t1 = [t348 * qJD(5) - t317 * qJ(3) - t321 * qJD(3) - t350 * t369 - t349 * t318 - t367 * t307 + (-t340 * pkin(1) - t366 * t368) * qJD(1), t315 * t349 - t316 * t342 - t341 * t344, -t315, qJD(5) * t347 - t311 * t350 + t312 * t367, t311, 0; -t346 * qJD(5) - t315 * qJ(3) + t323 * qJD(3) + t350 * t312 - t349 * t316 + t367 * t311 + (-t337 * pkin(1) + t364 * t368) * qJD(1), -t317 * t349 + t318 * t342 + t341 * t345, t317 -(-t321 * t335 + t354) * qJD(5) + t367 * t369 - t350 * t307, t307, 0; 0 (t342 * t339 * qJD(2) + (-qJD(2) * t349 + t341) * t336) * t332, t351, -t370 * qJD(5) - t367 * (-t335 * t351 + (t334 * t335 + t338 * t365) * qJD(4)) - t350 * t320, t320, 0;];
JaD_transl  = t1;
