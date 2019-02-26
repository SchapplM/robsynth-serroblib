% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:49
% EndTime: 2019-02-26 21:26:50
% DurationCPUTime: 0.48s
% Computational Cost: add. (368->86), mult. (1171->142), div. (0->0), fcn. (1132->8), ass. (0->61)
t336 = sin(qJ(2));
t333 = sin(pkin(9));
t334 = cos(pkin(9));
t385 = pkin(4) + pkin(3);
t388 = qJ(4) * t333 + t334 * t385 + pkin(2);
t339 = cos(qJ(2));
t370 = pkin(8) + r_i_i_C(2) - qJ(3);
t391 = t370 * t339;
t394 = t388 * t336 + t391;
t371 = t336 * qJD(3);
t393 = (t336 * pkin(2) + t391) * qJD(2) - t371;
t335 = sin(qJ(5));
t338 = cos(qJ(5));
t353 = t333 * t335 + t334 * t338;
t390 = qJD(5) * t353;
t381 = t334 * t335;
t354 = t333 * t338 - t381;
t389 = qJD(5) * t354;
t337 = sin(qJ(1));
t375 = qJD(2) * t336;
t367 = t337 * t375;
t340 = cos(qJ(1));
t379 = t340 * t333;
t369 = t339 * t379;
t377 = qJD(1) * t337;
t322 = qJD(1) * t369 - t333 * t367 - t334 * t377;
t378 = t340 * t334;
t327 = t337 * t333 + t339 * t378;
t323 = qJD(1) * t327 - t334 * t367;
t380 = t337 * t339;
t324 = t333 * t380 + t378;
t325 = t334 * t380 - t379;
t357 = t324 * t335 + t325 * t338;
t306 = qJD(5) * t357 - t322 * t338 + t323 * t335;
t358 = t324 * t338 - t325 * t335;
t387 = qJD(5) * t358 + t322 * t335 + t323 * t338;
t386 = -t333 * qJD(4) + qJD(6) * t354;
t384 = -r_i_i_C(1) - pkin(5);
t382 = r_i_i_C(3) + qJ(6);
t376 = qJD(1) * t340;
t374 = qJD(2) * t339;
t373 = qJD(2) * t340;
t368 = t333 * t374;
t366 = t336 * t373;
t365 = t339 * t373;
t362 = t370 * t336;
t326 = -t337 * t334 + t369;
t356 = t326 * t338 - t327 * t335;
t355 = t326 * t335 + t327 * t338;
t352 = qJD(1) * t354;
t351 = qJD(2) * t354;
t350 = qJD(2) * t353;
t344 = -pkin(2) * t339 - pkin(1) + t362;
t343 = t339 * t350;
t341 = t339 * qJD(3) + t386 * t336 + (-t339 * t388 + t362) * qJD(2);
t321 = -qJD(1) * t325 - t334 * t366;
t320 = qJD(1) * t324 + t333 * t366;
t316 = t336 * t390 - t338 * t368 + t374 * t381;
t305 = qJD(5) * t356 - t320 * t335 + t321 * t338;
t304 = qJD(5) * t355 + t320 * t338 + t321 * t335;
t1 = [t358 * qJD(6) - t322 * qJ(4) - t324 * qJD(4) - t385 * t323 + t384 * t387 - t382 * t306 + t393 * t337 + (-t337 * pkin(7) + t340 * t344) * qJD(1), -t384 * (-t340 * t343 + (-t340 * t389 + t353 * t377) * t336) + t382 * (t354 * t365 + (-t337 * t352 - t340 * t390) * t336) + t394 * t377 + t341 * t340, -t336 * t377 + t365, -t320, qJD(6) * t355 + t304 * t384 + t305 * t382, t304; -t356 * qJD(6) - t320 * qJ(4) + t326 * qJD(4) + t385 * t321 - t384 * t305 + t382 * t304 - t393 * t340 + (t340 * pkin(7) + t337 * t344) * qJD(1), -t384 * (-t337 * t343 + (-t337 * t389 - t353 * t376) * t336) - t382 * (-t351 * t380 + (t337 * t390 - t340 * t352) * t336) - t394 * t376 + t341 * t337, t336 * t376 + t337 * t374, t322, t357 * qJD(6) + t384 * t306 + t382 * t387, t306; 0, t371 - t384 * (-t336 * t350 + t339 * t389) - t382 * (-t336 * t351 - t339 * t390) - t386 * t339 - t394 * qJD(2), t375, t368, t353 * qJD(6) * t336 + t382 * (t336 * t389 + t353 * t374) + t384 * t316, t316;];
JaD_transl  = t1;
