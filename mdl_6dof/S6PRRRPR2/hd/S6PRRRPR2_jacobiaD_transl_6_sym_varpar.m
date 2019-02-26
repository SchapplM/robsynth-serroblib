% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:09
% EndTime: 2019-02-26 20:11:09
% DurationCPUTime: 0.54s
% Computational Cost: add. (761->87), mult. (1153->149), div. (0->0), fcn. (1154->14), ass. (0->62)
t370 = pkin(12) + qJ(6);
t366 = sin(t370);
t367 = cos(t370);
t395 = t366 * r_i_i_C(1) + t367 * r_i_i_C(2);
t423 = t395 * qJD(6);
t417 = r_i_i_C(3) + pkin(10) + qJ(5);
t396 = t367 * r_i_i_C(1) - t366 * r_i_i_C(2);
t392 = cos(pkin(12)) * pkin(5) + pkin(4) + t396;
t372 = qJ(3) + qJ(4);
t368 = sin(t372);
t369 = cos(t372);
t371 = qJD(3) + qJD(4);
t377 = sin(qJ(3));
t382 = (t392 * t368 - t417 * t369) * t371 + qJD(3) * t377 * pkin(3) + t369 * t423 - t368 * qJD(5);
t378 = sin(qJ(2));
t380 = cos(qJ(2));
t374 = sin(pkin(11));
t416 = cos(pkin(6));
t401 = t374 * t416;
t415 = cos(pkin(11));
t358 = -t378 * t401 + t415 * t380;
t414 = t368 * t371;
t413 = t369 * t371;
t375 = sin(pkin(6));
t412 = t374 * t375;
t411 = t375 * t378;
t410 = t375 * t380;
t409 = qJD(2) * t378;
t405 = t368 * t411;
t404 = t369 * t411;
t403 = t375 * t409;
t402 = qJD(2) * t410;
t400 = t375 * t415;
t397 = t368 * t400;
t357 = t415 * t378 + t380 * t401;
t353 = t357 * qJD(2);
t394 = t371 * t412 - t353;
t393 = t416 * t415;
t391 = t380 * t393;
t389 = t396 * qJD(6);
t388 = t416 * t371 + t402;
t387 = sin(pkin(12)) * pkin(5) + pkin(9) + pkin(8) + t395;
t356 = t374 * t380 + t378 * t393;
t379 = cos(qJ(3));
t386 = -t379 * pkin(3) - t417 * t368 - t392 * t369 - pkin(2);
t351 = -qJD(2) * t391 + t374 * t409;
t335 = -t351 * t368 + t356 * t413 - t371 * t397;
t336 = -t356 * t414 + (-t371 * t400 - t351) * t369;
t346 = t356 * t369 - t397;
t385 = t346 * qJD(5) - t423 * (-t356 * t368 - t369 * t400) + t417 * t336 - t392 * t335;
t337 = t358 * t413 + t394 * t368;
t338 = -t358 * t414 + t394 * t369;
t348 = t358 * t369 + t368 * t412;
t384 = t348 * qJD(5) - t423 * (-t358 * t368 + t369 * t412) + t417 * t338 - t392 * t337;
t343 = t388 * t368 + t371 * t404;
t344 = t388 * t369 - t371 * t405;
t350 = t416 * t368 + t404;
t383 = t350 * qJD(5) - t423 * (t416 * t369 - t405) + t417 * t344 - t392 * t343;
t355 = t374 * t378 - t391;
t354 = t358 * qJD(2);
t352 = t356 * qJD(2);
t1 = [0, -t387 * t353 + t386 * t354 + t382 * t357 + t358 * t389 (t353 * t377 + (-t358 * t379 - t377 * t412) * qJD(3)) * pkin(3) + t384, t384, t337 (-t338 * t366 + t354 * t367) * r_i_i_C(1) + (-t338 * t367 - t354 * t366) * r_i_i_C(2) + ((-t348 * t367 - t357 * t366) * r_i_i_C(1) + (t348 * t366 - t357 * t367) * r_i_i_C(2)) * qJD(6); 0, -t387 * t351 + t386 * t352 + t382 * t355 + t356 * t389 (t351 * t377 + (-t356 * t379 + t377 * t400) * qJD(3)) * pkin(3) + t385, t385, t335 (-t336 * t366 + t352 * t367) * r_i_i_C(1) + (-t336 * t367 - t352 * t366) * r_i_i_C(2) + ((-t346 * t367 - t355 * t366) * r_i_i_C(1) + (t346 * t366 - t355 * t367) * r_i_i_C(2)) * qJD(6); 0 ((t386 * qJD(2) + t389) * t378 + (t387 * qJD(2) - t382) * t380) * t375 (-t377 * t402 + (-t416 * t377 - t379 * t411) * qJD(3)) * pkin(3) + t383, t383, t343 (-t344 * t366 + t367 * t403) * r_i_i_C(1) + (-t344 * t367 - t366 * t403) * r_i_i_C(2) + ((-t350 * t367 + t366 * t410) * r_i_i_C(1) + (t350 * t366 + t367 * t410) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
