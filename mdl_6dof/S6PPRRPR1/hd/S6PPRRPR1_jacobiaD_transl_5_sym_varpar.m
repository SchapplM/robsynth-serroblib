% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:30
% EndTime: 2019-02-26 19:40:31
% DurationCPUTime: 0.27s
% Computational Cost: add. (411->56), mult. (1327->109), div. (0->0), fcn. (1510->14), ass. (0->54)
t360 = sin(pkin(12));
t363 = sin(pkin(6));
t370 = sin(qJ(3));
t372 = cos(qJ(3));
t365 = cos(pkin(12));
t367 = cos(pkin(7));
t390 = t365 * t367;
t362 = sin(pkin(7));
t368 = cos(pkin(6));
t393 = t362 * t368;
t346 = (t360 * t372 + t370 * t390) * t363 + t370 * t393;
t366 = cos(pkin(11));
t361 = sin(pkin(11));
t395 = t361 * t368;
t355 = -t360 * t395 + t366 * t365;
t354 = -t366 * t360 - t365 * t395;
t394 = t362 * t363;
t376 = t354 * t367 + t361 * t394;
t342 = t355 * t372 + t376 * t370;
t389 = t366 * t368;
t353 = t360 * t389 + t361 * t365;
t352 = -t361 * t360 + t365 * t389;
t386 = t366 * t394;
t377 = -t352 * t367 + t386;
t400 = -t353 * t372 + t377 * t370;
t399 = r_i_i_C(3) + qJ(5);
t398 = t353 * t370;
t392 = t363 * t365;
t391 = t363 * t367;
t388 = qJD(3) * t370;
t387 = qJD(3) * t372;
t384 = t362 * t387;
t383 = t367 * t387;
t347 = -t352 * t362 - t366 * t391;
t369 = sin(qJ(4));
t371 = cos(qJ(4));
t382 = t347 * t369 - t371 * t400;
t348 = -t354 * t362 + t361 * t391;
t381 = t342 * t371 + t348 * t369;
t351 = -t362 * t392 + t368 * t367;
t380 = t346 * t371 + t351 * t369;
t359 = sin(pkin(13));
t364 = cos(pkin(13));
t379 = t364 * r_i_i_C(1) - t359 * r_i_i_C(2) + pkin(4);
t378 = -t359 * r_i_i_C(1) - t364 * r_i_i_C(2) - pkin(9);
t374 = qJD(3) * (t399 * t369 + t379 * t371 + pkin(3));
t373 = t369 * qJD(5) + (-t379 * t369 + t399 * t371) * qJD(4);
t343 = t363 * t360 * t388 - t368 * t384 - t383 * t392;
t337 = -t361 * t363 * t384 - t354 * t383 + t355 * t388;
t335 = -t352 * t383 + (t372 * t386 + t398) * qJD(3);
t333 = t380 * qJD(4) - t343 * t369;
t331 = t381 * qJD(4) - t337 * t369;
t329 = t382 * qJD(4) - t335 * t369;
t1 = [0, 0, t378 * t337 + t373 * (-t355 * t370 + t376 * t372) - t342 * t374, t381 * qJD(5) + t399 * (-t337 * t371 + (-t342 * t369 + t348 * t371) * qJD(4)) - t379 * t331, t331, 0; 0, 0, t378 * t335 + t373 * (-t377 * t372 - t398) + t400 * t374, t382 * qJD(5) + t399 * (-t335 * t371 + (t347 * t371 + t369 * t400) * qJD(4)) - t379 * t329, t329, 0; 0, 0, t378 * t343 + t373 * (t372 * t393 + (-t360 * t370 + t372 * t390) * t363) - t346 * t374, t380 * qJD(5) + t399 * (-t343 * t371 + (-t346 * t369 + t351 * t371) * qJD(4)) - t379 * t333, t333, 0;];
JaD_transl  = t1;
