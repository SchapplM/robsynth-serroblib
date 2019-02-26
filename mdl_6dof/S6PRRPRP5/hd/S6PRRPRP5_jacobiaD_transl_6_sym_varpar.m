% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:33
% EndTime: 2019-02-26 20:03:33
% DurationCPUTime: 0.42s
% Computational Cost: add. (488->92), mult. (1522->154), div. (0->0), fcn. (1573->10), ass. (0->62)
t370 = sin(qJ(3));
t372 = cos(qJ(5));
t394 = -qJD(6) * t372 + qJD(4);
t404 = pkin(3) + pkin(9) + r_i_i_C(2);
t373 = cos(qJ(3));
t407 = qJD(3) * t373;
t375 = -qJ(4) * t407 + (t404 * qJD(3) - t394) * t370;
t371 = sin(qJ(2));
t374 = cos(qJ(2));
t367 = sin(pkin(10));
t413 = cos(pkin(6));
t398 = t367 * t413;
t412 = cos(pkin(10));
t359 = -t371 * t398 + t412 * t374;
t418 = (qJD(2) * t370 + qJD(5)) * t371 - t374 * t407;
t416 = pkin(4) + pkin(8);
t415 = -r_i_i_C(1) - pkin(5);
t414 = r_i_i_C(3) + qJ(6);
t368 = sin(pkin(6));
t411 = t368 * t370;
t410 = t368 * t373;
t409 = t368 * t374;
t408 = qJD(2) * t371;
t406 = qJD(5) * t370;
t369 = sin(qJ(5));
t405 = t369 * qJD(6);
t403 = t372 * t409;
t402 = t368 * t408;
t401 = qJD(2) * t409;
t397 = t368 * t412;
t392 = t370 * t397;
t389 = t413 * t412;
t385 = t374 * t389;
t352 = -qJD(2) * t385 + t367 * t408;
t356 = t367 * t371 - t385;
t391 = -t356 * t406 - t352;
t358 = t412 * t371 + t374 * t398;
t354 = t358 * qJD(2);
t390 = -t358 * t406 - t354;
t357 = t367 * t374 + t371 * t389;
t381 = -t357 * t370 - t373 * t397;
t388 = t356 * t372 - t369 * t381;
t384 = -t359 * t370 + t367 * t410;
t387 = t358 * t372 - t369 * t384;
t386 = (qJD(2) + t406) * t374;
t383 = t359 * t373 + t367 * t411;
t360 = t371 * t411 - t413 * t373;
t382 = t413 * t370 + t371 * t410;
t380 = -t370 * qJ(4) - t404 * t373 - pkin(2);
t379 = -t415 * t369 - t414 * t372 + qJ(4);
t353 = t357 * qJD(2);
t378 = qJD(5) * t357 + t353 * t370 + t356 * t407;
t355 = t359 * qJD(2);
t377 = qJD(5) * t359 + t355 * t370 + t358 * t407;
t376 = (t414 * t369 - t415 * t372) * qJD(5) + t394;
t350 = t382 * qJD(3) + t370 * t401;
t343 = t383 * qJD(3) - t354 * t370;
t341 = -qJD(3) * t392 - t352 * t370 + t357 * t407;
t335 = -t350 * t372 - qJD(5) * t403 + (qJD(5) * t360 + t402) * t369;
t329 = t387 * qJD(5) - t343 * t372 + t355 * t369;
t327 = t388 * qJD(5) - t341 * t372 + t353 * t369;
t1 = [0, t359 * t405 - t416 * t354 - t415 * (-t377 * t369 + t390 * t372) + t414 * (t390 * t369 + t377 * t372) + t375 * t358 + t380 * t355, -t404 * t343 + t379 * (t384 * qJD(3) - t354 * t373) + t376 * t383, t343, t387 * qJD(6) + t414 * (t343 * t369 + t355 * t372 + (-t358 * t369 - t372 * t384) * qJD(5)) + t415 * t329, t329; 0, t357 * t405 - t416 * t352 - t415 * (-t378 * t369 + t391 * t372) + t414 * (t391 * t369 + t378 * t372) + t375 * t356 + t380 * t353, -t404 * t341 + t379 * (t381 * qJD(3) - t352 * t373) + t376 * (t357 * t373 - t392) t341, t388 * qJD(6) + t414 * (t341 * t369 + t353 * t372 + (-t356 * t369 - t372 * t381) * qJD(5)) + t415 * t327, t327; 0 (-t415 * (-t418 * t369 + t372 * t386) + t414 * (t369 * t386 + t418 * t372) + t371 * t405 - t375 * t374 + (t380 * t371 + t416 * t374) * qJD(2)) * t368, -t404 * t350 + t379 * (-t360 * qJD(3) + t373 * t401) + t376 * t382, t350 -(-t360 * t369 + t403) * qJD(6) + t414 * (t372 * t402 + t350 * t369 + (t360 * t372 + t369 * t409) * qJD(5)) + t415 * t335, t335;];
JaD_transl  = t1;
