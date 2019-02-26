% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:40
% EndTime: 2019-02-26 21:40:41
% DurationCPUTime: 0.75s
% Computational Cost: add. (1031->116), mult. (3020->197), div. (0->0), fcn. (3324->12), ass. (0->80)
t453 = sin(qJ(1));
t456 = cos(qJ(1));
t449 = cos(pkin(6));
t447 = sin(pkin(11));
t452 = sin(qJ(2));
t501 = cos(pkin(11));
t503 = cos(qJ(2));
t476 = t503 * t501;
t463 = -t452 * t447 + t476;
t461 = t449 * t463;
t478 = t452 * t501;
t464 = t503 * t447 + t478;
t412 = -t453 * t464 + t456 * t461;
t450 = sin(qJ(6));
t454 = cos(qJ(6));
t430 = t464 * t449;
t413 = t430 * t456 + t453 * t463;
t451 = sin(qJ(4));
t455 = cos(qJ(4));
t448 = sin(pkin(6));
t495 = t448 * t456;
t511 = t413 * t451 + t455 * t495;
t515 = -t412 * t454 + t450 * t511;
t514 = -t412 * t450 - t454 * t511;
t475 = t454 * r_i_i_C(1) - t450 * r_i_i_C(2);
t462 = t475 * qJD(6) + qJD(5);
t474 = -r_i_i_C(1) * t450 - r_i_i_C(2) * t454;
t472 = qJ(5) - t474;
t487 = r_i_i_C(3) + pkin(10) + pkin(4);
t513 = (t487 * t451 - t472 * t455) * qJD(4) - t462 * t451;
t510 = qJD(1) * t461 + qJD(2) * t463;
t489 = qJD(2) * t452;
t509 = -qJD(2) * t476 + t447 * t489;
t425 = t509 * t449;
t479 = t503 * qJD(2);
t432 = -qJD(2) * t478 - t447 * t479;
t491 = qJD(1) * t453;
t508 = t430 * t491 - t453 * t432 + (-qJD(1) * t463 + t425) * t456;
t482 = t448 * t491;
t507 = qJD(4) * t511 - t451 * t482 + t455 * t508;
t505 = t472 * t451 + t487 * t455 + pkin(3);
t504 = pkin(5) + pkin(9);
t502 = pkin(2) * t449;
t416 = -t453 * t430 + t456 * t463;
t497 = t416 * t451;
t496 = t448 * t453;
t493 = t452 * t453;
t492 = t452 * t456;
t490 = qJD(1) * t456;
t488 = qJD(4) * t455;
t486 = pkin(2) * t489;
t485 = t451 * t495;
t484 = t503 * t453;
t483 = t503 * t456;
t481 = t448 * t490;
t429 = t464 * t448;
t473 = t429 * t455 + t449 * t451;
t417 = t429 * t451 - t449 * t455;
t470 = t475 + t504;
t468 = t416 * t455 + t451 * t496;
t467 = qJD(6) * t474;
t393 = -qJD(4) * t485 + t413 * t488 - t451 * t508 - t455 * t482;
t399 = -t413 * qJD(1) + t453 * t425 + t432 * t456;
t446 = t503 * pkin(2) + pkin(1);
t433 = -t448 * qJD(3) + t479 * t502;
t431 = t452 * t502 + (-pkin(8) - qJ(3)) * t448;
t428 = t463 * t448;
t426 = qJD(2) * t430;
t424 = qJD(2) * t429;
t423 = t509 * t448;
t415 = -t453 * t461 - t456 * t464;
t408 = -t455 * t496 + t497;
t403 = t473 * qJD(4) - t423 * t451;
t401 = -t456 * t426 - t510 * t453 - t464 * t490;
t398 = -t453 * t426 + t510 * t456 - t464 * t491;
t392 = t399 * t455 - qJD(4) * t497 + (t451 * t490 + t453 * t488) * t448;
t391 = t468 * qJD(4) + t399 * t451 - t455 * t481;
t390 = t391 * t450 + t398 * t454 + (t408 * t454 + t415 * t450) * qJD(6);
t389 = t391 * t454 - t398 * t450 + (-t408 * t450 + t415 * t454) * qJD(6);
t1 = [t453 * t486 + t508 * pkin(3) - t511 * qJD(5) - t456 * t433 - t472 * t393 + t470 * t401 + (t514 * r_i_i_C(1) + t515 * r_i_i_C(2)) * qJD(6) + (t431 * t453 - t446 * t456) * qJD(1) + t487 * t507, t416 * t467 + t470 * t399 + ((t449 * t493 - t483) * qJD(2) + (-t449 * t483 + t493) * qJD(1)) * pkin(2) - t505 * t398 - t513 * t415, t481, -t487 * t391 + t472 * t392 + t462 * t468, t391, r_i_i_C(1) * t389 - t390 * r_i_i_C(2); -t456 * t486 + t399 * pkin(3) + t390 * r_i_i_C(1) + t389 * r_i_i_C(2) + t391 * qJ(5) + t408 * qJD(5) - t453 * t433 + t504 * t398 + (-t431 * t456 - t446 * t453) * qJD(1) + t487 * t392, t413 * t467 - t470 * t508 + ((-t449 * t492 - t484) * qJD(2) + (-t449 * t484 - t492) * qJD(1)) * pkin(2) + t505 * t401 - t513 * t412, t482, t462 * (t413 * t455 - t485) - t472 * t507 - t487 * t393, t393 (t393 * t454 + t401 * t450) * r_i_i_C(1) + (-t393 * t450 + t401 * t454) * r_i_i_C(2) + (-t515 * r_i_i_C(1) + t514 * r_i_i_C(2)) * qJD(6); 0, -t470 * t423 - t424 * t505 - t513 * t428 + t429 * t467 - t448 * t486, 0, t462 * t473 + t472 * (-t417 * qJD(4) - t423 * t455) - t487 * t403, t403 (t403 * t454 - t424 * t450) * r_i_i_C(1) + (-t403 * t450 - t424 * t454) * r_i_i_C(2) + ((-t417 * t450 + t428 * t454) * r_i_i_C(1) + (-t417 * t454 - t428 * t450) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
