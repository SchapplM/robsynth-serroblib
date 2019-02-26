% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:26
% EndTime: 2019-02-26 20:53:26
% DurationCPUTime: 0.73s
% Computational Cost: add. (657->123), mult. (2123->219), div. (0->0), fcn. (2305->14), ass. (0->84)
t454 = sin(pkin(7));
t452 = sin(pkin(13));
t463 = cos(qJ(3));
t500 = cos(pkin(13));
t478 = qJD(3) * t500;
t460 = sin(qJ(3));
t487 = qJD(3) * t460;
t504 = t452 * t487 - t463 * t478;
t419 = t504 * t454;
t457 = cos(pkin(7));
t421 = t504 * t457;
t443 = -t463 * t452 - t460 * t500;
t424 = t443 * t454;
t426 = t443 * t457;
t458 = cos(pkin(6));
t453 = sin(pkin(12));
t464 = cos(qJ(1));
t492 = t464 * t453;
t456 = cos(pkin(12));
t461 = sin(qJ(1));
t493 = t461 * t456;
t473 = t458 * t493 + t492;
t429 = t473 * qJD(1);
t494 = t461 * t453;
t482 = t458 * t494;
t488 = qJD(1) * t464;
t430 = -qJD(1) * t482 + t456 * t488;
t491 = t464 * t456;
t434 = -t458 * t491 + t494;
t435 = t458 * t492 + t493;
t486 = qJD(3) * t463;
t439 = -t452 * t486 - t460 * t478;
t455 = sin(pkin(6));
t469 = -t460 * t452 + t463 * t500;
t489 = qJD(1) * t461;
t399 = t434 * t421 + t429 * t426 + t430 * t469 + t435 * t439 + (t419 * t464 - t424 * t489) * t455;
t481 = t455 * t489;
t414 = t429 * t454 + t457 * t481;
t459 = sin(qJ(5));
t462 = cos(qJ(5));
t510 = t399 * t459 - t414 * t462;
t509 = -t399 * t462 - t414 * t459;
t497 = t455 * t464;
t404 = -t424 * t497 - t434 * t426 - t435 * t469;
t415 = -t434 * t454 + t457 * t497;
t506 = t404 * t462 + t415 * t459;
t505 = -t404 * t459 + t415 * t462;
t427 = t434 * qJD(1);
t428 = t435 * qJD(1);
t437 = -t482 + t491;
t396 = t473 * t421 - t427 * t426 - t428 * t469 + t437 * t439 + (-t419 * t461 - t424 * t488) * t455;
t503 = t458 * t419 + (t421 * t456 - t439 * t453) * t455;
t502 = -r_i_i_C(3) - pkin(10);
t501 = pkin(9) + qJ(4);
t499 = t454 * t460;
t498 = t455 * t461;
t496 = t457 * t460;
t490 = pkin(3) * t499 + t501 * t457 + qJ(2);
t483 = pkin(3) * t486;
t485 = t457 * qJD(4) + t454 * t483 + qJD(2);
t484 = pkin(3) * t487;
t480 = t455 * t488;
t476 = -t459 * r_i_i_C(1) - t462 * r_i_i_C(2);
t474 = t462 * r_i_i_C(1) - t459 * r_i_i_C(2) + pkin(4);
t472 = qJD(5) * t476;
t468 = qJD(3) * t443;
t420 = t454 * t468;
t422 = t457 * t468;
t423 = t469 * t454;
t425 = t469 * t457;
t438 = t469 * qJD(3);
t465 = t420 * t497 + t434 * t422 - t423 * t481 + t429 * t425 - t430 * t443 + t435 * t438;
t451 = t463 * pkin(3) + pkin(2);
t441 = -t454 * qJD(4) + t457 * t483;
t433 = -t455 * t456 * t454 + t458 * t457;
t432 = pkin(3) * t496 - t501 * t454;
t417 = t454 * t473 + t457 * t498;
t412 = -t427 * t454 + t457 * t480;
t411 = -t458 * t424 + (-t426 * t456 + t453 * t469) * t455;
t406 = -t424 * t498 + t426 * t473 + t437 * t469;
t395 = -t473 * t422 + t427 * t425 - t428 * t443 - t437 * t438 + (t420 * t461 + t423 * t488) * t455;
t393 = t396 * t462 + t412 * t459 + (-t406 * t459 + t417 * t462) * qJD(5);
t392 = -t396 * t459 + t412 * t462 + (-t406 * t462 - t417 * t459) * qJD(5);
t1 = [t509 * r_i_i_C(1) + t510 * r_i_i_C(2) - t399 * pkin(4) - t430 * t451 + t435 * t484 + t429 * t432 + t434 * t441 - pkin(1) * t488 + t502 * t465 + (t505 * r_i_i_C(1) - t506 * r_i_i_C(2)) * qJD(5) + (t485 * t464 - t490 * t489) * t455, t480, -t502 * t396 + (t423 * t498 - t425 * t473 + t437 * t443) * t472 + t474 * t395 + (t428 * t460 + (t427 * t457 + t454 * t480) * t463 + (-t437 * t463 + (-t454 * t498 + t457 * t473) * t460) * qJD(3)) * pkin(3), t412, t392 * r_i_i_C(1) - t393 * r_i_i_C(2), 0; -t437 * t484 - pkin(1) * t489 + t396 * pkin(4) + t393 * r_i_i_C(1) + t392 * r_i_i_C(2) + t427 * t432 - t428 * t451 - t473 * t441 + t502 * t395 + (t485 * t461 + t490 * t488) * t455, t481, -t502 * t399 + (-t423 * t497 - t434 * t425 + t435 * t443) * t472 - t474 * t465 + (-t430 * t460 + (-t429 * t457 + t454 * t481) * t463 + (-t435 * t463 + (t434 * t457 + t454 * t497) * t460) * qJD(3)) * pkin(3), t414, -t510 * r_i_i_C(1) + t509 * r_i_i_C(2) + (t506 * r_i_i_C(1) + t505 * r_i_i_C(2)) * qJD(5), 0; 0, 0, t502 * t503 + (t458 * t423 + (t425 * t456 + t443 * t453) * t455) * t472 + t474 * (t458 * t420 + (t422 * t456 - t438 * t453) * t455) + (-t458 * t499 + (-t453 * t463 - t456 * t496) * t455) * qJD(3) * pkin(3), 0, -t476 * t503 + ((-t411 * t462 - t433 * t459) * r_i_i_C(1) + (t411 * t459 - t433 * t462) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
