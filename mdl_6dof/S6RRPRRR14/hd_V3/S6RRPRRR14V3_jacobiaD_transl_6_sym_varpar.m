% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14V3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14V3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_6_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:12
% EndTime: 2019-04-12 15:12:13
% DurationCPUTime: 1.00s
% Computational Cost: add. (479->149), mult. (1504->263), div. (0->0), fcn. (1518->10), ass. (0->94)
t459 = sin(qJ(4));
t460 = sin(qJ(2));
t458 = sin(qJ(5));
t464 = cos(qJ(4));
t537 = -qJD(5) * t464 + qJD(2);
t483 = t537 * t458;
t463 = cos(qJ(5));
t511 = qJD(4) * t463;
t545 = (t459 * t511 - t483) * t460;
t461 = sin(qJ(1));
t465 = cos(qJ(2));
t466 = cos(qJ(1));
t519 = t466 * t464;
t443 = t461 * t459 + t465 * t519;
t508 = qJD(4) * t466;
t490 = t464 * t508;
t512 = qJD(4) * t459;
t493 = t461 * t512;
t515 = qJD(2) * t460;
t496 = t461 * t515;
t429 = t443 * qJD(1) - t464 * t496 - t465 * t493 - t490;
t520 = t466 * t459;
t524 = t461 * t465;
t440 = t464 * t524 - t520;
t514 = qJD(2) * t465;
t517 = qJD(1) * t466;
t479 = t460 * t517 + t461 * t514;
t506 = qJD(5) * t460;
t419 = (-qJD(5) * t440 + t479) * t458 + (t461 * t506 + t429) * t463;
t510 = qJD(4) * t464;
t518 = qJD(1) * t461;
t428 = t510 * t524 - t464 * t518 + (t517 * t465 - t496 - t508) * t459;
t457 = sin(qJ(6));
t462 = cos(qJ(6));
t544 = t419 * t457 - t428 * t462;
t543 = -t419 * t462 - t428 * t457;
t521 = t465 * t463;
t526 = t460 * t464;
t480 = t458 * t526 + t521;
t481 = (t457 * r_i_i_C(1) + t462 * r_i_i_C(2)) * t459;
t541 = -t480 * r_i_i_C(3) + t465 * qJ(3) - t460 * t481;
t528 = t460 * t458;
t431 = t440 * t463 + t461 * t528;
t439 = t459 * t524 + t519;
t540 = t431 * t457 - t439 * t462;
t539 = t431 * t462 + t439 * t457;
t513 = qJD(2) * t466;
t478 = -t460 * t518 + t465 * t513;
t502 = qJD(6) * t462;
t503 = qJD(6) * t459;
t535 = (t457 * t510 + t459 * t502) * r_i_i_C(1) + (-t457 * t503 + t462 * t510) * r_i_i_C(2) - (t458 * t512 + t537 * t463) * r_i_i_C(3) + qJD(2) * qJ(3);
t534 = t458 * r_i_i_C(3);
t529 = t457 * t463;
t527 = t460 * t463;
t525 = t460 * t466;
t523 = t462 * t463;
t509 = qJD(4) * t465;
t507 = qJD(5) * t458;
t504 = qJD(6) * t457;
t501 = qJD(6) * t463;
t500 = t460 * qJD(3);
t499 = qJD(5) * t463 * r_i_i_C(3);
t492 = t460 * t513;
t488 = qJD(2) * t464 - qJD(5);
t425 = t488 * t521 - t545;
t484 = -t460 * t503 - t425;
t482 = t488 * t465;
t434 = t443 * t463 + t458 * t525;
t438 = -t465 * t458 + t463 * t526;
t475 = t457 * t501 + t462 * t507;
t474 = -t457 * t507 + t462 * t501;
t436 = t438 * t466;
t471 = -qJD(6) * t438 + t459 * t514 + t460 * t510;
t418 = -t431 * qJD(5) - t429 * t458 + t479 * t463;
t470 = r_i_i_C(3) * t507 + qJD(3) + (-t464 * t534 - t481) * qJD(2);
t469 = -t463 * t482 + t545;
t468 = t475 * r_i_i_C(1) + t474 * r_i_i_C(2) - t499;
t467 = -t535 * t460 + t470 * t465;
t442 = -t461 * t464 + t465 * t520;
t441 = t464 * t521 + t528;
t435 = t438 * t461;
t433 = -t443 * t458 + t463 * t525;
t430 = -t440 * t458 + t461 * t527;
t427 = (qJD(1) - t509) * t520 + (-t492 + (-qJD(1) * t465 + qJD(4)) * t461) * t464;
t426 = t439 * qJD(1) + t459 * t492 - t465 * t490 - t493;
t424 = t537 * t527 + (t460 * t512 - t482) * t458;
t423 = t465 * t483 + (-t459 * t509 - t488 * t460) * t463;
t422 = -qJD(1) * t436 + t469 * t461;
t421 = t438 * t518 + t469 * t466;
t417 = (t466 * t506 + t427) * t463 + (-qJD(5) * t443 + t478) * t458;
t416 = t434 * qJD(5) + t427 * t458 - t478 * t463;
t415 = t417 * t462 - t426 * t457 + (-t434 * t457 + t442 * t462) * qJD(6);
t414 = -t417 * t457 - t426 * t462 + (-t434 * t462 - t442 * t457) * qJD(6);
t1 = [t543 * r_i_i_C(1) + t544 * r_i_i_C(2) + t418 * r_i_i_C(3) - t461 * t500 + (t540 * r_i_i_C(1) + t539 * r_i_i_C(2)) * qJD(6) - t479 * qJ(3) (t421 * t462 + t436 * t504) * r_i_i_C(1) + (-t421 * t457 + t436 * t502) * r_i_i_C(2) - t541 * t518 + t467 * t466, t478 (t426 * t523 + t427 * t457 + t443 * t502) * r_i_i_C(1) + (-t426 * t529 + t427 * t462 - t443 * t504) * r_i_i_C(2) + t426 * t534 + t468 * t442, t417 * r_i_i_C(3) + (t416 * t457 - t433 * t502) * r_i_i_C(2) + (-t416 * t462 - t433 * t504) * r_i_i_C(1), t414 * r_i_i_C(1) - t415 * r_i_i_C(2); t415 * r_i_i_C(1) + t414 * r_i_i_C(2) + t416 * r_i_i_C(3) + t478 * qJ(3) + t466 * t500 (t422 * t462 + t435 * t504) * r_i_i_C(1) + (-t422 * t457 + t435 * t502) * r_i_i_C(2) + t541 * t517 + t467 * t461, t479 (-t428 * t523 + t429 * t457 + t440 * t502) * r_i_i_C(1) + (t428 * t529 + t429 * t462 - t440 * t504) * r_i_i_C(2) - t428 * t534 + t468 * t439, t419 * r_i_i_C(3) + (-t418 * t457 - t430 * t502) * r_i_i_C(2) + (t418 * t462 - t430 * t504) * r_i_i_C(1), -t544 * r_i_i_C(1) + t543 * r_i_i_C(2) + (-t539 * r_i_i_C(1) + t540 * r_i_i_C(2)) * qJD(6); 0 (t423 * t462 - t441 * t504) * r_i_i_C(1) + (-t423 * t457 - t441 * t502) * r_i_i_C(2) + t470 * t460 + t535 * t465, t515 ((t457 * t464 - t459 * t523) * r_i_i_C(1) + (t459 * t529 + t462 * t464) * r_i_i_C(2) - t459 * t534) * t514 + ((-qJD(4) * t534 + (-r_i_i_C(1) * t462 + r_i_i_C(2) * t457) * (-qJD(6) + t511)) * t464 + ((-qJD(4) * t457 + t475) * r_i_i_C(1) + (-qJD(4) * t462 + t474) * r_i_i_C(2) - t499) * t459) * t460, t425 * r_i_i_C(3) + (-t424 * t457 + t480 * t502) * r_i_i_C(2) + (t424 * t462 + t480 * t504) * r_i_i_C(1) (t471 * r_i_i_C(1) + t484 * r_i_i_C(2)) * t462 + (t484 * r_i_i_C(1) - t471 * r_i_i_C(2)) * t457;];
JaD_transl  = t1;
