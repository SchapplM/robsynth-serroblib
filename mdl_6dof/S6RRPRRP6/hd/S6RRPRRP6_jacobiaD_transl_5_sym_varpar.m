% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:52
% EndTime: 2019-02-26 21:48:53
% DurationCPUTime: 0.77s
% Computational Cost: add. (845->121), mult. (2486->210), div. (0->0), fcn. (2730->12), ass. (0->77)
t453 = cos(pkin(6));
t451 = sin(pkin(11));
t506 = cos(pkin(11));
t508 = cos(qJ(2));
t480 = t508 * t506;
t456 = sin(qJ(2));
t494 = qJD(2) * t456;
t512 = -qJD(2) * t480 + t451 * t494;
t431 = t512 * t453;
t483 = t456 * t506;
t469 = t508 * t451 + t483;
t436 = t469 * t453;
t484 = t508 * qJD(2);
t438 = -qJD(2) * t483 - t451 * t484;
t457 = sin(qJ(1));
t460 = cos(qJ(1));
t468 = -t456 * t451 + t480;
t496 = qJD(1) * t457;
t407 = t436 * t496 - t457 * t438 + (-qJD(1) * t468 + t431) * t460;
t455 = sin(qJ(4));
t459 = cos(qJ(4));
t477 = t460 * t436 + t457 * t468;
t452 = sin(pkin(6));
t487 = t452 * t496;
t500 = t452 * t460;
t490 = t459 * t500;
t401 = (-qJD(4) * t477 + t487) * t455 - qJD(4) * t490 - t407 * t459;
t432 = qJD(2) * t436;
t495 = qJD(1) * t460;
t467 = t468 * t453;
t513 = qJD(1) * t467 + t468 * qJD(2);
t408 = -t460 * t432 - t513 * t457 - t469 * t495;
t454 = sin(qJ(5));
t458 = cos(qJ(5));
t519 = t401 * t454 + t408 * t458;
t518 = -t401 * t458 + t408 * t454;
t414 = t455 * t500 - t459 * t477;
t418 = -t457 * t469 + t460 * t467;
t517 = -t414 * t454 + t418 * t458;
t516 = t414 * t458 + t418 * t454;
t472 = qJD(5) * (t454 * r_i_i_C(1) + t458 * r_i_i_C(2));
t475 = t458 * r_i_i_C(1) - t454 * r_i_i_C(2) + pkin(4);
t509 = r_i_i_C(3) + pkin(10);
t515 = (t475 * t455 - t509 * t459) * qJD(4) + t459 * t472;
t510 = t509 * t455 + t475 * t459 + pkin(3);
t507 = pkin(2) * t453;
t501 = t452 * t457;
t498 = t456 * t457;
t497 = t456 * t460;
t493 = qJD(5) * t454;
t492 = qJD(5) * t458;
t491 = pkin(2) * t494;
t489 = t508 * t457;
t488 = t508 * t460;
t486 = t452 * t495;
t435 = t469 * t452;
t424 = t435 * t459 + t453 * t455;
t478 = -t435 * t455 + t453 * t459;
t476 = -t457 * t436 + t460 * t468;
t473 = -t455 * t476 + t459 * t501;
t416 = t455 * t501 + t459 * t476;
t463 = -t477 * qJD(1) + t457 * t431 + t460 * t438;
t462 = t414 * qJD(4) + t407 * t455 + t459 * t487;
t450 = t508 * pkin(2) + pkin(1);
t439 = -t452 * qJD(3) + t484 * t507;
t437 = t456 * t507 + (-pkin(8) - qJ(3)) * t452;
t434 = t468 * t452;
t430 = qJD(2) * t435;
t429 = t512 * t452;
t421 = -t457 * t467 - t460 * t469;
t411 = t478 * qJD(4) - t429 * t459;
t405 = -t457 * t432 + t513 * t460 - t469 * t496;
t399 = t473 * qJD(4) + t455 * t486 + t459 * t463;
t398 = t416 * qJD(4) + t455 * t463 - t459 * t486;
t397 = t399 * t458 + t405 * t454 + (-t416 * t454 - t421 * t458) * qJD(5);
t396 = -t399 * t454 + t405 * t458 + (-t416 * t458 + t421 * t454) * qJD(5);
t1 = [t518 * r_i_i_C(1) + t519 * r_i_i_C(2) - t401 * pkin(4) + t407 * pkin(3) + t408 * pkin(9) + t457 * t491 - t460 * t439 + t509 * t462 + (t517 * r_i_i_C(1) - t516 * r_i_i_C(2)) * qJD(5) + (t457 * t437 - t460 * t450) * qJD(1) (t454 * t463 + t476 * t492) * r_i_i_C(1) + (t458 * t463 - t476 * t493) * r_i_i_C(2) + t463 * pkin(9) - t510 * t405 + ((t453 * t498 - t488) * qJD(2) + (-t453 * t488 + t498) * qJD(1)) * pkin(2) - t515 * t421, t486, -t475 * t398 + t509 * t399 - t473 * t472, t396 * r_i_i_C(1) - t397 * r_i_i_C(2), 0; -t460 * t491 + t463 * pkin(3) + t399 * pkin(4) + t405 * pkin(9) + t397 * r_i_i_C(1) + t396 * r_i_i_C(2) - t457 * t439 + t509 * t398 + (-t437 * t460 - t450 * t457) * qJD(1) (-t407 * t454 + t477 * t492) * r_i_i_C(1) + (-t407 * t458 - t477 * t493) * r_i_i_C(2) - t407 * pkin(9) + t510 * t408 + ((-t453 * t497 - t489) * qJD(2) + (-t453 * t489 - t497) * qJD(1)) * pkin(2) - t515 * t418, t487, t509 * t401 - (-t455 * t477 - t490) * t472 + t475 * t462, -t519 * r_i_i_C(1) + t518 * r_i_i_C(2) + (t516 * r_i_i_C(1) + t517 * r_i_i_C(2)) * qJD(5), 0; 0 (-t429 * t454 + t435 * t492) * r_i_i_C(1) + (-t429 * t458 - t435 * t493) * r_i_i_C(2) - t429 * pkin(9) - t452 * t491 - t510 * t430 - t515 * t434, 0, t509 * t411 - t478 * t472 + t475 * (-t424 * qJD(4) + t429 * t455) (-t411 * t454 + t430 * t458) * r_i_i_C(1) + (-t411 * t458 - t430 * t454) * r_i_i_C(2) + ((-t424 * t458 + t434 * t454) * r_i_i_C(1) + (t424 * t454 + t434 * t458) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
