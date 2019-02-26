% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:52
% EndTime: 2019-02-26 21:19:52
% DurationCPUTime: 0.61s
% Computational Cost: add. (695->95), mult. (1898->163), div. (0->0), fcn. (2026->14), ass. (0->74)
t448 = sin(qJ(3));
t442 = sin(pkin(7));
t443 = sin(pkin(6));
t452 = cos(qJ(1));
t494 = t443 * t452;
t478 = t442 * t494;
t445 = cos(pkin(7));
t446 = cos(pkin(6));
t444 = cos(pkin(13));
t489 = t452 * t444;
t441 = sin(pkin(13));
t449 = sin(qJ(1));
t492 = t449 * t441;
t462 = t446 * t489 - t492;
t506 = t462 * t445;
t465 = -t506 + t478;
t510 = t465 * t448;
t479 = t446 * t492;
t485 = qJD(1) * t452;
t424 = -qJD(1) * t479 + t444 * t485;
t490 = t452 * t441;
t491 = t449 * t444;
t427 = t446 * t490 + t491;
t451 = cos(qJ(3));
t461 = t446 * t491 + t490;
t423 = t461 * qJD(1);
t495 = t443 * t449;
t476 = qJD(1) * t495;
t460 = -t423 * t445 + t442 * t476;
t468 = t451 * t478;
t401 = -qJD(3) * t468 + (-qJD(3) * t427 + t460) * t448 - (-qJD(3) * t506 - t424) * t451;
t439 = qJD(4) + qJD(5);
t505 = t442 * t462 + t445 * t494;
t509 = t505 * t439 - t401;
t408 = -t427 * t451 + t510;
t500 = t423 * t442;
t415 = t445 * t476 + t500;
t508 = -t408 * t439 - t415;
t493 = t444 * t445;
t496 = t442 * t446;
t504 = (-t441 * t448 + t451 * t493) * t443 + t451 * t496;
t417 = t443 * (t441 * t451 + t448 * t493) + t448 * t496;
t502 = r_i_i_C(3) + pkin(11) + pkin(10);
t501 = t417 * t439;
t440 = qJ(4) + qJ(5);
t437 = sin(t440);
t438 = cos(t440);
t429 = -t479 + t489;
t464 = t442 * t495 - t445 * t461;
t410 = t429 * t451 + t448 * t464;
t413 = t505 * qJD(1);
t470 = -t410 * t439 + t413;
t422 = t427 * qJD(1);
t456 = -t429 * t448 + t451 * t464;
t399 = qJD(1) * t510 + t456 * qJD(3) - t422 * t451;
t420 = t442 * t461 + t445 * t495;
t475 = t420 * t439 + t399;
t396 = -t475 * t437 + t470 * t438;
t397 = t470 * t437 + t475 * t438;
t488 = t396 * r_i_i_C(1) - t397 * r_i_i_C(2);
t487 = (t437 * t509 - t438 * t508) * r_i_i_C(1) + (t437 * t508 + t438 * t509) * r_i_i_C(2);
t411 = t504 * qJD(3);
t425 = -t443 * t444 * t442 + t446 * t445;
t469 = -t425 * t439 - t411;
t486 = (t469 * t437 - t438 * t501) * r_i_i_C(1) + (t437 * t501 + t469 * t438) * r_i_i_C(2);
t483 = t443 * qJD(2);
t450 = cos(qJ(4));
t436 = t450 * pkin(4) + pkin(3);
t466 = t438 * r_i_i_C(1) - t437 * r_i_i_C(2) + t436;
t447 = sin(qJ(4));
t457 = -qJD(4) * t447 * pkin(4) + (-t437 * r_i_i_C(1) - t438 * r_i_i_C(2)) * t439;
t454 = qJD(3) * t408 - t424 * t448 + t451 * t460;
t398 = qJD(3) * t410 - t422 * t448 + (t451 * t506 - t468) * qJD(1);
t1 = [-pkin(9) * t500 + t452 * t483 - t424 * pkin(2) - t401 * t436 + (r_i_i_C(1) * t509 + r_i_i_C(2) * t508) * t438 + (r_i_i_C(1) * t508 - r_i_i_C(2) * t509) * t437 + t502 * t454 + (-t452 * pkin(1) + (-pkin(9) * t445 - qJ(2)) * t495) * qJD(1) + (-t415 * t447 + (-t408 * t447 + t450 * t505) * qJD(4)) * pkin(4), t443 * t485, -t466 * t398 + t502 * t399 + t457 * t456 (-t399 * t447 + t413 * t450 + (-t410 * t450 - t420 * t447) * qJD(4)) * pkin(4) + t488, t488, 0; t449 * t483 - t422 * pkin(2) + t397 * r_i_i_C(1) + t396 * r_i_i_C(2) + t399 * t436 + t502 * t398 + (t413 * t447 + (-t410 * t447 + t420 * t450) * qJD(4)) * pkin(4) + (-t449 * pkin(1) + pkin(9) * t505 + qJ(2) * t494) * qJD(1), t476, t502 * t401 + t457 * (-t427 * t448 - t465 * t451) + t466 * t454 (-t401 * t447 + t415 * t450 + (t408 * t450 + t447 * t505) * qJD(4)) * pkin(4) + t487, t487, 0; 0, 0, -t466 * t417 * qJD(3) + t502 * t411 + t457 * t504 (-t411 * t447 + (-t417 * t450 - t425 * t447) * qJD(4)) * pkin(4) + t486, t486, 0;];
JaD_transl  = t1;
